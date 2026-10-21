import argparse
import json
import shlex
import subprocess
import sys
import tomllib
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import NoReturn

import diff

REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST_PATH = REPO_ROOT / "corpus" / "manifest.json"
DATA_DIR = REPO_ROOT / "corpus" / "data"
OUT_DIR = REPO_ROOT / "harness" / "out"
EXPLAINED_PATH = REPO_ROOT / "harness" / "explained.toml"
GEMMI_DUMP_COMMAND = [sys.executable, "harness/gemmi_dump.py"]
DEFAULT_OURS_COMMAND = ["cargo", "run", "-q", "--release", "-p", "mmcif-dump", "--"]
DEFAULT_DUMPER_CRATE = REPO_ROOT / "tools" / "mmcif-dump"
EXPLAINED_KEYS = {"entry", "category", "match", "reason"}
STDERR_TAIL_LINES = 5
USAGE_ERROR_STATUS = 2


@dataclass(frozen=True)
class Rule:
    entry: str
    category: str
    match: str | None
    reason: str


@dataclass(frozen=True)
class EntryResult:
    entry: str
    unexplained: dict[str, list[str]]
    explained_count: int
    used_rules: frozenset[int]


def load_rules(path: Path) -> list[Rule]:
    with open(path, "rb") as handle:
        document = tomllib.load(handle)
    rules = []
    for table in document.get("explained", []):
        unknown = set(table) - EXPLAINED_KEYS
        missing = {"entry", "category", "reason"} - set(table)
        if unknown or missing:
            raise ValueError(f"{path}: bad [[explained]] table {table} (unknown {sorted(unknown)}, missing {sorted(missing)})")
        rules.append(Rule(table["entry"].lower(), table["category"], table.get("match"), table["reason"]))
    return rules


def rule_covers(rule: Rule, entry: str, category: str, message: str) -> bool:
    if rule.entry != entry or rule.category != category:
        return False
    return rule.match is None or rule.match in message


def apply_rules(entry: str, differences: diff.Differences, rules: list[Rule]) -> EntryResult:
    unexplained: dict[str, list[str]] = {}
    explained_count = 0
    used: set[int] = set()
    for category, messages in differences.items():
        for message in messages:
            covering = [i for i, rule in enumerate(rules) if rule_covers(rule, entry, category, message)]
            if covering:
                explained_count += 1
                used.update(covering)
            else:
                unexplained.setdefault(category, []).append(message)
    return EntryResult(entry, unexplained, explained_count, frozenset(used))


def fail(message: str) -> NoReturn:
    print(message, file=sys.stderr)
    raise SystemExit(USAGE_ERROR_STATUS)


def tail(text: str) -> str:
    return " | ".join(text.strip().splitlines()[-STDERR_TAIL_LINES:])


def run_dump(command: list[str], source: Path, destination: Path) -> str | None:
    relative_source = source.relative_to(REPO_ROOT)
    with open(destination, "wb") as output:
        completed = subprocess.run(
            command + [str(relative_source)], cwd=REPO_ROOT, stdout=output, stderr=subprocess.PIPE, text=True
        )
    if completed.returncode != 0:
        return f"exit {completed.returncode}: {tail(completed.stderr)}"
    return None


def failure_result(entry: str, category: str, message: str) -> EntryResult:
    return EntryResult(entry, {category: [message]}, 0, frozenset())


def process_entry(entry: str, ours_command: list[str], tolerance: float, rules: list[Rule]) -> EntryResult:
    source = DATA_DIR / f"{entry}.cif"
    gemmi_path = OUT_DIR / f"{entry}.gemmi.json"
    ours_path = OUT_DIR / f"{entry}.ours.json"
    gemmi_error = run_dump(GEMMI_DUMP_COMMAND, source, gemmi_path)
    if gemmi_error:
        return failure_result(entry, "gemmi_failed", gemmi_error)
    ours_error = run_dump(ours_command, source, ours_path)
    if ours_error:
        return failure_result(entry, "ours_failed", ours_error)
    try:
        differences = diff.compare(diff.load_dump(str(ours_path)), diff.load_dump(str(gemmi_path)), tolerance)
    except (json.JSONDecodeError, KeyError, TypeError) as error:
        return failure_result(entry, "ours_invalid_dump", f"{type(error).__name__}: {error}")
    return apply_rules(entry, differences, rules)


def select_entries(requested: list[str] | None) -> list[str]:
    manifest = json.loads(MANIFEST_PATH.read_text())
    if not requested:
        return sorted(manifest)
    wanted = [identifier.lower() for identifier in requested]
    unknown = [identifier for identifier in wanted if identifier not in manifest]
    if unknown:
        fail(f"ids not in {MANIFEST_PATH.relative_to(REPO_ROOT)}: {' '.join(unknown)}")
    return sorted(set(wanted))


def check_inputs(entries: list[str], ours_command: list[str]) -> None:
    missing = [entry for entry in entries if not (DATA_DIR / f"{entry}.cif").exists()]
    if missing:
        fail(f"missing corpus data for {' '.join(missing)}; run corpus/fetch.py first")
    if ours_command == DEFAULT_OURS_COMMAND and not DEFAULT_DUMPER_CRATE.is_dir():
        fail(
            f"Rust dumper crate not found at {DEFAULT_DUMPER_CRATE.relative_to(REPO_ROOT)}; "
            "create the mmcif-dump crate or pass --ours-cmd '<command prefix>' (the cif path is appended)"
        )


def category_summary(unexplained: dict[str, list[str]]) -> str:
    return ", ".join(f"{category}={len(messages)}" for category, messages in sorted(unexplained.items()))


def print_table(failing: list[EntryResult]) -> None:
    width = max(len(result.entry) for result in failing)
    print(f"{'entry'.ljust(width)}  unexplained  categories")
    for result in failing:
        total = sum(len(messages) for messages in result.unexplained.values())
        print(f"{result.entry.ljust(width)}  {total:11d}  {category_summary(result.unexplained)}")


def print_details(failing: list[EntryResult]) -> None:
    for result in failing:
        print(f"\n== {result.entry}")
        print(diff.format_report(result.unexplained))


def warn_unused_rules(results: list[EntryResult], rules: list[Rule]) -> None:
    used = set().union(*(result.used_rules for result in results)) if results else set()
    for index, rule in enumerate(rules):
        if index not in used:
            print(f"note: explained rule never matched: {rule.entry} {rule.category} {rule.match!r}", file=sys.stderr)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ids", nargs="*")
    parser.add_argument("--ours-cmd", help="command prefix that prints our JSON dump for the cif path appended to it")
    parser.add_argument("--tolerance", type=float, default=diff.DEFAULT_TOLERANCE)
    parser.add_argument("--explained", type=Path, default=EXPLAINED_PATH)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--details", action="store_true", help="print example differences per failing entry")
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    ours_command = shlex.split(arguments.ours_cmd) if arguments.ours_cmd else DEFAULT_OURS_COMMAND
    entries = select_entries(arguments.ids)
    check_inputs(entries, ours_command)
    rules = load_rules(arguments.explained)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=arguments.jobs) as pool:
        results = list(pool.map(lambda entry: process_entry(entry, ours_command, arguments.tolerance, rules), entries))
    failing = [result for result in results if result.unexplained]
    explained_total = sum(result.explained_count for result in results)
    print(f"{len(entries)} entries, {len(entries) - len(failing)} clean, {len(failing)} with unexplained differences, {explained_total} differences explained")
    if failing:
        print_table(failing)
        if arguments.details:
            print_details(failing)
    warn_unused_rules(results, rules)
    return 1 if failing else 0


if __name__ == "__main__":
    sys.exit(main())

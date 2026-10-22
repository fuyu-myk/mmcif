import argparse
import importlib.metadata
import json
import shlex
import subprocess
import sys
import tomllib
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any, NoReturn

import diff

REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST_PATH = REPO_ROOT / "corpus" / "manifest.json"
CORPUS_PATH = REPO_ROOT / "corpus" / "full.json"
DATA_DIR = REPO_ROOT / "corpus" / "data"
OUT_DIR = REPO_ROOT / "harness" / "out"
EXPLAINED_PATH = REPO_ROOT / "harness" / "explained.toml"
GEMMI_DUMP_COMMAND = [sys.executable, "harness/gemmi_dump.py"]
DUMPER_BUILD_COMMAND = ["cargo", "build", "--release", "-p", "mmcif-dump", "--message-format=json-render-diagnostics"]
DUMPER_CRATE = REPO_ROOT / "tools" / "mmcif-dump"
EXPLAINED_KEYS = {"entry", "example", "category", "match", "reason"}
CORPORA = ("seed", "full")
DEFAULT_TIMEOUT_SECONDS = 300.0
STDERR_TAIL_LINES = 5
USAGE_ERROR_STATUS = 2


@dataclass(frozen=True)
class Rule:
    entry: str | None
    example: str | None
    category: str
    match: str | None
    reason: str


@dataclass(frozen=True)
class EntryResult:
    entry: str
    unexplained: dict[str, list[str]]
    explained_count: int
    explained_by_rule: dict[int, int]


def parse_rule(path: Path, table: dict[str, str]) -> Rule:
    unknown = set(table) - EXPLAINED_KEYS
    required = {"category", "reason"} | ({"match"} if "example" in table else set())
    missing = required - set(table)
    empty = sorted(key for key, value in table.items() if not value)
    scoped = ("entry" in table) != ("example" in table)
    if unknown or missing or empty or not scoped:
        raise ValueError(
            f"{path}: bad [[explained]] table {table} (unknown {sorted(unknown)}, missing {sorted(missing)}, "
            f"empty {empty}; name exactly one of entry, for one entry, or example, for a rule over every entry "
            "that also needs a match)"
        )
    entry = table.get("entry")
    example = table.get("example")
    return Rule(
        entry.lower() if entry is not None else None,
        example.lower() if example is not None else None,
        table["category"],
        table.get("match"),
        table["reason"],
    )


def load_rules(path: Path) -> list[Rule]:
    with open(path, "rb") as handle:
        document = tomllib.load(handle)
    return [parse_rule(path, table) for table in document.get("explained", [])]


def rule_covers(rule: Rule, entry: str, category: str, message: str) -> bool:
    if rule.entry not in (None, entry) or rule.category != category:
        return False
    return rule.match is None or rule.match in message


def apply_rules(entry: str, differences: diff.Differences, rules: list[Rule]) -> EntryResult:
    unexplained: dict[str, list[str]] = {}
    explained_count = 0
    by_rule: dict[int, int] = defaultdict(int)
    for category, messages in differences.items():
        for message in messages:
            covering = [i for i, rule in enumerate(rules) if rule_covers(rule, entry, category, message)]
            if covering:
                explained_count += 1
                for index in covering:
                    by_rule[index] += 1
            else:
                unexplained.setdefault(category, []).append(message)
    return EntryResult(entry, unexplained, explained_count, dict(by_rule))


def as_bcif_differences(differences: diff.Differences) -> diff.Differences:
    return {
        "bcif_" + category.replace("_ours", "_bcif").replace("_gemmi", "_cif"): [
            message.replace("only in ours", "only in bcif").replace("only in gemmi", "only in cif") for message in messages
        ]
        for category, messages in differences.items()
    }


def fail(message: str) -> NoReturn:
    print(message, file=sys.stderr)
    raise SystemExit(USAGE_ERROR_STATUS)


def tail(text: str | bytes | None) -> str:
    if isinstance(text, bytes):
        text = text.decode(errors="replace")
    return " | ".join((text or "").strip().splitlines()[-STDERR_TAIL_LINES:])


def run_dump(command: list[str], source: Path, destination: Path, timeout: float) -> tuple[str, str] | None:
    relative_source = source.relative_to(REPO_ROOT)
    try:
        with open(destination, "wb") as output:
            completed = subprocess.run(
                command + [str(relative_source)],
                cwd=REPO_ROOT,
                stdout=output,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout,
            )
    except subprocess.TimeoutExpired as expired:
        stderr = tail(expired.stderr)
        return "timeout", f"no dump within {timeout:g} s" + (f": {stderr}" if stderr else "")
    except OSError as problem:
        return "failed", f"{type(problem).__name__}: {problem}"
    if completed.returncode != 0:
        return "failed", f"exit {completed.returncode}: {tail(completed.stderr)}"
    return None


def load_side(side: str, command: list[str], source: Path, timeout: float, differences: diff.Differences) -> dict[str, Any] | None:
    destination = OUT_DIR / f"{source.stem}.{side}.json"
    error = run_dump(command, source, destination, timeout)
    if error:
        differences[f"{side}_{error[0]}"] = [error[1]]
        return None
    try:
        return diff.load_dump(str(destination))
    except (ValueError, OSError) as problem:
        differences[f"{side}_invalid_dump"] = [f"{type(problem).__name__}: {problem}"]
        return None


def compare_into(differences: diff.Differences, side: str, ours: dict[str, Any], theirs: dict[str, Any], tolerance: float, binary: bool) -> None:
    try:
        found = diff.compare(ours, theirs, tolerance)
    except (KeyError, TypeError, IndexError, ValueError) as problem:
        differences[f"{side}_invalid_dump"] = [f"{type(problem).__name__}: {problem}"]
        return
    differences.update(as_bcif_differences(found) if binary else found)


def process_entry(entry: str, ours_command: list[str], tolerance: float, rules: list[Rule], timeout: float) -> EntryResult:
    cif = DATA_DIR / f"{entry}.cif"
    differences: diff.Differences = {}
    gemmi = load_side("gemmi", GEMMI_DUMP_COMMAND, cif, timeout, differences)
    ours = load_side("ours", ours_command, cif, timeout, differences)
    ours_bcif = load_side("ours_bcif", ours_command, DATA_DIR / f"{entry}.bcif", timeout, differences)
    if ours is not None and gemmi is not None:
        compare_into(differences, "ours", ours, gemmi, tolerance, binary=False)
    if ours is not None and ours_bcif is not None:
        compare_into(differences, "ours_bcif", ours_bcif, ours, tolerance, binary=True)
    if not differences:
        for side in ("gemmi", "ours", "ours_bcif"):
            (OUT_DIR / f"{entry}.{side}.json").unlink(missing_ok=True)
    return apply_rules(entry, differences, rules)


def corpus_entries(corpus: str) -> list[str]:
    if corpus == "full":
        return sorted(json.loads(CORPUS_PATH.read_text())["pins"])
    return sorted(json.loads(MANIFEST_PATH.read_text()))


def select_entries(corpus: str, requested: list[str] | None) -> list[str]:
    available = corpus_entries(corpus)
    if not requested:
        return available
    wanted = sorted({identifier.lower() for identifier in requested})
    unknown = [identifier for identifier in wanted if identifier not in set(available)]
    if unknown:
        fail(f"ids not in the {corpus} corpus: {' '.join(unknown)}")
    return wanted


def check_inputs(entries: list[str]) -> None:
    missing = [f"{entry}.{extension}" for entry in entries for extension in ("cif", "bcif") if not (DATA_DIR / f"{entry}.{extension}").exists()]
    if missing:
        fail(f"missing corpus data for {' '.join(missing)}; run corpus/fetch.py first")


def build_dumper() -> list[str]:
    if not DUMPER_CRATE.is_dir():
        fail(
            f"Rust dumper crate not found at {DUMPER_CRATE.relative_to(REPO_ROOT)}; "
            "pass --ours-cmd '<command prefix>' (the structure path is appended)"
        )
    completed = subprocess.run(DUMPER_BUILD_COMMAND, cwd=REPO_ROOT, stdout=subprocess.PIPE, text=True)
    if completed.returncode != 0:
        fail(f"{' '.join(DUMPER_BUILD_COMMAND)} failed")
    executables = [
        message["executable"]
        for message in map(json.loads, completed.stdout.splitlines())
        if message.get("reason") == "compiler-artifact" and message.get("executable") and message["target"]["name"] == "mmcif-dump"
    ]
    if not executables or not Path(executables[-1]).is_file():
        fail(f"{' '.join(DUMPER_BUILD_COMMAND)} reported no mmcif-dump executable")
    return [executables[-1]]


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


def describe_rule(rule: Rule) -> str:
    scope = rule.entry or f"every entry (example {rule.example})"
    return f"{scope} {rule.category} {rule.match!r}"


def warn_unused_rules(results: list[EntryResult], rules: list[Rule]) -> None:
    used = {index for result in results for index in result.explained_by_rule}
    for index, rule in enumerate(rules):
        if index not in used:
            print(f"note: explained rule never matched: {describe_rule(rule)}", file=sys.stderr)


def entry_record(result: EntryResult) -> dict[str, Any]:
    record: dict[str, Any] = {}
    if result.explained_count:
        record["explained"] = result.explained_count
    if result.unexplained:
        record["unexplained"] = {
            category: {"count": len(messages), "examples": messages[: diff.EXAMPLES_PER_CATEGORY]}
            for category, messages in result.unexplained.items()
        }
    return record


def rule_record(index: int, rule: Rule, results: list[EntryResult]) -> dict[str, Any]:
    fields = {"entry": rule.entry, "example": rule.example, "category": rule.category, "match": rule.match, "reason": rule.reason}
    record = {key: value for key, value in fields.items() if value is not None}
    record["explained"] = sum(result.explained_by_rule.get(index, 0) for result in results)
    record["entries"] = sorted(result.entry for result in results if index in result.explained_by_rule)
    return record


def results_document(results: list[EntryResult], rules: list[Rule], header: dict[str, Any]) -> dict[str, Any]:
    return header | {
        "entries": {result.entry: entry_record(result) for result in sorted(results, key=lambda result: result.entry)},
        "rules": [rule_record(index, rule, results) for index, rule in enumerate(rules)],
    }


def source_version() -> str:
    completed = subprocess.run(
        ["git", "describe", "--always", "--dirty", "--abbrev=12"], cwd=REPO_ROOT, capture_output=True, text=True
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unknown"


def default_results_path(corpus: str, ids: list[str] | None) -> Path:
    return OUT_DIR / f"results-{corpus}{'-subset' if ids else ''}.json"


def gemmi_version() -> str:
    try:
        return importlib.metadata.version("gemmi")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", choices=CORPORA, default="seed", help="seed: corpus/manifest.json; full: every entry pinned in corpus/full.json")
    parser.add_argument("--ids", nargs="*")
    parser.add_argument("--ours-cmd", help="command prefix that prints our JSON dump for the structure path appended to it")
    parser.add_argument("--tolerance", type=float, default=diff.DEFAULT_TOLERANCE)
    parser.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT_SECONDS, help="seconds each dump may take")
    parser.add_argument("--explained", type=Path, default=EXPLAINED_PATH)
    parser.add_argument("--results", type=Path, help="where to write the results file (default harness/out/results-<corpus>.json, or results-<corpus>-subset.json with --ids)")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--details", action="store_true", help="print example differences per failing entry")
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    entries = select_entries(arguments.corpus, arguments.ids)
    check_inputs(entries)
    ours_command = shlex.split(arguments.ours_cmd) if arguments.ours_cmd else build_dumper()
    rules = load_rules(arguments.explained)
    results_path = arguments.results or default_results_path(arguments.corpus, arguments.ids)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results_path.parent.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=arguments.jobs) as pool:
        results = list(
            pool.map(lambda entry: process_entry(entry, ours_command, arguments.tolerance, rules, arguments.timeout), entries)
        )
    header = {
        "corpus": arguments.corpus,
        "subset": bool(arguments.ids),
        "source": source_version(),
        "gemmi": gemmi_version(),
        "tolerance": arguments.tolerance,
        "timeout_seconds": arguments.timeout,
    }
    results_path.write_text(json.dumps(results_document(results, rules, header), indent=1, sort_keys=True) + "\n")
    failing = [result for result in results if result.unexplained]
    explained_total = sum(result.explained_count for result in results)
    print(f"{len(entries)} entries, {len(entries) - len(failing)} clean, {len(failing)} with unexplained differences, {explained_total} differences explained")
    if failing:
        print_table(failing)
        if arguments.details:
            print_details(failing)
    warn_unused_rules(results, rules)
    print(f"results written to {results_path}")
    return 1 if failing else 0


if __name__ == "__main__":
    sys.exit(main())

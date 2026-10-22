import argparse
import json
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_PATH = REPO_ROOT / "harness" / "out" / "results-full.json"
CORPUS_PATH = REPO_ROOT / "corpus" / "full.json"
MANIFEST_PATH = REPO_ROOT / "corpus" / "manifest.json"
REPORT_PATH = REPO_ROOT / "corpus" / "REPORT.md"
CLASS_LABELS = {
    "insertion_codes": "Insertion codes",
    "altlocs": "Alternate conformations",
    "microheterogeneity": "Microheterogeneity",
    "numbering_divergence": "Numbering divergence",
    "assembly_operators": "Assembly operators",
    "branched": "Branched entities",
    "multichar_chain_ids": "Multi-character chain ids",
    "multiple_models": "Multiple models",
}
OUTCOMES = ("agree", "explained", "unexplained")
FAILURE_SUFFIXES = ("_failed", "_timeout", "_invalid_dump")
BINARY_PREFIX = "bcif_"
USAGE_ERROR_STATUS = 2


def check_results(results: dict[str, Any], corpus: dict[str, Any]) -> None:
    if results.get("corpus") != "full" or results.get("subset"):
        raise ValueError("the report needs a complete run of harness/run.py --corpus full")
    if set(results["entries"]) != set(corpus["pins"]):
        raise ValueError("the results do not cover the entries pinned in corpus/full.json; run harness/run.py --corpus full again")


def plural(count: int, noun: str, nouns: str) -> str:
    return f"{count:,} {noun if count == 1 else nouns}"


def outcome(record: dict[str, Any]) -> str:
    if "unexplained" in record:
        return "unexplained"
    return "explained" if record else "agree"


def tally(results: dict[str, Any], ids: list[str]) -> dict[str, int]:
    counts = {"entries": len(ids)} | {name: 0 for name in OUTCOMES}
    for identifier in ids:
        counts[outcome(results["entries"][identifier])] += 1
    return counts


def entry_classes(corpus: dict[str, Any], seed: dict[str, list[str]]) -> dict[str, set[str]]:
    return {
        identifier: set(corpus["entries"].get(identifier, {})) | (set(seed.get(identifier, [])) & set(CLASS_LABELS))
        for identifier in sorted(corpus["pins"])
    }


def class_rows(results: dict[str, Any], corpus: dict[str, Any], seed: dict[str, list[str]]) -> list[tuple[str, dict[str, int]]]:
    classes = entry_classes(corpus, seed)
    groups = [("Random baseline", corpus["baseline"]), ("Hand-picked seed", sorted(seed))]
    groups += [(label, [identifier for identifier, found in classes.items() if name in found]) for name, label in CLASS_LABELS.items()]
    groups.append(("All entries", sorted(corpus["pins"])))
    return [(label, tally(results, ids)) for label, ids in groups]


def entries_with(results: dict[str, Any], test: Any) -> set[str]:
    unexplained = {
        identifier
        for identifier, record in results["entries"].items()
        if any(test(category) for category in record.get("unexplained", {}))
    }
    explained = {identifier for rule in results["rules"] if test(rule["category"]) for identifier in rule["entries"]}
    return unexplained | explained


def binary_sentence(results: dict[str, Any]) -> str:
    differing = entries_with(results, lambda category: category.startswith(BINARY_PREFIX))
    if not differing:
        return "No entry differs between the two."
    return f"{plural(len(differing), 'entry differs', 'entries differ')} between the two; see below."


def failure_sentence(results: dict[str, Any]) -> str:
    failed = entries_with(results, lambda category: category.endswith(FAILURE_SUFFIXES))
    if not failed:
        return "No reading failed or ran out of time."
    return f"{plural(len(failed), 'entry', 'entries')} could not be read by one side or ran out of time; see below."


def summary(results: dict[str, Any], total: dict[str, int]) -> list[str]:
    return [
        "## Summary",
        "",
        f"{plural(total['entries'], 'entry was', 'entries were')} read by mmcif-model and by gemmi {results['gemmi']}, "
        "and their atoms, residues, entities and assemblies were compared. "
        f"{total['agree']:,} agree exactly, {total['explained']:,} differ only in ways explained below, "
        f"and {total['unexplained']:,} have differences that are not explained yet.",
        "",
        f"Each entry was also read from its BinaryCIF file by mmcif-model and compared with its CIF reading. {binary_sentence(results)}",
        "",
        failure_sentence(results),
    ]


def corpus_section(corpus: dict[str, Any], seed: dict[str, list[str]]) -> list[str]:
    lines = [
        "## The corpus",
        "",
        f"- {len(corpus['baseline']):,} entries drawn at random from the {corpus['holdings']['count']:,} current PDB entries, "
        f"with seed {corpus['random_seed']}. The list they were drawn from has SHA-256 `{corpus['holdings']['sha256']}`.",
        f"- The {len(seed)} hand-picked entries in `corpus/manifest.json`.",
        f"- For each pathology class, up to {corpus['class_size']} entries: the hand-picked members first, "
        "then the hardest candidates by a measure for that class.",
    ]
    for name, pool in sorted(corpus["pools"].items()):
        lines.append(f"- Candidates for {CLASS_LABELS[name].lower()} come from {pool['source']} with `{pool['query']}`.")
    lines.append("- Every file is pinned in `corpus/full.json` by SHA-256, with its PDB revision.")
    return lines


def run_section(results: dict[str, Any]) -> list[str]:
    return [
        "## The run",
        "",
        "| Setting | Value |",
        "| --- | --- |",
        f"| Source | `{results['source']}` |",
        f"| gemmi | {results['gemmi']} |",
        f"| Coordinate tolerance | {results['tolerance']:g} Å |",
        f"| Time limit per reading | {results['timeout_seconds']:g} s |",
    ]


def class_section(rows: list[tuple[str, dict[str, int]]]) -> list[str]:
    lines = [
        "## Results by class",
        "",
        "An entry counts under every class it belongs to.",
        "",
        "| Class | Entries | Agree | Explained | Unexplained |",
        "| --- | --- | --- | --- | --- |",
    ]
    for label, counts in rows:
        lines.append(f"| {label} | {counts['entries']:,} | {counts['agree']:,} | {counts['explained']:,} | {counts['unexplained']:,} |")
    return lines


def rule_lines(rule: dict[str, Any]) -> list[str]:
    scope = rule.get("entry") or f"any entry (checked by hand on {rule['example']})"
    lines = [f"### {rule['category']} in {scope}", "", rule["reason"], ""]
    if "match" in rule:
        lines.append(f"- Covers messages containing `{rule['match']}`.")
    lines.append(
        f"- {plural(rule['explained'], 'difference', 'differences')} in {plural(len(rule['entries']), 'entry', 'entries')}: "
        + ", ".join(rule["entries"])
    )
    return lines


def explained_section(results: dict[str, Any]) -> list[str]:
    lines = ["## Explained differences", "", "Each rule in `harness/explained.toml` that matched in this run."]
    for rule in results["rules"]:
        if rule["explained"]:
            lines += [""] + rule_lines(rule)
    return lines


def unexplained_section(results: dict[str, Any]) -> list[str]:
    lines = ["## Unexplained differences", ""]
    open_entries = [(identifier, record["unexplained"]) for identifier, record in sorted(results["entries"].items()) if "unexplained" in record]
    if not open_entries:
        return lines + ["No differences are unexplained."]
    for identifier, categories in open_entries:
        described = "; ".join(
            f"{category} ({found['count']:,}), for example `{found['examples'][0]}`" for category, found in sorted(categories.items())
        )
        lines.append(f"- **{identifier}**: {described}")
    return lines


def render(results: dict[str, Any], corpus: dict[str, Any], seed: dict[str, list[str]]) -> str:
    rows = class_rows(results, corpus, seed)
    sections = [
        [
            "# Comparison with gemmi",
            "",
            "This file is generated by `harness/report.py` from a run of `harness/run.py --corpus full`. Do not edit it by hand.",
        ],
        summary(results, rows[-1][1]),
        corpus_section(corpus, seed),
        run_section(results),
        class_section(rows),
        explained_section(results),
        unexplained_section(results),
    ]
    return "\n\n".join("\n".join(section) for section in sections) + "\n"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=RESULTS_PATH)
    parser.add_argument("--output", type=Path, default=REPORT_PATH)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    if not arguments.results.exists():
        print(f"{arguments.results} not found; run harness/run.py --corpus full first", file=sys.stderr)
        return USAGE_ERROR_STATUS
    results = json.loads(arguments.results.read_text())
    corpus = json.loads(CORPUS_PATH.read_text())
    seed = json.loads(MANIFEST_PATH.read_text())
    try:
        check_results(results, corpus)
    except ValueError as problem:
        print(f"{arguments.results}: {problem}", file=sys.stderr)
        return USAGE_ERROR_STATUS
    arguments.output.write_text(render(results, corpus, seed))
    print(f"report written to {arguments.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

import argparse
import hashlib
import json
import random
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Callable

from gemmi import cif

import classify
import fetch
from classify import column

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = classify.DATA_DIR
CORPUS_PATH = REPO_ROOT / "corpus" / "full.json"
HOLDINGS_SNAPSHOT = DATA_DIR / "holdings.json"
POOL_SNAPSHOT = DATA_DIR / "microheterogeneity.json"
HOLDINGS_URL = "https://data.rcsb.org/rest/v1/holdings/current/entry_ids"
PDBJ_SQL_URL = "https://pdbj.org/rest/newweb/search/sql"
MICROHETEROGENEITY_QUERY = "SELECT DISTINCT pdbid FROM pdbj.entity_poly_seq WHERE hetero = 'y' ORDER BY pdbid"
USER_AGENT = "mmcif-corpus-select/1.0 (differential test corpus)"
REQUEST_TIMEOUT_SECONDS = 120
DEFAULT_SEED = 20261021
BASELINE_SIZE = 3000
CLASS_SIZE = 50


def normalized(ids: list[str]) -> list[str]:
    return sorted({identifier.lower() for identifier in ids})


def fingerprint(ids: list[str]) -> str:
    return hashlib.sha256("\n".join(normalized(ids)).encode()).hexdigest()


def draw_baseline(holdings: list[str], seed: int, size: int) -> list[str]:
    return sorted(random.Random(seed).sample(normalized(holdings), size))


def pdbj_ids(response: dict) -> list[str]:
    return normalized([row[0] for row in response["results"]])


def get_json(url: str) -> object:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
        return json.load(response)


def snapshot(path: Path, url: str, refresh: bool) -> object:
    if refresh or not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        partial = path.with_suffix(path.suffix + ".part")
        partial.write_text(json.dumps(get_json(url)))
        partial.replace(path)
        print(f"downloaded {path.relative_to(REPO_ROOT)}")
    else:
        print(f"using {path.relative_to(REPO_ROOT)}")
    return json.loads(path.read_text())


def distinct_residues(block: cif.Block, flag_tag: str) -> int:
    rows = zip(
        column(block, "_atom_site.label_asym_id"),
        column(block, "_atom_site.auth_seq_id"),
        column(block, "_atom_site.pdbx_PDB_ins_code"),
        column(block, flag_tag),
    )
    return len({row[:3] for row in rows if row[3]})


def insertion_code_severity(block: cif.Block) -> int:
    return distinct_residues(block, "_atom_site.pdbx_PDB_ins_code")


def altloc_severity(block: cif.Block) -> int:
    return distinct_residues(block, "_atom_site.label_alt_id")


def microheterogeneity_severity(block: cif.Block) -> int:
    rows = zip(
        column(block, "_entity_poly_seq.entity_id"),
        column(block, "_entity_poly_seq.num"),
        column(block, "_entity_poly_seq.hetero"),
    )
    declared = {(entity, number) for entity, number, hetero in rows if hetero.lower().startswith("y")}
    return max(len(declared), classify.mixed_residue_positions(block))


def numbering_shift(label: str, auth: str) -> int | str:
    try:
        return int(auth) - int(label)
    except ValueError:
        return f"{label}:{auth}"


def numbering_severity(block: cif.Block) -> int:
    rows = zip(
        column(block, "_atom_site.label_asym_id"),
        column(block, "_atom_site.label_seq_id"),
        column(block, "_atom_site.auth_seq_id"),
    )
    return len({(asym, numbering_shift(label, auth)) for asym, label, auth in rows if label and label != auth})


def assembly_severity(block: cif.Block) -> int:
    expressions = column(block, "_pdbx_struct_assembly_gen.oper_expression")
    return max((classify.operator_count(expression) for expression in expressions), default=0)


def branched_severity(block: cif.Block) -> int:
    residues = len(column(block, "_pdbx_branch_scheme.mon_id"))
    entities = len(column(block, "_pdbx_entity_branch.entity_id"))
    return max(residues, entities)


def multichar_severity(block: cif.Block) -> int:
    identifiers = set(column(block, "_atom_site.label_asym_id") + column(block, "_atom_site.auth_asym_id"))
    return sum(1 for identifier in identifiers if len(identifier) > 1)


def model_severity(block: cif.Block) -> int:
    return len(set(column(block, "_atom_site.pdbx_PDB_model_num")))


SEVERITY: dict[str, Callable[[cif.Block], int]] = {
    "insertion_codes": insertion_code_severity,
    "altlocs": altloc_severity,
    "microheterogeneity": microheterogeneity_severity,
    "numbering_divergence": numbering_severity,
    "assembly_operators": assembly_severity,
    "branched": branched_severity,
    "multichar_chain_ids": multichar_severity,
    "multiple_models": model_severity,
}


def severities(block: cif.Block) -> dict[str, int]:
    return {
        name: SEVERITY[name](block)
        for name, predicate in classify.PATHOLOGY_PREDICATES.items()
        if predicate(block)
    }


def classify_file(identifier: str) -> tuple[str, dict[str, int] | str]:
    path = DATA_DIR / f"{identifier}.cif"
    if not path.exists():
        return identifier, "not downloaded"
    try:
        return identifier, severities(cif.read(str(path)).sole_block())
    except (RuntimeError, ValueError, IndexError) as error:
        return identifier, f"{type(error).__name__}: {error}"


def revision_key(row: tuple[str, str, str, str]) -> tuple[int, int]:
    return int(row[1]), int(row[2])


def latest_revision(block: cif.Block) -> dict[str, str]:
    rows = zip(
        column(block, "_pdbx_audit_revision_history.data_content_type"),
        column(block, "_pdbx_audit_revision_history.major_revision"),
        column(block, "_pdbx_audit_revision_history.minor_revision"),
        column(block, "_pdbx_audit_revision_history.revision_date"),
    )
    models = [row for row in rows if row[0] == "Structure model"]
    if not models:
        return {}
    _, major, minor, date = max(models, key=revision_key)
    return {"revision": f"{major}.{minor}", "revision_date": date}


def pin_entry(identifier: str) -> tuple[str, dict[str, str] | str]:
    paths = {extension: DATA_DIR / f"{identifier}.{extension}" for extension, _, _ in fetch.SOURCES}
    for path in paths.values():
        if not path.exists():
            return identifier, f"{path.name} not downloaded"
    try:
        revision = latest_revision(cif.read(str(paths["cif"])).sole_block())
    except (RuntimeError, ValueError, IndexError) as error:
        return identifier, f"{type(error).__name__}: {error}"
    return identifier, {extension: fetch.sha256_of(path) for extension, path in paths.items()} | revision


def run_ids(corpus: dict, seed: dict[str, list[str]]) -> list[str]:
    chosen = {identifier for members in corpus["classes"].values() for identifier in members}
    return sorted(set(seed) | set(corpus["baseline"]) | chosen)


def fill_classes(seed: dict[str, list[str]], candidates: dict[str, dict[str, int]], size: int) -> dict[str, list[str]]:
    filled = {}
    for name in classify.PATHOLOGY_PREDICATES:
        chosen = sorted(identifier for identifier, classes in seed.items() if name in classes)
        ranked = sorted(
            (identifier for identifier, found in candidates.items() if name in found and identifier not in seed),
            key=lambda identifier: (-candidates[identifier][name], identifier),
        )
        filled[name] = chosen + ranked[: max(size - len(chosen), 0)]
    return filled


def write_corpus(corpus: dict) -> None:
    CORPUS_PATH.write_text(json.dumps(corpus, indent=1, sort_keys=True) + "\n")


def draw(arguments: argparse.Namespace) -> int:
    holdings = snapshot(HOLDINGS_SNAPSHOT, HOLDINGS_URL, arguments.refresh)
    query = urllib.parse.urlencode({"format": "json", "q": MICROHETEROGENEITY_QUERY})
    pool = snapshot(POOL_SNAPSHOT, f"{PDBJ_SQL_URL}?{query}", arguments.refresh)
    current = set(normalized(holdings))
    corpus = {
        "random_seed": arguments.seed,
        "holdings": {"count": len(current), "sha256": fingerprint(holdings)},
        "baseline": draw_baseline(holdings, arguments.seed, arguments.size),
        "pools": {
            "microheterogeneity": {
                "source": PDBJ_SQL_URL,
                "query": MICROHETEROGENEITY_QUERY,
                "ids": [identifier for identifier in pdbj_ids(pool) if identifier in current],
            }
        },
    }
    write_corpus(corpus)
    print(f"{len(corpus['baseline'])} baseline and {len(corpus['pools']['microheterogeneity']['ids'])} pooled ids written to {CORPUS_PATH.relative_to(REPO_ROOT)}")
    return 0


def report_failures(results: dict[str, object], action: str) -> bool:
    failures = {identifier: result for identifier, result in results.items() if isinstance(result, str)}
    for identifier, reason in sorted(failures.items()):
        print(f"FAILED {identifier}: {reason}")
    if failures:
        print(f"{len(failures)} of {len(results)} entries could not be {action}; run corpus/fetch.py --corpus corpus/full.json first")
    return bool(failures)


def classify_corpus(arguments: argparse.Namespace) -> int:
    corpus = json.loads(CORPUS_PATH.read_text())
    seed = json.loads(classify.MANIFEST_PATH.read_text())
    with ProcessPoolExecutor(arguments.jobs) as pool:
        results = dict(pool.map(classify_file, fetch.candidate_ids(corpus), chunksize=16))
        if report_failures(results, "classified"):
            return 1
        corpus["entries"] = results
        corpus["class_size"] = arguments.size
        corpus["classes"] = fill_classes(seed, results, arguments.size)
        pins = dict(pool.map(pin_entry, run_ids(corpus, seed), chunksize=16))
        if report_failures(pins, "pinned"):
            return 1
    corpus["pins"] = pins
    write_corpus(corpus)
    members: dict[str, int] = defaultdict(int)
    for found in results.values():
        for name in found:
            members[name] += 1
    for name, chosen in corpus["classes"].items():
        print(f"{name:22s} {len(chosen):3d} chosen from {members[name]:5d} candidates")
    print(f"{len(pins)} entries pinned")
    return 0


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    draw_parser = commands.add_parser("draw", help="draw the baseline and pools into corpus/full.json")
    draw_parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    draw_parser.add_argument("--size", type=int, default=BASELINE_SIZE)
    draw_parser.add_argument("--refresh", action="store_true", help="download new holdings and pool snapshots")
    draw_parser.set_defaults(run=draw)
    classify_parser = commands.add_parser("classify", help="classify the fetched candidates and fill each class")
    classify_parser.add_argument("--size", type=int, default=CLASS_SIZE)
    classify_parser.add_argument("--jobs", type=int, default=None)
    classify_parser.set_defaults(run=classify_corpus)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    return arguments.run(arguments)


if __name__ == "__main__":
    sys.exit(main())

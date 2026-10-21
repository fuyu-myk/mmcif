import json
import re
import sys
import tomllib
from collections import defaultdict
from pathlib import Path
from typing import Callable

from gemmi import cif

REPO_ROOT = Path(__file__).resolve().parent.parent
CANDIDATES_PATH = REPO_ROOT / "corpus" / "candidates.toml"
DATA_DIR = REPO_ROOT / "corpus" / "data"
MANIFEST_PATH = REPO_ROOT / "corpus" / "manifest.json"

OPERATOR_GROUP = re.compile(r"\(([^)]*)\)")
NUMERIC_RANGE = re.compile(r"^(\d+)-(\d+)$")


def column(block: cif.Block, tag: str) -> list[str]:
    return [cif.as_string(value) for value in block.find_values(tag)]


def has_any_value(block: cif.Block, tag: str) -> bool:
    return any(column(block, tag))


def range_size(part: str) -> int:
    match = NUMERIC_RANGE.match(part.strip())
    if match is None:
        return 1
    return int(match.group(2)) - int(match.group(1)) + 1


def operator_count(expression: str) -> int:
    groups = OPERATOR_GROUP.findall(expression) or [expression]
    total = 1
    for group in groups:
        total *= sum(range_size(part) for part in group.split(",") if part.strip())
    return total


def has_insertion_codes(block: cif.Block) -> bool:
    return has_any_value(block, "_atom_site.pdbx_PDB_ins_code")


def has_altlocs(block: cif.Block) -> bool:
    return has_any_value(block, "_atom_site.label_alt_id")


def declares_sequence_heterogeneity(block: cif.Block) -> bool:
    return any(value.lower().startswith("y") for value in column(block, "_entity_poly_seq.hetero"))


def has_mixed_residue_types(block: cif.Block) -> bool:
    asym = column(block, "_atom_site.label_asym_id")
    seq = column(block, "_atom_site.label_seq_id")
    comp = column(block, "_atom_site.label_comp_id")
    names_by_position: dict[tuple[str, str], set[str]] = defaultdict(set)
    for asym_id, seq_id, comp_id in zip(asym, seq, comp):
        if seq_id:
            names_by_position[(asym_id, seq_id)].add(comp_id)
    return any(len(names) > 1 for names in names_by_position.values())


def has_microheterogeneity(block: cif.Block) -> bool:
    return declares_sequence_heterogeneity(block) or has_mixed_residue_types(block)


def has_numbering_divergence(block: cif.Block) -> bool:
    label_seq = column(block, "_atom_site.label_seq_id")
    auth_seq = column(block, "_atom_site.auth_seq_id")
    return any(label and label != auth for label, auth in zip(label_seq, auth_seq))


def has_multi_operator_assembly(block: cif.Block) -> bool:
    expressions = column(block, "_pdbx_struct_assembly_gen.oper_expression")
    return any(operator_count(expression) > 1 for expression in expressions)


def has_branched_entities(block: cif.Block) -> bool:
    return has_any_value(block, "_pdbx_entity_branch.entity_id")


def has_multichar_chain_ids(block: cif.Block) -> bool:
    identifiers = column(block, "_atom_site.label_asym_id") + column(block, "_atom_site.auth_asym_id")
    return any(len(identifier) > 1 for identifier in set(identifiers))


def has_multiple_models(block: cif.Block) -> bool:
    return len(set(column(block, "_atom_site.pdbx_PDB_model_num"))) > 1


def has_atoms(block: cif.Block) -> bool:
    return has_any_value(block, "_atom_site.id")


PATHOLOGY_PREDICATES: dict[str, Callable[[cif.Block], bool]] = {
    "insertion_codes": has_insertion_codes,
    "altlocs": has_altlocs,
    "microheterogeneity": has_microheterogeneity,
    "numbering_divergence": has_numbering_divergence,
    "assembly_operators": has_multi_operator_assembly,
    "branched": has_branched_entities,
    "multichar_chain_ids": has_multichar_chain_ids,
    "multiple_models": has_multiple_models,
}
LISTED_CLASSES = ("baseline",)


def load_candidates(path: Path) -> dict[str, list[str]]:
    with open(path, "rb") as handle:
        raw = tomllib.load(handle)
    return {name: [identifier.lower() for identifier in entry["ids"]] for name, entry in raw.items()}


def classify_block(block: cif.Block, listed_as: set[str]) -> list[str]:
    classes = [name for name, predicate in PATHOLOGY_PREDICATES.items() if predicate(block)]
    if "baseline" in listed_as and has_atoms(block):
        classes.append("baseline")
    return sorted(classes)


def classify_entry(identifier: str, listed_as: set[str]) -> list[str] | None:
    path = DATA_DIR / f"{identifier}.cif"
    if not path.exists():
        return None
    return classify_block(cif.read(str(path)).sole_block(), listed_as)


def listed_classes_by_id(candidates: dict[str, list[str]]) -> dict[str, set[str]]:
    listed: dict[str, set[str]] = defaultdict(set)
    for class_name, ids in candidates.items():
        for identifier in ids:
            listed[identifier].add(class_name)
    return listed


def classify_all(candidates: dict[str, list[str]]) -> dict[str, list[str] | None]:
    listed = listed_classes_by_id(candidates)
    return {identifier: classify_entry(identifier, listed[identifier]) for identifier in sorted(listed)}


def find_failures(candidates: dict[str, list[str]], results: dict[str, list[str] | None]) -> list[tuple[str, str, str]]:
    failures = []
    for class_name, ids in candidates.items():
        for identifier in ids:
            classes = results[identifier]
            if classes is None:
                failures.append((class_name, identifier, "not downloaded"))
            elif class_name not in classes:
                failures.append((class_name, identifier, "predicate false"))
    return failures


def build_manifest(results: dict[str, list[str] | None]) -> dict[str, list[str]]:
    return {identifier: classes for identifier, classes in sorted(results.items()) if classes}


def summarize(manifest: dict[str, list[str]], candidates: dict[str, list[str]]) -> dict[str, list[str]]:
    members: dict[str, list[str]] = {name: [] for name in candidates}
    for identifier, classes in manifest.items():
        for class_name in classes:
            members.setdefault(class_name, []).append(identifier)
    return members


def main() -> int:
    candidates = load_candidates(CANDIDATES_PATH)
    results = classify_all(candidates)
    manifest = build_manifest(results)
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    for class_name, ids in sorted(summarize(manifest, candidates).items()):
        print(f"{class_name:22s} {len(ids):3d}  {' '.join(ids)}")
    failures = find_failures(candidates, results)
    for class_name, identifier, reason in failures:
        print(f"FAILED {identifier} listed under {class_name}: {reason}")
    print(f"{len(manifest)} entries in manifest, {len(failures)} candidate failures")
    return 0


if __name__ == "__main__":
    sys.exit(main())

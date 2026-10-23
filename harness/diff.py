import argparse
import json
import sys
from collections import defaultdict
from typing import Any, Callable

DEFAULT_TOLERANCE = 1e-3
CENTROID_TOLERANCE = 1e-2
EXAMPLES_PER_CATEGORY = 10
FLOAT_SLACK = 1e-9

AUTH, LABEL, SEQ, ICODE, COMP, NAME, ALTLOC, ELEMENT, X, Y, Z, OCC, BISO = range(13)
RESIDUE_LABEL_SEQ = 5

Differences = dict[str, list[str]]


def atom_key(atom: list[Any]) -> tuple:
    return tuple(atom[:ALTLOC + 1])


def describe_atom_key(key: tuple) -> str:
    auth, label, seq, icode, comp, name, altloc = key
    suffix = f".{altloc}" if altloc else ""
    return f"{auth}/{label} {comp} {seq}{icode} {name}{suffix}"


def residue_key(residue: list[Any]) -> tuple:
    return tuple(residue[:RESIDUE_LABEL_SEQ])


def describe_residue_key(key: tuple) -> str:
    auth, label, seq, icode, comp = key
    return f"{auth}/{label} {comp} {seq}{icode}"


def group_by(records: list[list[Any]], key_of: Callable[[list[Any]], tuple]) -> dict[tuple, list[list[Any]]]:
    groups: dict[tuple, list[list[Any]]] = defaultdict(list)
    for record in records:
        groups[key_of(record)].append(record)
    return groups


def record(differences: Differences, category: str, message: str) -> None:
    differences[category].append(message)


def compare_atom_pair(differences: Differences, prefix: str, key: tuple, ours: list[Any], theirs: list[Any], tolerance: float) -> None:
    label = f"{prefix} {describe_atom_key(key)}"
    if ours[ELEMENT] != theirs[ELEMENT]:
        record(differences, "element", f"{label}: {ours[ELEMENT]} vs {theirs[ELEMENT]}")
    deviation = max(abs(ours[axis] - theirs[axis]) for axis in (X, Y, Z))
    if deviation > tolerance + FLOAT_SLACK:
        record(differences, "coordinate", f"{label}: {ours[X:Z + 1]} vs {theirs[X:Z + 1]}")
    for column, category in ((OCC, "occupancy"), (BISO, "b_iso")):
        mine, other = ours[column], theirs[column]
        if mine is None or other is None:
            differs = mine is not other
        else:
            differs = abs(mine - other) > tolerance + FLOAT_SLACK
        if differs:
            record(differences, category, f"{label}: {mine} vs {other}")


def compare_atoms(differences: Differences, number: int, ours: list[list[Any]], theirs: list[list[Any]], tolerance: float) -> None:
    prefix = f"model {number}:"
    if len(ours) != len(theirs):
        record(differences, "atom_count", f"{prefix} {len(ours)} vs {len(theirs)}")
    ours_by_key = group_by(ours, atom_key)
    theirs_by_key = group_by(theirs, atom_key)
    for key in sorted(set(ours_by_key) | set(theirs_by_key), key=str):
        mine = ours_by_key.get(key, [])
        other = theirs_by_key.get(key, [])
        for left, right in zip(mine, other):
            compare_atom_pair(differences, prefix, key, left, right, tolerance)
        for _ in mine[len(other):]:
            record(differences, "atoms_only_ours", f"{prefix} {describe_atom_key(key)}")
        for _ in other[len(mine):]:
            record(differences, "atoms_only_gemmi", f"{prefix} {describe_atom_key(key)}")


def compare_residues(differences: Differences, number: int, ours: list[list[Any]], theirs: list[list[Any]]) -> None:
    prefix = f"model {number}:"
    ours_by_key = group_by(ours, residue_key)
    theirs_by_key = group_by(theirs, residue_key)
    for key in sorted(set(ours_by_key) | set(theirs_by_key), key=str):
        mine = ours_by_key.get(key, [])
        other = theirs_by_key.get(key, [])
        for left, right in zip(mine, other):
            if left[RESIDUE_LABEL_SEQ] != right[RESIDUE_LABEL_SEQ]:
                record(differences, "residue_label_seq", f"{prefix} {describe_residue_key(key)}: {left[RESIDUE_LABEL_SEQ]} vs {right[RESIDUE_LABEL_SEQ]}")
        for _ in mine[len(other):]:
            record(differences, "residues_only_ours", f"{prefix} {describe_residue_key(key)}")
        for _ in other[len(mine):]:
            record(differences, "residues_only_gemmi", f"{prefix} {describe_residue_key(key)}")


def compare_models(differences: Differences, ours: dict[str, Any], theirs: dict[str, Any], tolerance: float) -> None:
    ours_by_number = {model["number"]: model for model in ours["models"]}
    theirs_by_number = {model["number"]: model for model in theirs["models"]}
    for number in sorted(set(ours_by_number) - set(theirs_by_number)):
        record(differences, "models", f"model {number} only in ours")
    for number in sorted(set(theirs_by_number) - set(ours_by_number)):
        record(differences, "models", f"model {number} only in gemmi")
    for number in sorted(set(ours_by_number) & set(theirs_by_number)):
        mine, other = ours_by_number[number], theirs_by_number[number]
        compare_atoms(differences, number, mine["atoms"], other["atoms"], tolerance)
        compare_residues(differences, number, mine["residues"], other["residues"])


def compare_entities(differences: Differences, ours: dict[str, Any], theirs: dict[str, Any]) -> None:
    ours_by_id = {entity["id"]: entity for entity in ours["entities"]}
    theirs_by_id = {entity["id"]: entity for entity in theirs["entities"]}
    for identifier in sorted(set(ours_by_id) - set(theirs_by_id)):
        record(differences, "entity_set", f"entity {identifier} only in ours")
    for identifier in sorted(set(theirs_by_id) - set(ours_by_id)):
        record(differences, "entity_set", f"entity {identifier} only in gemmi")
    for identifier in sorted(set(ours_by_id) & set(theirs_by_id)):
        mine, other = ours_by_id[identifier], theirs_by_id[identifier]
        for field, category in (("type", "entity_type"), ("polymer_type", "entity_polymer_type")):
            if mine[field] != other[field]:
                record(differences, category, f"entity {identifier}: {mine[field]} vs {other[field]}")
        if sorted(mine["subchains"]) != sorted(other["subchains"]):
            record(differences, "entity_subchains", f"entity {identifier}: {sorted(mine['subchains'])} vs {sorted(other['subchains'])}")


def copies_match(ours: dict[str, Any], theirs: dict[str, Any]) -> bool:
    if ours["atoms"] != theirs["atoms"]:
        return False
    return all(abs(a - b) <= CENTROID_TOLERANCE + FLOAT_SLACK for a, b in zip(ours["centroid"], theirs["centroid"]))


def unmatched_copies(ours: list[dict[str, Any]], theirs: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    remaining = list(theirs)
    only_ours = []
    for copy in ours:
        index = next((i for i, candidate in enumerate(remaining) if copies_match(copy, candidate)), None)
        if index is None:
            only_ours.append(copy)
        else:
            remaining.pop(index)
    return only_ours, remaining


def describe_copy(copy: dict[str, Any]) -> str:
    return f"{copy['atoms']} atoms at {copy['centroid']}"


def compare_assembly_copies(differences: Differences, identifier: str, ours: list[dict[str, Any]], theirs: list[dict[str, Any]]) -> None:
    only_ours, only_theirs = unmatched_copies(ours, theirs)
    for copy in only_ours:
        record(differences, "assembly_copies", f"assembly {identifier}: only in ours: {describe_copy(copy)}")
    for copy in only_theirs:
        record(differences, "assembly_copies", f"assembly {identifier}: only in gemmi: {describe_copy(copy)}")


def compare_assemblies(differences: Differences, ours: dict[str, Any], theirs: dict[str, Any]) -> None:
    ours_by_id = {assembly["id"]: assembly for assembly in ours["assemblies"]}
    theirs_by_id = {assembly["id"]: assembly for assembly in theirs["assemblies"]}
    for identifier in sorted(set(ours_by_id) - set(theirs_by_id)):
        record(differences, "assembly_set", f"assembly {identifier} only in ours")
    for identifier in sorted(set(theirs_by_id) - set(ours_by_id)):
        record(differences, "assembly_set", f"assembly {identifier} only in gemmi")
    for identifier in sorted(set(ours_by_id) & set(theirs_by_id)):
        mine, other = ours_by_id[identifier], theirs_by_id[identifier]
        if mine["atom_count"] != other["atom_count"]:
            record(differences, "assembly_atom_count", f"assembly {identifier}: {mine['atom_count']} vs {other['atom_count']}")
        compare_assembly_copies(differences, identifier, mine["copies"], other["copies"])


def compare(ours: dict[str, Any], theirs: dict[str, Any], tolerance: float = DEFAULT_TOLERANCE) -> Differences:
    differences: Differences = defaultdict(list)
    if ours["entry"] != theirs["entry"]:
        record(differences, "entry", f"{ours['entry']} vs {theirs['entry']}")
    compare_models(differences, ours, theirs, tolerance)
    compare_entities(differences, ours, theirs)
    compare_assemblies(differences, ours, theirs)
    return dict(differences)


def format_report(differences: Differences, limit: int = EXAMPLES_PER_CATEGORY) -> str:
    lines = []
    for category in sorted(differences):
        messages = differences[category]
        lines.append(f"{category}: {len(messages)}")
        lines.extend(f"  {message}" for message in messages[:limit])
        if len(messages) > limit:
            lines.append(f"  ... {len(messages) - limit} more")
    return "\n".join(lines)


def load_dump(path: str) -> dict[str, Any]:
    with open(path) as handle:
        return json.load(handle)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ours", required=True)
    parser.add_argument("--gemmi", required=True)
    parser.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    differences = compare(load_dump(arguments.ours), load_dump(arguments.gemmi), arguments.tolerance)
    if not differences:
        print("identical")
        return 0
    print(format_report(differences))
    return 1


if __name__ == "__main__":
    sys.exit(main())

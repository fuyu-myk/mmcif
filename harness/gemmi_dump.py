import json
import sys
from typing import Any

import gemmi

COORDINATE_DECIMALS = 3

POLYMER_TYPE_NAMES = {
    gemmi.PolymerType.PeptideL: "polypeptide(L)",
    gemmi.PolymerType.PeptideD: "polypeptide(D)",
    gemmi.PolymerType.Dna: "polydeoxyribonucleotide",
    gemmi.PolymerType.Rna: "polyribonucleotide",
    gemmi.PolymerType.DnaRnaHybrid: "polydeoxyribonucleotide/polyribonucleotide hybrid",
    gemmi.PolymerType.SaccharideD: "polysaccharide(D)",
    gemmi.PolymerType.SaccharideL: "polysaccharide(L)",
    gemmi.PolymerType.Pna: "peptide nucleic acid",
    gemmi.PolymerType.CyclicPseudoPeptide: "cyclic-pseudo-peptide",
    gemmi.PolymerType.Other: "other",
}

ENTITY_TYPE_NAMES = {
    gemmi.EntityType.Polymer: "polymer",
    gemmi.EntityType.NonPolymer: "non-polymer",
    gemmi.EntityType.Branched: "branched",
    gemmi.EntityType.Water: "water",
    gemmi.EntityType.Unknown: "unknown",
}


def rounded(value: float) -> float:
    return round(value, COORDINATE_DECIMALS) + 0.0


def insertion_code(residue: gemmi.Residue) -> str:
    code = residue.seqid.icode
    return "" if code == " " else code


def alternate_location(atom: gemmi.Atom) -> str:
    return "" if atom.altloc == "\0" else atom.altloc


def atom_record(chain: gemmi.Chain, residue: gemmi.Residue, atom: gemmi.Atom) -> list[Any]:
    position = atom.pos
    return [
        chain.name,
        residue.subchain,
        residue.seqid.num,
        insertion_code(residue),
        residue.name,
        atom.name,
        alternate_location(atom),
        atom.element.name.upper(),
        rounded(position.x),
        rounded(position.y),
        rounded(position.z),
        rounded(atom.occ),
        rounded(atom.b_iso),
    ]


def residue_record(chain: gemmi.Chain, residue: gemmi.Residue) -> list[Any]:
    return [
        chain.name,
        residue.subchain,
        residue.seqid.num,
        insertion_code(residue),
        residue.name,
        residue.label_seq,
    ]


def dump_model(model: gemmi.Model) -> dict[str, Any]:
    atoms = []
    residues = []
    for chain in model:
        for residue in chain:
            residues.append(residue_record(chain, residue))
            atoms.extend(atom_record(chain, residue, atom) for atom in residue)
    return {"number": model.num, "atoms": atoms, "residues": residues}


def dump_entity(entity: gemmi.Entity) -> dict[str, Any]:
    is_polymer = entity.entity_type == gemmi.EntityType.Polymer
    polymer_type = POLYMER_TYPE_NAMES.get(entity.polymer_type) if is_polymer else None
    return {
        "id": entity.name,
        "type": ENTITY_TYPE_NAMES[entity.entity_type],
        "polymer_type": polymer_type,
        "subchains": sorted(entity.subchains),
    }


def subchain_copy(residues: list[gemmi.Residue]) -> dict[str, Any]:
    positions = [atom.pos for residue in residues for atom in residue]
    count = len(positions)
    centroid = [sum(axis) / count for axis in zip(*((p.x, p.y, p.z) for p in positions))]
    return {"atoms": count, "centroid": [rounded(value) for value in centroid]}


def chain_copies(chain: gemmi.Chain) -> list[dict[str, Any]]:
    residues_by_subchain: dict[str, list[gemmi.Residue]] = {}
    for residue in chain:
        residues_by_subchain.setdefault(residue.subchain, []).append(residue)
    return [subchain_copy(residues) for residues in residues_by_subchain.values() if any(len(r) for r in residues)]


def dump_assembly(assembly: gemmi.Assembly, model: gemmi.Model) -> dict[str, Any]:
    expanded = gemmi.make_assembly(assembly, model, gemmi.HowToNameCopiedChain.Dup)
    copies = [copy for chain in expanded for copy in chain_copies(chain)]
    return {"id": assembly.name, "atom_count": sum(copy["atoms"] for copy in copies), "copies": copies}


def entry_identifier(structure: gemmi.Structure) -> str | None:
    return structure.info["_entry.id"] if "_entry.id" in structure.info else None


def dump_structure(path: str) -> dict[str, Any]:
    structure = gemmi.read_structure(path)
    first_model = structure[0] if len(structure) else None
    assemblies = [dump_assembly(a, first_model) for a in structure.assemblies] if first_model is not None else []
    return {
        "entry": entry_identifier(structure),
        "models": [dump_model(model) for model in structure],
        "entities": [dump_entity(entity) for entity in structure.entities],
        "assemblies": assemblies,
    }


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: gemmi_dump.py <file.cif>", file=sys.stderr)
        return 2
    json.dump(dump_structure(sys.argv[1]), sys.stdout, separators=(",", ":"))
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

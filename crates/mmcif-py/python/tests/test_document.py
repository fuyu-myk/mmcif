import pytest

import mmcif_io as m

TEXT = """data_demo
_entry.id DEMO
loop_
_atom_site.id
_atom_site.label_alt_id
_atom_site.pdbx_PDB_ins_code
1 . ?
2 A 'B'
"""


def test_table_layer():
    document = m.Document.parse(TEXT)
    assert len(document) == 1
    block = document.blocks[0]
    assert block.name == "demo" and document.block("DEMO").name == "demo"
    assert block.category_names == ["entry", "atom_site"]
    atoms = block["atom_site"]
    assert atoms.items == ["id", "label_alt_id", "pdbx_PDB_ins_code"]
    assert len(atoms) == atoms.row_count == 2
    assert atoms["label_alt_id"] == [None, "A"]
    assert atoms.raw_column("pdbx_PDB_ins_code") == ["?", "B"]
    assert "id" in atoms and "missing" not in atoms
    assert atoms.column("missing") is None
    with pytest.raises(KeyError):
        block["missing"]
    assert block.category("entry")["id"] == ["DEMO"]


def test_round_trip_and_handles_outlive_document():
    document = m.Document.parse(TEXT)
    category = document.blocks[0].category("atom_site")
    again = m.Document.parse(document.to_cif())
    assert again.blocks[0]["atom_site"]["label_alt_id"] == [None, "A"]
    del document
    assert category["id"] == ["1", "2"]


def test_document_errors(tmp_path):
    with pytest.raises(m.ReadError, match="line 2"):
        m.Document.parse("data_x\n_a.b\n")
    with pytest.raises(m.ReadError):
        m.Document.open(tmp_path / "missing.cif")
    with pytest.raises(m.ReadError):
        m.Document.from_bcif(b"\x00\x01")


def test_structure_from_document(pathologies_path):
    document = m.Document.open(pathologies_path)
    structure = m.Structure.from_document(document)
    assert structure.models[0].atom_count == 16


def test_binary_document(crambin_bcif_path):
    document = m.Document.open_bcif(crambin_bcif_path)
    assert document.blocks[0]["atom_site"].row_count == 327

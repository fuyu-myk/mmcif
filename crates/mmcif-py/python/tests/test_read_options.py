import pytest

import mmcif_io as m

BROKEN_ENTITIES = ("1 2 SER y\n", "1 x SER y\n")
BROKEN_ASSEMBLIES = ("1 1 AA,B,C,D\n", "1 ? AA,B,C,D\n")
BROKEN_REMARK_350 = (
    "REMARK 350 BIOMOLECULE: 1\n"
    "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A\n"
    "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000\n"
    "REMARK 350   BIOMT2   1  0.000000  x.000000  0.000000        0.00000\n"
    "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000\n"
)


def broken_cif(pathologies_path, tmp_path, *breaks):
    text = pathologies_path.read_text()
    for before, after in breaks:
        assert before in text
        text = text.replace(before, after)
    path = tmp_path / "broken.cif"
    path.write_text(text)
    return path


def test_first_model_only_reads_one_model(pathologies_path, legacy_pdb_path):
    assert len(m.read(pathologies_path).models) > 1
    for structure in [
        m.read(pathologies_path, first_model_only=True),
        m.read_mmcif(pathologies_path, first_model_only=True),
    ]:
        assert [model.number for model in structure.models] == [1]
        assert structure.models[0].atom_count == 16
    assert len(m.read_pdb(legacy_pdb_path).models) == 2
    assert len(m.read_pdb(legacy_pdb_path, first_model_only=True).models) == 1
    assert len(m.read(legacy_pdb_path, first_model_only=True).models) == 1


def test_a_lenient_read_skips_malformed_optional_categories(pathologies_path, tmp_path):
    path = broken_cif(pathologies_path, tmp_path, BROKEN_ENTITIES, BROKEN_ASSEMBLIES)
    with pytest.raises(m.ReadError):
        m.read(path)
    with pytest.raises(m.ReadError):
        m.read_mmcif(path, first_model_only=True)
    for structure in [m.read(path, lenient=True), m.read_mmcif(path, lenient=True)]:
        assert structure.entities == [] and structure.assemblies == []
        assert len(structure.skipped) == 2
        assert structure.skipped[0].startswith("_entity_poly_seq.")
        assert structure.skipped[1].startswith("_pdbx_struct_assembly_gen.")
        assert structure.models[0].atom_count == 16


def test_a_lenient_pdb_read_skips_malformed_remark_350(legacy_pdb_path, tmp_path):
    header, rest = legacy_pdb_path.read_text().split("\n", 1)
    path = tmp_path / "broken.pdb"
    path.write_text(f"{header}\n{BROKEN_REMARK_350}{rest}")
    with pytest.raises(m.ReadError, match="PDB line"):
        m.read_pdb(path)
    for structure in [m.read_pdb(path, lenient=True), m.read(path, lenient=True)]:
        assert structure.assemblies == []
        assert len(structure.skipped) == 1
        assert structure.skipped[0].startswith("PDB line 5:")
        assert len(structure.models) == 2


def test_nothing_is_skipped_without_a_malformed_category(pathologies_path, crambin_bcif_path):
    assert m.read(pathologies_path).skipped == []
    assert m.read(pathologies_path, lenient=True).skipped == []
    crambin = m.read_bcif(crambin_bcif_path, first_model_only=True, lenient=True)
    assert crambin.skipped == [] and crambin.models[0].atom_count == 327
    document = m.Document.open(pathologies_path)
    assert m.Structure.from_document(document).skipped == []


def test_the_options_are_keywords(pathologies_path):
    with pytest.raises(TypeError):
        m.read(pathologies_path, True)

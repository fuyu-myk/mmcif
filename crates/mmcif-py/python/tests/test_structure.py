import gzip
import threading

import numpy as np
import pytest

import mmcif_io as m


@pytest.fixture
def structure(pathologies_path):
    return m.read_mmcif(pathologies_path)


def test_hierarchy(structure):
    assert structure.entry_id == "TEST"
    assert [model.number for model in structure.models] == [1, 2]
    model = structure.models[0]
    assert model.atom_count == 16 and len(model) == 4
    assert [chain.label_asym for chain in model.chains] == ["AA", "B", "C", "D"]
    chain = model.chain("AA")
    assert chain.auth_asym == "H1"
    assert chain.entity_id == "1"
    assert [residue.comp for residue in chain.residues] == ["SER", "THR", "ALA", "LYS", "VAL"]
    assert model.chain("missing") is None


def test_insertion_codes_and_microheterogeneity(structure):
    chain = structure.models[0].chain("AA")
    assert [r.comp for r in chain.residues_at_auth(100)] == ["SER", "THR"]
    inserted = chain.residues_at_auth(100, "A")
    assert [r.comp for r in inserted] == ["ALA"]
    assert inserted[0].ins_code == "A"
    assert inserted[0].key == ("AA", 3, 100, "A", "ALA", ())
    assert inserted[0].auth_key == ("H1", 100, "A")
    assert [r.comp for r in chain.residues_at_label(2)] == ["SER", "THR"]


def test_atoms_and_altlocs(structure):
    model = structure.models[0]
    alanine = model.chain("AA").residues_at_auth(100, "A")[0]
    assert alanine.altlocs == ["A", "B"]
    zinc = model.chain("C").residues[0].atom("ZN")
    assert (zinc.element, zinc.charge, zinc.serial) == ("Zn", 2, 14)
    assert zinc.position == (0.0, 0.0, 5.0)
    assert zinc.residue.chain.label_asym == "C"
    single = model.single_conformer("highest_occupancy")
    assert single.atom_count == 13
    assert [r.comp for r in model.single_conformer("B").chain("AA").residues][0] == "THR"
    with pytest.raises(ValueError):
        model.single_conformer("nonsense")


def test_struct_of_arrays_columns(structure):
    model = structure.models[0]
    coordinates = model.coordinates
    assert coordinates.shape == (16, 3) and coordinates.dtype == np.float64
    np.testing.assert_allclose(coordinates[0], [1.0, 2.0, 3.0])
    assert model.elements[13] == "Zn"
    assert model.atom_names[:2] == ["N", "CA"]
    assert model.altlocs[:4] == ["A", "A", "B", "B"]
    assert model.occupancies.shape == (16,)
    residues = model.residue_indices
    for residue in model.residues():
        start, end = residue.atom_range
        assert (residues[start:end] == residue.index).all()
    chain = model.chain("D")
    start, end = chain.atom_range
    np.testing.assert_array_equal(chain.coordinates, coordinates[start:end])
    assert (model.chain_indices[start:end] == chain.index).all()


def test_numbering(structure):
    assert structure.to_auth("AA", 1) == ("H1", 99, None)
    assert structure.to_auth("AA", 3) == ("H1", 100, "A")
    assert structure.to_label("H1", 100, "A") == ("AA", 3)
    with pytest.raises(m.NumberingError, match="not a polymer residue"):
        structure.to_label("H1", 201)
    with pytest.raises(m.NumberingError):
        structure.to_auth("AA", 0)
    assert issubclass(m.NumberingError, ValueError)


def test_entities_and_assemblies(structure):
    kinds = [entity.kind for entity in structure.entities]
    assert kinds == ["polymer", "branched", "non-polymer", "water"]
    polymer = structure.entity("1")
    assert polymer.polymer_type == "polypeptide(L)"
    assert polymer.canonical_sequence == "GSAKVT"
    assert polymer.sequence[1] == ["SER", "THR"]
    assert [a.operator_count for a in structure.assemblies] == [1, 2, 1]
    dimer = structure.assemble("2")
    assert [c.operators for c in dimer.chains] == [["1"], ["1"], ["2"], ["2"]]
    product = structure.assemble("3")
    np.testing.assert_allclose(product.coordinates[0], [-11.0, -2.0, 3.0])
    with pytest.raises(m.AssemblyError):
        structure.assemble("9")
    with pytest.raises(IndexError):
        structure.assemble("1", model_index=5)


def test_assembly_oligomeric_state_and_deposited_check(structure, crambin_bcif_path):
    crambin = m.read(crambin_bcif_path).assemblies[0]
    assert (crambin.oligomeric_details, crambin.oligomeric_count) == ("monomeric", 1)
    assert [(a.oligomeric_details, a.oligomeric_count) for a in structure.assemblies] == [(None, None)] * 3
    assert [structure.assembly_is_deposited(i) for i in ("1", "2", "3")] == [True, False, True]
    assert structure.assembly_is_deposited("1", tolerance=0.0)
    with pytest.raises(m.AssemblyError):
        structure.assembly_is_deposited("9")


def test_handles_outlive_their_structure(pathologies_path):
    structure = m.read(pathologies_path)
    atom = structure.models[0].chain("C").residues[0].atoms[0]
    del structure
    assert atom.name == "ZN"
    assert atom.residue.chain.model.atom_count == 16


def test_handle_identity(structure):
    model = structure.models[0]
    first = model.chain("AA").residues[0]
    again = model.residues()[0]
    assert first == again and hash(first) == hash(again)
    assert first != model.residues()[1]
    assert len({first, again}) == 1
    assert structure.models[0] == structure.models[0]
    assert hash(model.chain("AA").model) == hash(model)
    assert structure.models[0] != structure.models[1]


def test_missing_files_raise_os_errors(tmp_path):
    missing = tmp_path / "missing.cif"
    with pytest.raises(FileNotFoundError) as raised:
        m.read(missing)
    assert raised.value.filename == str(missing)
    with pytest.raises(FileNotFoundError):
        m.Document.open(missing)
    with pytest.raises(FileNotFoundError):
        m.CcdIndex.open(missing)
    directory = tmp_path / "directory.cif"
    directory.mkdir()
    with pytest.raises(OSError):
        m.read_mmcif(directory)


def test_read_errors(tmp_path):
    bad = tmp_path / "bad.cif"
    bad.write_text("data_x\nloop_\n_atom_site.id\n")
    with pytest.raises(m.ReadError, match="atom_site"):
        m.read_mmcif(bad)
    with pytest.raises(m.ReadError):
        m.parse_mmcif("data_x\n_a.b 'unterminated\n")
    with pytest.raises(m.ReadError):
        m.read(tmp_path / "structure.xyz")


def test_reads_release_the_gil(pathologies_path):
    results = []
    threads = [threading.Thread(target=lambda: results.append(m.read(pathologies_path).models[0].atom_count)) for _ in range(8)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert results == [16] * 8


def test_legacy_pdb(legacy_pdb_path):
    structure = m.read_pdb(legacy_pdb_path)
    assert structure.entry_id == "1ABC"
    model = structure.models[0]
    assert [c.label_asym for c in model.chains] == ["A", "B", "C", "D", "E", "F"]
    with pytest.raises(m.NumberingError, match="no label numbering"):
        structure.to_label("A", 1)


def test_binary_cif_matches_text(crambin_bcif_path):
    binary = m.read_bcif(crambin_bcif_path)
    text = m.parse_bcif(crambin_bcif_path.read_bytes())
    assert binary.models[0].atom_count == text.models[0].atom_count == 327
    np.testing.assert_array_equal(binary.models[0].coordinates, text.models[0].coordinates)


def test_reads_gzip_compressed_files(tmp_path, crambin_bcif_path, legacy_pdb_path):
    for source, name in ((crambin_bcif_path, "1crn.bcif.gz"), (legacy_pdb_path, "legacy.pdb1.gz")):
        target = tmp_path / name
        target.write_bytes(gzip.compress(source.read_bytes()))
        compressed, plain = m.read(target), m.read(source)
        assert compressed.models[0].atom_count == plain.models[0].atom_count
        np.testing.assert_array_equal(compressed.models[0].coordinates, plain.models[0].coordinates)


def test_a_missing_b_factor_is_none_and_nan():
    text = (
        "data_x\nloop_\n_atom_site.label_asym_id\n_atom_site.label_comp_id\n_atom_site.label_atom_id\n"
        "_atom_site.type_symbol\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n_atom_site.B_iso_or_equiv\nA GLY N N 1 0 0 0 .\nA GLY CA C 1 1.5 0 0 12.5\n"
    )
    model = m.parse_mmcif(text).models[0]
    assert [atom.b_factor for atom in model.atoms()] == [None, 12.5]
    assert np.isnan(model.b_factors[0]) and model.b_factors[1] == 12.5


def test_assembled_residue_keys_name_their_operators(structure):
    dimer = structure.assemble("2")
    keys = [residue.key for residue in dimer.residues()]
    assert len(set(keys)) == len(keys) == 12
    copy = next(r for r in dimer.residues() if r.chain.operators == ["2"])
    assert copy.key == ("AA", 2, 100, None, "SER", ("2",))
    assert repr(copy) == "<Residue AA/2/100/SER@2>"

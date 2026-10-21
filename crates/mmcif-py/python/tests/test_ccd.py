import pytest

import mmcif_io as m


def test_build_save_and_query(components_path, tmp_path):
    out = tmp_path / "ccd.idx"
    built = m.CcdIndex.build(components_path, "test-release", out)
    index = m.CcdIndex.open(out)
    assert index.release == built.release == "test-release"
    assert len(index) == 8
    assert index.ids() == sorted(index.ids())
    alanine = index["ALA"]
    assert (alanine.one_letter_code, alanine.kind) == ("A", "L-PEPTIDE LINKING")
    assert len(alanine.atoms) == 13 and len(alanine.bonds) == 12
    name, element, charge, leaving, aromatic, stereo = alanine.atoms[1]
    assert (name, element, stereo) == ("CA", "C", "S")
    assert index["MSE"].parent == "MET"
    assert "HEM" in index and "NOPE" not in index
    assert index.get("NOPE") is None
    with pytest.raises(KeyError):
        index["NOPE"]


def test_corrupt_index():
    with pytest.raises(m.ReadError):
        m.CcdIndex.from_bytes(b"not an index")

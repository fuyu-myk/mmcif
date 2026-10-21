from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[4]
MODEL_FIXTURES = REPO / "crates" / "mmcif-model" / "tests" / "fixtures"


@pytest.fixture
def pathologies_path():
    return MODEL_FIXTURES / "pathologies.cif"


@pytest.fixture
def legacy_pdb_path():
    return MODEL_FIXTURES / "legacy.pdb"


@pytest.fixture
def crambin_bcif_path():
    return REPO / "fuzz" / "seeds" / "bcif" / "1crn.bcif"


@pytest.fixture
def components_path():
    return REPO / "crates" / "mmcif-ccd" / "tests" / "fixtures" / "components-sample.cif"

import doctest
import inspect
import io
import shutil

import pytest

import mmcif_io


def public_objects():
    yield "mmcif_io", mmcif_io
    for name in sorted(dir(mmcif_io)):
        if name.startswith("_") or name == "mmcif_io":
            continue
        obj = getattr(mmcif_io, name)
        yield name, obj
        if inspect.isclass(obj) and not issubclass(obj, BaseException):
            for attr, member in vars(obj).items():
                if not attr.startswith("_"):
                    yield f"{name}.{attr}", member


def test_every_public_name_has_a_docstring():
    missing = [
        name
        for name, obj in public_objects()
        if not (getattr(obj, "__doc__", None) or "").strip()
    ]
    assert missing == []


@pytest.fixture
def example_files(
    tmp_path, monkeypatch, crambin_bcif_path, pathologies_path, legacy_pdb_path, components_path
):
    shutil.copy(crambin_bcif_path, tmp_path / "1crn.bcif")
    crambin = mmcif_io.Document.open_bcif(crambin_bcif_path)
    (tmp_path / "1crn.cif").write_text(crambin.to_cif())
    shutil.copy(pathologies_path, tmp_path / "example.cif")
    shutil.copy(legacy_pdb_path, tmp_path / "example.pdb")
    shutil.copy(components_path, tmp_path / "components.cif")
    monkeypatch.chdir(tmp_path)


def test_docstring_examples_pass(example_files):
    parser = doctest.DocTestParser()
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
    report = io.StringIO()
    examples = 0
    for name, obj in public_objects():
        text = getattr(obj, "__doc__", None) or ""
        test = parser.get_doctest(text, {"mmcif_io": mmcif_io}, name, None, 0)
        examples += len(test.examples)
        runner.run(test, out=report.write)
    assert runner.failures == 0, report.getvalue()
    assert examples > 0

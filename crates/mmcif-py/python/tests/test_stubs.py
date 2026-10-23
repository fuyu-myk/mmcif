import ast
import importlib.resources
import subprocess
import sys
import types
import typing
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

import mmcif_io as m

ALLOWLIST = Path(__file__).with_name("stubtest-allowlist.txt")


def test_the_package_ships_its_stubs_and_is_marked_typed():
    package = importlib.resources.files("mmcif_io")
    assert package.joinpath("__init__.pyi").is_file()
    assert package.joinpath("py.typed").is_file()


def test_the_stubs_match_the_binding():
    result = subprocess.run(
        [sys.executable, "-m", "mypy.stubtest", "mmcif_io", "--allowlist", str(ALLOWLIST)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def stub_properties():
    tree = ast.parse(importlib.resources.files("mmcif_io").joinpath("__init__.pyi").read_text())
    namespace = {"np": np, "NDArray": NDArray, **vars(m)}
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom, ast.Assign)):
            exec(compile(ast.Module([node], []), "stub", "exec"), namespace)
    properties = {}
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            properties[node.name] = {
                item.name: eval(ast.unparse(item.returns), namespace)
                for item in node.body
                if isinstance(item, ast.FunctionDef)
                and any(ast.unparse(d) == "property" for d in item.decorator_list)
            }
    return properties


def conforms(value, annotation):
    origin, args = typing.get_origin(annotation), typing.get_args(annotation)
    if annotation is None:
        return value is None
    if origin in (typing.Union, types.UnionType):
        return any(conforms(value, arg) for arg in args)
    if origin is list:
        return isinstance(value, list) and all(conforms(item, args[0]) for item in value)
    if origin is tuple:
        if len(args) == 2 and args[1] is Ellipsis:
            return isinstance(value, tuple) and all(conforms(item, args[0]) for item in value)
        return isinstance(value, tuple) and len(value) == len(args) and all(map(conforms, value, args))
    if origin is NDArray:
        return isinstance(value, np.ndarray) and value.dtype == args[0]
    if origin is np.ndarray:
        return isinstance(value, np.ndarray) and value.dtype == typing.get_args(args[1])[0]
    if annotation is int:
        return isinstance(value, int) and not isinstance(value, bool)
    return isinstance(value, annotation)


def samples(pathologies_path, components_path, tmp_path):
    structure = m.read(pathologies_path)
    model = structure.assemble("2")
    document = m.Document.open(pathologies_path)
    index = m.CcdIndex.build(components_path, "test-release", tmp_path / "components.idx")
    objects = [
        structure,
        m.read(pathologies_path, lenient=True),
        *structure.entities,
        *structure.assemblies,
        *structure.models,
        model,
        *model.chains,
        *model.residues(),
        *model.atoms(),
        m.NeighbourGrid.for_model(model),
        document,
        *document.blocks,
        *[document.blocks[0].category(name) for name in document.blocks[0].category_names],
        index,
        *[index[id] for id in index.ids()],
    ]
    return objects


def test_every_property_returns_the_type_its_stub_gives(pathologies_path, components_path, tmp_path):
    properties = stub_properties()
    checked = set()
    for sample in samples(pathologies_path, components_path, tmp_path):
        name = type(sample).__name__
        for attribute, annotation in properties[name].items():
            value = getattr(sample, attribute)
            assert conforms(value, annotation), f"{name}.{attribute} = {value!r} is not {annotation}"
            checked.add((name, attribute))
    expected = {(name, attribute) for name, attributes in properties.items() for attribute in attributes}
    assert checked == expected

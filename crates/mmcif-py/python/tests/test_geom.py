import numpy as np
import pytest

import mmcif_io as m


def brute_force_pairs(points, radius):
    distances = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=-1)
    first, second = np.nonzero(np.triu(distances <= radius, k=1))
    return np.stack([first, second], axis=1)


def test_queries_match_brute_force():
    rng = np.random.default_rng(7)
    points = rng.uniform(-20.0, 20.0, size=(400, 3))
    grid = m.NeighbourGrid(points, cell_size=3.0)
    assert len(grid) == 400 and grid.cell_size == 3.0
    query = (1.0, -2.0, 0.5)
    expected = np.nonzero(np.linalg.norm(points - query, axis=1) <= 6.0)[0]
    np.testing.assert_array_equal(grid.within(query, 6.0), expected)
    np.testing.assert_array_equal(grid.pairs_within(4.0), brute_force_pairs(points, 4.0))
    index, distance = grid.nearest(query)
    assert index == np.argmin(np.linalg.norm(points - query, axis=1))
    assert distance == pytest.approx(np.linalg.norm(points[index] - query))
    others = rng.uniform(-20.0, 20.0, size=(30, 3))
    contacts = grid.contacts(others, 5.0)
    assert contacts.shape[1] == 2
    for other, own in contacts:
        assert np.linalg.norm(others[other] - points[own]) <= 5.0


def test_grid_over_model_coordinates(pathologies_path):
    model = m.read(pathologies_path).models[0]
    grid = m.NeighbourGrid.for_model(model, cell_size=2.0)
    assert len(grid) == model.atom_count
    zinc = model.chain("C").residues[0].atoms[0]
    near = grid.within(zinc.position, 0.5)
    assert list(near) == [zinc.index]
    assert grid.nearest((100.0, 100.0, 100.0), 1.0) is None


def test_invalid_input():
    with pytest.raises(ValueError, match="shape"):
        m.NeighbourGrid(np.zeros((4, 2)))
    with pytest.raises(ValueError):
        m.NeighbourGrid(np.array([[0.0, np.nan, 0.0]]))
    grid = m.NeighbourGrid(np.zeros((1, 3)))
    with pytest.raises(ValueError):
        grid.within((0.0, 0.0, 0.0), -1.0)

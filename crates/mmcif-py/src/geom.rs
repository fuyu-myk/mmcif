use mmcif_geom::{GeomError, NeighbourGrid, Radius};
use numpy::{PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::hierarchy::PyModel;

fn geom_error(error: GeomError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn checked_radius(radius: f64) -> PyResult<Radius> {
    Radius::new(radius).map_err(geom_error)
}

fn is_finite(point: &[f64; 3]) -> bool {
    point.iter().all(|value| value.is_finite())
}

fn query_point(point: [f64; 3]) -> PyResult<[f64; 3]> {
    if is_finite(&point) {
        Ok(point)
    } else {
        Err(PyValueError::new_err(
            "query point has a non-finite coordinate",
        ))
    }
}

fn rows(array: &PyReadonlyArray2<'_, f64>) -> PyResult<Vec<[f64; 3]>> {
    let view = array.as_array();
    if view.ncols() != 3 {
        return Err(PyValueError::new_err(format!(
            "expected an (n, 3) array, got shape {:?}",
            view.shape()
        )));
    }
    Ok(view
        .rows()
        .into_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect())
}

fn pair_array<'py>(py: Python<'py>, pairs: Vec<(u32, u32)>) -> PyResult<Bound<'py, PyArray2<u32>>> {
    let count = pairs.len();
    let flat: Vec<u32> = pairs.into_iter().flat_map(|(a, b)| [a, b]).collect();
    PyArray1::from_vec(py, flat).reshape([count, 2])
}

/// A spatial index over points for distance queries.
///
/// ``NeighbourGrid(points, cell_size=4.0)`` builds a grid over ``points``, an
/// ``(n, 3)`` float64 NumPy array, with cubic cells of edge length ``cell_size``.
/// The grid keeps its own copy of the points, and queries return indices into
/// them. ``len(grid)`` is the number of points.
///
/// Raises
/// ------
/// ValueError
///     If ``points`` is not an ``(n, 3)`` array, a point has a NaN or infinite
///     coordinate, or ``cell_size`` is not positive and finite.
///
/// Examples
/// --------
/// >>> import numpy as np
/// >>> grid = mmcif_io.NeighbourGrid(np.array([[0.0, 0, 0], [1, 0, 0], [5, 0, 0]]))
/// >>> grid
/// <NeighbourGrid with 3 points, cell size 4>
/// >>> grid.within([0, 0, 0], 1.0).tolist()
/// [0, 1]
#[pyclass(name = "NeighbourGrid", module = "mmcif_io", frozen)]
pub struct PyNeighbourGrid {
    inner: NeighbourGrid,
}

#[pymethods]
impl PyNeighbourGrid {
    #[new]
    #[pyo3(signature = (points, cell_size = 4.0))]
    fn new(py: Python<'_>, points: PyReadonlyArray2<'_, f64>, cell_size: f64) -> PyResult<Self> {
        let points = rows(&points)?;
        let inner = py
            .detach(|| NeighbourGrid::new(&points, cell_size))
            .map_err(geom_error)?;
        Ok(PyNeighbourGrid { inner })
    }

    /// Build a grid over the atom positions of ``model``, with cubic cells of edge
    /// length ``cell_size``.
    ///
    /// Queries return atom indices, as in ``Model.atoms()``.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If ``cell_size`` is not positive and finite.
    ///
    /// Examples
    /// --------
    /// >>> model = mmcif_io.read("1crn.bcif").models[0]
    /// >>> mmcif_io.NeighbourGrid.for_model(model)
    /// <NeighbourGrid with 327 points, cell size 4>
    #[staticmethod]
    #[pyo3(signature = (model, cell_size = 4.0))]
    fn for_model(py: Python<'_>, model: &PyModel, cell_size: f64) -> PyResult<Self> {
        let model = std::sync::Arc::clone(&model.inner);
        let inner = py
            .detach(move || NeighbourGrid::new(model.positions(), cell_size))
            .map_err(geom_error)?;
        Ok(PyNeighbourGrid { inner })
    }

    /// The edge length of a grid cell.
    #[getter]
    fn cell_size(&self) -> f64 {
        self.inner.cell_size()
    }

    /// Return the indices of the points within ``radius`` of ``point``, in ascending
    /// order, as a uint32 NumPy array.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If ``radius`` is negative or NaN, or ``point`` has a NaN or infinite
    ///     coordinate.
    ///
    /// Examples
    /// --------
    /// >>> model = mmcif_io.read("1crn.bcif").models[0]
    /// >>> grid = mmcif_io.NeighbourGrid.for_model(model)
    /// >>> atoms = model.atoms()
    /// >>> [atoms[index].name for index in grid.within(atoms[1].position, 1.6)]
    /// ['N', 'CA', 'C', 'CB']
    fn within<'py>(
        &self,
        py: Python<'py>,
        point: [f64; 3],
        radius: f64,
    ) -> PyResult<Bound<'py, PyArray1<u32>>> {
        let radius = checked_radius(radius)?;
        let point = query_point(point)?;
        let found = py.detach(|| self.inner.within(point, radius));
        Ok(PyArray1::from_vec(py, found))
    }

    /// Return ``(index, distance)`` for the point closest to ``point``, or ``None`` if
    /// no point lies within ``max_radius``.
    ///
    /// When several points are equally close, the one with the lowest index is
    /// returned.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If ``max_radius`` is negative or NaN, or ``point`` has a NaN or infinite
    ///     coordinate.
    ///
    /// Examples
    /// --------
    /// >>> model = mmcif_io.read("1crn.bcif").models[0]
    /// >>> grid = mmcif_io.NeighbourGrid.for_model(model)
    /// >>> grid.nearest(model.atoms()[1].position)
    /// (1, 0.0)
    /// >>> grid.nearest([0, 0, 0], max_radius=1.0) is None
    /// True
    #[pyo3(signature = (point, max_radius = f64::INFINITY))]
    fn nearest(
        &self,
        py: Python<'_>,
        point: [f64; 3],
        max_radius: f64,
    ) -> PyResult<Option<(u32, f64)>> {
        let max_radius = checked_radius(max_radius)?;
        let point = query_point(point)?;
        Ok(py.detach(|| self.inner.nearest(point, max_radius)))
    }

    /// Return every pair of points within ``radius`` of each other, as a uint32 NumPy
    /// array of shape ``(n, 2)``.
    ///
    /// Each row is ``(lower, higher)`` by index, and the rows are in ascending order.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If ``radius`` is negative or NaN.
    ///
    /// Examples
    /// --------
    /// >>> model = mmcif_io.read("1crn.bcif").models[0]
    /// >>> grid = mmcif_io.NeighbourGrid.for_model(model)
    /// >>> grid.pairs_within(1.6)[:3].tolist()
    /// [[0, 1], [1, 2], [1, 4]]
    fn pairs_within<'py>(
        &self,
        py: Python<'py>,
        radius: f64,
    ) -> PyResult<Bound<'py, PyArray2<u32>>> {
        let radius = checked_radius(radius)?;
        let pairs = py.detach(|| self.inner.pairs_within(radius));
        pair_array(py, pairs)
    }

    /// Return an ``(other, point)`` index pair for each point in ``others`` that lies
    /// within ``radius`` of a point in the grid, as a uint32 NumPy array of shape
    /// ``(n, 2)``.
    ///
    /// ``others`` is an ``(m, 3)`` float64 NumPy array. Rows are ordered by ``other``
    /// and then by ``point``.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If ``others`` is not an ``(m, 3)`` array, a point in ``others`` has a NaN
    ///     or infinite coordinate, or ``radius`` is negative or NaN.
    ///
    /// Examples
    /// --------
    /// >>> import numpy as np
    /// >>> grid = mmcif_io.NeighbourGrid(np.array([[0.0, 0, 0], [1, 0, 0], [5, 0, 0]]))
    /// >>> grid.contacts(np.array([[4.5, 0, 0], [0.5, 0, 0]]), 1.0).tolist()
    /// [[0, 2], [1, 0], [1, 1]]
    fn contacts<'py>(
        &self,
        py: Python<'py>,
        others: PyReadonlyArray2<'_, f64>,
        radius: f64,
    ) -> PyResult<Bound<'py, PyArray2<u32>>> {
        let others = rows(&others)?;
        let radius = checked_radius(radius)?;
        if let Some(index) = others.iter().position(|point| !is_finite(point)) {
            return Err(geom_error(GeomError::NonFiniteCoordinate { index }));
        }
        let pairs = py
            .detach(|| self.inner.contacts(&others, radius))
            .map_err(geom_error)?;
        pair_array(py, pairs)
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "<NeighbourGrid with {} points, cell size {}>",
            self.inner.len(),
            self.inner.cell_size()
        )
    }
}

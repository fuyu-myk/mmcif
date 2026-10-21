use mmcif_geom::NeighbourGrid;
use numpy::{PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::hierarchy::PyModel;

fn geom_error(error: mmcif_geom::GeomError) -> PyErr {
    PyValueError::new_err(error.to_string())
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

    #[staticmethod]
    #[pyo3(signature = (model, cell_size = 4.0))]
    fn for_model(py: Python<'_>, model: &PyModel, cell_size: f64) -> PyResult<Self> {
        let model = std::sync::Arc::clone(&model.inner);
        let inner = py
            .detach(move || NeighbourGrid::new(model.positions(), cell_size))
            .map_err(geom_error)?;
        Ok(PyNeighbourGrid { inner })
    }

    #[getter]
    fn cell_size(&self) -> f64 {
        self.inner.cell_size()
    }

    fn within<'py>(
        &self,
        py: Python<'py>,
        point: [f64; 3],
        radius: f64,
    ) -> PyResult<Bound<'py, PyArray1<u32>>> {
        let found = py
            .detach(|| self.inner.within(point, radius))
            .map_err(geom_error)?;
        Ok(PyArray1::from_vec(py, found))
    }

    #[pyo3(signature = (point, max_radius = f64::INFINITY))]
    fn nearest(
        &self,
        py: Python<'_>,
        point: [f64; 3],
        max_radius: f64,
    ) -> PyResult<Option<(u32, f64)>> {
        py.detach(|| self.inner.nearest(point, max_radius))
            .map_err(geom_error)
    }

    fn pairs_within<'py>(
        &self,
        py: Python<'py>,
        radius: f64,
    ) -> PyResult<Bound<'py, PyArray2<u32>>> {
        let pairs = py
            .detach(|| self.inner.pairs_within(radius))
            .map_err(geom_error)?;
        pair_array(py, pairs)
    }

    fn contacts<'py>(
        &self,
        py: Python<'py>,
        others: PyReadonlyArray2<'_, f64>,
        radius: f64,
    ) -> PyResult<Bound<'py, PyArray2<u32>>> {
        let others = rows(&others)?;
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

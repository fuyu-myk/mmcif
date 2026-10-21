use std::path::PathBuf;
use std::sync::Arc;

use mmcif_cif::{BlockHandle, CategoryHandle, Document};
use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;

use crate::errors::read_error;

#[pyclass(name = "Document", module = "mmcif_io", frozen)]
pub struct PyDocument {
    pub(crate) inner: Arc<Document>,
}

fn wrap(document: Document) -> PyDocument {
    PyDocument {
        inner: Arc::new(document),
    }
}

#[pymethods]
impl PyDocument {
    #[staticmethod]
    fn parse(py: Python<'_>, text: String) -> PyResult<PyDocument> {
        py.detach(|| Document::parse(text))
            .map(wrap)
            .map_err(read_error)
    }

    #[staticmethod]
    fn open(py: Python<'_>, path: PathBuf) -> PyResult<PyDocument> {
        py.detach(|| Document::open(path))
            .map(wrap)
            .map_err(read_error)
    }

    #[staticmethod]
    fn from_bcif(py: Python<'_>, data: Vec<u8>) -> PyResult<PyDocument> {
        py.detach(|| Document::from_bcif(&data))
            .map(wrap)
            .map_err(read_error)
    }

    #[staticmethod]
    fn open_bcif(py: Python<'_>, path: PathBuf) -> PyResult<PyDocument> {
        py.detach(|| Document::open_bcif(path))
            .map(wrap)
            .map_err(read_error)
    }

    #[getter]
    fn blocks(&self) -> Vec<PyBlock> {
        BlockHandle::all(&self.inner)
            .into_iter()
            .map(|inner| PyBlock { inner })
            .collect()
    }

    fn block(&self, name: &str) -> Option<PyBlock> {
        BlockHandle::all(&self.inner)
            .into_iter()
            .find(|handle| handle.view().name().eq_ignore_ascii_case(name))
            .map(|inner| PyBlock { inner })
    }

    fn to_cif(&self) -> String {
        self.inner.to_cif_string()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __repr__(&self) -> String {
        format!("<Document with {} blocks>", self.inner.len())
    }
}

#[pyclass(name = "Block", module = "mmcif_io", frozen)]
pub struct PyBlock {
    inner: BlockHandle,
}

#[pymethods]
impl PyBlock {
    #[getter]
    fn name(&self) -> String {
        self.inner.view().name().to_owned()
    }

    #[getter]
    fn category_names(&self) -> Vec<String> {
        self.inner
            .view()
            .category_names()
            .into_iter()
            .map(str::to_owned)
            .collect()
    }

    fn category(&self, name: &str) -> Option<PyCategory> {
        self.inner.category(name).map(|inner| PyCategory { inner })
    }

    fn __getitem__(&self, name: &str) -> PyResult<PyCategory> {
        self.category(name)
            .ok_or_else(|| PyKeyError::new_err(name.to_owned()))
    }

    fn __contains__(&self, name: &str) -> bool {
        self.inner.view().category(name).is_some()
    }

    fn __repr__(&self) -> String {
        format!("<Block {}>", self.inner.view().name())
    }
}

#[pyclass(name = "Category", module = "mmcif_io", frozen)]
pub struct PyCategory {
    inner: CategoryHandle,
}

#[pymethods]
impl PyCategory {
    #[getter]
    fn name(&self) -> String {
        self.inner.name().to_owned()
    }

    #[getter]
    fn items(&self) -> Vec<String> {
        self.inner
            .view()
            .items()
            .into_iter()
            .map(str::to_owned)
            .collect()
    }

    #[getter]
    fn row_count(&self) -> usize {
        self.inner.view().row_count()
    }

    fn column(&self, item: &str) -> Option<Vec<Option<String>>> {
        let column = self.inner.view().column(item)?;
        Some(
            column
                .iter()
                .map(|value| value.as_str().map(str::to_owned))
                .collect(),
        )
    }

    fn raw_column(&self, item: &str) -> Option<Vec<String>> {
        let column = self.inner.view().column(item)?;
        Some(column.iter().map(|value| value.raw().to_owned()).collect())
    }

    fn __getitem__(&self, item: &str) -> PyResult<Vec<Option<String>>> {
        self.column(item)
            .ok_or_else(|| PyKeyError::new_err(item.to_owned()))
    }

    fn __contains__(&self, item: &str) -> bool {
        self.inner.view().has(item)
    }

    fn __len__(&self) -> usize {
        self.inner.view().row_count()
    }

    fn __repr__(&self) -> String {
        format!(
            "<Category {} with {} rows>",
            self.inner.name(),
            self.inner.view().row_count()
        )
    }
}

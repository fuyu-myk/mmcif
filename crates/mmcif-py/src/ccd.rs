use std::path::PathBuf;
use std::sync::Arc;

use mmcif_ccd::{AtomStereo, BondStereo, CcdIndex, Component, build_index_from_file};
use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;

use crate::errors::{ccd_error, os_error, read_error};

#[pyclass(name = "CcdIndex", module = "mmcif_io", frozen)]
pub struct PyCcdIndex {
    inner: Arc<CcdIndex>,
}

fn wrap(index: CcdIndex) -> PyCcdIndex {
    PyCcdIndex {
        inner: Arc::new(index),
    }
}

#[pymethods]
impl PyCcdIndex {
    #[staticmethod]
    fn open(py: Python<'_>, path: PathBuf) -> PyResult<PyCcdIndex> {
        py.detach(|| CcdIndex::open(&path))
            .map(wrap)
            .map_err(|error| ccd_error(error, &path))
    }

    #[staticmethod]
    #[pyo3(signature = (source, release, out = None))]
    fn build(
        py: Python<'_>,
        source: PathBuf,
        release: String,
        out: Option<PathBuf>,
    ) -> PyResult<PyCcdIndex> {
        let bytes = py
            .detach(|| {
                let mut bytes = Vec::new();
                build_index_from_file(&source, &release, &mut bytes).map(|_| bytes)
            })
            .map_err(|error| ccd_error(error, &source))?;
        if let Some(out) = out {
            py.detach(|| std::fs::write(&out, &bytes))
                .map_err(|error| os_error(&error, &out))?;
        }
        CcdIndex::from_bytes(bytes).map(wrap).map_err(read_error)
    }

    #[staticmethod]
    fn from_bytes(data: Vec<u8>) -> PyResult<PyCcdIndex> {
        CcdIndex::from_bytes(data).map(wrap).map_err(read_error)
    }

    #[getter]
    fn release(&self) -> String {
        self.inner.release().to_owned()
    }

    fn ids(&self) -> Vec<String> {
        self.inner.ids().map(str::to_owned).collect()
    }

    fn get(&self, id: &str) -> PyResult<Option<PyComponent>> {
        Ok(self
            .inner
            .get(id)
            .map_err(read_error)?
            .map(PyComponent::from))
    }

    fn __getitem__(&self, id: &str) -> PyResult<PyComponent> {
        self.get(id)?
            .ok_or_else(|| PyKeyError::new_err(id.to_owned()))
    }

    fn __contains__(&self, id: &str) -> PyResult<bool> {
        Ok(self.get(id)?.is_some())
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "<CcdIndex {} with {} components>",
            self.inner.release(),
            self.inner.len()
        )
    }
}

type AtomRow = (String, String, i8, bool, bool, &'static str);
type BondRow = (u16, u16, &'static str, bool, &'static str);

#[pyclass(name = "Component", module = "mmcif_io", frozen, get_all)]
pub struct PyComponent {
    id: String,
    name: String,
    kind: String,
    one_letter_code: Option<String>,
    parent: Option<String>,
    formula: Option<String>,
    formal_charge: i32,
    release_status: String,
    atoms: Vec<AtomRow>,
    bonds: Vec<BondRow>,
}

fn atom_stereo(stereo: AtomStereo) -> &'static str {
    match stereo {
        AtomStereo::None => "N",
        AtomStereo::R => "R",
        AtomStereo::S => "S",
    }
}

fn bond_stereo(stereo: BondStereo) -> &'static str {
    match stereo {
        BondStereo::None => "N",
        BondStereo::E => "E",
        BondStereo::Z => "Z",
    }
}

impl From<Component> for PyComponent {
    fn from(component: Component) -> PyComponent {
        PyComponent {
            id: component.id().to_owned(),
            name: component.name().to_owned(),
            kind: component.kind().to_owned(),
            one_letter_code: component.one_letter_code().map(str::to_owned),
            parent: component.parent().map(str::to_owned),
            formula: component.formula().map(str::to_owned),
            formal_charge: component.formal_charge(),
            release_status: component.release_status().to_owned(),
            atoms: component
                .atoms()
                .iter()
                .map(|atom| {
                    (
                        atom.name().to_owned(),
                        atom.element().to_owned(),
                        atom.charge(),
                        atom.is_leaving(),
                        atom.is_aromatic(),
                        atom_stereo(atom.stereo()),
                    )
                })
                .collect(),
            bonds: component
                .bonds()
                .iter()
                .map(|bond| {
                    let [first, second] = bond.atom_indices();
                    (
                        first,
                        second,
                        bond.order().code(),
                        bond.is_aromatic(),
                        bond_stereo(bond.stereo()),
                    )
                })
                .collect(),
        }
    }
}

#[pymethods]
impl PyComponent {
    fn __repr__(&self) -> String {
        format!("<Component {} {}>", self.id, self.name)
    }
}

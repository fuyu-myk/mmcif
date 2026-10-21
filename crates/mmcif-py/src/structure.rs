use std::path::PathBuf;
use std::sync::Arc;

use mmcif_model::{
    AssemblyId, AuthResidueKey, AuthSeqId, Entity, InsCode, LabelResidueKey, LabelSeqId,
    ModelNumber, ReadError, Structure,
};
use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

use crate::document::PyDocument;
use crate::errors::{assembly_error, numbering_error, read_error};
use crate::hierarchy::PyModel;

#[pyclass(name = "Structure", module = "mmcif_io", frozen)]
pub struct PyStructure {
    inner: Arc<Structure>,
}

fn wrap(result: Result<Structure, ReadError>) -> PyResult<PyStructure> {
    result
        .map(|structure| PyStructure {
            inner: Arc::new(structure),
        })
        .map_err(read_error)
}

#[pyfunction]
pub fn read(py: Python<'_>, path: PathBuf) -> PyResult<PyStructure> {
    wrap(py.detach(|| mmcif_model::read(path)))
}

#[pyfunction]
pub fn read_mmcif(py: Python<'_>, path: PathBuf) -> PyResult<PyStructure> {
    wrap(py.detach(|| mmcif_model::read_mmcif(path)))
}

#[pyfunction]
pub fn read_bcif(py: Python<'_>, path: PathBuf) -> PyResult<PyStructure> {
    wrap(py.detach(|| mmcif_model::read_bcif(path)))
}

#[pyfunction]
pub fn read_pdb(py: Python<'_>, path: PathBuf) -> PyResult<PyStructure> {
    wrap(py.detach(|| mmcif_model::read_pdb(path)))
}

#[pyfunction]
pub fn parse_mmcif(py: Python<'_>, text: String) -> PyResult<PyStructure> {
    wrap(py.detach(|| mmcif_model::parse_mmcif(&text)))
}

#[pyfunction]
pub fn parse_bcif(py: Python<'_>, data: Vec<u8>) -> PyResult<PyStructure> {
    wrap(py.detach(|| mmcif_model::parse_bcif(&data)))
}

#[pyfunction]
pub fn parse_pdb(py: Python<'_>, text: String) -> PyResult<PyStructure> {
    wrap(py.detach(|| mmcif_model::parse_pdb(&text)))
}

#[pymethods]
impl PyStructure {
    #[staticmethod]
    fn from_document(py: Python<'_>, document: &PyDocument) -> PyResult<PyStructure> {
        let document = Arc::clone(&document.inner);
        wrap(py.detach(move || Structure::from_document(&document)))
    }

    #[getter]
    fn entry_id(&self) -> Option<String> {
        self.inner.entry_id().map(str::to_owned)
    }

    #[getter]
    fn models(&self) -> Vec<PyModel> {
        self.inner
            .models()
            .iter()
            .map(|model| PyModel {
                inner: Arc::clone(model),
            })
            .collect()
    }

    fn model(&self, number: i32) -> Option<PyModel> {
        self.inner
            .model(ModelNumber::new(number))
            .map(|model| PyModel {
                inner: Arc::clone(model),
            })
    }

    #[getter]
    fn entities(&self) -> Vec<PyEntity> {
        self.inner.entities().iter().map(PyEntity::from).collect()
    }

    fn entity(&self, id: &str) -> Option<PyEntity> {
        self.inner.entity(id).map(PyEntity::from)
    }

    #[getter]
    fn assemblies(&self) -> Vec<PyAssembly> {
        self.inner
            .assemblies()
            .iter()
            .map(|assembly| PyAssembly {
                id: assembly.id().to_string(),
                details: assembly.details().map(str::to_owned),
                operator_count: assembly.operator_count(),
                generators: assembly
                    .generators()
                    .iter()
                    .map(|generator| {
                        (
                            generator.expression().to_owned(),
                            generator
                                .asym_ids()
                                .iter()
                                .map(|asym| asym.to_string())
                                .collect(),
                        )
                    })
                    .collect(),
            })
            .collect()
    }

    #[pyo3(signature = (assembly_id, model_index = 0))]
    fn assemble(&self, py: Python<'_>, assembly_id: &str, model_index: usize) -> PyResult<PyModel> {
        let model = self
            .inner
            .models()
            .get(model_index)
            .cloned()
            .ok_or_else(|| PyIndexError::new_err(format!("no model at index {model_index}")))?;
        let structure = Arc::clone(&self.inner);
        let id = AssemblyId::new(assembly_id);
        let assembled = py
            .detach(move || structure.assemble(&id, &model))
            .map_err(assembly_error)?;
        Ok(PyModel {
            inner: Arc::new(assembled),
        })
    }

    fn to_auth(&self, label_asym: &str, label_seq: u32) -> PyResult<(String, i32, Option<String>)> {
        let seq = LabelSeqId::new(label_seq)
            .ok_or_else(|| numbering_error("label_seq_id starts at 1"))?;
        let key = LabelResidueKey {
            asym: label_asym.into(),
            seq,
        };
        let auth = self
            .inner
            .numbering()
            .to_auth(&key)
            .map_err(numbering_error)?;
        Ok((
            auth.chain.to_string(),
            auth.seq.get(),
            auth.ins_code.get().map(String::from),
        ))
    }

    #[pyo3(signature = (chain, seq, ins_code = None))]
    fn to_label(&self, chain: &str, seq: i32, ins_code: Option<char>) -> PyResult<(String, u32)> {
        let key = AuthResidueKey {
            chain: chain.into(),
            seq: AuthSeqId::new(seq),
            ins_code: ins_code.map_or(InsCode::NONE, InsCode::new),
        };
        let label = self
            .inner
            .numbering()
            .to_label(&key)
            .map_err(numbering_error)?;
        Ok((label.asym.to_string(), label.seq.get()))
    }

    fn __repr__(&self) -> String {
        format!(
            "<Structure {} with {} models>",
            self.inner.entry_id().unwrap_or("?"),
            self.inner.models().len()
        )
    }
}

#[pyclass(name = "Entity", module = "mmcif_io", frozen, get_all)]
pub struct PyEntity {
    id: String,
    kind: String,
    description: Option<String>,
    polymer_type: Option<String>,
    canonical_sequence: Option<String>,
    sequence: Vec<Vec<String>>,
    asym_ids: Vec<String>,
}

impl From<&Entity> for PyEntity {
    fn from(entity: &Entity) -> PyEntity {
        let polymer = entity.polymer();
        PyEntity {
            id: entity.id().to_string(),
            kind: entity.kind().as_mmcif_str().to_owned(),
            description: entity.description().map(str::to_owned),
            polymer_type: polymer.map(|polymer| polymer.kind().as_mmcif_str().to_owned()),
            canonical_sequence: polymer
                .and_then(|polymer| polymer.canonical_sequence().map(str::to_owned)),
            sequence: polymer
                .map(|polymer| {
                    polymer
                        .sequence()
                        .iter()
                        .map(|position| {
                            position
                                .comps()
                                .iter()
                                .map(|comp| comp.to_string())
                                .collect()
                        })
                        .collect()
                })
                .unwrap_or_default(),
            asym_ids: entity
                .asym_ids()
                .iter()
                .map(|asym| asym.to_string())
                .collect(),
        }
    }
}

#[pymethods]
impl PyEntity {
    fn __repr__(&self) -> String {
        format!("<Entity {} {}>", self.id, self.kind)
    }
}

#[pyclass(name = "Assembly", module = "mmcif_io", frozen, get_all)]
pub struct PyAssembly {
    id: String,
    details: Option<String>,
    operator_count: usize,
    generators: Vec<(String, Vec<String>)>,
}

#[pymethods]
impl PyAssembly {
    fn __repr__(&self) -> String {
        format!(
            "<Assembly {} with {} operators>",
            self.id, self.operator_count
        )
    }
}

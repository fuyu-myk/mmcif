use std::path::{Path, PathBuf};
use std::sync::Arc;

use mmcif_model::{
    AssemblyId, AuthResidueKey, AuthSeqId, Entity, InsCode, LabelResidueKey, LabelSeqId,
    ModelNumber, ReadError, Structure,
};
use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

use crate::document::PyDocument;
use crate::errors::{assembly_error, model_error, numbering_error, read_error};
use crate::hierarchy::PyModel;

/// A macromolecular structure: its models, entities, biological assemblies and
/// residue numbering.
///
/// Read one with ``read`` or one of the other readers.
///
/// Examples
/// --------
/// >>> structure = mmcif_io.read("1crn.bcif")
/// >>> structure.entry_id
/// '1CRN'
/// >>> structure.models
/// [<Model 1 with 1 chains, 327 atoms>]
/// >>> structure.entities
/// [<Entity 1 polymer>]
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

fn wrap_at(result: Result<Structure, ReadError>, path: &Path) -> PyResult<PyStructure> {
    result
        .map(|structure| PyStructure {
            inner: Arc::new(structure),
        })
        .map_err(|error| model_error(error, path))
}

/// Read the structure in the file at ``path``.
///
/// The format follows the file extension: ``.cif`` and ``.mmcif`` are read as
/// mmCIF, ``.bcif`` as BinaryCIF, and ``.pdb``, ``.ent`` and assembly files such
/// as ``.pdb1`` as PDB. A gzip-compressed file, such as ``1crn.cif.gz``, is
/// decompressed. ``path`` is a string or a path-like object.
///
/// Raises
/// ------
/// OSError
///     If the file cannot be read.
/// ReadError
///     If the extension is not one of these, or the contents are not valid for
///     the format.
///
/// Examples
/// --------
/// >>> mmcif_io.read("1crn.bcif")
/// <Structure 1CRN with 1 models>
/// >>> mmcif_io.read("example.pdb").entry_id
/// '1ABC'
#[pyfunction]
pub fn read(py: Python<'_>, path: PathBuf) -> PyResult<PyStructure> {
    wrap_at(py.detach(|| mmcif_model::read(&path)), &path)
}

/// Read the structure in the mmCIF file at ``path``.
///
/// Raises
/// ------
/// OSError
///     If the file cannot be read.
/// ReadError
///     If the contents are not valid for the format, or lack the categories a
///     structure needs.
///
/// Examples
/// --------
/// >>> structure = mmcif_io.read_mmcif("1crn.cif")
/// >>> structure.entry_id
/// '1CRN'
#[pyfunction]
pub fn read_mmcif(py: Python<'_>, path: PathBuf) -> PyResult<PyStructure> {
    wrap_at(py.detach(|| mmcif_model::read_mmcif(&path)), &path)
}

/// Read the structure in the BinaryCIF file at ``path``.
///
/// Raises
/// ------
/// OSError
///     If the file cannot be read.
/// ReadError
///     If the contents are not valid for the format, or lack the categories a
///     structure needs.
///
/// Examples
/// --------
/// >>> structure = mmcif_io.read_bcif("1crn.bcif")
/// >>> structure.models[0].atom_count
/// 327
#[pyfunction]
pub fn read_bcif(py: Python<'_>, path: PathBuf) -> PyResult<PyStructure> {
    wrap_at(py.detach(|| mmcif_model::read_bcif(&path)), &path)
}

/// Read the structure in the PDB file at ``path``.
///
/// Raises
/// ------
/// OSError
///     If the file cannot be read.
/// ReadError
///     If the contents are not valid for the format, or lack the categories a
///     structure needs.
///
/// Examples
/// --------
/// >>> structure = mmcif_io.read_pdb("example.pdb")
/// >>> len(structure.models)
/// 2
#[pyfunction]
pub fn read_pdb(py: Python<'_>, path: PathBuf) -> PyResult<PyStructure> {
    wrap_at(py.detach(|| mmcif_model::read_pdb(&path)), &path)
}

/// Parse a structure from mmCIF text.
///
/// Raises
/// ------
/// ReadError
///     If ``text`` is not valid mmCIF, or lacks the categories a structure
///     needs.
///
/// Examples
/// --------
/// >>> with open("1crn.cif") as file:
/// ...     structure = mmcif_io.parse_mmcif(file.read())
/// >>> structure.entry_id
/// '1CRN'
#[pyfunction]
pub fn parse_mmcif(py: Python<'_>, text: String) -> PyResult<PyStructure> {
    wrap(py.detach(|| mmcif_model::parse_mmcif(&text)))
}

/// Parse a structure from BinaryCIF bytes.
///
/// Raises
/// ------
/// ReadError
///     If ``data`` is not valid BinaryCIF, or lacks the categories a structure
///     needs.
///
/// Examples
/// --------
/// >>> with open("1crn.bcif", "rb") as file:
/// ...     structure = mmcif_io.parse_bcif(file.read())
/// >>> structure.entry_id
/// '1CRN'
#[pyfunction]
pub fn parse_bcif(py: Python<'_>, data: Vec<u8>) -> PyResult<PyStructure> {
    wrap(py.detach(|| mmcif_model::parse_bcif(&data)))
}

/// Parse a structure from PDB text.
///
/// Raises
/// ------
/// ReadError
///     If ``text`` is not valid PDB.
///
/// Examples
/// --------
/// >>> with open("example.pdb") as file:
/// ...     structure = mmcif_io.parse_pdb(file.read())
/// >>> structure.entry_id
/// '1ABC'
#[pyfunction]
pub fn parse_pdb(py: Python<'_>, text: String) -> PyResult<PyStructure> {
    wrap(py.detach(|| mmcif_model::parse_pdb(&text)))
}

#[pymethods]
impl PyStructure {
    /// Build a structure from a parsed CIF document.
    ///
    /// The first data block of ``document`` is read.
    ///
    /// Raises
    /// ------
    /// ReadError
    ///     If the document has no data block, or lacks the categories a structure
    ///     needs.
    ///
    /// Examples
    /// --------
    /// >>> document = mmcif_io.Document.open("1crn.cif")
    /// >>> mmcif_io.Structure.from_document(document).entry_id
    /// '1CRN'
    #[staticmethod]
    fn from_document(py: Python<'_>, document: &PyDocument) -> PyResult<PyStructure> {
        let document = Arc::clone(&document.inner);
        wrap(py.detach(move || Structure::from_document(&document)))
    }

    /// The ``_entry.id`` of the structure, or ``None`` if the file gives none.
    #[getter]
    fn entry_id(&self) -> Option<String> {
        self.inner.entry_id().map(str::to_owned)
    }

    /// The models of the structure, in file order.
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

    /// Return the model with model number ``number``, or ``None``.
    ///
    /// Examples
    /// --------
    /// >>> structure = mmcif_io.read("example.pdb")
    /// >>> structure.model(2).number
    /// 2
    /// >>> structure.model(3) is None
    /// True
    fn model(&self, number: i32) -> Option<PyModel> {
        self.inner
            .model(ModelNumber::new(number))
            .map(|model| PyModel {
                inner: Arc::clone(model),
            })
    }

    /// The entities of the structure, in file order.
    #[getter]
    fn entities(&self) -> Vec<PyEntity> {
        self.inner.entities().iter().map(PyEntity::from).collect()
    }

    /// Return the entity with id ``id``, or ``None``.
    ///
    /// Examples
    /// --------
    /// >>> structure = mmcif_io.read("1crn.bcif")
    /// >>> structure.entity("1").description
    /// 'CRAMBIN'
    /// >>> structure.entity("9") is None
    /// True
    fn entity(&self, id: &str) -> Option<PyEntity> {
        self.inner.entity(id).map(PyEntity::from)
    }

    /// The biological assemblies of the structure, in file order.
    #[getter]
    fn assemblies(&self) -> Vec<PyAssembly> {
        self.inner
            .assemblies()
            .iter()
            .map(|assembly| PyAssembly {
                id: assembly.id().to_string(),
                details: assembly.details().map(str::to_owned),
                oligomeric_details: assembly.oligomeric_details().map(str::to_owned),
                oligomeric_count: assembly.oligomeric_count(),
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

    /// Return a new model holding the assembly ``assembly_id``, built from the model
    /// at index ``model_index`` in ``models``.
    ///
    /// Each chain of the new model records the operators that placed it in
    /// ``Chain.operators``.
    ///
    /// Raises
    /// ------
    /// IndexError
    ///     If there is no model at ``model_index``.
    /// AssemblyError
    ///     If the structure has no assembly ``assembly_id``, or it cannot be built.
    ///
    /// Examples
    /// --------
    /// >>> structure = mmcif_io.read("1crn.bcif")
    /// >>> assembled = structure.assemble("1")
    /// >>> assembled.chains[0].operators
    /// ['1']
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

    /// Return ``True`` if assembly ``assembly_id`` shows the deposited chains in
    /// their deposited arrangement.
    ///
    /// That is, it places no chain more than once, and places every chain by the
    /// same transform, to within ``tolerance`` in each rotation and translation
    /// element. Such an assembly may hold fewer chains than the deposited model, or
    /// move all of them together, but it shows nothing the deposited model does
    /// not.
    ///
    /// Raises
    /// ------
    /// AssemblyError
    ///     If the structure has no assembly ``assembly_id``, or it names an operator
    ///     the structure does not define.
    ///
    /// Examples
    /// --------
    /// >>> structure = mmcif_io.read("example.cif")
    /// >>> structure.assembly_is_deposited("1")
    /// True
    /// >>> structure.assembly_is_deposited("2")
    /// False
    #[pyo3(signature = (assembly_id, tolerance = 1e-3))]
    fn assembly_is_deposited(&self, assembly_id: &str, tolerance: f64) -> PyResult<bool> {
        self.inner
            .assembly_is_deposited(&AssemblyId::new(assembly_id), tolerance)
            .map_err(assembly_error)
    }

    /// Return the author address ``(chain, seq, ins_code)`` of the polymer residue
    /// at label address ``(label_asym, label_seq)``.
    ///
    /// ``ins_code`` is ``None`` when the residue has no insertion code.
    ///
    /// Raises
    /// ------
    /// NumberingError
    ///     If ``label_seq`` is 0, or no polymer residue has that label address.
    ///
    /// Examples
    /// --------
    /// >>> structure = mmcif_io.read("example.cif")
    /// >>> structure.to_auth("AA", 3)
    /// ('H1', 100, 'A')
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

    /// Return the label address ``(label_asym, label_seq)`` of the polymer residue
    /// at author address ``(chain, seq, ins_code)``.
    ///
    /// Raises
    /// ------
    /// NumberingError
    ///     If no polymer residue has that author address.
    ///
    /// Examples
    /// --------
    /// >>> structure = mmcif_io.read("example.cif")
    /// >>> structure.to_label("H1", 100, "A")
    /// ('AA', 3)
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

/// An entity of a structure, from ``_entity``.
///
/// Examples
/// --------
/// >>> entity = mmcif_io.read("1crn.bcif").entities[0]
/// >>> entity.kind, entity.polymer_type, entity.description
/// ('polymer', 'polypeptide(L)', 'CRAMBIN')
/// >>> entity.sequence[:3]
/// [['THR'], ['THR'], ['CYS']]
#[pyclass(name = "Entity", module = "mmcif_io", frozen, get_all)]
pub struct PyEntity {
    /// The ``_entity.id``.
    id: String,
    /// The entity type: ``"polymer"``, ``"non-polymer"``, ``"branched"``,
    /// ``"water"`` or ``"macrolide"``.
    kind: String,
    /// The ``_entity.pdbx_description``, or ``None``.
    description: Option<String>,
    /// The ``_entity_poly.type``, such as ``"polypeptide(L)"``, or ``None`` if
    /// the entity is not a polymer.
    polymer_type: Option<String>,
    /// The canonical one-letter sequence, from
    /// ``_entity_poly.pdbx_seq_one_letter_code_can``, or ``None``.
    canonical_sequence: Option<String>,
    /// The component ids at each position of the polymer sequence. A position
    /// holds more than one id where the sequence is heterogeneous. Empty if
    /// the entity is not a polymer.
    sequence: Vec<Vec<String>>,
    /// The label asym ids of the chains that belong to the entity.
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

/// A biological assembly, from ``_pdbx_struct_assembly`` and
/// ``_pdbx_struct_assembly_gen``.
///
/// Build its coordinates with ``Structure.assemble``.
///
/// Examples
/// --------
/// >>> assembly = mmcif_io.read("1crn.bcif").assemblies[0]
/// >>> assembly.id, assembly.details
/// ('1', 'author_defined_assembly')
/// >>> assembly.oligomeric_details, assembly.oligomeric_count
/// ('monomeric', 1)
/// >>> assembly.generators
/// [('1', ['A'])]
#[pyclass(name = "Assembly", module = "mmcif_io", frozen, get_all)]
pub struct PyAssembly {
    /// The ``_pdbx_struct_assembly.id``.
    id: String,
    /// The ``_pdbx_struct_assembly.details``, or ``None``.
    details: Option<String>,
    /// The oligomeric state from ``_pdbx_struct_assembly.oligomeric_details``,
    /// such as ``"dimeric"``, or ``None``.
    oligomeric_details: Option<String>,
    /// The number of polymer chains in the assembly from
    /// ``_pdbx_struct_assembly.oligomeric_count``, or ``None``.
    oligomeric_count: Option<u32>,
    /// The number of operator combinations across all generators.
    operator_count: usize,
    /// The ``(operator expression, label asym ids)`` pairs from
    /// ``_pdbx_struct_assembly_gen``, in file order.
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

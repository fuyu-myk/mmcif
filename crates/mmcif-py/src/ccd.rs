use std::path::PathBuf;
use std::sync::Arc;

use mmcif_ccd::{AtomStereo, BondStereo, CcdIndex, Component, build_index_from_file};
use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;

use crate::errors::{ccd_error, os_error, read_error};

/// A compiled Chemical Component Dictionary index that returns components by id.
///
/// The index records the name of the dictionary release it was built from.
/// ``index[id]`` returns the component ``id`` as ``get`` does and raises
/// ``KeyError`` if there is none, ``id in index`` tests for one, and
/// ``len(index)`` is the number of components.
///
/// Examples
/// --------
/// >>> index = mmcif_io.CcdIndex.build("components.cif", "2026-10-20")
/// >>> index
/// <CcdIndex 2026-10-20 with 8 components>
/// >>> index["ALA"].name
/// 'ALANINE'
/// >>> "XYZ" in index
/// False
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
    /// Read the index file at ``path``, as written by ``build``.
    ///
    /// Raises
    /// ------
    /// OSError
    ///     If the file cannot be read.
    /// ReadError
    ///     If the contents are not a component index.
    ///
    /// Examples
    /// --------
    /// >>> mmcif_io.CcdIndex.build("components.cif", "2026-10-20", out="ccd.idx")
    /// <CcdIndex 2026-10-20 with 8 components>
    /// >>> mmcif_io.CcdIndex.open("ccd.idx")
    /// <CcdIndex 2026-10-20 with 8 components>
    #[staticmethod]
    fn open(py: Python<'_>, path: PathBuf) -> PyResult<PyCcdIndex> {
        py.detach(|| CcdIndex::open(&path))
            .map(wrap)
            .map_err(|error| ccd_error(error, &path))
    }

    /// Compile the Chemical Component Dictionary file at ``source`` into an index for
    /// the release named ``release``.
    ///
    /// ``source`` is CIF with one data block per component, as in ``components.cif``.
    /// If ``out`` is given, the index is also written to that path for ``open`` to
    /// read.
    ///
    /// Raises
    /// ------
    /// OSError
    ///     If ``source`` cannot be read or is not valid UTF-8, or ``out`` cannot be
    ///     written.
    /// ReadError
    ///     If ``source`` is not valid CIF, a block cannot be read as a component, a
    ///     component id is not 1 to 8 ASCII characters without whitespace, or two
    ///     blocks have the same component id.
    ///
    /// Examples
    /// --------
    /// >>> index = mmcif_io.CcdIndex.build("components.cif", "2026-10-20")
    /// >>> index.release
    /// '2026-10-20'
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

    /// Read an index from ``data``, the contents of a file written by ``build``.
    ///
    /// Raises
    /// ------
    /// ReadError
    ///     If ``data`` is not a component index or was written in a format version
    ///     this release does not read.
    ///
    /// Examples
    /// --------
    /// >>> mmcif_io.CcdIndex.build("components.cif", "2026-10-20", out="ccd.idx")
    /// <CcdIndex 2026-10-20 with 8 components>
    /// >>> with open("ccd.idx", "rb") as file:
    /// ...     index = mmcif_io.CcdIndex.from_bytes(file.read())
    /// >>> len(index)
    /// 8
    #[staticmethod]
    fn from_bytes(data: Vec<u8>) -> PyResult<PyCcdIndex> {
        CcdIndex::from_bytes(data).map(wrap).map_err(read_error)
    }

    /// The name of the dictionary release given when the index was built.
    #[getter]
    fn release(&self) -> String {
        self.inner.release().to_owned()
    }

    /// Return the component ids in the index, in ascending byte order.
    ///
    /// Examples
    /// --------
    /// >>> index = mmcif_io.CcdIndex.build("components.cif", "2026-10-20")
    /// >>> index.ids()[:3]
    /// ['A1AAA', 'ALA', 'ATP']
    fn ids(&self) -> Vec<String> {
        self.inner.ids().map(str::to_owned).collect()
    }

    /// Return the component with id ``id``, or ``None`` if the index does not hold
    /// it.
    ///
    /// Ids are matched exactly, including case.
    ///
    /// Raises
    /// ------
    /// ReadError
    ///     If the index record for ``id`` is corrupt.
    ///
    /// Examples
    /// --------
    /// >>> index = mmcif_io.CcdIndex.build("components.cif", "2026-10-20")
    /// >>> index.get("MSE").parent
    /// 'MET'
    /// >>> index.get("mse") is None
    /// True
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

/// A chemical component from the dictionary, with its atoms and the bonds between
/// them.
///
/// Examples
/// --------
/// >>> component = mmcif_io.CcdIndex.build("components.cif", "2026-10-20")["ALA"]
/// >>> component
/// <Component ALA ALANINE>
/// >>> component.kind, component.formula
/// ('L-PEPTIDE LINKING', 'C3 H7 N O2')
/// >>> component.atoms[1]
/// ('CA', 'C', 0, False, False, 'S')
/// >>> component.bonds[0]
/// (0, 1, 'SING', False, 'N')
#[pyclass(name = "Component", module = "mmcif_io", frozen, get_all)]
pub struct PyComponent {
    /// The component id, such as ``ALA``.
    id: String,
    /// The component name, such as ``ALANINE``, or an empty string if the
    /// dictionary gives none.
    name: String,
    /// The component type, such as ``L-PEPTIDE LINKING``, or an empty string if the
    /// dictionary gives none.
    kind: String,
    /// The one-letter code, such as ``A`` for alanine, or ``None``.
    one_letter_code: Option<String>,
    /// The id of the parent component, such as ``MET`` for selenomethionine, or
    /// ``None``.
    parent: Option<String>,
    /// The chemical formula, such as ``C3 H7 N O2``, or ``None``.
    formula: Option<String>,
    /// The formal charge of the component, or 0 if the dictionary gives none.
    formal_charge: i32,
    /// The release status, such as ``REL``, or an empty string if the dictionary
    /// gives none.
    release_status: String,
    /// The atoms, in the order the dictionary lists them, as
    /// ``(name, element, charge, is_leaving, is_aromatic, stereo)`` tuples.
    ///
    /// ``element`` is as the dictionary writes it, such as ``C`` or ``FE``.
    /// ``charge`` is 0 if the dictionary gives none, and ``stereo`` is ``"R"``,
    /// ``"S"`` or ``"N"`` for none.
    atoms: Vec<AtomRow>,
    /// The bonds, in the order the dictionary lists them, as
    /// ``(first, second, order, is_aromatic, stereo)`` tuples.
    ///
    /// ``first`` and ``second`` are positions in ``atoms``. ``order`` is the
    /// ``_chem_comp_bond.value_order`` code: ``SING``, ``DOUB``, ``TRIP``, ``QUAD``,
    /// ``AROM``, ``DELO``, ``PI`` or ``POLY``. ``stereo`` is ``"E"``, ``"Z"`` or
    /// ``"N"`` for none.
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

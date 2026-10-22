use std::path::PathBuf;
use std::sync::Arc;

use mmcif_cif::{BlockHandle, CategoryHandle, Document};
use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;

use crate::errors::{cif_error, read_error};

/// A CIF document: data blocks of categories, each a table of named items.
///
/// ``len(document)`` is the number of data blocks.
///
/// Examples
/// --------
/// >>> document = mmcif_io.Document.open("1crn.cif")
/// >>> document
/// <Document with 1 blocks>
/// >>> block = document.blocks[0]
/// >>> block["atom_site"].row_count
/// 327
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
    /// Parse a CIF document from text.
    ///
    /// Raises
    /// ------
    /// ReadError
    ///     If ``text`` is not valid CIF.
    ///
    /// Examples
    /// --------
    /// >>> document = mmcif_io.Document.parse("data_x\n_cell.length_a 10.0\n")
    /// >>> document.blocks[0]["cell"]["length_a"]
    /// ['10.0']
    #[staticmethod]
    fn parse(py: Python<'_>, text: String) -> PyResult<PyDocument> {
        py.detach(|| Document::parse(text))
            .map(wrap)
            .map_err(read_error)
    }

    /// Read the CIF document in the file at ``path``.
    ///
    /// Raises
    /// ------
    /// OSError
    ///     If the file cannot be read.
    /// ReadError
    ///     If the contents are not valid CIF.
    ///
    /// Examples
    /// --------
    /// >>> mmcif_io.Document.open("1crn.cif").blocks[0].name
    /// '1CRN'
    #[staticmethod]
    fn open(py: Python<'_>, path: PathBuf) -> PyResult<PyDocument> {
        py.detach(|| Document::open(&path))
            .map(wrap)
            .map_err(|error| cif_error(error, &path))
    }

    /// Decode a CIF document from BinaryCIF bytes.
    ///
    /// Raises
    /// ------
    /// ReadError
    ///     If ``data`` is not valid BinaryCIF.
    ///
    /// Examples
    /// --------
    /// >>> with open("1crn.bcif", "rb") as file:
    /// ...     document = mmcif_io.Document.from_bcif(file.read())
    /// >>> document.blocks[0].name
    /// '1CRN'
    #[staticmethod]
    fn from_bcif(py: Python<'_>, data: Vec<u8>) -> PyResult<PyDocument> {
        py.detach(|| Document::from_bcif(&data))
            .map(wrap)
            .map_err(read_error)
    }

    /// Read the BinaryCIF document in the file at ``path``.
    ///
    /// Raises
    /// ------
    /// OSError
    ///     If the file cannot be read.
    /// ReadError
    ///     If the contents are not valid BinaryCIF.
    ///
    /// Examples
    /// --------
    /// >>> mmcif_io.Document.open_bcif("1crn.bcif").blocks[0].name
    /// '1CRN'
    #[staticmethod]
    fn open_bcif(py: Python<'_>, path: PathBuf) -> PyResult<PyDocument> {
        py.detach(|| Document::open_bcif(&path))
            .map(wrap)
            .map_err(|error| cif_error(error, &path))
    }

    /// The data blocks of the document, in file order.
    #[getter]
    fn blocks(&self) -> Vec<PyBlock> {
        BlockHandle::all(&self.inner)
            .into_iter()
            .map(|inner| PyBlock { inner })
            .collect()
    }

    /// Return the data block named ``name``, or ``None``.
    ///
    /// Names are compared without regard to ASCII case.
    ///
    /// Examples
    /// --------
    /// >>> document = mmcif_io.Document.open("1crn.cif")
    /// >>> document.block("1crn").name
    /// '1CRN'
    /// >>> document.block("2abc") is None
    /// True
    fn block(&self, name: &str) -> Option<PyBlock> {
        BlockHandle::all(&self.inner)
            .into_iter()
            .find(|handle| handle.view().name().eq_ignore_ascii_case(name))
            .map(|inner| PyBlock { inner })
    }

    /// Return the document written as CIF text.
    ///
    /// Parsing the text gives a document with the same blocks, categories and values.
    ///
    /// Examples
    /// --------
    /// >>> document = mmcif_io.Document.parse("data_x\n_cell.length_a 10.0\n")
    /// >>> print(document.to_cif())
    /// data_x
    /// #
    /// _cell.length_a 10.0
    /// #
    /// <BLANKLINE>
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

/// A data block of a CIF document.
///
/// ``block[name]`` returns the category ``name`` and raises ``KeyError`` if there
/// is none, and ``name in block`` tests for one. Category names may be given
/// with or without the leading underscore.
///
/// Examples
/// --------
/// >>> block = mmcif_io.Document.open("1crn.cif").blocks[0]
/// >>> "atom_site" in block
/// True
/// >>> block["_entry"]["id"]
/// ['1CRN']
#[pyclass(name = "Block", module = "mmcif_io", frozen)]
pub struct PyBlock {
    inner: BlockHandle,
}

#[pymethods]
impl PyBlock {
    /// The name of the block, without the ``data_`` prefix.
    #[getter]
    fn name(&self) -> String {
        self.inner.view().name().to_owned()
    }

    /// The names of the categories in the block, in file order, without the leading
    /// underscore.
    #[getter]
    fn category_names(&self) -> Vec<String> {
        self.inner
            .view()
            .category_names()
            .into_iter()
            .map(str::to_owned)
            .collect()
    }

    /// Return the category ``name``, or ``None``.
    ///
    /// ``name`` may be given with or without the leading underscore.
    ///
    /// Examples
    /// --------
    /// >>> block = mmcif_io.Document.open("1crn.cif").blocks[0]
    /// >>> block.category("cell").name
    /// 'cell'
    /// >>> block.category("nothing") is None
    /// True
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

/// A category of a data block: a table whose columns are items.
///
/// ``category[item]`` returns the column ``item`` as ``column`` does and raises
/// ``KeyError`` if there is none, ``item in category`` tests for one, and
/// ``len(category)`` is the number of rows.
///
/// Examples
/// --------
/// >>> category = mmcif_io.Document.open("1crn.cif").blocks[0]["atom_site"]
/// >>> category
/// <Category atom_site with 327 rows>
/// >>> category["label_comp_id"][:3]
/// ['THR', 'THR', 'THR']
#[pyclass(name = "Category", module = "mmcif_io", frozen)]
pub struct PyCategory {
    inner: CategoryHandle,
}

#[pymethods]
impl PyCategory {
    /// The name of the category, without the leading underscore.
    #[getter]
    fn name(&self) -> String {
        self.inner.name().to_owned()
    }

    /// The names of the items in the category, in file order, without the category
    /// prefix.
    #[getter]
    fn items(&self) -> Vec<String> {
        self.inner
            .view()
            .items()
            .into_iter()
            .map(str::to_owned)
            .collect()
    }

    /// The number of rows in the category.
    #[getter]
    fn row_count(&self) -> usize {
        self.inner.view().row_count()
    }

    /// Return the values of ``item`` in row order, or ``None`` if the category has no
    /// such item.
    ///
    /// Values that are ``?`` or ``.`` in the file are returned as ``None``.
    ///
    /// Examples
    /// --------
    /// >>> category = mmcif_io.Document.open("1crn.cif").blocks[0]["atom_site"]
    /// >>> category.column("label_seq_id")[:3]
    /// ['1', '1', '1']
    /// >>> category.column("pdbx_PDB_ins_code")[:3]
    /// [None, None, None]
    fn column(&self, item: &str) -> Option<Vec<Option<String>>> {
        let column = self.inner.view().column(item)?;
        Some(
            column
                .iter()
                .map(|value| value.as_str().map(str::to_owned))
                .collect(),
        )
    }

    /// Return the values of ``item`` in row order as written in the file, or
    /// ``None`` if the category has no such item.
    ///
    /// ``?`` and ``.`` are returned as the strings ``"?"`` and ``"."``.
    ///
    /// Examples
    /// --------
    /// >>> category = mmcif_io.Document.open("1crn.cif").blocks[0]["atom_site"]
    /// >>> category.raw_column("pdbx_PDB_ins_code")[:3]
    /// ['?', '?', '?']
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

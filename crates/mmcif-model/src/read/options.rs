use std::path::Path;

use mmcif_cif::Document;

use crate::error::ReadError;
use crate::pdb;
use crate::structure::Structure;

/// How [`read_with`] and [`Structure::from_block_with`] read a structure.
///
/// The default reads every model and fails on any malformed category, as
/// [`read`](crate::read) does.
///
/// # Examples
///
/// ```
/// use mmcif_model::{ReadOptions, read_with};
///
/// let options = ReadOptions::new().first_model_only(true).lenient(true);
/// let loaded = read_with("tests/fixtures/pathologies.cif", &options)?;
///
/// assert_eq!(loaded.structure.models().len(), 1);
/// assert!(loaded.skipped.is_empty());
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ReadOptions {
    first_model_only: bool,
    lenient: bool,
}

impl ReadOptions {
    /// Creates the default options.
    pub fn new() -> ReadOptions {
        ReadOptions::default()
    }

    /// Sets whether to read only the first model, skipping the atoms of every
    /// later one.
    pub fn first_model_only(self, first_model_only: bool) -> ReadOptions {
        ReadOptions {
            first_model_only,
            ..self
        }
    }

    /// Sets whether a malformed optional group of categories is skipped
    /// rather than failing the read.
    ///
    /// The optional groups are the entities, the numbering schemes
    /// `_pdbx_poly_seq_scheme` and `_pdbx_branch_scheme`, the operators and
    /// the assemblies, and in a PDB file the `REMARK 350` assemblies. A
    /// skipped group is left empty, a skipped numbering scheme leaves the
    /// numbering that the atoms give, and its error is reported in
    /// [`Loaded::skipped`]. The atoms are never skipped.
    pub fn lenient(self, lenient: bool) -> ReadOptions {
        ReadOptions { lenient, ..self }
    }

    pub(crate) fn reads_first_model_only(&self) -> bool {
        self.first_model_only
    }

    pub(crate) fn skipping(&self) -> Skipping {
        Skipping {
            lenient: self.lenient,
            errors: Vec::new(),
        }
    }
}

/// A structure read with [`ReadOptions`].
#[derive(Debug)]
#[non_exhaustive]
pub struct Loaded {
    /// The structure.
    pub structure: Structure,
    /// The document the structure was read from, when [`read_with`] read an
    /// mmCIF or BinaryCIF file.
    pub document: Option<Document>,
    /// The errors of the optional groups a lenient read skipped, in the order
    /// they were read.
    pub skipped: Vec<ReadError>,
}

impl Loaded {
    pub(crate) fn new(structure: Structure, skipped: Vec<ReadError>) -> Loaded {
        Loaded {
            structure,
            document: None,
            skipped,
        }
    }
}

/// A structure file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Format {
    /// mmCIF text.
    Mmcif,
    /// BinaryCIF.
    Bcif,
    /// The legacy PDB format.
    Pdb,
}

impl Format {
    /// Returns the format a file's extension names, ignoring case, or `None`.
    ///
    /// `cif` and `mmcif` are mmCIF, `bcif` is BinaryCIF, and `pdb`, `ent` and
    /// the assembly extensions `pdb1`, `pdb2` and so on are PDB. A trailing
    /// `.gz` is passed over, since the readers decompress gzip-compressed
    /// files when the `gzip` feature is on.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::Format;
    ///
    /// assert_eq!(Format::from_path("1ABC.CIF"), Some(Format::Mmcif));
    /// assert_eq!(Format::from_path("1abc.bcif.gz"), Some(Format::Bcif));
    /// assert_eq!(Format::from_path("1abc.pdb1"), Some(Format::Pdb));
    /// assert_eq!(Format::from_path("1abc.xyz"), None);
    /// ```
    pub fn from_path(path: impl AsRef<Path>) -> Option<Format> {
        let name = path.as_ref().file_name()?.to_str()?.to_ascii_lowercase();
        let name = name.strip_suffix(".gz").unwrap_or(&name);
        let (_, extension) = name.rsplit_once('.')?;
        match extension {
            "cif" | "mmcif" => Some(Format::Mmcif),
            "bcif" => Some(Format::Bcif),
            "pdb" | "ent" => Some(Format::Pdb),
            _ if is_assembly_extension(extension) => Some(Format::Pdb),
            _ => None,
        }
    }
}

fn is_assembly_extension(extension: &str) -> bool {
    extension.strip_prefix("pdb").is_some_and(|number| {
        !number.is_empty() && number.bytes().all(|byte| byte.is_ascii_digit())
    })
}

/// Reads a structure from a file with `options`, choosing the format with
/// [`Format::from_path`].
///
/// An mmCIF or BinaryCIF file is read from its first data block, and the
/// parsed document is kept in [`Loaded::document`].
///
/// # Errors
///
/// Returns [`ReadError::UnknownFormat`] if the extension names no format, and
/// otherwise the errors of [`read`](crate::read), except those of the
/// optional groups that a lenient read skips.
///
/// # Examples
///
/// ```
/// use mmcif_model::{ReadOptions, read_with};
///
/// let loaded = read_with("tests/fixtures/pathologies.cif", &ReadOptions::new())?;
/// let document = loaded.document.unwrap();
///
/// assert_eq!(document.block(0).unwrap().name(), "TEST");
/// assert_eq!(loaded.structure.entry_id(), Some("TEST"));
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
pub fn read_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<Loaded, ReadError> {
    let path = path.as_ref();
    match Format::from_path(path) {
        Some(Format::Mmcif) => from_document_with(Document::open(path)?, options),
        Some(Format::Bcif) => from_document_with(Document::open_bcif(path)?, options),
        Some(Format::Pdb) => pdb::read_pdb_with(path, options),
        None => Err(ReadError::UnknownFormat(path.display().to_string())),
    }
}

fn from_document_with(document: Document, options: &ReadOptions) -> Result<Loaded, ReadError> {
    let block = document.block(0).ok_or(ReadError::NoDataBlock)?;
    let loaded = Structure::from_block_with(&block, options)?;
    Ok(Loaded {
        document: Some(document),
        ..loaded
    })
}

pub(crate) struct Skipping {
    lenient: bool,
    errors: Vec<ReadError>,
}

impl Skipping {
    pub(crate) fn or_default<T: Default>(
        &mut self,
        result: Result<T, ReadError>,
    ) -> Result<T, ReadError> {
        match result {
            Err(error) if self.lenient => {
                self.errors.push(error);
                Ok(T::default())
            }
            result => result,
        }
    }

    pub(crate) fn into_errors(self) -> Vec<ReadError> {
        self.errors
    }
}

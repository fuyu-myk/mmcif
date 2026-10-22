use crate::expression::InvalidExpression;
use crate::ids::{AuthAsymId, LabelAsymId};

/// An error from reading a structure.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ReadError {
    /// The file could not be read, or its contents are not valid CIF or
    /// BinaryCIF.
    #[error(transparent)]
    Cif(#[from] mmcif_cif::Error),
    /// The text is not valid CIF.
    #[error(transparent)]
    Parse(#[from] mmcif_cif::ParseError),
    /// The document has no data block.
    #[error("document has no data block")]
    NoDataBlock,
    /// A required category is absent. Holds the category name without the
    /// leading underscore.
    #[error("missing category _{0}")]
    MissingCategory(&'static str),
    /// A category lacks a required item.
    #[error("missing item _{category}.{item}")]
    MissingItem {
        /// The category name, without the leading underscore.
        category: &'static str,
        /// The item name.
        item: &'static str,
    },
    /// A value is missing where one is required, or cannot be read.
    #[error("_{category}.{item} row {row}: {problem}: {value:?}")]
    InvalidValue {
        /// The category name, without the leading underscore.
        category: &'static str,
        /// The item name.
        item: &'static str,
        /// The zero-based row within the category.
        row: usize,
        /// The value as text.
        value: String,
        /// A short description of the problem.
        problem: &'static str,
    },
    /// The atoms of one label chain have different author chain ids.
    #[error("chain {label_asym} has atoms in author chains {first} and {second}")]
    InconsistentAuthChain {
        /// The label chain id.
        label_asym: LabelAsymId,
        /// The author chain id of the chain's first atom.
        first: AuthAsymId,
        /// The differing author chain id.
        second: AuthAsymId,
    },
    /// An operator expression in `_pdbx_struct_assembly_gen` cannot be
    /// parsed.
    #[error(transparent)]
    OperatorExpression(#[from] InvalidExpression),
    /// A PDB file has a malformed record or is not valid UTF-8.
    #[error("PDB line {line}: {problem}")]
    Pdb {
        /// The one-based line number.
        line: usize,
        /// A short description of the problem.
        problem: &'static str,
    },
    /// The file extension is not one that [`read`](crate::read) recognises.
    /// Holds the path.
    #[error("unrecognised structure file extension in {0}")]
    UnknownFormat(String),
}

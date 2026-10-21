use crate::expression::InvalidExpression;
use crate::ids::{AuthAsymId, LabelAsymId};

#[derive(Debug, thiserror::Error)]
pub enum ReadError {
    #[error(transparent)]
    Cif(#[from] mmcif_cif::Error),
    #[error(transparent)]
    Parse(#[from] mmcif_cif::ParseError),
    #[error("document has no data block")]
    NoDataBlock,
    #[error("missing category _{0}")]
    MissingCategory(&'static str),
    #[error("missing item _{category}.{item}")]
    MissingItem {
        category: &'static str,
        item: &'static str,
    },
    #[error("_{category}.{item} row {row}: {problem}: {value:?}")]
    InvalidValue {
        category: &'static str,
        item: &'static str,
        row: usize,
        value: String,
        problem: &'static str,
    },
    #[error("chain {label_asym} has atoms in author chains {first} and {second}")]
    InconsistentAuthChain {
        label_asym: LabelAsymId,
        first: AuthAsymId,
        second: AuthAsymId,
    },
    #[error(transparent)]
    OperatorExpression(#[from] InvalidExpression),
    #[error("PDB line {line}: {problem}")]
    Pdb { line: usize, problem: &'static str },
    #[error("unrecognised structure file extension in {0}")]
    UnknownFormat(String),
}

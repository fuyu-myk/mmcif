use std::io;

use mmcif_cif::ParseError;

use crate::format::{FORMAT_VERSION, MAX_ID_LEN};

#[derive(Debug, thiserror::Error)]
pub enum CcdError {
    #[error(transparent)]
    Io(#[from] io::Error),
    #[error(transparent)]
    Parse(#[from] ParseError),
    #[error("block {block}: missing required item {item}")]
    MissingItem { block: String, item: String },
    #[error("block {block}: {item} has invalid value {value:?}")]
    InvalidValue {
        block: String,
        item: String,
        value: String,
    },
    #[error("block {block}: unknown bond order {order:?}")]
    UnknownBondOrder { block: String, order: String },
    #[error("block {block}: bond references unknown atom {atom:?}")]
    UnknownBondAtom { block: String, atom: String },
    #[error("block {block}: {count} atoms exceed the supported maximum")]
    TooManyAtoms { block: String, count: usize },
    #[error("component id {id:?} is longer than {max} bytes", max = MAX_ID_LEN)]
    IdTooLong { id: String },
    #[error("component id {id:?} must be non-empty printable ASCII")]
    InvalidId { id: String },
    #[error("duplicate component id {id:?}")]
    DuplicateId { id: String },
    #[error("not a component index: bad magic")]
    BadMagic,
    #[error("unsupported index format version {found}, expected {expected}", expected = FORMAT_VERSION)]
    UnsupportedVersion { found: u32 },
    #[error("corrupt index: {reason}")]
    CorruptIndex { reason: &'static str },
}

pub(crate) fn corrupt(reason: &'static str) -> CcdError {
    CcdError::CorruptIndex { reason }
}

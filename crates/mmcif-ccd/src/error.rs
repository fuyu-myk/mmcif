use std::io;

use mmcif_cif::ParseError;

use crate::format::{FORMAT_VERSION, MAX_ID_LEN};

/// An error from building, reading or querying a [`CcdIndex`](crate::CcdIndex).
///
/// # Examples
///
/// ```
/// use mmcif_ccd::{CcdError, CcdIndex};
///
/// let error = CcdIndex::from_bytes(b"not an index".to_vec()).unwrap_err();
/// assert!(matches!(error, CcdError::BadMagic));
/// ```
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum CcdError {
    /// Reading or writing a file or stream failed, or a source file is not
    /// valid UTF-8.
    #[error(transparent)]
    Io(#[from] io::Error),
    /// The dictionary source is not valid CIF.
    #[error(transparent)]
    Parse(#[from] ParseError),
    /// A required item is missing from a data block, or holds `?` or `.`.
    #[error("block {block}: missing required item {item}")]
    MissingItem {
        /// The name of the data block.
        block: String,
        /// The tag of the missing item, such as `_chem_comp.id`.
        item: String,
    },
    /// An item holds a value that is not valid for that item.
    #[error("block {block}: {item} has invalid value {value:?}")]
    InvalidValue {
        /// The name of the data block.
        block: String,
        /// The tag of the item, such as `_chem_comp_atom.charge`.
        item: String,
        /// The value as written in the block.
        value: String,
    },
    /// A bond has an order that is not a [`BondOrder`](crate::BondOrder) code.
    #[error("block {block}: unknown bond order {order:?}")]
    UnknownBondOrder {
        /// The name of the data block.
        block: String,
        /// The order as written in the block.
        order: String,
    },
    /// A bond names an atom that is not in the component's atom list.
    #[error("block {block}: bond references unknown atom {atom:?}")]
    UnknownBondAtom {
        /// The name of the data block.
        block: String,
        /// The atom name as written in the block.
        atom: String,
    },
    /// A component has more than 65,536 atoms.
    #[error("block {block}: {count} atoms exceed the supported maximum")]
    TooManyAtoms {
        /// The name of the data block.
        block: String,
        /// The number of atoms in the block.
        count: usize,
    },
    /// A component id is longer than 8 bytes.
    #[error("component id {id:?} is longer than {max} bytes", max = MAX_ID_LEN)]
    IdTooLong {
        /// The component id.
        id: String,
    },
    /// A component id is empty or contains a whitespace, control or non-ASCII
    /// character.
    #[error("component id {id:?} must be non-empty printable ASCII")]
    InvalidId {
        /// The component id.
        id: String,
    },
    /// Two components added to one index have the same id.
    #[error("duplicate component id {id:?}")]
    DuplicateId {
        /// The repeated component id.
        id: String,
    },
    /// The data does not begin with the signature of a component index.
    #[error("not a component index: bad magic")]
    BadMagic,
    /// The index was written in a format version this crate does not read.
    #[error("unsupported index format version {found}, expected {expected}", expected = FORMAT_VERSION)]
    UnsupportedVersion {
        /// The format version recorded in the index.
        found: u32,
    },
    /// The index is truncated or malformed.
    #[error("corrupt index: {reason}")]
    CorruptIndex {
        /// A description of the defect found.
        reason: &'static str,
    },
}

pub(crate) fn corrupt(reason: &'static str) -> CcdError {
    CcdError::CorruptIndex { reason }
}

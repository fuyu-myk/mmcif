use crate::builder::BuildError;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum BcifError {
    #[error("MessagePack data ends unexpectedly at byte {offset}")]
    Truncated { offset: usize },
    #[error("unsupported MessagePack marker 0x{marker:02x} at byte {offset}")]
    UnsupportedMarker { marker: u8, offset: usize },
    #[error("MessagePack nesting is deeper than {limit} levels")]
    TooDeep { limit: usize },
    #[error("{count} unexpected bytes after the MessagePack document")]
    TrailingBytes { count: usize },
    #[error("missing field {0:?}")]
    MissingField(&'static str),
    #[error("field {field:?} is not {expected}")]
    UnexpectedType {
        field: &'static str,
        expected: &'static str,
    },
    #[error("unknown encoding kind {0:?}")]
    UnknownEncoding(String),
    #[error("unsupported data type code {0}")]
    UnsupportedDataType(i64),
    #[error("{encoding} cannot decode {found} data")]
    UnexpectedArray {
        encoding: &'static str,
        found: &'static str,
    },
    #[error("invalid {encoding} parameter {parameter:?}")]
    InvalidParameter {
        encoding: &'static str,
        parameter: &'static str,
    },
    #[error("{len} bytes do not split into {element_size}-byte elements")]
    MisalignedBytes { len: usize, element_size: usize },
    #[error("{what} has {actual} elements, expected {expected}")]
    LengthMismatch {
        what: &'static str,
        expected: usize,
        actual: usize,
    },
    #[error("decoding needs more than {limit} array elements")]
    LimitExceeded { limit: usize },
    #[error("decoded text exceeds the limit of {limit} bytes")]
    TextLimitExceeded { limit: usize },
    #[error("integer overflow while decoding {0}")]
    Overflow(&'static str),
    #[error("RunLength data has an odd number of elements ({0})")]
    OddRunLength(usize),
    #[error("IntegerPacking data ends inside a packed value")]
    UnterminatedPacking,
    #[error("string is not valid UTF-8")]
    InvalidUtf8,
    #[error("string offsets {start}..{end} do not fit the string data")]
    InvalidStringOffsets { start: i64, end: i64 },
    #[error("string index {0} is out of range")]
    InvalidStringIndex(i64),
    #[error("mask value {0} is not 0, 1 or 2")]
    InvalidMask(i64),
    #[error("value of {tag} cannot be written as CIF 1.1 text")]
    UnrepresentableValue { tag: String },
    #[error(transparent)]
    Build(#[from] BuildError),
}

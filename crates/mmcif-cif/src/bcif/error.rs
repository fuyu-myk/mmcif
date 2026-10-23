/// An error from decoding BinaryCIF data.
///
/// # Examples
///
/// ```
/// use mmcif_cif::{BcifError, Document};
///
/// assert_eq!(
///     Document::from_bcif(&[0x81]).unwrap_err(),
///     BcifError::Truncated { offset: 1 }
/// );
/// ```
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum BcifError {
    /// The input ends before a MessagePack value is complete.
    #[error("MessagePack data ends unexpectedly at byte {offset}")]
    Truncated {
        /// The byte offset at which more data was needed.
        offset: usize,
    },
    /// The input has a MessagePack type marker the decoder does not read, such
    /// as an extension type.
    #[error("unsupported MessagePack marker 0x{marker:02x} at byte {offset}")]
    UnsupportedMarker {
        /// The marker byte.
        marker: u8,
        /// The byte offset of the marker.
        offset: usize,
    },
    /// MessagePack arrays and maps are nested more deeply than the decoder
    /// allows.
    #[error("MessagePack nesting is deeper than {limit} levels")]
    TooDeep {
        /// The maximum nesting depth.
        limit: usize,
    },
    /// Bytes follow the end of the MessagePack value.
    #[error("{count} unexpected bytes after the MessagePack document")]
    TrailingBytes {
        /// The number of bytes after the value.
        count: usize,
    },
    /// A field the format requires is missing.
    #[error("missing field {0:?}")]
    MissingField(&'static str),
    /// A field or record does not have the MessagePack type the format requires.
    #[error("field {field:?} is not {expected}")]
    UnexpectedType {
        /// The name of the field or record.
        field: &'static str,
        /// A description of the expected type.
        expected: &'static str,
    },
    /// An encoding has a `kind` the decoder does not know.
    #[error("unknown encoding kind {0:?}")]
    UnknownEncoding(String),
    /// A data type code is unknown, or is not a floating-point type where one
    /// is required.
    #[error("unsupported data type code {0}")]
    UnsupportedDataType(i64),
    /// Data reached an encoding or column part that cannot take data of its
    /// kind.
    #[error("{encoding} cannot decode {found} data")]
    UnexpectedArray {
        /// The encoding or column part that received the data.
        encoding: &'static str,
        /// The kind of data it received.
        found: &'static str,
    },
    /// An encoding parameter has an invalid value.
    #[error("invalid {encoding} parameter {parameter:?}")]
    InvalidParameter {
        /// The name of the encoding.
        encoding: &'static str,
        /// The name of the parameter.
        parameter: &'static str,
    },
    /// The length of `ByteArray` data is not a multiple of its element size.
    #[error("{len} bytes do not split into {element_size}-byte elements")]
    MisalignedBytes {
        /// The number of bytes.
        len: usize,
        /// The size of one element in bytes.
        element_size: usize,
    },
    /// Decoded data has a different number of elements than declared.
    #[error("{what} has {actual} elements, expected {expected}")]
    LengthMismatch {
        /// The data whose number of elements differs from the declared number.
        what: &'static str,
        /// The declared number of elements.
        expected: usize,
        /// The decoded number of elements.
        actual: usize,
    },
    /// Decoding needs more elements than
    /// [`BcifLimits::max_elements`](crate::BcifLimits::max_elements) allows.
    #[error("decoding needs more than {limit} array elements")]
    LimitExceeded {
        /// The element limit.
        limit: usize,
    },
    /// The decoded document needs more text than
    /// [`BcifLimits::max_text_bytes`](crate::BcifLimits::max_text_bytes) allows.
    #[error("decoded text exceeds the limit of {limit} bytes")]
    TextLimitExceeded {
        /// The text limit in bytes.
        limit: usize,
    },
    /// An integer overflows while decoding the named data.
    #[error("integer overflow while decoding {0}")]
    Overflow(&'static str),
    /// `RunLength` data has an odd number of elements.
    #[error("RunLength data has an odd number of elements ({0})")]
    OddRunLength(usize),
    /// `IntegerPacking` data ends inside a packed value.
    #[error("IntegerPacking data ends inside a packed value")]
    UnterminatedPacking,
    /// A string is not valid UTF-8.
    #[error("string is not valid UTF-8")]
    InvalidUtf8,
    /// A pair of `StringArray` offsets does not select text within its string
    /// data.
    #[error("string offsets {start}..{end} do not fit the string data")]
    InvalidStringOffsets {
        /// The start offset.
        start: i64,
        /// The end offset.
        end: i64,
    },
    /// A `StringArray` index does not refer to one of its strings.
    #[error("string index {0} is out of range")]
    InvalidStringIndex(i64),
    /// A mask value is not 0, 1 or 2.
    #[error("mask value {0} is not 0, 1 or 2")]
    InvalidMask(i64),
    /// A string contains a line break followed by `;`.
    #[error("value of {tag} cannot be written as CIF 1.1 text")]
    UnrepresentableValue {
        /// The tag of the column holding the string.
        tag: String,
    },
    /// A block header is empty or contains whitespace.
    #[error("invalid block name {0:?}")]
    InvalidBlockName(String),
    /// A category name and a column name do not form a tag without whitespace.
    #[error("invalid tag {0:?}")]
    InvalidTag(String),
    /// A tag appears more than once in a block, ignoring ASCII case.
    #[error("tag {0} appears more than once in the block")]
    DuplicateTag(String),
    /// The text of the decoded document exceeds `u32::MAX` bytes, or a value
    /// is longer than `2^29 - 1` bytes.
    #[error("document exceeds the supported size")]
    TooLarge,
}

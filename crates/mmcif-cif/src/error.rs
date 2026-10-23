use std::io;

use crate::bcif::BcifError;

/// An error from parsing CIF text, with the position where it was found.
///
/// # Examples
///
/// ```
/// use mmcif_cif::{Document, ParseErrorKind};
///
/// let error = Document::parse("data_x\n_entry.id\n").unwrap_err();
///
/// assert_eq!(error.kind, ParseErrorKind::MissingValue);
/// assert_eq!((error.line, error.column, error.offset), (2, 1, 7));
/// assert_eq!(error.to_string(), "line 2, column 1: tag has no value");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("line {line}, column {column}: {kind}")]
pub struct ParseError {
    /// The kind of problem found.
    pub kind: ParseErrorKind,
    /// The byte offset of the error from the start of the input.
    pub offset: usize,
    /// The line of the error, counting from 1.
    pub line: usize,
    /// The byte position of the error within its line, counting from 1.
    pub column: usize,
}

/// The kind of problem a [`ParseError`] reports.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum ParseErrorKind {
    /// A quoted value has no closing quote before the end of its line.
    #[error("unterminated quoted value")]
    UnterminatedQuote,
    /// A text field has no closing `;` at the start of a later line.
    #[error("unterminated text field")]
    UnterminatedTextField,
    /// A data block header is `data_` with no name.
    #[error("data block header has no name")]
    EmptyBlockName,
    /// Something other than whitespace or a comment comes before the first data block header.
    #[error("content before the first data block")]
    OutsideBlock,
    /// A tag is not followed by a value.
    #[error("tag has no value")]
    MissingValue,
    /// A value outside a loop does not follow a tag.
    #[error("value has no tag")]
    UnexpectedValue,
    /// A `loop_` is not followed by any tags.
    #[error("loop_ has no tags")]
    EmptyLoop,
    /// A tag appears more than once in a block or save frame, ignoring ASCII case.
    #[error("tag {0} appears more than once")]
    DuplicateTag(String),
    /// The number of values in a loop is not a multiple of its number of tags.
    #[error("loop has {values} values for {tags} tags")]
    RaggedLoop {
        /// The number of values in the loop.
        values: usize,
        /// The number of tags in the loop.
        tags: usize,
    },
    /// A save frame starts inside another save frame.
    #[error("save frame inside a save frame")]
    NestedSaveFrame,
    /// A `save_` terminator appears outside a save frame.
    #[error("save frame terminator outside a save frame")]
    UnmatchedSaveEnd,
    /// A save frame is not closed by `save_` before its block ends.
    #[error("save frame is not terminated")]
    UnterminatedSaveFrame,
    /// The reserved word `global_` or `stop_` appears.
    #[error("{0} is not allowed in CIF")]
    Reserved(&'static str),
    /// The input is not valid UTF-8.
    #[error("input is not valid UTF-8")]
    InvalidUtf8,
    /// The input is longer than `u32::MAX` bytes, or a value is longer than `2^29 - 1` bytes.
    #[error("input or value exceeds the supported size")]
    TooLarge,
}

impl ParseError {
    pub(crate) fn at(src: &[u8], offset: usize, kind: ParseErrorKind) -> Self {
        let offset = offset.min(src.len());
        let before = &src[..offset];
        let line = memchr::memchr_iter(b'\n', before).count() + 1;
        let line_start = memchr::memrchr(b'\n', before).map_or(0, |i| i + 1);
        ParseError {
            kind,
            offset,
            line,
            column: offset - line_start + 1,
        }
    }
}

/// An error from reading a [`Document`](crate::Document) from a file.
///
/// # Examples
///
/// ```
/// use mmcif_cif::{Document, Error};
///
/// let error = Document::open("missing.cif").unwrap_err();
/// assert!(matches!(error, Error::Io(_)));
/// ```
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// The file could not be opened or read.
    #[error(transparent)]
    Io(#[from] io::Error),
    /// The file is not valid CIF text.
    #[error(transparent)]
    Parse(#[from] ParseError),
    /// The file is not valid BinaryCIF.
    #[error(transparent)]
    Bcif(#[from] BcifError),
    /// The file is gzip-compressed, and this build cannot decompress it
    /// because the `gzip` feature is off.
    #[error("the file is gzip-compressed and the gzip feature is off")]
    Compressed,
}

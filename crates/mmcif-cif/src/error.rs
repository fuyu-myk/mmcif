use std::io;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("line {line}, column {column}: {kind}")]
pub struct ParseError {
    pub kind: ParseErrorKind,
    pub offset: usize,
    pub line: usize,
    pub column: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ParseErrorKind {
    #[error("unterminated quoted value")]
    UnterminatedQuote,
    #[error("unterminated text field")]
    UnterminatedTextField,
    #[error("data block header has no name")]
    EmptyBlockName,
    #[error("content before the first data block")]
    OutsideBlock,
    #[error("tag has no value")]
    MissingValue,
    #[error("value has no tag")]
    UnexpectedValue,
    #[error("loop_ has no tags")]
    EmptyLoop,
    #[error("loop has {values} values for {tags} tags")]
    RaggedLoop { values: usize, tags: usize },
    #[error("save frame inside a save frame")]
    NestedSaveFrame,
    #[error("save frame terminator outside a save frame")]
    UnmatchedSaveEnd,
    #[error("save frame is not terminated")]
    UnterminatedSaveFrame,
    #[error("{0} is not allowed in CIF")]
    Reserved(&'static str),
    #[error("input is not valid UTF-8")]
    InvalidUtf8,
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

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error(transparent)]
    Io(#[from] io::Error),
    #[error(transparent)]
    Parse(#[from] ParseError),
}

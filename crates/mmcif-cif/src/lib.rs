mod bcif;
mod builder;
mod document;
mod error;
mod handle;
mod parse;
mod source;
mod token;
mod value;
mod write;

pub use bcif::{BcifError, BcifLimits};
pub use builder::{BuildError, DocumentBuilder};
pub use document::{Block, Category, Column, Document, Entry, Frame, Loop, Pair, StreamedBlock};
pub use error::{Error, ParseError, ParseErrorKind};
pub use handle::{BlockHandle, CategoryHandle};
pub use parse::BlockStream;
pub use value::{InvalidValue, Value, ValueKind};

pub fn stream_blocks(text: &str) -> BlockStream<'_> {
    BlockStream::new(text)
}

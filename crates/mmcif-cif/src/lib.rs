//! Lossless reading and writing of CIF and mmCIF text, with a BinaryCIF decoder.
//!
//! A [`Document`] holds data blocks whose categories read as tables of
//! [`Value`]s. Each value keeps the quoting it was read with, and a document
//! writes back to CIF text.
//!
//! # Examples
//!
//! ```
//! use mmcif_cif::Document;
//!
//! let text = "\
//! data_1ABC
//! _entry.id 1ABC
//! loop_
//! _atom_site.id
//! _atom_site.label_atom_id
//! 1 N
//! 2 CA
//! ";
//! let document = Document::parse(text)?;
//! let block = document.block_named("1ABC").unwrap();
//!
//! let atoms = block.category("atom_site").unwrap();
//! assert_eq!(atoms.row_count(), 2);
//! assert_eq!(atoms.column("label_atom_id").unwrap().get(1).unwrap().as_str(), Some("CA"));
//!
//! let written = document.to_cif_string();
//! assert_eq!(Document::parse(written)?, document);
//! # Ok::<(), mmcif_cif::ParseError>(())
//! ```

#![warn(missing_docs)]

mod bcif;
mod builder;
mod columns;
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
pub use columns::TypedValue;
pub use document::{Block, Category, Column, Document, Entry, Frame, Loop, Pair, StreamedBlock};
pub use error::{Error, ParseError, ParseErrorKind};
pub use handle::{BlockHandle, CategoryHandle};
pub use parse::BlockStream;
pub use value::{InvalidValue, Value, ValueKind};

/// Returns an iterator that parses the data blocks of `text` one at a time.
///
/// Each item is a [`StreamedBlock`] or the [`ParseError`] that stopped
/// parsing, and no items follow an error. Errors are those of
/// [`Document::parse`].
///
/// # Examples
///
/// ```
/// use mmcif_cif::stream_blocks;
///
/// let text = "data_a\n_x.y 1\ndata_b\n_x.y 2\n";
///
/// let mut names = Vec::new();
/// for streamed in stream_blocks(text) {
///     names.push(streamed?.block().name().to_owned());
/// }
/// assert_eq!(names, ["a", "b"]);
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
pub fn stream_blocks(text: &str) -> BlockStream<'_> {
    BlockStream::new(text)
}

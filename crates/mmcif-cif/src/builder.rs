use std::collections::HashSet;

use crate::columns::{RawColumn, RawTable};
use crate::document::{Cell, Document, RawBlock, RawEntry, RawLoop, Span};
use crate::value::{InvalidValue, Value};

/// An error from building a document with a [`DocumentBuilder`].
///
/// # Examples
///
/// ```
/// use mmcif_cif::{BuildError, DocumentBuilder, Value};
///
/// let mut builder = DocumentBuilder::new();
/// assert_eq!(builder.pair("_entry.id", Value::text("1ABC")), Err(BuildError::NoBlock));
/// ```
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum BuildError {
    /// An entry was added before any block was started.
    #[error("no data block has been started")]
    NoBlock,
    /// A value was pushed while no loop was open.
    #[error("no loop has been started")]
    NoLoop,
    /// A block name is empty or contains whitespace.
    #[error("invalid block name {0:?}")]
    InvalidBlockName(String),
    /// A tag does not start with `_` or contains whitespace.
    #[error("invalid tag {0:?}")]
    InvalidTag(String),
    /// The kind of a value cannot hold its text.
    #[error(transparent)]
    InvalidValue(#[from] InvalidValue),
    /// A loop was started with no tags.
    #[error("loop_ has no tags")]
    EmptyLoop,
    /// A tag appears more than once in a block, ignoring ASCII case.
    #[error("tag {0} appears more than once in the block")]
    DuplicateTag(String),
    /// A loop was closed with no values.
    #[error("loop_ has no rows")]
    LoopWithoutRows,
    /// The number of values in a loop is not a multiple of its number of tags.
    #[error("loop has {values} values for {tags} tags")]
    RaggedLoop {
        /// The number of values in the loop.
        values: usize,
        /// The number of tags in the loop.
        tags: usize,
    },
    /// The text of the document exceeds `u32::MAX` bytes, or a value is
    /// longer than `2^29 - 1` bytes.
    #[error("document exceeds the supported size")]
    TooLarge,
}

/// Builds a [`Document`] from block names, tags and values.
///
/// Entries are added to the most recently started block. A loop takes values
/// from [`push`](DocumentBuilder::push) until another block, pair or loop is
/// started, or the builder is finished.
///
/// # Examples
///
/// ```
/// use mmcif_cif::{DocumentBuilder, Value};
///
/// let mut builder = DocumentBuilder::new();
/// builder.start_block("1ABC")?;
/// builder.pair("_entry.id", Value::text("1ABC"))?;
/// builder.start_loop(["_atom.id", "_atom.name"])?;
/// for text in ["1", "N", "2", "CA"] {
///     builder.push(Value::text(text))?;
/// }
/// let document = builder.finish()?;
///
/// assert_eq!(
///     document.to_cif_string(),
///     "data_1ABC\n#\n_entry.id 1ABC\n#\nloop_\n_atom.id\n_atom.name\n1 N\n2 CA\n#\n"
/// );
/// # Ok::<(), mmcif_cif::BuildError>(())
/// ```
#[derive(Default)]
pub struct DocumentBuilder {
    text: String,
    blocks: Vec<RawBlock>,
    open_loop: Option<RawLoop>,
    block_tags: HashSet<String>,
}

fn is_name_token(name: &str) -> bool {
    !name.is_empty() && !name.bytes().any(|b| b.is_ascii_whitespace())
}

impl DocumentBuilder {
    /// Creates a builder with no blocks.
    pub fn new() -> Self {
        DocumentBuilder::default()
    }

    pub(crate) fn text_len(&self) -> usize {
        self.text.len()
    }

    fn store(&mut self, text: &str) -> Result<(usize, usize), BuildError> {
        let start = self.text.len();
        self.text.push_str(text);
        if u32::try_from(self.text.len()).is_err() {
            return Err(BuildError::TooLarge);
        }
        Ok((start, self.text.len()))
    }

    fn store_span(&mut self, text: &str) -> Result<Span, BuildError> {
        let (start, end) = self.store(text)?;
        Ok(Span {
            start: start as u32,
            len: (end - start) as u32,
        })
    }

    pub(crate) fn text_cell(&mut self, value: Value<'_>) -> Result<Cell, BuildError> {
        self.store_value(value)
    }

    pub(crate) fn push_table(
        &mut self,
        tags: &[&str],
        columns: Vec<RawColumn>,
    ) -> Result<(), BuildError> {
        self.close_loop()?;
        if self.blocks.is_empty() {
            return Err(BuildError::NoBlock);
        }
        if tags.is_empty() {
            return Err(BuildError::EmptyLoop);
        }
        let rows = columns.first().map_or(0, RawColumn::len);
        if rows == 0 {
            return Err(BuildError::LoopWithoutRows);
        }
        if columns.len() != tags.len() || columns.iter().any(|column| column.len() != rows) {
            return Err(BuildError::RaggedLoop {
                values: columns.iter().map(RawColumn::len).sum(),
                tags: tags.len(),
            });
        }
        let tags = tags
            .iter()
            .map(|tag| self.tag_span(tag))
            .collect::<Result<Vec<_>, _>>()?;
        let block = self.blocks.last_mut().ok_or(BuildError::NoBlock)?;
        block.entries.push(RawEntry::Table(RawTable {
            tags,
            columns,
            rows,
        }));
        Ok(())
    }

    fn store_value(&mut self, value: Value<'_>) -> Result<Cell, BuildError> {
        if !value.kind().can_hold(value.raw()) {
            return Err(InvalidValue {
                kind: value.kind(),
                text: value.raw().to_owned(),
            }
            .into());
        }
        let (start, end) = self.store(value.raw())?;
        Cell::new(start, end, value.kind()).ok_or(BuildError::TooLarge)
    }

    fn close_loop(&mut self) -> Result<(), BuildError> {
        let Some(raw_loop) = self.open_loop.take() else {
            return Ok(());
        };
        if raw_loop.cells.is_empty() {
            return Err(BuildError::LoopWithoutRows);
        }
        if raw_loop.cells.len() % raw_loop.tags.len() != 0 {
            return Err(BuildError::RaggedLoop {
                values: raw_loop.cells.len(),
                tags: raw_loop.tags.len(),
            });
        }
        let block = self.blocks.last_mut().ok_or(BuildError::NoBlock)?;
        block.entries.push(RawEntry::Loop(raw_loop));
        Ok(())
    }

    /// Starts a data block called `name`, closing any open loop.
    ///
    /// `name` is given without the `data_` prefix.
    ///
    /// # Errors
    ///
    /// Returns [`BuildError::InvalidBlockName`] if `name` is empty or contains
    /// whitespace, and [`BuildError::TooLarge`] if the document text would
    /// exceed `u32::MAX` bytes. Closing an open loop returns
    /// [`BuildError::LoopWithoutRows`] if it has no values and
    /// [`BuildError::RaggedLoop`] if its number of values is not a multiple of
    /// its number of tags.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{BuildError, DocumentBuilder};
    ///
    /// let mut builder = DocumentBuilder::new();
    /// builder.start_block("1ABC")?;
    ///
    /// assert_eq!(
    ///     builder.start_block("two words"),
    ///     Err(BuildError::InvalidBlockName("two words".into()))
    /// );
    /// # Ok::<(), BuildError>(())
    /// ```
    pub fn start_block(&mut self, name: &str) -> Result<(), BuildError> {
        self.close_loop()?;
        if !is_name_token(name) {
            return Err(BuildError::InvalidBlockName(name.to_owned()));
        }
        let name = self.store_span(name)?;
        self.block_tags.clear();
        self.blocks.push(RawBlock {
            name,
            entries: Vec::new(),
        });
        Ok(())
    }

    fn tag_span(&mut self, tag: &str) -> Result<Span, BuildError> {
        if !tag.starts_with('_') || !is_name_token(tag) {
            return Err(BuildError::InvalidTag(tag.to_owned()));
        }
        if !self.block_tags.insert(tag.to_ascii_lowercase()) {
            return Err(BuildError::DuplicateTag(tag.to_owned()));
        }
        self.store_span(tag)
    }

    /// Adds `tag` with `value` to the current block, closing any open loop.
    ///
    /// # Errors
    ///
    /// Returns [`BuildError::NoBlock`] if no block has been started,
    /// [`BuildError::InvalidTag`] if `tag` does not start with `_` or contains
    /// whitespace, [`BuildError::DuplicateTag`] if the block already has `tag`
    /// in any ASCII case, [`BuildError::InvalidValue`] if the kind of `value`
    /// cannot hold its text, and [`BuildError::TooLarge`] if the document text
    /// would exceed `u32::MAX` bytes. Closing an open loop returns
    /// [`BuildError::LoopWithoutRows`] if it has no values and
    /// [`BuildError::RaggedLoop`] if its number of values is not a multiple of
    /// its number of tags.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{BuildError, DocumentBuilder, Value};
    ///
    /// let mut builder = DocumentBuilder::new();
    /// builder.start_block("1ABC")?;
    /// builder.pair("_struct.title", Value::text("two words"))?;
    ///
    /// assert_eq!(
    ///     builder.pair("_STRUCT.TITLE", Value::UNKNOWN),
    ///     Err(BuildError::DuplicateTag("_STRUCT.TITLE".into()))
    /// );
    /// assert_eq!(
    ///     builder.finish()?.to_cif_string(),
    ///     "data_1ABC\n#\n_struct.title 'two words'\n#\n"
    /// );
    /// # Ok::<(), BuildError>(())
    /// ```
    pub fn pair(&mut self, tag: &str, value: Value<'_>) -> Result<(), BuildError> {
        self.close_loop()?;
        if self.blocks.is_empty() {
            return Err(BuildError::NoBlock);
        }
        let tag = self.tag_span(tag)?;
        let value = self.store_value(value)?;
        let block = self.blocks.last_mut().ok_or(BuildError::NoBlock)?;
        block.entries.push(RawEntry::Pair { tag, value });
        Ok(())
    }

    /// Starts a loop with `tags` in the current block, closing any open loop.
    ///
    /// The loop then takes values from [`push`](DocumentBuilder::push), row by
    /// row.
    ///
    /// # Errors
    ///
    /// Returns [`BuildError::NoBlock`] if no block has been started,
    /// [`BuildError::EmptyLoop`] if `tags` is empty, [`BuildError::InvalidTag`]
    /// if a tag does not start with `_` or contains whitespace,
    /// [`BuildError::DuplicateTag`] if the block already has a tag in any ASCII
    /// case, and [`BuildError::TooLarge`] if the document text would exceed
    /// `u32::MAX` bytes. Closing an open loop returns
    /// [`BuildError::LoopWithoutRows`] if it has no values and
    /// [`BuildError::RaggedLoop`] if its number of values is not a multiple of
    /// its number of tags.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{BuildError, DocumentBuilder};
    ///
    /// let mut builder = DocumentBuilder::new();
    /// builder.start_block("1ABC")?;
    /// builder.start_loop(["_atom.id"])?;
    ///
    /// assert_eq!(builder.start_loop(["_atom.name"]), Err(BuildError::LoopWithoutRows));
    /// # Ok::<(), BuildError>(())
    /// ```
    pub fn start_loop<'t>(
        &mut self,
        tags: impl IntoIterator<Item = &'t str>,
    ) -> Result<(), BuildError> {
        self.close_loop()?;
        if self.blocks.is_empty() {
            return Err(BuildError::NoBlock);
        }
        let tags = tags
            .into_iter()
            .map(|tag| self.tag_span(tag))
            .collect::<Result<Vec<_>, _>>()?;
        if tags.is_empty() {
            return Err(BuildError::EmptyLoop);
        }
        self.open_loop = Some(RawLoop {
            tags,
            cells: Vec::new(),
        });
        Ok(())
    }

    /// Adds `value` to the open loop.
    ///
    /// Values fill the loop row by row, one for each tag.
    ///
    /// # Errors
    ///
    /// Returns [`BuildError::NoLoop`] if no loop is open,
    /// [`BuildError::InvalidValue`] if the kind of `value` cannot hold its
    /// text, and [`BuildError::TooLarge`] if the document text would exceed
    /// `u32::MAX` bytes.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{BuildError, DocumentBuilder, Value};
    ///
    /// let mut builder = DocumentBuilder::new();
    /// builder.start_block("1ABC")?;
    /// assert_eq!(builder.push(Value::UNKNOWN), Err(BuildError::NoLoop));
    ///
    /// builder.start_loop(["_atom.id", "_atom.alt_id"])?;
    /// builder.push(Value::text("1"))?;
    /// builder.push(Value::INAPPLICABLE)?;
    ///
    /// assert_eq!(
    ///     builder.finish()?.to_cif_string(),
    ///     "data_1ABC\n#\nloop_\n_atom.id\n_atom.alt_id\n1 .\n#\n"
    /// );
    /// # Ok::<(), BuildError>(())
    /// ```
    pub fn push(&mut self, value: Value<'_>) -> Result<(), BuildError> {
        if self.open_loop.is_none() {
            return Err(BuildError::NoLoop);
        }
        let cell = self.store_value(value)?;
        self.open_loop
            .as_mut()
            .ok_or(BuildError::NoLoop)?
            .cells
            .push(cell);
        Ok(())
    }

    /// Closes any open loop and returns the built document.
    ///
    /// # Errors
    ///
    /// Returns [`BuildError::LoopWithoutRows`] if the open loop has no values,
    /// and [`BuildError::RaggedLoop`] if its number of values is not a
    /// multiple of its number of tags.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{BuildError, DocumentBuilder, Value};
    ///
    /// let mut builder = DocumentBuilder::new();
    /// builder.start_block("1ABC")?;
    /// builder.start_loop(["_atom.id", "_atom.name"])?;
    /// builder.push(Value::text("1"))?;
    ///
    /// assert_eq!(builder.finish().unwrap_err(), BuildError::RaggedLoop { values: 1, tags: 2 });
    /// # Ok::<(), BuildError>(())
    /// ```
    pub fn finish(mut self) -> Result<Document, BuildError> {
        self.close_loop()?;
        Ok(Document::from_parts(self.text, self.blocks))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::value::ValueKind;

    #[test]
    fn builds_pairs_and_loops() {
        let mut builder = DocumentBuilder::new();
        builder.start_block("test").unwrap();
        builder.pair("_entry.id", Value::text("1ABC")).unwrap();
        builder.start_loop(["_atom.id", "_atom.name"]).unwrap();
        for value in ["1", "N", "2", "two words"] {
            builder.push(Value::text(value)).unwrap();
        }
        let document = builder.finish().unwrap();
        let block = document.block(0).unwrap();
        let atoms = block.category("atom").unwrap();
        assert_eq!(atoms.row_count(), 2);
        let names = atoms.column("name").unwrap();
        assert_eq!(names.get(1).unwrap().as_str(), Some("two words"));
        assert_eq!(names.get(1).unwrap().kind(), ValueKind::SingleQuoted);
    }

    #[test]
    fn rejects_ragged_loops() {
        let mut builder = DocumentBuilder::new();
        builder.start_block("test").unwrap();
        builder.start_loop(["_a.x", "_a.y"]).unwrap();
        builder.push(Value::text("1")).unwrap();
        assert_eq!(
            builder.finish().unwrap_err(),
            BuildError::RaggedLoop { values: 1, tags: 2 }
        );
    }

    #[test]
    fn rejects_duplicate_tags_within_a_block() {
        let mut builder = DocumentBuilder::new();
        builder.start_block("a").unwrap();
        builder.pair("_x.y", Value::UNKNOWN).unwrap();
        assert_eq!(
            builder.start_loop(["_X.Y"]),
            Err(BuildError::DuplicateTag("_X.Y".into()))
        );
        builder.start_block("b").unwrap();
        builder.pair("_x.y", Value::UNKNOWN).unwrap();
    }

    #[test]
    fn rejects_loops_without_rows() {
        let mut builder = DocumentBuilder::new();
        builder.start_block("test").unwrap();
        builder.start_loop(["_a.x"]).unwrap();
        assert_eq!(
            builder.pair("_b.y", Value::UNKNOWN),
            Err(BuildError::LoopWithoutRows)
        );
    }

    #[test]
    fn rejects_unrepresentable_values() {
        let mut builder = DocumentBuilder::new();
        builder.start_block("test").unwrap();
        assert!(matches!(
            builder.pair("_a.x", Value::text("a\n;b")),
            Err(BuildError::InvalidValue(_))
        ));
    }

    #[test]
    fn requires_a_block_and_valid_tags() {
        let mut builder = DocumentBuilder::new();
        assert_eq!(
            builder.pair("_a.x", Value::UNKNOWN),
            Err(BuildError::NoBlock)
        );
        builder.start_block("b").unwrap();
        assert_eq!(
            builder.pair("a.x", Value::UNKNOWN),
            Err(BuildError::InvalidTag("a.x".into()))
        );
        assert_eq!(builder.push(Value::UNKNOWN), Err(BuildError::NoLoop));
    }
}

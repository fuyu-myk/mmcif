use std::collections::HashSet;

use crate::columns::{RawColumn, RawTable};
use crate::document::{Cell, Document, RawBlock, RawEntry, RawLoop, Span};
use crate::value::{InvalidValue, Value};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum BuildError {
    #[error("no data block has been started")]
    NoBlock,
    #[error("no loop has been started")]
    NoLoop,
    #[error("invalid block name {0:?}")]
    InvalidBlockName(String),
    #[error("invalid tag {0:?}")]
    InvalidTag(String),
    #[error(transparent)]
    InvalidValue(#[from] InvalidValue),
    #[error("loop_ has no tags")]
    EmptyLoop,
    #[error("tag {0} appears more than once in the block")]
    DuplicateTag(String),
    #[error("loop_ has no rows")]
    LoopWithoutRows,
    #[error("loop has {values} values for {tags} tags")]
    RaggedLoop { values: usize, tags: usize },
    #[error("document exceeds the supported size")]
    TooLarge,
}

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

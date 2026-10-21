use std::fmt;
use std::path::Path;

use crate::error::{Error, ParseError};
use crate::parse::Parser;
use crate::source::{MappedText, Source, read_text};
use crate::value::{Value, ValueKind};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Span {
    pub start: u32,
    pub len: u32,
}

impl Span {
    pub fn resolve(self, src: &str) -> &str {
        let start = self.start as usize;
        &src[start..start + self.len as usize]
    }
}

const KIND_SHIFT: u32 = 29;
const LEN_MASK: u32 = (1 << KIND_SHIFT) - 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Cell {
    start: u32,
    len_kind: u32,
}

impl Cell {
    pub fn new(start: usize, end: usize, kind: ValueKind) -> Option<Cell> {
        let len = u32::try_from(end - start)
            .ok()
            .filter(|&len| len <= LEN_MASK)?;
        Some(Cell {
            start: u32::try_from(start).ok()?,
            len_kind: len | (kind_code(kind) << KIND_SHIFT),
        })
    }

    pub fn resolve(self, src: &str) -> Value<'_> {
        let start = self.start as usize;
        let len = (self.len_kind & LEN_MASK) as usize;
        Value::from_parts(
            kind_from_code(self.len_kind >> KIND_SHIFT),
            &src[start..start + len],
        )
    }
}

fn kind_code(kind: ValueKind) -> u32 {
    match kind {
        ValueKind::Bare => 0,
        ValueKind::SingleQuoted => 1,
        ValueKind::DoubleQuoted => 2,
        ValueKind::TextField => 3,
        ValueKind::Inapplicable => 4,
        ValueKind::Unknown => 5,
    }
}

fn kind_from_code(code: u32) -> ValueKind {
    match code {
        0 => ValueKind::Bare,
        1 => ValueKind::SingleQuoted,
        2 => ValueKind::DoubleQuoted,
        3 => ValueKind::TextField,
        4 => ValueKind::Inapplicable,
        _ => ValueKind::Unknown,
    }
}

#[derive(Debug, Clone)]
pub(crate) struct RawBlock {
    pub name: Span,
    pub entries: Vec<RawEntry>,
}

#[derive(Debug, Clone)]
pub(crate) enum RawEntry {
    Pair { tag: Span, value: Cell },
    Loop(RawLoop),
    Frame(RawFrame),
}

#[derive(Debug, Clone)]
pub(crate) struct RawLoop {
    pub tags: Vec<Span>,
    pub cells: Vec<Cell>,
}

#[derive(Debug, Clone)]
pub(crate) struct RawFrame {
    pub name: Span,
    pub entries: Vec<RawEntry>,
}

pub struct Document {
    source: Source,
    blocks: Vec<RawBlock>,
}

impl Document {
    pub fn parse(text: impl Into<String>) -> Result<Document, ParseError> {
        Document::from_source(Source::Owned(text.into()))
    }

    pub fn from_bytes(bytes: Vec<u8>) -> Result<Document, ParseError> {
        Document::from_source(Source::Owned(read_text(bytes)?))
    }

    pub fn open(path: impl AsRef<Path>) -> Result<Document, Error> {
        let bytes = std::fs::read(path)?;
        Ok(Document::from_bytes(bytes)?)
    }

    /// # Safety
    /// The file must not be modified or truncated while the document is alive.
    #[allow(unsafe_code)]
    pub unsafe fn open_mapped(path: impl AsRef<Path>) -> Result<Document, Error> {
        // SAFETY: forwarded from this function's contract.
        let source = unsafe { MappedText::open(path.as_ref())? };
        Ok(Document::from_source(source)?)
    }

    fn from_source(source: Source) -> Result<Document, ParseError> {
        let blocks = Parser::new(source.text())?.collect::<Result<Vec<_>, _>>()?;
        Ok(Document { source, blocks })
    }

    pub(crate) fn from_parts(text: String, blocks: Vec<RawBlock>) -> Document {
        Document {
            source: Source::Owned(text),
            blocks,
        }
    }

    pub(crate) fn text(&self) -> &str {
        self.source.text()
    }

    pub fn blocks(&self) -> impl ExactSizeIterator<Item = Block<'_>> {
        let src = self.text();
        self.blocks.iter().map(move |raw| Block { src, raw })
    }

    pub fn block(&self, index: usize) -> Option<Block<'_>> {
        self.blocks.get(index).map(|raw| Block {
            src: self.text(),
            raw,
        })
    }

    pub fn block_named(&self, name: &str) -> Option<Block<'_>> {
        self.blocks()
            .find(|block| block.name().eq_ignore_ascii_case(name))
    }

    pub fn len(&self) -> usize {
        self.blocks.len()
    }

    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }
}

impl PartialEq for Document {
    fn eq(&self, other: &Document) -> bool {
        self.blocks().eq(other.blocks())
    }
}

impl fmt::Debug for Document {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.blocks()).finish()
    }
}

pub struct StreamedBlock<'a> {
    src: &'a str,
    raw: RawBlock,
}

impl<'a> StreamedBlock<'a> {
    pub(crate) fn new(src: &'a str, raw: RawBlock) -> Self {
        StreamedBlock { src, raw }
    }

    pub fn block(&self) -> Block<'_> {
        Block {
            src: self.src,
            raw: &self.raw,
        }
    }
}

#[derive(Clone, Copy)]
pub struct Block<'a> {
    src: &'a str,
    raw: &'a RawBlock,
}

impl<'a> Block<'a> {
    pub fn name(&self) -> &'a str {
        self.raw.name.resolve(self.src)
    }

    pub fn entries(&self) -> impl ExactSizeIterator<Item = Entry<'a>> + use<'a> {
        let src = self.src;
        self.raw.entries.iter().map(move |raw| Entry::new(src, raw))
    }

    pub fn category(&self, name: &str) -> Option<Category<'a>> {
        let name = name.strip_prefix('_').unwrap_or(name);
        let src = self.src;
        let matching =
            |tag: Span| category_of(tag.resolve(src)).filter(|c| c.eq_ignore_ascii_case(name));
        for raw in &self.raw.entries {
            if let RawEntry::Loop(lp) = raw
                && let Some(found) = lp.tags.iter().find_map(|&tag| matching(tag))
            {
                return Some(Category {
                    src,
                    name: found,
                    layout: Layout::Loop(lp),
                });
            }
        }
        let mut found = None;
        let mut pairs = Vec::new();
        for raw in &self.raw.entries {
            if let RawEntry::Pair { tag, value } = raw
                && let Some(category) = matching(*tag)
            {
                found.get_or_insert(category);
                pairs.push((*tag, *value));
            }
        }
        Some(Category {
            src,
            name: found?,
            layout: Layout::Pairs(pairs),
        })
    }

    pub fn category_names(&self) -> Vec<&'a str> {
        let mut names: Vec<&'a str> = Vec::new();
        let tags = self.entries().flat_map(|entry| match entry {
            Entry::Pair(pair) => vec![pair.tag()],
            Entry::Loop(lp) => lp.tags().collect(),
            Entry::Frame(_) => Vec::new(),
        });
        for tag in tags {
            if let Some(name) = category_of(tag)
                && !names.iter().any(|n| n.eq_ignore_ascii_case(name))
            {
                names.push(name);
            }
        }
        names
    }
}

impl PartialEq for Block<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.name() == other.name() && self.entries().eq(other.entries())
    }
}

impl fmt::Debug for Block<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Block")
            .field("name", &self.name())
            .field("entries", &self.entries().collect::<Vec<_>>())
            .finish()
    }
}

pub(crate) fn category_of(tag: &str) -> Option<&str> {
    let body = tag.strip_prefix('_')?;
    body.split_once('.').map(|(category, _)| category)
}

fn item_of(tag: &str) -> &str {
    tag.split_once('.').map_or(tag, |(_, item)| item)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Entry<'a> {
    Pair(Pair<'a>),
    Loop(Loop<'a>),
    Frame(Frame<'a>),
}

impl<'a> Entry<'a> {
    fn new(src: &'a str, raw: &'a RawEntry) -> Self {
        match raw {
            RawEntry::Pair { tag, value } => Entry::Pair(Pair {
                src,
                tag: *tag,
                value: *value,
            }),
            RawEntry::Loop(raw) => Entry::Loop(Loop { src, raw }),
            RawEntry::Frame(raw) => Entry::Frame(Frame { src, raw }),
        }
    }
}

#[derive(Clone, Copy)]
pub struct Pair<'a> {
    src: &'a str,
    tag: Span,
    value: Cell,
}

impl<'a> Pair<'a> {
    pub fn tag(&self) -> &'a str {
        self.tag.resolve(self.src)
    }

    pub fn value(&self) -> Value<'a> {
        self.value.resolve(self.src)
    }
}

impl PartialEq for Pair<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.tag() == other.tag() && self.value() == other.value()
    }
}

impl fmt::Debug for Pair<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("Pair")
            .field(&self.tag())
            .field(&self.value())
            .finish()
    }
}

#[derive(Clone, Copy)]
pub struct Loop<'a> {
    src: &'a str,
    raw: &'a RawLoop,
}

impl<'a> Loop<'a> {
    pub fn tags(&self) -> impl ExactSizeIterator<Item = &'a str> + use<'a> {
        let src = self.src;
        self.raw.tags.iter().map(move |tag| tag.resolve(src))
    }

    pub fn row_count(&self) -> usize {
        self.raw.cells.len() / self.raw.tags.len()
    }

    pub fn values(&self) -> impl ExactSizeIterator<Item = Value<'a>> + use<'a> {
        let src = self.src;
        self.raw.cells.iter().map(move |cell| cell.resolve(src))
    }

    pub fn value(&self, row: usize, column: usize) -> Option<Value<'a>> {
        if column >= self.raw.tags.len() {
            return None;
        }
        let index = row.checked_mul(self.raw.tags.len())?.checked_add(column)?;
        self.raw.cells.get(index).map(|cell| cell.resolve(self.src))
    }
}

impl PartialEq for Loop<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.tags().eq(other.tags()) && self.values().eq(other.values())
    }
}

impl fmt::Debug for Loop<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Loop")
            .field("tags", &self.tags().collect::<Vec<_>>())
            .field("rows", &self.row_count())
            .finish()
    }
}

#[derive(Clone, Copy)]
pub struct Frame<'a> {
    src: &'a str,
    raw: &'a RawFrame,
}

impl<'a> Frame<'a> {
    pub fn name(&self) -> &'a str {
        self.raw.name.resolve(self.src)
    }

    pub fn entries(&self) -> impl ExactSizeIterator<Item = Entry<'a>> + use<'a> {
        let src = self.src;
        self.raw.entries.iter().map(move |raw| Entry::new(src, raw))
    }
}

impl PartialEq for Frame<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.name() == other.name() && self.entries().eq(other.entries())
    }
}

impl fmt::Debug for Frame<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Frame")
            .field("name", &self.name())
            .field("entries", &self.entries().collect::<Vec<_>>())
            .finish()
    }
}

#[derive(Debug, Clone)]
enum Layout<'a> {
    Loop(&'a RawLoop),
    Pairs(Vec<(Span, Cell)>),
}

#[derive(Debug, Clone)]
pub struct Category<'a> {
    src: &'a str,
    name: &'a str,
    layout: Layout<'a>,
}

impl<'a> Category<'a> {
    pub fn name(&self) -> &'a str {
        self.name
    }

    pub fn row_count(&self) -> usize {
        match &self.layout {
            Layout::Loop(lp) => lp.cells.len() / lp.tags.len(),
            Layout::Pairs(_) => 1,
        }
    }

    fn owns(&self, tag: &str) -> bool {
        category_of(tag).is_some_and(|category| category.eq_ignore_ascii_case(self.name))
    }

    pub fn items(&self) -> Vec<&'a str> {
        let src = self.src;
        match &self.layout {
            Layout::Loop(lp) => lp
                .tags
                .iter()
                .map(|t| t.resolve(src))
                .filter(|tag| self.owns(tag))
                .map(item_of)
                .collect(),
            Layout::Pairs(pairs) => pairs.iter().map(|(t, _)| item_of(t.resolve(src))).collect(),
        }
    }

    pub fn has(&self, item: &str) -> bool {
        self.items().iter().any(|i| i.eq_ignore_ascii_case(item))
    }

    pub fn column(&self, item: &str) -> Option<Column<'a>> {
        let src = self.src;
        let matches = |tag: &Span| {
            let tag = tag.resolve(src);
            self.owns(tag) && item_of(tag).eq_ignore_ascii_case(item)
        };
        let cells = match &self.layout {
            Layout::Loop(lp) => ColumnCells::Strided {
                cells: &lp.cells,
                offset: lp.tags.iter().position(matches)?,
                stride: lp.tags.len(),
            },
            Layout::Pairs(pairs) => ColumnCells::Single(pairs.iter().find(|(t, _)| matches(t))?.1),
        };
        Some(Column { src, cells })
    }
}

#[derive(Debug, Clone, Copy)]
enum ColumnCells<'a> {
    Strided {
        cells: &'a [Cell],
        offset: usize,
        stride: usize,
    },
    Single(Cell),
}

#[derive(Debug, Clone, Copy)]
pub struct Column<'a> {
    src: &'a str,
    cells: ColumnCells<'a>,
}

impl<'a> Column<'a> {
    pub fn len(&self) -> usize {
        match self.cells {
            ColumnCells::Strided { cells, stride, .. } => cells.len() / stride,
            ColumnCells::Single(_) => 1,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn get(&self, row: usize) -> Option<Value<'a>> {
        match self.cells {
            ColumnCells::Strided {
                cells,
                offset,
                stride,
            } => {
                let index = row.checked_mul(stride)?.checked_add(offset)?;
                cells.get(index).map(|cell| cell.resolve(self.src))
            }
            ColumnCells::Single(cell) => (row == 0).then(|| cell.resolve(self.src)),
        }
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = Value<'a>> + use<'a> {
        let column = *self;
        (0..column.len()).map(move |row| column.get(row).expect("row is within the column"))
    }
}

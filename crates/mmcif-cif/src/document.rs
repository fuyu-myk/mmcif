use std::fmt;
use std::path::Path;

use crate::columns::{NumericColumn, RawColumn, RawTable, TypedValue};
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
    Table(RawTable),
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

/// A CIF document: its data blocks and the text they were read from.
///
/// Blocks, categories and values are views that borrow from the document. Two
/// documents are equal when their blocks have equal names and entries.
///
/// # Examples
///
/// ```
/// use mmcif_cif::Document;
///
/// let document = Document::parse("data_a\n_x.y 1\ndata_b\n_x.y 2\n")?;
/// let names: Vec<_> = document.blocks().map(|block| block.name()).collect();
///
/// assert_eq!(names, ["a", "b"]);
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
pub struct Document {
    source: Source,
    blocks: Vec<RawBlock>,
}

impl Document {
    /// Parses CIF text into a document.
    ///
    /// A leading UTF-8 byte order mark is skipped. Text with no data block,
    /// such as an empty string or only comments, gives an empty document.
    ///
    /// # Errors
    ///
    /// Returns a [`ParseError`] for the first problem in `text`, with one of
    /// these kinds:
    ///
    /// - [`UnterminatedQuote`] or [`UnterminatedTextField`] if a quoted value
    ///   or text field is not closed.
    /// - [`OutsideBlock`] if anything other than whitespace and comments comes
    ///   before the first `data_` header, and [`EmptyBlockName`] if a header
    ///   has no name.
    /// - [`MissingValue`] if a tag has no value, and [`UnexpectedValue`] if a
    ///   value outside a loop has no tag.
    /// - [`EmptyLoop`] if a `loop_` has no tags, and [`RaggedLoop`] if the
    ///   number of values in a loop is not a multiple of its number of tags.
    /// - [`DuplicateTag`] if a tag appears twice in a block or save frame,
    ///   ignoring ASCII case.
    /// - [`NestedSaveFrame`], [`UnmatchedSaveEnd`] or [`UnterminatedSaveFrame`]
    ///   if a save frame opens inside another, a `save_` terminator has no
    ///   frame to close, or a frame is not closed before its block ends.
    /// - [`Reserved`] if `global_` or `stop_` appears.
    /// - [`TooLarge`] if `text` is longer than `u32::MAX` bytes or holds a
    ///   value longer than `2^29 - 1` bytes.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{Document, ParseErrorKind};
    ///
    /// let document = Document::parse("data_1ABC\n_entry.id 1ABC\n")?;
    /// assert_eq!(document.len(), 1);
    ///
    /// let error = Document::parse("_entry.id 1ABC\n").unwrap_err();
    /// assert_eq!(error.kind, ParseErrorKind::OutsideBlock);
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    ///
    /// [`UnterminatedQuote`]: crate::ParseErrorKind::UnterminatedQuote
    /// [`UnterminatedTextField`]: crate::ParseErrorKind::UnterminatedTextField
    /// [`OutsideBlock`]: crate::ParseErrorKind::OutsideBlock
    /// [`EmptyBlockName`]: crate::ParseErrorKind::EmptyBlockName
    /// [`MissingValue`]: crate::ParseErrorKind::MissingValue
    /// [`UnexpectedValue`]: crate::ParseErrorKind::UnexpectedValue
    /// [`EmptyLoop`]: crate::ParseErrorKind::EmptyLoop
    /// [`RaggedLoop`]: crate::ParseErrorKind::RaggedLoop
    /// [`DuplicateTag`]: crate::ParseErrorKind::DuplicateTag
    /// [`NestedSaveFrame`]: crate::ParseErrorKind::NestedSaveFrame
    /// [`UnmatchedSaveEnd`]: crate::ParseErrorKind::UnmatchedSaveEnd
    /// [`UnterminatedSaveFrame`]: crate::ParseErrorKind::UnterminatedSaveFrame
    /// [`Reserved`]: crate::ParseErrorKind::Reserved
    /// [`TooLarge`]: crate::ParseErrorKind::TooLarge
    pub fn parse(text: impl Into<String>) -> Result<Document, ParseError> {
        Document::from_source(Source::Owned(text.into()))
    }

    /// Parses CIF text from UTF-8 bytes.
    ///
    /// # Errors
    ///
    /// Returns a [`ParseError`] of kind
    /// [`InvalidUtf8`](crate::ParseErrorKind::InvalidUtf8), positioned at the
    /// first invalid byte, if `bytes` is not valid UTF-8. Otherwise returns the
    /// errors of [`Document::parse`].
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{Document, ParseErrorKind};
    ///
    /// let document = Document::from_bytes(b"data_x\n_a.b 1\n".to_vec())?;
    /// assert_eq!(document.len(), 1);
    ///
    /// let error = Document::from_bytes(b"data_x\n_a.b \xff\n".to_vec()).unwrap_err();
    /// assert_eq!(error.kind, ParseErrorKind::InvalidUtf8);
    /// assert_eq!((error.line, error.column), (2, 6));
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn from_bytes(bytes: Vec<u8>) -> Result<Document, ParseError> {
        Document::from_source(Source::Owned(read_text(bytes)?))
    }

    /// Reads and parses the CIF file at `path`.
    ///
    /// The file is read with [`read_file`](crate::read_file), so it may be
    /// gzip-compressed when the `gzip` feature is on.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] if the file cannot be read, [`Error::Compressed`]
    /// if it is compressed and the `gzip` feature is off, and [`Error::Parse`]
    /// with the errors of [`Document::from_bytes`] if its contents are not
    /// valid CIF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use mmcif_cif::Document;
    ///
    /// let document = Document::open("1abc.cif")?;
    /// for block in document.blocks() {
    ///     println!("{}", block.name());
    /// }
    /// # Ok::<(), mmcif_cif::Error>(())
    /// ```
    pub fn open(path: impl AsRef<Path>) -> Result<Document, Error> {
        Ok(Document::from_bytes(crate::read_file(path)?)?)
    }

    /// Parses the CIF file at `path` through a memory map of its contents.
    ///
    /// The file must hold uncompressed CIF text. A file that is empty or is
    /// not a regular file is read into memory, not mapped.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] if the file cannot be opened, read or mapped, and
    /// [`Error::Parse`] with the errors of [`Document::from_bytes`] if its
    /// contents are not valid CIF.
    ///
    /// # Safety
    ///
    /// The file must not be modified or truncated while the document is alive.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use mmcif_cif::Document;
    ///
    /// let document = unsafe { Document::open_mapped("1abc.cif")? };
    /// println!("{} blocks", document.len());
    /// # Ok::<(), mmcif_cif::Error>(())
    /// ```
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

    /// Returns an iterator over the blocks, in the order they appear.
    pub fn blocks(&self) -> impl ExactSizeIterator<Item = Block<'_>> {
        let src = self.text();
        self.blocks.iter().map(move |raw| Block { src, raw })
    }

    /// Returns the block at `index`, or `None` if `index` is out of range.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    ///
    /// let document = Document::parse("data_a\n_x.y 1\n")?;
    ///
    /// assert_eq!(document.block(0).unwrap().name(), "a");
    /// assert!(document.block(1).is_none());
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn block(&self, index: usize) -> Option<Block<'_>> {
        self.blocks.get(index).map(|raw| Block {
            src: self.text(),
            raw,
        })
    }

    /// Returns the first block called `name`, ignoring ASCII case, or `None`
    /// if there is no such block.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    ///
    /// let document = Document::parse("data_1ABC\n_entry.id 1ABC\n")?;
    ///
    /// assert_eq!(document.block_named("1abc").unwrap().name(), "1ABC");
    /// assert!(document.block_named("2XYZ").is_none());
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn block_named(&self, name: &str) -> Option<Block<'_>> {
        self.blocks()
            .find(|block| block.name().eq_ignore_ascii_case(name))
    }

    /// Returns the number of blocks.
    pub fn len(&self) -> usize {
        self.blocks.len()
    }

    /// Returns `true` if the document has no blocks.
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

/// A data block parsed by a [`BlockStream`](crate::BlockStream).
///
/// It holds the entries of the block and borrows the text it was parsed from.
///
/// # Examples
///
/// ```
/// use mmcif_cif::stream_blocks;
///
/// for streamed in stream_blocks("data_a\nloop_\n_x.y\n1\n2\n") {
///     let streamed = streamed?;
///     assert_eq!(streamed.block().category("x").unwrap().row_count(), 2);
/// }
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
pub struct StreamedBlock<'a> {
    src: &'a str,
    raw: RawBlock,
}

impl<'a> StreamedBlock<'a> {
    pub(crate) fn new(src: &'a str, raw: RawBlock) -> Self {
        StreamedBlock { src, raw }
    }

    /// Returns a view of the block.
    pub fn block(&self) -> Block<'_> {
        Block {
            src: self.src,
            raw: &self.raw,
        }
    }
}

/// A view of one data block.
///
/// A block gives its entries in file order and its categories by name. Two
/// blocks are equal when their names and entries are equal.
///
/// # Examples
///
/// ```
/// use mmcif_cif::Document;
///
/// let document = Document::parse("data_1ABC\n_entry.id 1ABC\n_cell.length_a 10.0\n")?;
/// let block = document.block(0).unwrap();
///
/// assert_eq!(block.name(), "1ABC");
/// assert_eq!(block.entries().len(), 2);
/// assert_eq!(block.category_names(), ["entry", "cell"]);
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
#[derive(Clone, Copy)]
pub struct Block<'a> {
    src: &'a str,
    raw: &'a RawBlock,
}

impl<'a> Block<'a> {
    /// Returns the name of the block, without the `data_` prefix.
    pub fn name(&self) -> &'a str {
        self.raw.name.resolve(self.src)
    }

    /// Returns an iterator over the pairs, loops and save frames of the block,
    /// in file order.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{Document, Entry};
    ///
    /// let text = "data_x\n_entry.id 1ABC\nloop_\n_atom.id\n1\n2\nsave_f\n_item.name x\nsave_\n";
    /// let document = Document::parse(text)?;
    ///
    /// let kinds: Vec<_> = document
    ///     .block(0)
    ///     .unwrap()
    ///     .entries()
    ///     .map(|entry| match entry {
    ///         Entry::Pair(_) => "pair",
    ///         Entry::Loop(_) => "loop",
    ///         Entry::Frame(_) => "frame",
    ///     })
    ///     .collect();
    /// assert_eq!(kinds, ["pair", "loop", "frame"]);
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn entries(&self) -> impl ExactSizeIterator<Item = Entry<'a>> + use<'a> {
        let src = self.src;
        self.raw.entries.iter().map(move |raw| Entry::new(src, raw))
    }

    /// Returns the category called `name`, or `None` if no tag in the block
    /// belongs to it.
    ///
    /// A tag belongs to the category named by the text between its leading
    /// `_` and its first `.`, and a tag without a `.` belongs to none. `name`
    /// is matched ignoring ASCII case, with or without a leading `_`. If a
    /// loop holds a tag of the category, the category is read from the first
    /// such loop. Otherwise it is read as a single row from the pairs of the
    /// block whose tags belong to it. Save frames are not searched.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    ///
    /// let document = Document::parse("data_x\n_cell.length_a 10.0\n_cell.length_b 20.0\n")?;
    /// let block = document.block(0).unwrap();
    ///
    /// let cell = block.category("_CELL").unwrap();
    /// assert_eq!(cell.row_count(), 1);
    /// assert_eq!(cell.items(), ["length_a", "length_b"]);
    /// assert!(block.category("atom_site").is_none());
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn category(&self, name: &str) -> Option<Category<'a>> {
        self.category_named(name)
            .or_else(|| self.category_named(name.strip_prefix('_')?))
    }

    fn category_named(&self, name: &str) -> Option<Category<'a>> {
        let src = self.src;
        let matching =
            |tag: Span| category_of(tag.resolve(src)).filter(|c| c.eq_ignore_ascii_case(name));
        for raw in &self.raw.entries {
            let (tags, layout) = match raw {
                RawEntry::Loop(lp) => (&lp.tags, Layout::Loop(lp)),
                RawEntry::Table(table) => (&table.tags, Layout::Table(table)),
                _ => continue,
            };
            if let Some(found) = tags.iter().find_map(|&tag| matching(tag)) {
                return Some(Category {
                    src,
                    name: found,
                    layout,
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

    /// Returns the name of each category in the block, in order of first
    /// appearance.
    ///
    /// Each name is written as in the first tag of its category, and names
    /// that differ only in ASCII case are listed once. Tags in save frames are
    /// not included.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    ///
    /// let text = "data_x\n_entry.id 1ABC\nloop_\n_atom.id\n_ATOM.name\n1 N\n";
    /// let document = Document::parse(text)?;
    ///
    /// assert_eq!(document.block(0).unwrap().category_names(), ["entry", "atom"]);
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
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

/// An entry of a block or save frame.
///
/// # Examples
///
/// ```
/// use mmcif_cif::{Document, Entry};
///
/// let document = Document::parse("data_x\n_entry.id 1ABC\nloop_\n_atom.id\n1\n2\n")?;
///
/// for entry in document.block(0).unwrap().entries() {
///     match entry {
///         Entry::Pair(pair) => assert_eq!(pair.tag(), "_entry.id"),
///         Entry::Loop(atoms) => assert_eq!(atoms.row_count(), 2),
///         Entry::Frame(frame) => println!("save frame {}", frame.name()),
///     }
/// }
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Entry<'a> {
    /// A tag and its value.
    Pair(Pair<'a>),
    /// A `loop_` with its tags and rows.
    Loop(Loop<'a>),
    /// A save frame.
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
            RawEntry::Loop(raw) => Entry::Loop(Loop {
                src,
                data: LoopData::Rows(raw),
            }),
            RawEntry::Table(raw) => Entry::Loop(Loop {
                src,
                data: LoopData::Columns(raw),
            }),
            RawEntry::Frame(raw) => Entry::Frame(Frame { src, raw }),
        }
    }
}

/// A tag and its value outside a loop.
///
/// # Examples
///
/// ```
/// use mmcif_cif::{Document, Entry};
///
/// let document = Document::parse("data_x\n_entry.id 1ABC\n")?;
///
/// let Some(Entry::Pair(pair)) = document.block(0).unwrap().entries().next() else {
///     panic!("expected a pair");
/// };
/// assert_eq!(pair.tag(), "_entry.id");
/// assert_eq!(pair.value().as_str(), Some("1ABC"));
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
#[derive(Clone, Copy)]
pub struct Pair<'a> {
    src: &'a str,
    tag: Span,
    value: Cell,
}

impl<'a> Pair<'a> {
    /// Returns the tag, including its leading `_`.
    pub fn tag(&self) -> &'a str {
        self.tag.resolve(self.src)
    }

    /// Returns the value.
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
enum LoopData<'a> {
    Rows(&'a RawLoop),
    Columns(&'a RawTable),
}

/// A `loop_`: a list of tags and rows holding one value for each tag.
///
/// Rows and columns are numbered from 0, and columns follow the order of the
/// tags.
///
/// # Examples
///
/// ```
/// use mmcif_cif::{Document, Entry};
///
/// let document = Document::parse("data_x\nloop_\n_atom.id\n_atom.name\n1 N\n2 CA\n")?;
///
/// let Some(Entry::Loop(atoms)) = document.block(0).unwrap().entries().next() else {
///     panic!("expected a loop");
/// };
/// assert_eq!(atoms.tags().collect::<Vec<_>>(), ["_atom.id", "_atom.name"]);
/// assert_eq!(atoms.row_count(), 2);
/// assert_eq!(atoms.value(1, 1).unwrap().as_str(), Some("CA"));
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
#[derive(Clone, Copy)]
pub struct Loop<'a> {
    src: &'a str,
    data: LoopData<'a>,
}

impl<'a> Loop<'a> {
    fn raw_tags(&self) -> &'a [Span] {
        match self.data {
            LoopData::Rows(raw) => &raw.tags,
            LoopData::Columns(raw) => &raw.tags,
        }
    }

    /// Returns an iterator over the tags of the loop, in order.
    pub fn tags(&self) -> impl ExactSizeIterator<Item = &'a str> + use<'a> {
        let src = self.src;
        self.raw_tags().iter().map(move |tag| tag.resolve(src))
    }

    /// Returns the number of rows.
    pub fn row_count(&self) -> usize {
        match self.data {
            LoopData::Rows(raw) => raw.cells.len() / raw.tags.len(),
            LoopData::Columns(raw) => raw.rows,
        }
    }

    /// Returns an iterator over every value of the loop, row by row.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{Document, Entry};
    ///
    /// let document = Document::parse("data_x\nloop_\n_atom.id\n_atom.name\n1 N\n2 CA\n")?;
    /// let Some(Entry::Loop(atoms)) = document.block(0).unwrap().entries().next() else {
    ///     panic!("expected a loop");
    /// };
    ///
    /// let values: Vec<_> = atoms.values().map(|value| value.raw()).collect();
    /// assert_eq!(values, ["1", "N", "2", "CA"]);
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn values(&self) -> impl ExactSizeIterator<Item = Value<'a>> + use<'a> {
        let src = self.src;
        let data = self.data;
        let width = self.raw_tags().len();
        (0..self.row_count() * width).map(move |index| match data {
            LoopData::Rows(raw) => raw.cells[index].resolve(src),
            LoopData::Columns(raw) => column_cells(&raw.columns[index % width])
                .get(src, index / width)
                .expect("index is within the loop"),
        })
    }

    /// Returns the value at `row` and `column`, or `None` if either is out of
    /// range.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{Document, Entry};
    ///
    /// let document = Document::parse("data_x\nloop_\n_atom.id\n_atom.name\n1 N\n2 CA\n")?;
    /// let Some(Entry::Loop(atoms)) = document.block(0).unwrap().entries().next() else {
    ///     panic!("expected a loop");
    /// };
    ///
    /// assert_eq!(atoms.value(0, 1).unwrap().as_str(), Some("N"));
    /// assert!(atoms.value(2, 0).is_none());
    /// assert!(atoms.value(0, 2).is_none());
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn value(&self, row: usize, column: usize) -> Option<Value<'a>> {
        match self.data {
            LoopData::Rows(raw) => {
                if column >= raw.tags.len() {
                    return None;
                }
                let index = row.checked_mul(raw.tags.len())?.checked_add(column)?;
                raw.cells.get(index).map(|cell| cell.resolve(self.src))
            }
            LoopData::Columns(raw) => column_cells(raw.columns.get(column)?).get(self.src, row),
        }
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

/// A save frame: a named group of entries inside a block.
///
/// # Examples
///
/// ```
/// use mmcif_cif::{Document, Entry};
///
/// let document = Document::parse("data_dict\nsave_atom_site\n_item.name x\nsave_\n")?;
///
/// let Some(Entry::Frame(frame)) = document.block(0).unwrap().entries().next() else {
///     panic!("expected a save frame");
/// };
/// assert_eq!(frame.name(), "atom_site");
/// assert_eq!(frame.entries().len(), 1);
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
#[derive(Clone, Copy)]
pub struct Frame<'a> {
    src: &'a str,
    raw: &'a RawFrame,
}

impl<'a> Frame<'a> {
    /// Returns the name of the frame, without the `save_` prefix.
    pub fn name(&self) -> &'a str {
        self.raw.name.resolve(self.src)
    }

    /// Returns an iterator over the pairs and loops of the frame, in file
    /// order.
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
    Table(&'a RawTable),
    Pairs(Vec<(Span, Cell)>),
}

/// A view of one category of a block, as a table of rows and items.
///
/// Each item is a column with one value per row. A category read from a loop
/// has the rows of the loop, and a category read from pairs has one row. Item
/// names are matched ignoring ASCII case.
///
/// # Examples
///
/// ```
/// use mmcif_cif::Document;
///
/// let text = "data_x\nloop_\n_atom_site.id\n_atom_site.label_atom_id\n1 N\n2 CA\n";
/// let document = Document::parse(text)?;
/// let atoms = document.block(0).unwrap().category("atom_site").unwrap();
///
/// assert_eq!(atoms.row_count(), 2);
/// assert_eq!(atoms.items(), ["id", "label_atom_id"]);
/// assert_eq!(atoms.column("LABEL_ATOM_ID").unwrap().get(0).unwrap().as_str(), Some("N"));
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
#[derive(Debug, Clone)]
pub struct Category<'a> {
    src: &'a str,
    name: &'a str,
    layout: Layout<'a>,
}

impl<'a> Category<'a> {
    /// Returns the name of the category as written in its first tag, without
    /// the leading `_`.
    pub fn name(&self) -> &'a str {
        self.name
    }

    /// Returns the number of rows, which is 1 for a category read from pairs.
    pub fn row_count(&self) -> usize {
        match &self.layout {
            Layout::Loop(lp) => lp.cells.len() / lp.tags.len(),
            Layout::Table(table) => table.rows,
            Layout::Pairs(_) => 1,
        }
    }

    fn owns(&self, tag: &str) -> bool {
        category_of(tag).is_some_and(|category| category.eq_ignore_ascii_case(self.name))
    }

    /// Returns the item names of the category, in order.
    ///
    /// An item name is the part of a tag after its first `.`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    ///
    /// let document = Document::parse("data_x\n_cell.length_a 10.0\n_cell.length_b 20.0\n")?;
    /// let cell = document.block(0).unwrap().category("cell").unwrap();
    ///
    /// assert_eq!(cell.items(), ["length_a", "length_b"]);
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn items(&self) -> Vec<&'a str> {
        let src = self.src;
        let owned = |tags: &'a [Span]| -> Vec<&'a str> {
            tags.iter()
                .map(|t| t.resolve(src))
                .filter(|tag| self.owns(tag))
                .map(item_of)
                .collect()
        };
        match &self.layout {
            Layout::Loop(lp) => owned(&lp.tags),
            Layout::Table(table) => owned(&table.tags),
            Layout::Pairs(pairs) => pairs.iter().map(|(t, _)| item_of(t.resolve(src))).collect(),
        }
    }

    /// Returns `true` if the category has an item called `item`, ignoring
    /// ASCII case.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    ///
    /// let document = Document::parse("data_x\n_cell.length_a 10.0\n")?;
    /// let cell = document.block(0).unwrap().category("cell").unwrap();
    ///
    /// assert!(cell.has("LENGTH_A"));
    /// assert!(!cell.has("length_b"));
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn has(&self, item: &str) -> bool {
        self.items().iter().any(|i| i.eq_ignore_ascii_case(item))
    }

    /// Returns the values of the item called `item`, ignoring ASCII case, or
    /// `None` if the category has no such item.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    ///
    /// let document = Document::parse("data_x\nloop_\n_atom.id\n_atom.name\n1 N\n2 CA\n")?;
    /// let atoms = document.block(0).unwrap().category("atom").unwrap();
    ///
    /// let names = atoms.column("name").unwrap();
    /// assert_eq!(names.get(1).unwrap().as_str(), Some("CA"));
    /// assert!(atoms.column("charge").is_none());
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
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
            Layout::Table(table) => {
                column_cells(&table.columns[table.tags.iter().position(matches)?])
            }
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
    Text(&'a [Cell]),
    Numeric(&'a NumericColumn),
}

fn column_cells(column: &RawColumn) -> ColumnCells<'_> {
    match column {
        RawColumn::Text(cells) => ColumnCells::Text(cells),
        RawColumn::Numeric(numeric) => ColumnCells::Numeric(numeric),
    }
}

impl<'a> ColumnCells<'a> {
    fn get(self, src: &'a str, row: usize) -> Option<Value<'a>> {
        match self {
            ColumnCells::Strided {
                cells,
                offset,
                stride,
            } => {
                let index = row.checked_mul(stride)?.checked_add(offset)?;
                cells.get(index).map(|cell| cell.resolve(src))
            }
            ColumnCells::Single(cell) => (row == 0).then(|| cell.resolve(src)),
            ColumnCells::Text(cells) => cells.get(row).map(|cell| cell.resolve(src)),
            ColumnCells::Numeric(column) => column.value(row),
        }
    }
}

/// The values of one item of a [`Category`], one for each row.
///
/// # Examples
///
/// ```
/// use mmcif_cif::Document;
///
/// let document = Document::parse("data_x\nloop_\n_atom.name\nN\nCA\nC\n")?;
/// let atoms = document.block(0).unwrap().category("atom").unwrap();
/// let names = atoms.column("name").unwrap();
///
/// assert_eq!(names.len(), 3);
/// assert_eq!(names.iter().map(|value| value.raw()).collect::<Vec<_>>(), ["N", "CA", "C"]);
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
#[derive(Debug, Clone, Copy)]
pub struct Column<'a> {
    src: &'a str,
    cells: ColumnCells<'a>,
}

impl<'a> Column<'a> {
    /// Returns the number of values, which is the row count of the category.
    pub fn len(&self) -> usize {
        match self.cells {
            ColumnCells::Strided { cells, stride, .. } => cells.len() / stride,
            ColumnCells::Single(_) => 1,
            ColumnCells::Text(cells) => cells.len(),
            ColumnCells::Numeric(column) => column.len(),
        }
    }

    /// Returns `true` if the column has no values.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the value at `row`, or `None` if `row` is out of range.
    ///
    /// A value of a numeric BinaryCIF column is its number written as bare
    /// text, or `.` or `?` where the column is masked.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    ///
    /// let document = Document::parse("data_x\nloop_\n_atom.alt_id\nA\n.\n")?;
    /// let atoms = document.block(0).unwrap().category("atom").unwrap();
    /// let alt_ids = atoms.column("alt_id").unwrap();
    ///
    /// assert_eq!(alt_ids.get(0).unwrap().as_str(), Some("A"));
    /// assert!(alt_ids.get(1).unwrap().is_inapplicable());
    /// assert!(alt_ids.get(2).is_none());
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn get(&self, row: usize) -> Option<Value<'a>> {
        self.cells.get(self.src, row)
    }

    /// Returns the value at `row` as a number, if the column holds numbers
    /// decoded from BinaryCIF.
    ///
    /// Returns `None` if `row` is out of range or the column holds text.
    /// Columns parsed from CIF text, BinaryCIF string columns and the columns
    /// of single-row BinaryCIF categories hold text. A
    /// [`TypedValue::Float`] equals the text of [`get`](Column::get) parsed
    /// as an `f64`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{Document, TypedValue};
    ///
    /// # let bytes = b"\x81\xaadataBlocks\x91\x82\xa6header\xa41ABC\xaacategories\x91\x83\
    /// #     \xa4name\xa5_atom\xa8rowCount\x02\xa7columns\x91\x82\xa4name\xa2id\xa4data\
    /// #     \x82\xa4data\xc4\x02\x01\x02\xa8encoding\x91\x82\xa4kind\xa9ByteArray\xa4type\x01";
    /// let document = Document::from_bcif(bytes)?;
    /// let atoms = document.block(0).unwrap().category("atom").unwrap();
    /// let ids = atoms.column("id").unwrap();
    ///
    /// assert_eq!(ids.get(1).unwrap().raw(), "2");
    /// assert_eq!(ids.typed(1), Some(TypedValue::Int(2)));
    /// assert_eq!(ids.typed(2), None);
    /// # Ok::<(), mmcif_cif::BcifError>(())
    /// ```
    pub fn typed(&self, row: usize) -> Option<TypedValue> {
        match self.cells {
            ColumnCells::Numeric(column) => column.typed(row),
            _ => None,
        }
    }

    /// Returns an iterator over the values, in row order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = Value<'a>> + use<'a> {
        let column = *self;
        (0..column.len()).map(move |row| column.get(row).expect("row is within the column"))
    }
}

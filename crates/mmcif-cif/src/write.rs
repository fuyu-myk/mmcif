use std::io::{self, Write};

use crate::document::{Block, Document, Entry};
use crate::value::{Value, ValueKind};

struct CifWriter<W> {
    out: W,
    at_line_start: bool,
}

impl<W: Write> CifWriter<W> {
    fn raw(&mut self, text: &str) -> io::Result<()> {
        if !text.is_empty() {
            self.out.write_all(text.as_bytes())?;
            self.at_line_start = text.ends_with('\n');
        }
        Ok(())
    }

    fn end_line(&mut self) -> io::Result<()> {
        if !self.at_line_start {
            self.raw("\n")?;
        }
        Ok(())
    }

    fn separator(&mut self) -> io::Result<()> {
        if !self.at_line_start {
            self.raw(" ")?;
        }
        Ok(())
    }

    fn value(&mut self, value: Value<'_>) -> io::Result<()> {
        let raw = value.raw();
        match value.kind() {
            ValueKind::TextField => {
                self.end_line()?;
                self.raw(";")?;
                self.raw(raw)?;
                if raw.ends_with('\r') {
                    self.raw("\r")?;
                }
                self.raw("\n;\n")
            }
            ValueKind::Bare if raw.starts_with(';') => {
                self.raw(" ")?;
                self.raw(raw)
            }
            ValueKind::Bare | ValueKind::Inapplicable | ValueKind::Unknown => {
                self.separator()?;
                self.raw(raw)
            }
            ValueKind::SingleQuoted => {
                self.separator()?;
                self.raw("'")?;
                self.raw(raw)?;
                self.raw("'")
            }
            ValueKind::DoubleQuoted => {
                self.separator()?;
                self.raw("\"")?;
                self.raw(raw)?;
                self.raw("\"")
            }
        }
    }

    fn entry(&mut self, entry: Entry<'_>) -> io::Result<()> {
        match entry {
            Entry::Pair(pair) => {
                self.raw(pair.tag())?;
                self.value(pair.value())?;
                self.end_line()
            }
            Entry::Loop(lp) => {
                self.raw("loop_\n")?;
                for tag in lp.tags() {
                    self.raw(tag)?;
                    self.raw("\n")?;
                }
                let width = lp.tags().len();
                for (index, value) in lp.values().enumerate() {
                    self.value(value)?;
                    if (index + 1) % width == 0 {
                        self.end_line()?;
                    }
                }
                Ok(())
            }
            Entry::Frame(frame) => {
                self.raw("save_")?;
                self.raw(frame.name())?;
                self.raw("\n")?;
                for inner in frame.entries() {
                    self.entry(inner)?;
                }
                self.raw("save_\n")
            }
        }
    }

    fn block(&mut self, block: Block<'_>) -> io::Result<()> {
        self.raw("data_")?;
        self.raw(block.name())?;
        self.raw("\n#\n")?;
        for entry in block.entries() {
            self.entry(entry)?;
            self.raw("#\n")?;
        }
        Ok(())
    }
}

impl Document {
    /// Writes the document to `out` as CIF text.
    ///
    /// Each block starts with its `data_` header line and a `#` line, and
    /// each entry of the block is followed by a `#` line. A pair is written as
    /// its tag followed by its value, a loop as `loop_`, its tags one per line
    /// and then its rows, each starting on a new line, and a save frame as its
    /// entries between `save_` lines. Values keep their kind, so quoted values
    /// are written in their quotes, and a text field starts on a new line and
    /// ends with a line holding only `;`. Parsing the output gives a document
    /// equal to this one.
    ///
    /// # Errors
    ///
    /// Returns any I/O error from writing to or flushing `out`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    ///
    /// let document = Document::parse("data_x\n_struct.title   'two words'\n")?;
    ///
    /// let mut out = Vec::new();
    /// document.write_to(&mut out)?;
    /// assert_eq!(out, b"data_x\n#\n_struct.title 'two words'\n#\n");
    /// # Ok::<(), mmcif_cif::Error>(())
    /// ```
    pub fn write_to(&self, out: impl Write) -> io::Result<()> {
        let mut writer = CifWriter {
            out: io::BufWriter::new(out),
            at_line_start: true,
        };
        for block in self.blocks() {
            writer.block(block)?;
        }
        writer.out.flush()
    }

    /// Returns the document as CIF text, in the layout of
    /// [`Document::write_to`].
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    ///
    /// let text = "data_x\n_a.note\n;line one\nline two\n;\nloop_\n_b.id _b.name\n1 N 2 'C A'\n";
    /// let document = Document::parse(text)?;
    ///
    /// assert_eq!(
    ///     document.to_cif_string(),
    ///     "data_x\n#\n_a.note\n;line one\nline two\n;\n#\nloop_\n_b.id\n_b.name\n1 N\n2 'C A'\n#\n"
    /// );
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn to_cif_string(&self) -> String {
        let mut out = Vec::new();
        self.write_to(&mut out)
            .expect("writing to a Vec<u8> cannot fail");
        String::from_utf8(out).expect("documents only hold UTF-8 text")
    }
}

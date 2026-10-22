/// A range of byte offsets into a source text.
///
/// The start is inclusive and the end exclusive. A span never ends before it
/// starts.
///
/// # Examples
///
/// ```
/// use mmcif_spec::Span;
///
/// let text = "id: B";
/// let span = Span::new(4, 5);
///
/// assert_eq!(&text[span.start()..span.end()], "B");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    start: usize,
    end: usize,
}

impl Span {
    /// Creates a span from byte offset `start` up to, but not including,
    /// `end`.
    ///
    /// If `end` is less than `start`, the span is empty and ends at `start`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::Span;
    ///
    /// assert_eq!(Span::new(2, 6).end(), 6);
    /// assert_eq!(Span::new(9, 4).end(), 9);
    /// ```
    pub fn new(start: usize, end: usize) -> Span {
        Span {
            start,
            end: end.max(start),
        }
    }

    /// Returns the offset of the first byte in the span.
    pub fn start(self) -> usize {
        self.start
    }

    /// Returns the offset one past the last byte in the span.
    pub fn end(self) -> usize {
        self.end
    }
}

/// The id of a file added to a [`Sources`].
///
/// # Examples
///
/// ```
/// use mmcif_spec::Sources;
///
/// let mut sources = Sources::new();
/// let id = sources.add("design.yaml", "entities: []\n");
///
/// assert_eq!(sources.file(id).map(|file| file.name()), Some("design.yaml"));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SourceId(u32);

/// A line and column in a source text, both counted from 1.
///
/// Columns count characters, not bytes.
///
/// # Examples
///
/// ```
/// use mmcif_spec::{Location, SourceFile};
///
/// let file = SourceFile::new("design.yaml", "entities:\n  - protein:\n");
///
/// assert_eq!(file.location(14), Location { line: 2, column: 5 });
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Location {
    /// The line number, counted from 1.
    pub line: usize,
    /// The column number in characters, counted from 1.
    pub column: usize,
}

/// A named source text, indexed by line.
///
/// # Examples
///
/// ```
/// use mmcif_spec::SourceFile;
///
/// let file = SourceFile::new("design.yaml", "entities:\n  - protein:\n");
///
/// assert_eq!(file.name(), "design.yaml");
/// assert_eq!(file.line(2), "  - protein:");
/// ```
#[derive(Debug, Clone)]
pub struct SourceFile {
    name: String,
    text: String,
    line_starts: Vec<usize>,
}

impl SourceFile {
    /// Creates a source file called `name` holding `text`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::SourceFile;
    ///
    /// let file = SourceFile::new("design.yaml", "a: 1\n");
    ///
    /// assert_eq!(file.text(), "a: 1\n");
    /// ```
    pub fn new(name: impl Into<String>, text: impl Into<String>) -> SourceFile {
        let text: String = text.into();
        let mut line_starts = vec![0];
        line_starts.extend(text.match_indices('\n').map(|(offset, _)| offset + 1));
        SourceFile {
            name: name.into(),
            text,
            line_starts,
        }
    }

    /// Returns the name the file was created with.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the full text of the file.
    pub fn text(&self) -> &str {
        &self.text
    }

    /// Returns the line and column of the byte at `offset`.
    ///
    /// An offset past the end of the text is taken as the end of the text,
    /// and an offset inside a multi-byte character as the start of that
    /// character.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Location, SourceFile};
    ///
    /// let file = SourceFile::new("design.yaml", "entities:\n  - protein:\n");
    ///
    /// assert_eq!(file.location(0), Location { line: 1, column: 1 });
    /// assert_eq!(file.location(14), Location { line: 2, column: 5 });
    /// assert_eq!(file.location(999), Location { line: 3, column: 1 });
    /// ```
    pub fn location(&self, offset: usize) -> Location {
        let mut offset = offset.min(self.text.len());
        while !self.text.is_char_boundary(offset) {
            offset -= 1;
        }
        let line = self
            .line_starts
            .partition_point(|start| *start <= offset)
            .max(1);
        let start = self.line_starts[line - 1];
        Location {
            line,
            column: self.text[start..offset].chars().count() + 1,
        }
    }

    /// Returns the byte offset at which line `line` starts, counting lines
    /// from 1.
    ///
    /// Returns the length of the text for line 0 and for lines past the end.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::SourceFile;
    ///
    /// let file = SourceFile::new("design.yaml", "a: 1\nb: 2");
    ///
    /// assert_eq!(file.line_start(1), 0);
    /// assert_eq!(file.line_start(2), 5);
    /// assert_eq!(file.line_start(3), 9);
    /// ```
    pub fn line_start(&self, line: usize) -> usize {
        line.checked_sub(1)
            .and_then(|index| self.line_starts.get(index))
            .copied()
            .unwrap_or(self.text.len())
    }

    /// Returns the text of line `line` without its line ending, counting lines
    /// from 1.
    ///
    /// Returns an empty string for line 0 and for lines past the end.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::SourceFile;
    ///
    /// let file = SourceFile::new("design.yaml", "a: 1\r\nb: 2\r\n");
    ///
    /// assert_eq!(file.line(1), "a: 1");
    /// assert_eq!(file.line(2), "b: 2");
    /// assert_eq!(file.line(7), "");
    /// ```
    pub fn line(&self, line: usize) -> &str {
        match line
            .checked_sub(1)
            .and_then(|index| self.line_starts.get(index))
        {
            None => "",
            Some(start) => {
                let end = self
                    .line_starts
                    .get(line)
                    .map(|next| next - 1)
                    .unwrap_or(self.text.len());
                self.text[*start..end].trim_end_matches('\r')
            }
        }
    }
}

/// A set of source files that diagnostics point into.
///
/// # Examples
///
/// ```
/// use mmcif_spec::Sources;
///
/// let mut sources = Sources::new();
/// let yaml = sources.add("design.yaml", "entities: []\n");
/// let cif = sources.add("target.cif", "data_target\n");
///
/// assert_eq!(sources.file(yaml).unwrap().name(), "design.yaml");
/// assert_eq!(sources.file(cif).unwrap().text(), "data_target\n");
/// ```
#[derive(Debug, Clone, Default)]
pub struct Sources {
    files: Vec<SourceFile>,
}

impl Sources {
    /// Creates an empty set of sources.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::Sources;
    ///
    /// let mut sources = Sources::new();
    /// let id = sources.add("design.yaml", "entities: []\n");
    ///
    /// assert!(sources.file(id).is_some());
    /// ```
    pub fn new() -> Sources {
        Sources::default()
    }

    /// Adds a file called `name` holding `text` and returns its id.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::Sources;
    ///
    /// let mut sources = Sources::new();
    /// let first = sources.add("design.yaml", "");
    /// let second = sources.add("target.cif", "");
    ///
    /// assert_ne!(first, second);
    /// assert_eq!(sources.file(second).unwrap().name(), "target.cif");
    /// ```
    pub fn add(&mut self, name: impl Into<String>, text: impl Into<String>) -> SourceId {
        self.files.push(SourceFile::new(name, text));
        SourceId((self.files.len() - 1) as u32)
    }

    /// Returns the file with id `id`, or `None` if this set holds no such
    /// file.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::Sources;
    ///
    /// let mut sources = Sources::new();
    /// let id = sources.add("design.yaml", "entities: []\n");
    ///
    /// assert_eq!(sources.file(id).unwrap().text(), "entities: []\n");
    /// assert!(Sources::new().file(id).is_none());
    /// ```
    pub fn file(&self, id: SourceId) -> Option<&SourceFile> {
        self.files.get(id.0 as usize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locations_are_one_based_in_lines_and_characters() {
        let file = SourceFile::new("design.yaml", "entities:\n  - protein:\n      id: B\n");
        assert_eq!(file.location(0), Location { line: 1, column: 1 });
        assert_eq!(file.location(10), Location { line: 2, column: 1 });
        assert_eq!(file.location(14), Location { line: 2, column: 5 });
        assert_eq!(file.line(2), "  - protein:");
    }

    #[test]
    fn columns_count_characters_not_bytes() {
        let file = SourceFile::new("design.yaml", "name: \u{00e5}\u{00e5}x\n");
        let offset = file.text().find('x').unwrap();
        assert_eq!(file.location(offset), Location { line: 1, column: 9 });
    }

    #[test]
    fn an_offset_past_the_end_lands_on_the_last_line() {
        let file = SourceFile::new("design.yaml", "a: 1");
        assert_eq!(file.location(999), Location { line: 1, column: 5 });
        assert_eq!(file.line(7), "");
        assert_eq!(file.line(0), "");
    }

    #[test]
    fn carriage_returns_are_trimmed_from_reported_lines() {
        let file = SourceFile::new("design.yaml", "a: 1\r\nb: 2\r\n");
        assert_eq!(file.line(1), "a: 1");
        assert_eq!(file.location(6), Location { line: 2, column: 1 });
    }

    #[test]
    fn a_span_never_ends_before_it_starts() {
        assert_eq!(Span::new(9, 4).end(), 9);
    }

    #[test]
    fn sources_hand_back_what_they_were_given() {
        let mut sources = Sources::new();
        let yaml = sources.add("design.yaml", "a: 1");
        let cif = sources.add("8r3a.cif", "data_8r3a");
        assert_eq!(sources.file(yaml).unwrap().name(), "design.yaml");
        assert_eq!(sources.file(cif).unwrap().text(), "data_8r3a");
    }
}

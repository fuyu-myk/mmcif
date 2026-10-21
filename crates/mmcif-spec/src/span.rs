#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    start: usize,
    end: usize,
}

impl Span {
    pub fn new(start: usize, end: usize) -> Span {
        Span {
            start,
            end: end.max(start),
        }
    }

    pub fn start(self) -> usize {
        self.start
    }

    pub fn end(self) -> usize {
        self.end
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SourceId(u32);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Location {
    pub line: usize,
    pub column: usize,
}

#[derive(Debug, Clone)]
pub struct SourceFile {
    name: String,
    text: String,
    line_starts: Vec<usize>,
}

impl SourceFile {
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

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn text(&self) -> &str {
        &self.text
    }

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

    pub fn line_start(&self, line: usize) -> usize {
        line.checked_sub(1)
            .and_then(|index| self.line_starts.get(index))
            .copied()
            .unwrap_or(self.text.len())
    }

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

#[derive(Debug, Clone, Default)]
pub struct Sources {
    files: Vec<SourceFile>,
}

impl Sources {
    pub fn new() -> Sources {
        Sources::default()
    }

    pub fn add(&mut self, name: impl Into<String>, text: impl Into<String>) -> SourceId {
        self.files.push(SourceFile::new(name, text));
        SourceId((self.files.len() - 1) as u32)
    }

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

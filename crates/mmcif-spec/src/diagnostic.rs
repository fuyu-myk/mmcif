use std::fmt::Write as _;

use crate::span::{SourceFile, SourceId, Sources, Span};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    Warning,
    Error,
}

impl Severity {
    pub fn as_str(self) -> &'static str {
        match self {
            Severity::Warning => "warning",
            Severity::Error => "error",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Label {
    pub source: SourceId,
    pub span: Span,
    pub message: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    pub severity: Severity,
    pub code: &'static str,
    pub message: String,
    pub labels: Vec<Label>,
    pub help: Option<String>,
}

impl Diagnostic {
    pub fn error(code: &'static str, message: impl Into<String>) -> Diagnostic {
        Diagnostic::new(Severity::Error, code, message)
    }

    pub fn warning(code: &'static str, message: impl Into<String>) -> Diagnostic {
        Diagnostic::new(Severity::Warning, code, message)
    }

    fn new(severity: Severity, code: &'static str, message: impl Into<String>) -> Diagnostic {
        Diagnostic {
            severity,
            code,
            message: message.into(),
            labels: Vec::new(),
            help: None,
        }
    }

    #[must_use]
    pub fn at(mut self, source: SourceId, span: Span, message: impl Into<String>) -> Diagnostic {
        self.labels.push(Label {
            source,
            span,
            message: message.into(),
        });
        self
    }

    #[must_use]
    pub fn help(mut self, help: impl Into<String>) -> Diagnostic {
        self.help = Some(help.into());
        self
    }

    pub fn is_error(&self) -> bool {
        self.severity == Severity::Error
    }
}

#[derive(Debug, Clone, Default)]
pub struct Report {
    diagnostics: Vec<Diagnostic>,
}

impl Report {
    pub fn new() -> Report {
        Report::default()
    }

    pub fn push(&mut self, diagnostic: Diagnostic) {
        self.diagnostics.push(diagnostic);
    }

    pub fn extend(&mut self, diagnostics: impl IntoIterator<Item = Diagnostic>) {
        self.diagnostics.extend(diagnostics);
    }

    pub fn diagnostics(&self) -> &[Diagnostic] {
        &self.diagnostics
    }

    pub fn errors(&self) -> impl Iterator<Item = &Diagnostic> {
        self.diagnostics.iter().filter(|entry| entry.is_error())
    }

    pub fn has_errors(&self) -> bool {
        self.errors().next().is_some()
    }

    pub fn codes(&self) -> Vec<&'static str> {
        self.diagnostics.iter().map(|entry| entry.code).collect()
    }

    pub fn is_empty(&self) -> bool {
        self.diagnostics.is_empty()
    }

    pub fn len(&self) -> usize {
        self.diagnostics.len()
    }
}

impl Sources {
    pub fn render(&self, diagnostic: &Diagnostic) -> String {
        let mut out = String::new();
        let _ = writeln!(
            out,
            "{}[{}]: {}",
            diagnostic.severity.as_str(),
            diagnostic.code,
            diagnostic.message
        );
        let width = diagnostic
            .labels
            .iter()
            .filter_map(|label| {
                let file = self.file(label.source)?;
                Some(file.location(label.span.start()).line.to_string().len())
            })
            .max()
            .unwrap_or(1);
        for label in &diagnostic.labels {
            match self.file(label.source) {
                None => {
                    let _ = writeln!(out, "{:width$} = {}", "", label.message);
                }
                Some(file) => render_label(&mut out, file, label, width),
            }
        }
        if let Some(help) = &diagnostic.help {
            let _ = writeln!(out, "{:width$} = help: {help}", "");
        }
        out
    }

    pub fn render_all(&self, report: &Report) -> String {
        report
            .diagnostics()
            .iter()
            .map(|diagnostic| self.render(diagnostic))
            .collect::<Vec<String>>()
            .join("\n")
    }
}

fn render_label(out: &mut String, file: &SourceFile, label: &Label, width: usize) {
    let start = file.location(label.span.start());
    let text = file.line(start.line);
    let _ = writeln!(
        out,
        "{:width$}--> {}:{}:{}",
        "",
        file.name(),
        start.line,
        start.column
    );
    let _ = writeln!(out, "{:width$} |", "");
    let _ = writeln!(out, "{:width$} | {text}", start.line);
    let _ = writeln!(
        out,
        "{:width$} | {:column$}{} {}",
        "",
        "",
        "^".repeat(caret_width(file, label.span, start.line)),
        label.message,
        column = start.column - 1
    );
}

fn caret_width(file: &SourceFile, span: Span, line: usize) -> usize {
    let line_end = file.line_start(line) + file.line(line).len();
    let span_end = span.end().min(line_end);
    file.text()
        .get(span.start()..span_end)
        .map(|text| text.chars().count())
        .unwrap_or(1)
        .max(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sources() -> (Sources, SourceId) {
        let mut sources = Sources::new();
        let id = sources.add(
            "design.yaml",
            "entities:\n  - protein:\n      id: B\n      sequence: 80..140\n",
        );
        (sources, id)
    }

    #[test]
    fn a_rendered_diagnostic_points_at_its_span() {
        let (sources, id) = sources();
        let rendered = sources.render(
            &Diagnostic::error("unknown-chain", "no chain B in 8r3a.cif")
                .at(id, Span::new(33, 34), "not in the target")
                .help("available chains: A"),
        );
        assert_eq!(
            rendered,
            concat!(
                "error[unknown-chain]: no chain B in 8r3a.cif\n",
                " --> design.yaml:3:11\n",
                "  |\n",
                "3 |       id: B\n",
                "  |           ^ not in the target\n",
                "  = help: available chains: A\n",
            )
        );
    }

    #[test]
    fn a_caret_covers_the_whole_span_on_one_line() {
        let (sources, id) = sources();
        let rendered = sources.render(&Diagnostic::warning("length", "wide range").at(
            id,
            Span::new(51, 58),
            "80 to 140",
        ));
        assert!(rendered.starts_with("warning[length]: wide range\n"));
        assert!(rendered.contains("4 |       sequence: 80..140\n"));
        assert!(rendered.contains("  |                 ^^^^^^^ 80 to 140\n"));
    }

    #[test]
    fn a_span_crossing_lines_stops_at_the_first_line() {
        let (sources, id) = sources();
        let rendered = sources.render(&Diagnostic::error("entity", "bad entity").at(
            id,
            Span::new(14, 40),
            "here",
        ));
        let carets = rendered
            .lines()
            .find(|line| line.contains('^'))
            .expect("a caret line");
        assert_eq!(carets, "  |     ^^^^^^^^ here");
    }

    #[test]
    fn labels_without_a_known_file_still_print_their_message() {
        let (sources, _) = sources();
        let mut other = Sources::new();
        other.add("design.yaml", "");
        let unknown = other.add("other.cif", "");
        let rendered = sources.render(&Diagnostic::error("io", "unreadable").at(
            unknown,
            Span::new(0, 1),
            "no file",
        ));
        assert_eq!(rendered, "error[io]: unreadable\n  = no file\n");
    }

    #[test]
    fn a_report_separates_errors_from_warnings() {
        let mut report = Report::new();
        assert!(report.is_empty() && !report.has_errors());
        report.push(Diagnostic::warning("wide-range", "wide"));
        assert!(!report.has_errors());
        report.push(Diagnostic::error("unknown-chain", "missing"));
        assert!(report.has_errors());
        assert_eq!(report.len(), 2);
        assert_eq!(report.codes(), ["wide-range", "unknown-chain"]);
        assert_eq!(report.errors().count(), 1);
    }

    #[test]
    fn gutters_line_up_across_labels_with_wide_line_numbers() {
        let mut sources = Sources::new();
        let text = "a\n".repeat(12);
        let id = sources.add("design.yaml", text);
        let rendered = sources.render(
            &Diagnostic::error("dup", "twice")
                .at(id, Span::new(0, 1), "first")
                .at(id, Span::new(22, 23), "second"),
        );
        assert!(rendered.contains(" --> design.yaml:1:1\n"));
        assert!(rendered.contains(" --> design.yaml:12:1\n"));
        assert!(rendered.contains("12 | a\n"));
        assert!(rendered.contains(" 1 | a\n"));
    }
}

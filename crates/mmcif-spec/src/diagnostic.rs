use std::fmt::{self, Write as _};

use crate::span::{SourceFile, SourceId, Sources, Span};

/// How serious a [`Diagnostic`] is.
///
/// Warnings order before errors.
///
/// # Examples
///
/// ```
/// use mmcif_spec::{Diagnostic, Severity};
///
/// let diagnostic = Diagnostic::warning("unsupported-key", "`msa` is not checked");
///
/// assert_eq!(diagnostic.severity, Severity::Warning);
/// assert!(Severity::Warning < Severity::Error);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    /// A warning, which [`Report::has_errors`] does not count.
    Warning,
    /// An error, which [`Report::has_errors`] counts.
    Error,
}

impl Severity {
    /// Returns the lowercase name printed when a diagnostic is rendered.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::Severity;
    ///
    /// assert_eq!(Severity::Warning.as_str(), "warning");
    /// assert_eq!(Severity::Error.as_str(), "error");
    /// ```
    pub fn as_str(self) -> &'static str {
        match self {
            Severity::Warning => "warning",
            Severity::Error => "error",
        }
    }
}

/// A message attached to a span of a source file.
///
/// # Examples
///
/// ```
/// use mmcif_spec::{Diagnostic, Sources, Span};
///
/// let mut sources = Sources::new();
/// let id = sources.add("design.yaml", "id: B\n");
/// let diagnostic =
///     Diagnostic::error("unknown-chain", "no chain B").at(id, Span::new(4, 5), "not in the target");
///
/// let label = &diagnostic.labels[0];
/// assert_eq!(label.source, id);
/// assert_eq!(label.span, Span::new(4, 5));
/// assert_eq!(label.message, "not in the target");
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Label {
    /// The file the span points into.
    pub source: SourceId,
    /// The text the label points at.
    pub span: Span,
    /// The message printed after the carets.
    pub message: String,
}

/// A problem found in a specification, with labelled spans and an optional
/// help text.
///
/// # Examples
///
/// ```
/// use mmcif_spec::{Diagnostic, Sources, Span};
///
/// let mut sources = Sources::new();
/// let id = sources.add("design.yaml", "entities:\n  - protein:\n      id: B\n");
///
/// let diagnostic = Diagnostic::error("unknown-chain", "no chain B in target.cif")
///     .at(id, Span::new(33, 34), "not in the target")
///     .help("available chains: A");
///
/// assert_eq!(
///     sources.render(&diagnostic),
///     "\
/// error[unknown-chain]: no chain B in target.cif
///  --> design.yaml:3:11
///   |
/// 3 |       id: B
///   |           ^ not in the target
///   = help: available chains: A
/// "
/// );
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    /// Whether the diagnostic is a warning or an error.
    pub severity: Severity,
    /// A short name for the kind of problem, such as `unknown-chain`.
    pub code: &'static str,
    /// The description printed after the code.
    pub message: String,
    /// The spans the diagnostic points at, in the order they are printed.
    pub labels: Vec<Label>,
    /// A suggestion printed after the labels.
    pub help: Option<String>,
}

impl Diagnostic {
    /// Creates an error with no labels and no help text.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Diagnostic, Severity};
    ///
    /// let diagnostic = Diagnostic::error("unknown-chain", "chain Q is not declared");
    ///
    /// assert_eq!(diagnostic.severity, Severity::Error);
    /// assert_eq!(diagnostic.code, "unknown-chain");
    /// assert!(diagnostic.labels.is_empty());
    /// ```
    pub fn error(code: &'static str, message: impl Into<String>) -> Diagnostic {
        Diagnostic::new(Severity::Error, code, message)
    }

    /// Creates a warning with no labels and no help text.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Diagnostic, Severity};
    ///
    /// let diagnostic = Diagnostic::warning("unsupported-key", "`msa` is not checked");
    ///
    /// assert_eq!(diagnostic.severity, Severity::Warning);
    /// assert_eq!(diagnostic.help, None);
    /// ```
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

    /// Adds a label pointing at `span` of the file `source`, with `message`
    /// printed after the carets, and returns the diagnostic.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Diagnostic, Sources, Span};
    ///
    /// let mut sources = Sources::new();
    /// let id = sources.add("design.yaml", "atom2: [Q, 1, CK]\n");
    ///
    /// let diagnostic = Diagnostic::error("unknown-chain", "chain Q is not declared")
    ///     .at(id, Span::new(8, 9), "no entity declares this chain");
    ///
    /// assert_eq!(diagnostic.labels.len(), 1);
    /// assert_eq!(diagnostic.labels[0].span, Span::new(8, 9));
    /// ```
    #[must_use]
    pub fn at(mut self, source: SourceId, span: Span, message: impl Into<String>) -> Diagnostic {
        self.labels.push(Label {
            source,
            span,
            message: message.into(),
        });
        self
    }

    /// Sets the help text, replacing any earlier one, and returns the
    /// diagnostic.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::Diagnostic;
    ///
    /// let diagnostic = Diagnostic::error("unknown-chain", "chain Q is not declared")
    ///     .help("declared chains: S");
    ///
    /// assert_eq!(diagnostic.help.as_deref(), Some("declared chains: S"));
    /// ```
    #[must_use]
    pub fn help(mut self, help: impl Into<String>) -> Diagnostic {
        self.help = Some(help.into());
        self
    }

    /// Returns `true` if the diagnostic is an error.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::Diagnostic;
    ///
    /// assert!(Diagnostic::error("unknown-chain", "chain Q is not declared").is_error());
    /// assert!(!Diagnostic::warning("unsupported-key", "`msa` is not checked").is_error());
    /// ```
    pub fn is_error(&self) -> bool {
        self.severity == Severity::Error
    }
}

/// The diagnostics found in a specification, in the order they were added.
///
/// # Examples
///
/// ```
/// use mmcif_spec::{Diagnostic, Report};
///
/// let mut report = Report::new();
/// report.push(Diagnostic::warning("unsupported-key", "`msa` is not checked"));
/// report.push(Diagnostic::error("unknown-chain", "chain Q is not declared"));
///
/// assert_eq!(report.len(), 2);
/// assert_eq!(report.errors().count(), 1);
/// assert!(report.has_errors());
/// ```
#[derive(Debug, Clone, Default)]
pub struct Report {
    diagnostics: Vec<Diagnostic>,
}

impl Report {
    /// Creates an empty report.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::Report;
    ///
    /// let report = Report::new();
    ///
    /// assert!(report.is_empty());
    /// ```
    pub fn new() -> Report {
        Report::default()
    }

    /// Appends `diagnostic` to the report.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Diagnostic, Report};
    ///
    /// let mut report = Report::new();
    /// report.push(Diagnostic::error("unknown-chain", "chain Q is not declared"));
    ///
    /// assert_eq!(report.codes(), ["unknown-chain"]);
    /// ```
    pub fn push(&mut self, diagnostic: Diagnostic) {
        self.diagnostics.push(diagnostic);
    }

    /// Appends every diagnostic in `diagnostics` to the report, in order.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Diagnostic, Report};
    ///
    /// let mut report = Report::new();
    /// report.extend([
    ///     Diagnostic::error("unknown-chain", "chain Q is not declared"),
    ///     Diagnostic::error("unknown-residue", "chain S has no residue 99"),
    /// ]);
    ///
    /// assert_eq!(report.codes(), ["unknown-chain", "unknown-residue"]);
    /// ```
    pub fn extend(&mut self, diagnostics: impl IntoIterator<Item = Diagnostic>) {
        self.diagnostics.extend(diagnostics);
    }

    /// Returns every diagnostic in the order they were added.
    pub fn diagnostics(&self) -> &[Diagnostic] {
        &self.diagnostics
    }

    /// Returns an iterator over the errors in the report, in the order they
    /// were added.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Diagnostic, Report};
    ///
    /// let mut report = Report::new();
    /// report.push(Diagnostic::warning("unsupported-key", "`msa` is not checked"));
    /// report.push(Diagnostic::error("unknown-chain", "chain Q is not declared"));
    ///
    /// let codes: Vec<&str> = report.errors().map(|error| error.code).collect();
    /// assert_eq!(codes, ["unknown-chain"]);
    /// ```
    pub fn errors(&self) -> impl Iterator<Item = &Diagnostic> {
        self.diagnostics.iter().filter(|entry| entry.is_error())
    }

    /// Returns `true` if the report holds at least one error.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Diagnostic, Report};
    ///
    /// let mut report = Report::new();
    /// report.push(Diagnostic::warning("unsupported-key", "`msa` is not checked"));
    /// assert!(!report.has_errors());
    ///
    /// report.push(Diagnostic::error("unknown-chain", "chain Q is not declared"));
    /// assert!(report.has_errors());
    /// ```
    pub fn has_errors(&self) -> bool {
        self.errors().next().is_some()
    }

    /// Returns the code of every diagnostic, in the order they were added.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Diagnostic, Report};
    ///
    /// let mut report = Report::new();
    /// report.push(Diagnostic::warning("unsupported-key", "`msa` is not checked"));
    /// report.push(Diagnostic::error("unknown-chain", "chain Q is not declared"));
    ///
    /// assert_eq!(report.codes(), ["unsupported-key", "unknown-chain"]);
    /// ```
    pub fn codes(&self) -> Vec<&'static str> {
        self.diagnostics.iter().map(|entry| entry.code).collect()
    }

    /// Returns `true` if the report holds no diagnostics.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Diagnostic, Report};
    ///
    /// let mut report = Report::new();
    /// assert!(report.is_empty());
    ///
    /// report.push(Diagnostic::warning("unsupported-key", "`msa` is not checked"));
    /// assert!(!report.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.diagnostics.is_empty()
    }

    /// Returns the number of diagnostics in the report.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Diagnostic, Report};
    ///
    /// let mut report = Report::new();
    /// report.push(Diagnostic::warning("unsupported-key", "`msa` is not checked"));
    ///
    /// assert_eq!(report.len(), 1);
    /// ```
    pub fn len(&self) -> usize {
        self.diagnostics.len()
    }
}

impl Sources {
    /// Renders `diagnostic` as text.
    ///
    /// The first line holds the severity, code and message. Each label
    /// follows with its file name, line and column, the source line, and
    /// carets under the span, which stop at the end of the span's first line.
    /// A label whose file is not in this set prints only its message. The
    /// help text comes last. The result ends with a newline.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Diagnostic, Sources, Span};
    ///
    /// let mut sources = Sources::new();
    /// let yaml = "entities:\n  - protein:\n      id: S\n      sequence: 80..140\n";
    /// let id = sources.add("design.yaml", yaml);
    ///
    /// let diagnostic =
    ///     Diagnostic::warning("wide-range", "a wide range").at(id, Span::new(51, 58), "80 to 140");
    ///
    /// assert_eq!(
    ///     sources.render(&diagnostic),
    ///     "\
    /// warning[wide-range]: a wide range
    ///  --> design.yaml:4:17
    ///   |
    /// 4 |       sequence: 80..140
    ///   |                 ^^^^^^^ 80 to 140
    /// "
    /// );
    /// ```
    pub fn render(&self, diagnostic: &Diagnostic) -> String {
        let mut out = String::new();
        let _ = self.write(&mut out, diagnostic);
        out
    }

    fn write(&self, out: &mut String, diagnostic: &Diagnostic) -> fmt::Result {
        writeln!(
            out,
            "{}[{}]: {}",
            diagnostic.severity.as_str(),
            diagnostic.code,
            diagnostic.message
        )?;
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
                None => writeln!(out, "{:width$} = {}", "", label.message)?,
                Some(file) => write_label(out, file, label, width)?,
            }
        }
        if let Some(help) = &diagnostic.help {
            writeln!(out, "{:width$} = help: {help}", "")?;
        }
        Ok(())
    }

    /// Renders every diagnostic in `report` with [`Sources::render`],
    /// separated by blank lines.
    ///
    /// Returns an empty string for an empty report.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Diagnostic, Report, Sources};
    ///
    /// let sources = Sources::new();
    /// let mut report = Report::new();
    /// report.push(Diagnostic::error("unknown-chain", "chain Q is not declared"));
    /// report.push(Diagnostic::warning("unsupported-key", "`msa` is not checked"));
    ///
    /// assert_eq!(
    ///     sources.render_all(&report),
    ///     "\
    /// error[unknown-chain]: chain Q is not declared
    ///
    /// warning[unsupported-key]: `msa` is not checked
    /// "
    /// );
    /// assert_eq!(sources.render_all(&Report::new()), "");
    /// ```
    pub fn render_all(&self, report: &Report) -> String {
        report
            .diagnostics()
            .iter()
            .map(|diagnostic| self.render(diagnostic))
            .collect::<Vec<String>>()
            .join("\n")
    }
}

fn write_label(out: &mut String, file: &SourceFile, label: &Label, width: usize) -> fmt::Result {
    let start = file.location(label.span.start());
    let text = file.line(start.line);
    writeln!(
        out,
        "{:width$}--> {}:{}:{}",
        "",
        file.name(),
        start.line,
        start.column
    )?;
    writeln!(out, "{:width$} |", "")?;
    writeln!(out, "{:width$} | {text}", start.line)?;
    writeln!(
        out,
        "{:width$} | {:column$}{} {}",
        "",
        "",
        "^".repeat(caret_width(file, label.span, start.line)),
        label.message,
        column = start.column - 1
    )
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

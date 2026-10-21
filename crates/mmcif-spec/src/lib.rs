//! Parsing and checking of design specifications. Every type here may change until 1.0.

mod check;
mod design;
mod diagnostic;
mod document;
mod geometry;
mod lower;
mod read;
mod resolve;
mod sequence;
mod shape;
mod span;
mod yaml;

use std::path::{Path, PathBuf};

use mmcif_ir::Instantiated;

pub use design::{Design, ElementSource, SegmentKind, SegmentSource, TargetPath};
pub use diagnostic::{Diagnostic, Label, Report, Severity};
pub use lower::Addressing;
pub use sequence::{component, one_letter};
pub use span::{Location, SourceFile, SourceId, Sources, Span};

#[derive(Debug, thiserror::Error)]
pub enum SpecError {
    #[error("{path} could not be read: {cause}")]
    Unreadable {
        path: PathBuf,
        cause: std::io::Error,
    },
}

#[derive(Debug)]
pub struct Checked {
    pub sources: Sources,
    pub report: Report,
    pub source: SourceId,
    pub design: Option<Design>,
    pub shortest: Option<Instantiated>,
}

impl Checked {
    pub fn render(&self) -> String {
        self.sources.render_all(&self.report)
    }

    pub fn has_errors(&self) -> bool {
        self.report.has_errors()
    }
}

pub fn check_path(path: impl AsRef<Path>) -> Result<Checked, SpecError> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path).map_err(|cause| SpecError::Unreadable {
        path: path.to_path_buf(),
        cause,
    })?;
    let base = path.parent().unwrap_or(Path::new("."));
    Ok(check_text(resolve::display(path), &text, base))
}

pub fn check_text(name: impl Into<String>, text: &str, base: &Path) -> Checked {
    let mut sources = Sources::new();
    let source = sources.add(name, text);
    let mut report = Report::new();
    let root = match yaml::parse(source, text) {
        Ok(root) => root,
        Err(diagnostic) => {
            report.push(diagnostic);
            return Checked {
                sources,
                report,
                source,
                design: None,
                shortest: None,
            };
        }
    };
    let document = read::document(source, &root, &mut report);
    let lowered = lower::lower(source, base, &document, &mut report);
    let shortest = check::check(source, &lowered, &mut report);
    let mut design = lowered.design();
    design.note_skipped(
        report
            .diagnostics()
            .iter()
            .filter(|diagnostic| diagnostic.code.starts_with("unsupported-"))
            .cloned(),
    );
    Checked {
        sources,
        report,
        source,
        design: Some(design),
        shortest,
    }
}

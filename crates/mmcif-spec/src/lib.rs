//! Parsing and checking of design specifications.
//!
//! [`check_text`] and [`check_path`] read a BoltzGen design specification,
//! resolve its residue references against the structure files it names, and
//! return a [`Checked`] holding the lowered [`Design`] and a [`Report`] of
//! diagnostics that render with carets under the specification text. Every
//! type here may change until 1.0.
//!
//! # Examples
//!
//! ```
//! use mmcif_spec::check_text;
//! # let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
//!
//! let yaml = "\
//! entities:
//!   - file:
//!       path: target.cif
//!       include:
//!         - chain:
//!             id: C
//!   - protein:
//!       id: S
//!       sequence: 80..140
//! ";
//! let checked = check_text("design.yaml", yaml, &base);
//!
//! assert!(checked.has_errors());
//! assert_eq!(
//!     checked.render(),
//!     "\
//! error[unknown-chain]: target.cif has no chain C
//!  --> design.yaml:6:17
//!   |
//! 6 |             id: C
//!   |                 ^ not in the target
//!   = help: polymer chains in the file: A, B
//! "
//! );
//! ```

#![warn(missing_docs)]

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

/// An error from reading a specification file in [`check_path`].
///
/// # Examples
///
/// ```
/// use mmcif_spec::{SpecError, check_path};
///
/// let error = check_path("missing/design.yaml").unwrap_err();
///
/// assert!(matches!(error, SpecError::Unreadable { .. }));
/// ```
#[derive(Debug, thiserror::Error)]
pub enum SpecError {
    /// The specification file could not be read as UTF-8 text.
    #[error("{path} could not be read: {cause}")]
    Unreadable {
        /// The path passed to [`check_path`].
        path: PathBuf,
        /// The error from reading the file.
        cause: std::io::Error,
    },
}

/// The result of checking a specification.
///
/// # Examples
///
/// ```
/// use std::path::Path;
///
/// use mmcif_spec::check_text;
///
/// let yaml = "entities:\n  - protein:\n      id: S\n      sequence: 10..14C6C3\n";
/// let checked = check_text("design.yaml", yaml, Path::new("."));
///
/// assert!(checked.report.is_empty());
/// assert_eq!(checked.design.unwrap().segments.len(), 1);
/// assert_eq!(checked.shortest.unwrap().slot_count(), 21);
/// ```
#[derive(Debug)]
pub struct Checked {
    /// The source files the diagnostics point into, holding the
    /// specification text.
    pub sources: Sources,
    /// Every diagnostic found, in the order it was found.
    pub report: Report,
    /// The id of the specification text in `sources`.
    pub source: SourceId,
    /// The lowered design, or `None` if the text could not be parsed as YAML.
    ///
    /// The design is present even when `report` holds errors.
    pub design: Option<Design>,
    /// The design instantiated with every element at its minimum length, or
    /// `None` if there is no design or its lengths cannot be satisfied.
    pub shortest: Option<Instantiated>,
}

impl Checked {
    /// Renders every diagnostic in the report with [`Sources::render_all`].
    ///
    /// # Examples
    ///
    /// ```
    /// use std::path::Path;
    ///
    /// use mmcif_spec::check_text;
    ///
    /// let yaml = "\
    /// entities:
    ///   - protein:
    ///       id: S
    ///       sequence: 17
    /// constraints:
    ///   - bond:
    ///       atom1: [S, 1, SG]
    ///       atom2: [Q, 1, CK]
    /// ";
    /// let checked = check_text("design.yaml", yaml, Path::new("."));
    ///
    /// assert_eq!(
    ///     checked.render(),
    ///     "\
    /// error[unknown-chain]: chain Q is not declared
    ///  --> design.yaml:8:15
    ///   |
    /// 8 |       atom2: [Q, 1, CK]
    ///   |               ^ no entity declares this chain
    ///   = help: declared chains: S
    /// "
    /// );
    /// ```
    pub fn render(&self) -> String {
        self.sources.render_all(&self.report)
    }

    /// Returns `true` if the report holds at least one error.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::path::Path;
    ///
    /// use mmcif_spec::check_text;
    ///
    /// let yaml = "\
    /// entities:
    ///   - protein:
    ///       id: S
    ///       sequence: 20
    ///   - dna:
    ///       id: D
    ///       sequence: ATGC
    /// ";
    /// let checked = check_text("design.yaml", yaml, Path::new("."));
    ///
    /// assert_eq!(checked.report.codes(), ["unsupported-entity"]);
    /// assert!(!checked.has_errors());
    /// ```
    pub fn has_errors(&self) -> bool {
        self.report.has_errors()
    }
}

/// Reads and checks the specification file at `path`.
///
/// Structure files named in the specification are read relative to the
/// directory holding `path`. Diagnostics name the specification by its file
/// name.
///
/// # Errors
///
/// Returns [`SpecError::Unreadable`] if the file cannot be read as UTF-8
/// text. Problems in the specification itself are reported in
/// [`Checked::report`].
///
/// # Examples
///
/// ```
/// use mmcif_spec::check_path;
///
/// let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/peptide.yaml");
/// let checked = check_path(&path)?;
///
/// assert!(!checked.has_errors());
/// assert_eq!(checked.sources.file(checked.source).unwrap().name(), "peptide.yaml");
/// # Ok::<(), mmcif_spec::SpecError>(())
/// ```
pub fn check_path(path: impl AsRef<Path>) -> Result<Checked, SpecError> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path).map_err(|cause| SpecError::Unreadable {
        path: path.to_path_buf(),
        cause,
    })?;
    let base = path.parent().unwrap_or(Path::new("."));
    Ok(check_text(resolve::display(path), &text, base))
}

/// Checks the specification `text`, naming it `name` in diagnostics.
///
/// Structure files named in the specification are read relative to `base`.
/// Every problem found, including a structure file that cannot be read, is
/// reported in [`Checked::report`].
///
/// # Examples
///
/// ```
/// use std::path::Path;
///
/// use mmcif_spec::check_text;
///
/// let yaml = "entities:\n  - protein:\n      id: S\n      sequence: 10..14C6C3\n";
/// let checked = check_text("design.yaml", yaml, Path::new("."));
/// assert!(checked.report.is_empty());
///
/// let checked = check_text("design.yaml", "entities: [", Path::new("."));
/// assert_eq!(checked.report.codes(), ["yaml-syntax"]);
/// assert!(checked.design.is_none());
/// ```
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

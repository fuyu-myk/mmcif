//! Backend input written from a design.
//!
//! [`lower`] writes a checked [`Design`] as a BoltzGen design specification
//! in YAML or as an RFdiffusion3 input in JSON, and returns a [`Refusal`]
//! listing every part of the design the [`Backend`] cannot express. The
//! output format tracks each backend and may change until 1.0.
//!
//! # Examples
//!
//! ```
//! use mmcif_lower::{Backend, lower};
//! use mmcif_spec::check_text;
//! # let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
//!
//! let yaml = "\
//! entities:
//!   - file:
//!       path: target.cif
//!       include:
//!         - chain:
//!             id: A
//!       design:
//!         - chain:
//!             id: A
//!             res_index: 4..5
//! ";
//! let checked = check_text("design.yaml", yaml, &base);
//! let design = checked.design.as_ref().unwrap();
//!
//! assert_eq!(
//!     lower(Backend::Rfd3, checked.source, design)?,
//!     r#"{
//!   "design": {
//!     "contig": "A1-3,2,A6",
//!     "input": "target.cif"
//!   }
//! }
//! "#
//! );
//! # Ok::<(), mmcif_lower::Refusal>(())
//! ```

#![warn(missing_docs)]

mod boltzgen;
mod capability;
mod ranges;
mod rfd3;
mod slots;
mod text;

use mmcif_spec::{Design, Diagnostic, SourceId};

pub use capability::Capabilities;

const ALL: [Backend; 2] = [Backend::BoltzGen, Backend::Rfd3];

/// A protein design tool that [`lower`] writes input for.
///
/// # Examples
///
/// ```
/// use mmcif_lower::Backend;
///
/// let backend = Backend::named("rfd3").unwrap();
///
/// assert_eq!(backend, Backend::Rfd3);
/// assert!(!backend.capabilities().bonds);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Backend {
    /// BoltzGen, which takes a design specification in YAML.
    BoltzGen,
    /// RFdiffusion3, which takes a design input in JSON.
    Rfd3,
}

/// The problems that stop [`lower`] from writing a design for a backend.
///
/// # Examples
///
/// ```
/// use mmcif_lower::{Backend, lower};
/// use mmcif_spec::{Report, check_text};
/// # let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
///
/// let yaml = "\
/// entities:
///   - file:
///       path: target.cif
///   - file:
///       path: other.cif
/// ";
/// let checked = check_text("design.yaml", yaml, &base);
/// let design = checked.design.as_ref().unwrap();
///
/// let refusal = lower(Backend::Rfd3, checked.source, design).unwrap_err();
/// assert_eq!(refusal.backend, Backend::Rfd3);
/// assert_eq!(refusal.to_string(), "rfd3 cannot express this specification");
///
/// let mut report = Report::new();
/// report.extend(refusal.problems().iter().cloned());
/// assert_eq!(
///     checked.sources.render_all(&report),
///     "\
/// error[unsupported-by-backend]: rfd3 takes one input structure, not 2
///  --> design.yaml:5:13
///   |
/// 5 |       path: other.cif
///   |             ^^^^^^^^^ this file would be ignored
///   = help: merge the chains into one structure file
/// "
/// );
/// ```
#[derive(Debug, thiserror::Error)]
#[error("{} cannot express this specification", backend.name())]
pub struct Refusal {
    /// The backend that refused the design.
    pub backend: Backend,
    problems: Vec<Diagnostic>,
}

impl Refusal {
    /// Returns an error for each problem, with labels pointing into the
    /// specification text.
    pub fn problems(&self) -> &[Diagnostic] {
        &self.problems
    }
}

impl Backend {
    /// Returns the backend called `name`, or `None` if no backend has that
    /// name.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_lower::Backend;
    ///
    /// assert_eq!(Backend::named("boltzgen"), Some(Backend::BoltzGen));
    /// assert_eq!(Backend::named("rfd3"), Some(Backend::Rfd3));
    /// assert_eq!(Backend::named("BoltzGen"), None);
    /// ```
    pub fn named(name: &str) -> Option<Backend> {
        match name {
            "boltzgen" => Some(Backend::BoltzGen),
            "rfd3" => Some(Backend::Rfd3),
            _ => None,
        }
    }

    /// Returns the name of the backend, as accepted by [`Backend::named`].
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_lower::Backend;
    ///
    /// assert_eq!(Backend::BoltzGen.name(), "boltzgen");
    /// assert_eq!(Backend::named(Backend::Rfd3.name()), Some(Backend::Rfd3));
    /// ```
    pub fn name(self) -> &'static str {
        self.capabilities().name
    }

    /// Returns what the backend can express.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_lower::Backend;
    ///
    /// assert!(Backend::BoltzGen.capabilities().ligands);
    /// assert!(!Backend::Rfd3.capabilities().ligands);
    /// ```
    pub fn capabilities(self) -> &'static Capabilities {
        match self {
            Backend::BoltzGen => &boltzgen::CAPABILITIES,
            Backend::Rfd3 => &rfd3::CAPABILITIES,
        }
    }

    /// Returns the name of every backend, in the order of [`Backend::all`].
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_lower::Backend;
    ///
    /// assert_eq!(Backend::names(), ["boltzgen", "rfd3"]);
    /// ```
    pub fn names() -> Vec<&'static str> {
        ALL.iter().map(|backend| backend.name()).collect()
    }

    /// Returns every backend.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_lower::Backend;
    ///
    /// assert_eq!(Backend::all(), [Backend::BoltzGen, Backend::Rfd3]);
    /// ```
    pub fn all() -> &'static [Backend] {
        &ALL
    }

    /// Returns the version of the backend's input format that the output is
    /// written against.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_lower::Backend;
    ///
    /// for backend in Backend::all() {
    ///     println!("{}: {}", backend.name(), backend.schema());
    /// }
    /// ```
    pub fn schema(self) -> &'static str {
        self.capabilities().schema
    }
}

/// Writes `design` as input for `backend`.
///
/// [`Backend::BoltzGen`] output is a BoltzGen design specification in YAML
/// and [`Backend::Rfd3`] output is an RFdiffusion3 input in JSON. Residues
/// taken from structure files are written with their `label_seq_id` numbers.
/// `source` is the id of the specification text that the diagnostics of a
/// [`Refusal`] point into. Only `design` is read: errors in the report it
/// was checked with do not cause a refusal.
///
/// # Errors
///
/// Returns a [`Refusal`] if `backend` cannot express the design. It holds an
/// error for each problem found, with one of these codes:
///
/// - the code of an entry of [`Design::skipped`], for each part of the
///   specification that was parsed but is not represented in the design;
/// - `unsupported-by-backend`, for a ligand, a residue without fixed
///   coordinates that has a fixed or restricted identity, or a constraint
///   that the backend's [`Capabilities`] exclude, and for more distinct
///   structure files than it takes;
/// - `cannot-lower`, for a part with no form in the backend's input, such as
///   a gap in the residues of a chain taken from a file for
///   [`Backend::Rfd3`].
///
/// # Examples
///
/// ```
/// use std::path::Path;
///
/// use mmcif_lower::{Backend, lower};
/// use mmcif_spec::check_text;
///
/// let yaml = "\
/// entities:
///   - protein:
///       id: R
///       sequence: 3..5C6C3
///   - ligand:
///       id: Q
///       ccd: WHL
/// constraints:
///   - bond:
///       atom1: [R, 4, SG]
///       atom2: [Q, 1, CK]
/// ";
/// let checked = check_text("design.yaml", yaml, Path::new("."));
/// let design = checked.design.as_ref().unwrap();
///
/// assert_eq!(lower(Backend::BoltzGen, checked.source, design)?, yaml);
///
/// let refusal = lower(Backend::Rfd3, checked.source, design).unwrap_err();
/// let problems = refusal.problems();
/// assert!(problems.iter().all(|problem| problem.code == "unsupported-by-backend"));
/// assert!(problems.iter().any(|problem| problem.message == "rfd3 cannot take a ligand"));
/// # Ok::<(), mmcif_lower::Refusal>(())
/// ```
pub fn lower(backend: Backend, source: SourceId, design: &Design) -> Result<String, Refusal> {
    let emitted = match backend {
        Backend::BoltzGen => boltzgen::emit(source, design),
        Backend::Rfd3 => rfd3::emit(source, design),
    };
    emitted.map_err(|problems| Refusal { backend, problems })
}

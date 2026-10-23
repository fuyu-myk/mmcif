use std::path::PathBuf;

use mmcif_ir::{ElementIndex, SegmentId, Spec, SymbolicAddress};

use crate::diagnostic::Diagnostic;
use crate::lower::{Addressing, Lowered};
use crate::span::Span;

/// The kind of entity that declared a chain of a [`Design`].
///
/// # Examples
///
/// ```
/// use std::path::Path;
///
/// use mmcif_spec::{SegmentKind, check_text};
///
/// let yaml = "\
/// entities:
///   - protein:
///       id: S
///       sequence: 20
///   - ligand:
///       id: Q
///       ccd: WHL
/// ";
/// let checked = check_text("design.yaml", yaml, Path::new("."));
/// let design = checked.design.as_ref().unwrap();
///
/// let kinds: Vec<SegmentKind> = design.segments.iter().map(|segment| segment.kind).collect();
/// assert_eq!(kinds, [SegmentKind::Protein, SegmentKind::Ligand]);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SegmentKind {
    /// A designed chain declared by a `protein` entity.
    Protein,
    /// A ligand declared by a `ligand` entity.
    Ligand,
    /// A chain taken from a structure file by a `file` entity.
    File,
}

impl SegmentKind {
    /// Returns how residue indices in a specification address a chain of
    /// this kind.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::{Addressing, SegmentKind};
    ///
    /// assert_eq!(SegmentKind::Protein.addressing(), Addressing::Positional);
    /// assert_eq!(SegmentKind::Ligand.addressing(), Addressing::Positional);
    /// assert_eq!(SegmentKind::File.addressing(), Addressing::Label);
    /// ```
    pub fn addressing(self) -> Addressing {
        match self {
            SegmentKind::File => Addressing::Label,
            SegmentKind::Protein | SegmentKind::Ligand => Addressing::Positional,
        }
    }
}

/// The source of one element of a chain of a [`Design`].
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
/// let design = checked.design.as_ref().unwrap();
///
/// let cysteine = &design.segments[0].elements[1];
/// assert_eq!(&yaml[cysteine.span.start()..cysteine.span.end()], "C");
/// assert_eq!(cysteine.label, None);
/// assert_eq!(cysteine.target, None);
/// ```
#[derive(Debug, Clone)]
pub struct ElementSource {
    /// The span of the specification text that declared the element.
    pub span: Span,
    /// The `label_seq_id` of the residue in its structure file, for chains
    /// taken from a file.
    pub label: Option<u32>,
    /// The position in [`Design::targets`] of the structure file holding the
    /// residue, for chains taken from a file.
    pub target: Option<usize>,
}

/// A chain of a [`Design`] and the source of each of its elements.
///
/// `elements` holds one entry per element of the chain's segment in
/// [`Design::spec`], in the same order.
///
/// # Examples
///
/// ```
/// use mmcif_spec::{SegmentKind, check_text};
/// # let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
///
/// let yaml = "\
/// entities:
///   - file:
///       path: target.cif
///       include:
///         - chain:
///             id: A
///             res_index: 2..4
/// ";
/// let checked = check_text("design.yaml", yaml, &base);
/// let chain = &checked.design.as_ref().unwrap().segments[0];
///
/// assert_eq!(chain.id.as_str(), "A");
/// assert_eq!(chain.kind, SegmentKind::File);
/// assert_eq!(chain.labels().collect::<Vec<u32>>(), [2, 3, 4]);
/// ```
#[derive(Debug, Clone)]
pub struct SegmentSource {
    /// The chain id, which is also the id of its segment in [`Design::spec`].
    pub id: SegmentId,
    /// The kind of entity that declared the chain.
    pub kind: SegmentKind,
    /// The span of the entity or chain selection that declared the chain.
    pub span: Span,
    /// The source of each element of the chain, in order.
    pub elements: Vec<ElementSource>,
}

impl SegmentSource {
    /// Returns the position in [`Design::targets`] of the structure file the
    /// chain is taken from, or `None` if it is not taken from a file.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::check_text;
    /// # let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    ///
    /// let yaml = "\
    /// entities:
    ///   - protein:
    ///       id: S
    ///       sequence: 20
    ///   - file:
    ///       path: target.cif
    ///       include:
    ///         - chain:
    ///             id: A
    /// ";
    /// let checked = check_text("design.yaml", yaml, &base);
    /// let design = checked.design.as_ref().unwrap();
    ///
    /// assert_eq!(design.segments[0].target(), None);
    /// assert_eq!(design.segments[1].target(), Some(0));
    /// ```
    pub fn target(&self) -> Option<usize> {
        self.elements.iter().find_map(|element| element.target)
    }

    /// Returns the `label_seq_id` of each residue taken from a structure
    /// file, in chain order.
    ///
    /// The iterator is empty for chains not taken from a file.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::check_text;
    /// # let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    ///
    /// let yaml = "\
    /// entities:
    ///   - file:
    ///       path: target.cif
    ///       include:
    ///         - chain:
    ///             id: A
    ///             res_index: 1,4..5
    /// ";
    /// let checked = check_text("design.yaml", yaml, &base);
    /// let chain = &checked.design.as_ref().unwrap().segments[0];
    ///
    /// assert_eq!(chain.labels().collect::<Vec<u32>>(), [1, 4, 5]);
    /// ```
    pub fn labels(&self) -> impl Iterator<Item = u32> + use<'_> {
        self.elements.iter().filter_map(|element| element.label)
    }

    /// Returns the source of the element at `index`, or `None` if the chain
    /// has no such element.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::path::Path;
    ///
    /// use mmcif_ir::ElementIndex;
    /// use mmcif_spec::check_text;
    ///
    /// let yaml = "entities:\n  - protein:\n      id: S\n      sequence: 10..14C6C3\n";
    /// let checked = check_text("design.yaml", yaml, Path::new("."));
    /// let chain = &checked.design.as_ref().unwrap().segments[0];
    ///
    /// let run = chain.element(ElementIndex::new(0)).unwrap();
    /// assert_eq!(&yaml[run.span.start()..run.span.end()], "10..14");
    /// assert!(chain.element(ElementIndex::new(5)).is_none());
    /// ```
    pub fn element(&self, index: ElementIndex) -> Option<&ElementSource> {
        self.elements.get(index.as_usize())
    }
}

/// A structure file named by a `file` entity.
///
/// # Examples
///
/// ```
/// use mmcif_spec::check_text;
/// # let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
///
/// let yaml = "entities:\n  - file:\n      path: target.cif\n";
/// let checked = check_text("design.yaml", yaml, &base);
/// let target = &checked.design.as_ref().unwrap().targets[0];
///
/// assert_eq!(target.declared, "target.cif");
/// assert_eq!(target.resolved, base.join("target.cif"));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TargetPath {
    /// The path as written in the specification.
    pub declared: String,
    /// The path the file was read from, which is `declared` joined onto the
    /// base directory of the check.
    pub resolved: PathBuf,
}

/// A specification lowered to a [`Spec`], with the source of every chain.
///
/// # Examples
///
/// ```
/// use mmcif_spec::{SegmentKind, check_text};
/// # let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
///
/// let yaml = "\
/// entities:
///   - file:
///       path: target.cif
///       include:
///         - chain:
///             id: A
///   - protein:
///       id: S
///       sequence: 10..14C6C3
/// constraints:
///   - bond:
///       atom1: [S, 11, SG]
///       atom2: [S, 18, SG]
/// ";
/// let checked = check_text("design.yaml", yaml, &base);
/// let design = checked.design.as_ref().unwrap();
///
/// assert_eq!(design.spec.segments().len(), 2);
/// assert_eq!(design.spec.constraints().len(), 1);
/// assert_eq!(design.segments[0].kind, SegmentKind::File);
/// assert_eq!(design.targets[0].declared, "target.cif");
/// ```
#[derive(Debug, Clone)]
pub struct Design {
    /// The lowered specification.
    pub spec: Spec,
    /// The chains of the design, in the order they were declared.
    pub segments: Vec<SegmentSource>,
    /// The structure file of each `file` entity whose file could be read, in
    /// the order they were declared.
    pub targets: Vec<TargetPath>,
    constraint_spans: Vec<Span>,
    skipped: Vec<Diagnostic>,
}

impl Design {
    /// Returns the warnings for parts of the specification that were parsed
    /// but are not represented in the design.
    ///
    /// These are the diagnostics whose code starts with `unsupported-`.
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
    /// let design = checked.design.as_ref().unwrap();
    ///
    /// assert_eq!(design.skipped()[0].code, "unsupported-entity");
    /// ```
    pub fn skipped(&self) -> &[Diagnostic] {
        &self.skipped
    }

    pub(crate) fn note_skipped(&mut self, diagnostics: impl IntoIterator<Item = Diagnostic>) {
        self.skipped.extend(diagnostics);
    }

    /// Returns the structure file `segment` is taken from, or `None` if it
    /// is not taken from a file.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_spec::check_text;
    /// # let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    ///
    /// let yaml = "\
    /// entities:
    ///   - file:
    ///       path: target.cif
    ///       include:
    ///         - chain:
    ///             id: A
    ///   - protein:
    ///       id: S
    ///       sequence: 20
    /// ";
    /// let checked = check_text("design.yaml", yaml, &base);
    /// let design = checked.design.as_ref().unwrap();
    ///
    /// let target = design.target_of(&design.segments[0]).unwrap();
    /// assert_eq!(target.declared, "target.cif");
    /// assert!(design.target_of(&design.segments[1]).is_none());
    /// ```
    pub fn target_of(&self, segment: &SegmentSource) -> Option<&TargetPath> {
        self.targets.get(segment.target()?)
    }

    /// Returns the chain with id `id`, or `None` if the design has no such
    /// chain.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::path::Path;
    ///
    /// use mmcif_ir::SegmentId;
    /// use mmcif_spec::{SegmentKind, check_text};
    ///
    /// let yaml = "entities:\n  - protein:\n      id: S\n      sequence: 20\n";
    /// let checked = check_text("design.yaml", yaml, Path::new("."));
    /// let design = checked.design.as_ref().unwrap();
    ///
    /// assert_eq!(design.segment(&SegmentId::new("S")).unwrap().kind, SegmentKind::Protein);
    /// assert!(design.segment(&SegmentId::new("Q")).is_none());
    /// ```
    pub fn segment(&self, id: &SegmentId) -> Option<&SegmentSource> {
        self.segments.iter().find(|segment| &segment.id == id)
    }

    /// Returns the span of the specification text that declared the element
    /// `address` points into.
    ///
    /// Returns the span of the whole chain if the chain has no such element,
    /// and `None` if the design has no such chain.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::path::Path;
    ///
    /// use mmcif_ir::{ElementIndex, SegmentId, SymbolicAddress};
    /// use mmcif_spec::check_text;
    ///
    /// let yaml = "entities:\n  - protein:\n      id: S\n      sequence: 10..14C6C3\n";
    /// let checked = check_text("design.yaml", yaml, Path::new("."));
    /// let design = checked.design.as_ref().unwrap();
    ///
    /// let address = SymbolicAddress::at(SegmentId::new("S"), ElementIndex::new(2), 0);
    /// let span = design.span_of(&address).unwrap();
    /// assert_eq!(&yaml[span.start()..span.end()], "6");
    /// ```
    pub fn span_of(&self, address: &SymbolicAddress) -> Option<Span> {
        let segment = self.segment(&address.segment)?;
        segment
            .element(address.element)
            .map(|element| element.span)
            .or(Some(segment.span))
    }

    /// Returns the span of the specification text that declared the
    /// constraint at `index` in [`Spec::constraints`], or `None` if there is
    /// no such constraint.
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
    ///       sequence: 10..14C6C3
    /// constraints:
    ///   - bond:
    ///       atom1: [S, 11, SG]
    ///       atom2: [S, 18, SG]
    /// ";
    /// let checked = check_text("design.yaml", yaml, Path::new("."));
    /// let design = checked.design.as_ref().unwrap();
    ///
    /// let span = design.constraint_span(0).unwrap();
    /// assert!(yaml[span.start()..span.end()].starts_with("bond:"));
    /// assert!(design.constraint_span(1).is_none());
    /// ```
    pub fn constraint_span(&self, index: usize) -> Option<Span> {
        self.constraint_spans.get(index).copied()
    }

    /// Returns the residue index a specification would use for `address`.
    ///
    /// For a chain taken from a file this is the `label_seq_id` of the
    /// residue. For a designed chain or a ligand it is the position of the
    /// residue counted from 1, with every run of designed residues before it
    /// at its minimum length. Returns `None` if the design has no such chain
    /// or element, or if the index does not fit in a `u32`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{ElementIndex, SegmentId, SymbolicAddress};
    /// use mmcif_spec::check_text;
    /// # let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    ///
    /// let yaml = "\
    /// entities:
    ///   - file:
    ///       path: target.cif
    ///       include:
    ///         - chain:
    ///             id: A
    ///             res_index: 4..6
    ///   - protein:
    ///       id: S
    ///       sequence: 10..14C6C3
    /// ";
    /// let checked = check_text("design.yaml", yaml, &base);
    /// let design = checked.design.as_ref().unwrap();
    /// let at = |chain: &str, element: u32, offset: u32| {
    ///     let element = ElementIndex::new(element);
    ///     design.index_of(&SymbolicAddress::at(SegmentId::new(chain), element, offset))
    /// };
    ///
    /// assert_eq!(at("A", 0, 0), Some(4));
    /// assert_eq!(at("S", 1, 0), Some(11));
    /// assert_eq!(at("S", 2, 3), Some(15));
    /// assert_eq!(at("Q", 0, 0), None);
    /// ```
    pub fn index_of(&self, address: &SymbolicAddress) -> Option<u32> {
        let segment = self.segment(&address.segment)?;
        match segment.kind.addressing() {
            Addressing::Label => segment.element(address.element)?.label,
            Addressing::Positional => {
                let definition = self.spec.segment(&address.segment).ok()?;
                let before: u64 = definition
                    .elements()
                    .iter()
                    .take(address.element.as_usize())
                    .map(|element| u64::from(element.length().min()))
                    .sum();
                let index = before + u64::from(address.offset.unwrap_or(0)) + 1;
                u32::try_from(index).ok()
            }
        }
    }
}

impl Lowered {
    pub fn design(self) -> Design {
        Design {
            spec: self.spec,
            segments: self.segments,
            targets: self
                .targets
                .into_iter()
                .map(|target| TargetPath {
                    declared: target.declared,
                    resolved: target.path,
                })
                .collect(),
            constraint_spans: self.constraint_spans,
            skipped: Vec::new(),
        }
    }
}

/// ```compile_fail,E0624
/// use std::path::Path;
///
/// use mmcif_spec::{Diagnostic, check_text};
///
/// let yaml = "entities:\n  - protein:\n      id: S\n      sequence: 20\n";
/// let mut design = check_text("design.yaml", yaml, Path::new(".")).design.unwrap();
///
/// design.note_skipped([Diagnostic::warning("unsupported-key", "`msa` is not checked")]);
/// ```
#[cfg(doctest)]
struct NoteSkippedIsPrivate;

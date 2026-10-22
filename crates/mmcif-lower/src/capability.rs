use mmcif_ir::{Constraint, Geometry, Identity, Slot, Spec};
use mmcif_spec::{Design, Diagnostic, SegmentKind, Severity, SourceId, Span};

/// What a [`Backend`](crate::Backend) can express.
///
/// # Examples
///
/// ```
/// use mmcif_lower::Backend;
///
/// let rfd3 = Backend::Rfd3.capabilities();
///
/// assert_eq!(rfd3.name, "rfd3");
/// assert!(!rfd3.secondary_structure);
/// assert_eq!(rfd3.inputs, Some(1));
/// ```
#[non_exhaustive]
pub struct Capabilities {
    /// The name of the backend, as accepted by
    /// [`Backend::named`](crate::Backend::named).
    pub name: &'static str,
    /// The version of the backend's input format that the output is written
    /// against.
    pub schema: &'static str,
    /// Whether a residue without fixed coordinates can have a fixed or
    /// restricted identity.
    pub designed_identity: bool,
    /// Whether bond constraints can be expressed.
    pub bonds: bool,
    /// Whether distance constraints can be expressed.
    pub distances: bool,
    /// Whether secondary structure constraints can be expressed.
    pub secondary_structure: bool,
    /// Whether ligands can be expressed.
    pub ligands: bool,
    /// The largest number of distinct structure files the input can name, or
    /// `None` for no limit.
    pub inputs: Option<usize>,
}

pub fn unsupported(
    source: SourceId,
    design: &Design,
    capabilities: &Capabilities,
) -> Vec<Diagnostic> {
    let name = capabilities.name;
    let mut gaps: Vec<Diagnostic> = design
        .skipped()
        .iter()
        .map(|skipped| Diagnostic {
            severity: Severity::Error,
            help: Some(format!("{name} would be given a design without it")),
            ..skipped.clone()
        })
        .collect();
    for segment in &design.segments {
        let Ok(definition) = design.spec.segment(&segment.id) else {
            continue;
        };
        if segment.kind == SegmentKind::Ligand && !capabilities.ligands {
            gaps.push(
                refusal(
                    source,
                    segment.span,
                    format!("{name} cannot take a ligand"),
                    "this entity would be dropped",
                )
                .help(format!("{name} designs protein chains only")),
            );
        }
        if capabilities.designed_identity {
            continue;
        }
        for (position, element) in definition.elements().iter().enumerate() {
            if !restricts_designed_identity(element.slot()) {
                continue;
            }
            let span = segment
                .elements
                .get(position)
                .map(|element| element.span)
                .unwrap_or(segment.span);
            gaps.push(
                refusal(
                    source,
                    span,
                    format!("{name} cannot fix the identity of a designed position"),
                    "this residue would be designed freely",
                )
                .help(format!(
                    "{name} constrains sequence only where it also fixes coordinates"
                )),
            );
        }
    }
    for (index, constraint) in design.spec.constraints().iter().enumerate() {
        let span = design
            .constraint_span(index)
            .or_else(|| design.segments.first().map(|segment| segment.span))
            .unwrap_or(Span::new(0, 0));
        let gap = match constraint {
            Constraint::Bond { .. } if !capabilities.bonds => Some("a bond"),
            Constraint::Distance { .. } if !capabilities.distances => Some("a distance"),
            Constraint::SecondaryStructure { .. } if !capabilities.secondary_structure => {
                Some("secondary structure")
            }
            Constraint::Bond { .. } | Constraint::Distance { .. } => None,
            Constraint::SecondaryStructure { .. } => None,
            _ => Some("this"),
        };
        if let Some(what) = gap {
            gaps.push(refusal(
                source,
                span,
                format!("{name} cannot express {what} constraint"),
                "this constraint has no counterpart",
            ));
        }
    }
    let mut distinct: Vec<&std::path::Path> = design
        .targets
        .iter()
        .map(|target| target.resolved.as_path())
        .collect();
    distinct.sort_unstable();
    distinct.dedup();
    if capabilities
        .inputs
        .is_some_and(|allowed| distinct.len() > allowed)
    {
        let allowed = capabilities.inputs.unwrap_or_default();
        let span = second_target(design).unwrap_or(Span::new(0, 0));
        gaps.push(
            refusal(
                source,
                span,
                match allowed {
                    1 => format!("{name} takes one input structure, not {}", distinct.len()),
                    allowed => format!(
                        "{name} takes {allowed} input structures, not {}",
                        distinct.len()
                    ),
                },
                "this file would be ignored",
            )
            .help("merge the chains into one structure file"),
        );
    }
    gaps
}

fn restricts_designed_identity(slot: &Slot) -> bool {
    slot.geometry == Geometry::Free && slot.identity != Identity::Any
}

fn second_target(design: &Design) -> Option<Span> {
    design
        .segments
        .iter()
        .find(|segment| segment.target().is_some_and(|index| index > 0))
        .map(|segment| segment.span)
}

fn refusal(source: SourceId, span: Span, message: String, label: &str) -> Diagnostic {
    Diagnostic::error("unsupported-by-backend", message).at(source, span, label)
}

pub fn anchored_outside_a_file(source: SourceId, design: &Design, spec: &Spec) -> Vec<Diagnostic> {
    let mut gaps = Vec::new();
    for segment in &design.segments {
        if segment.kind == SegmentKind::File {
            continue;
        }
        let Ok(definition) = spec.segment(&segment.id) else {
            continue;
        };
        for (position, element) in definition.elements().iter().enumerate() {
            if element.slot().geometry.anchor().is_none() {
                continue;
            }
            let span = segment
                .elements
                .get(position)
                .map(|element| element.span)
                .unwrap_or(segment.span);
            gaps.push(
                Diagnostic::error(
                    "unsupported-by-backend",
                    format!("chain {} carries coordinates without a file", segment.id),
                )
                .at(source, span, "this residue is anchored")
                .help("declare anchored residues through a file entity"),
            );
        }
    }
    gaps
}

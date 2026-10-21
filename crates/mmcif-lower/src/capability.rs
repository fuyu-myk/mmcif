use mmcif_ir::{Constraint, Geometry, Identity, Slot, Spec};
use mmcif_spec::{Design, Diagnostic, SegmentKind, SourceId, Span};

pub struct Capabilities {
    pub name: &'static str,
    pub designed_identity: bool,
    pub bonds: bool,
    pub distances: bool,
    pub secondary_structure: bool,
    pub ligands: bool,
    pub inputs: usize,
}

pub fn unsupported(
    source: SourceId,
    design: &Design,
    capabilities: &Capabilities,
) -> Vec<Diagnostic> {
    let name = capabilities.name;
    let mut gaps = Vec::new();
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
            _ => None,
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
    if design.targets.len() > capabilities.inputs {
        let span = second_target(design).unwrap_or(Span::new(0, 0));
        gaps.push(
            refusal(
                source,
                span,
                format!(
                    "{name} takes {} input structure, not {}",
                    capabilities.inputs,
                    design.targets.len()
                ),
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
        .find(|segment| segment.target() > Some(0))
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

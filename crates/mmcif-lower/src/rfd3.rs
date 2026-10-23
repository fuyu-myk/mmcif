use mmcif_ir::{Element, Identity};
use mmcif_spec::{Design, Diagnostic, SegmentKind, SegmentSource, SourceId, Span};

use crate::capability::{self, Capabilities};
use crate::ranges;
use crate::slots::{self, Role};
use crate::text;

pub const CAPABILITIES: Capabilities = Capabilities {
    name: "rfd3",
    schema: "RosettaCommons/foundry v0.1.9 models/rfd3",
    designed_identity: false,
    bonds: false,
    distances: false,
    secondary_structure: false,
    ligands: false,
    inputs: Some(1),
};

pub fn emit(source: SourceId, design: &Design) -> Result<String, Vec<Diagnostic>> {
    let mut problems = capability::unsupported(source, design, &CAPABILITIES);
    problems.extend(capability::anchored_outside_a_file(
        source,
        design,
        &design.spec,
    ));
    if !problems.is_empty() {
        return Err(problems);
    }

    let mut pieces = Vec::new();
    let mut unfixed = Vec::new();
    let mut takes_input = false;
    for segment in &design.segments {
        match segment.kind {
            SegmentKind::File => match contig(source, design, segment) {
                Ok(piece) => {
                    takes_input |= !slots::motif_labels(design, segment).is_empty();
                    pieces.push(piece);
                    let labels = free_sequence(design, segment);
                    if !labels.is_empty() {
                        unfixed.push(ranges::contig(segment.id.as_str(), &labels));
                    }
                }
                Err(problem) => problems.push(problem),
            },
            SegmentKind::Protein => match designed(source, design, segment) {
                Ok(piece) => pieces.push(piece),
                Err(problem) => problems.push(problem),
            },
            SegmentKind::Ligand => {}
        }
    }
    if !problems.is_empty() {
        return Err(problems);
    }
    if pieces.is_empty() {
        return Err(vec![
            Diagnostic::error("cannot-lower", "the design has no chain to build").at(
                source,
                Span::new(0, 0),
                "nothing to diffuse",
            ),
        ]);
    }

    let mut fields: Vec<(&str, String)> = vec![("contig", pieces.join(",/0,"))];
    if let Some(path) = design.targets.first().filter(|_| takes_input) {
        fields.push(("input", path.declared.clone()));
    }
    if !unfixed.is_empty() {
        fields.push(("select_unfixed_sequence", unfixed.join(",")));
    }
    fields.sort_by_key(|(key, _)| *key);

    let mut out = String::from("{\n  \"design\": {\n");
    for (position, (key, value)) in fields.iter().enumerate() {
        let comma = match position + 1 == fields.len() {
            true => "",
            false => ",",
        };
        out.push_str(&format!(
            "    {}: {}{comma}\n",
            text::json(key),
            text::json(value)
        ));
    }
    out.push_str("  }\n}\n");
    Ok(out)
}

fn contig(
    source: SourceId,
    design: &Design,
    segment: &SegmentSource,
) -> Result<String, Diagnostic> {
    let mut pieces: Vec<String> = Vec::new();
    let mut motif: Vec<u32> = Vec::new();
    let mut designed = 0u32;
    for (role, label) in slots::roles(design, segment) {
        match role {
            Role::Designed => {
                flush(&mut pieces, segment, &mut motif, source)?;
                designed += 1;
            }
            _ => {
                if designed > 0 {
                    pieces.push(designed.to_string());
                    designed = 0;
                }
                let Some(label) = label else {
                    return Err(Diagnostic::error(
                        "cannot-lower",
                        format!("a residue of chain {} has no index", segment.id),
                    )
                    .at(source, segment.span, "this residue cannot be named"));
                };
                motif.push(label);
            }
        }
    }
    flush(&mut pieces, segment, &mut motif, source)?;
    if designed > 0 {
        pieces.push(designed.to_string());
    }
    Ok(pieces.join(","))
}

fn flush(
    pieces: &mut Vec<String>,
    segment: &SegmentSource,
    motif: &mut Vec<u32>,
    source: SourceId,
) -> Result<(), Diagnostic> {
    if motif.is_empty() {
        return Ok(());
    }
    let runs = ranges::runs(motif);
    motif.clear();
    match runs.as_slice() {
        [(start, end)] => {
            pieces.push(ranges::span(segment.id.as_str(), *start, *end));
            Ok(())
        }
        [(_, first), (second, _), ..] => Err(Diagnostic::error(
            "cannot-lower",
            format!("rfd3 cannot leave a gap inside chain {}", segment.id),
        )
        .at(
            source,
            segment.span,
            format!(
                "chain {} is missing residues {} to {}",
                segment.id,
                first + 1,
                second - 1
            ),
        )
        .help("select a range with no missing residues, or redesign the gap")),
        [] => Ok(()),
    }
}

fn free_sequence(design: &Design, segment: &SegmentSource) -> Vec<u32> {
    slots::roles(design, segment)
        .into_iter()
        .filter(|(role, _)| *role == Role::MotifFreeSequence)
        .filter_map(|(_, label)| label)
        .collect()
}

fn designed(
    source: SourceId,
    design: &Design,
    segment: &SegmentSource,
) -> Result<String, Diagnostic> {
    let Ok(definition) = design.spec.segment(&segment.id) else {
        return Err(Diagnostic::error(
            "cannot-lower",
            format!("chain {} has no elements", segment.id),
        )
        .at(source, segment.span, "nothing to diffuse"));
    };
    let mut pieces = Vec::new();
    for (position, element) in definition.elements().iter().enumerate() {
        let span = segment
            .elements
            .get(position)
            .map(|info| info.span)
            .unwrap_or(segment.span);
        let Some(text) = piece(element) else {
            return Err(Diagnostic::error(
                "cannot-lower",
                "this run has no contig form".to_string(),
            )
            .at(source, span, "a contig takes a length or a range"));
        };
        pieces.push(text);
    }
    Ok(pieces.join(","))
}

fn piece(element: &Element) -> Option<String> {
    match element.slot().identity {
        Identity::Any => {
            let length = element.length();
            Some(if length.is_symbolic() {
                format!("{}-{}", length.min(), length.max())
            } else {
                length.min().to_string()
            })
        }
        _ => None,
    }
}

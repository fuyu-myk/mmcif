use mmcif_ir::{Element, Identity, Length};
use mmcif_spec::{Design, Diagnostic, SegmentKind, SegmentSource, SourceId, Span};

use crate::capability::{self, Capabilities};
use crate::ranges;

pub const CAPABILITIES: Capabilities = Capabilities {
    name: "rfd3",
    designed_identity: false,
    bonds: false,
    distances: false,
    secondary_structure: false,
    ligands: false,
    inputs: 1,
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
    for segment in &design.segments {
        match segment.kind {
            SegmentKind::File => {
                let mut labels: Vec<u32> = segment.labels().collect();
                labels.sort_unstable();
                pieces.push(ranges::contig(segment.id.as_str(), &labels));
                let redesigned = redesigned(design, segment);
                if !redesigned.is_empty() {
                    unfixed.push(ranges::contig(segment.id.as_str(), &redesigned));
                }
            }
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
    if let Some(path) = design.targets.first() {
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
        out.push_str(&format!("    \"{key}\": {}{comma}\n", quoted(value)));
    }
    out.push_str("  }\n}\n");
    Ok(out)
}

fn redesigned(design: &Design, segment: &SegmentSource) -> Vec<u32> {
    let Ok(definition) = design.spec.segment(&segment.id) else {
        return Vec::new();
    };
    let mut labels: Vec<u32> = definition
        .elements()
        .iter()
        .enumerate()
        .filter(|(_, element)| element.slot().geometry.anchor().is_none())
        .filter_map(|(position, _)| segment.elements.get(position).and_then(|info| info.label))
        .collect();
    labels.sort_unstable();
    labels
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
        match piece(element) {
            Some(text) => pieces.push(text),
            None => {
                return Err(Diagnostic::error(
                    "cannot-lower",
                    "this run has no contig form".to_string(),
                )
                .at(source, span, "a contig takes a length or a range"));
            }
        }
    }
    Ok(pieces.join(","))
}

fn piece(element: &Element) -> Option<String> {
    match element.slot().identity {
        Identity::Any => match element.length() {
            Length::Exact(count) => Some(count.to_string()),
            Length::Range { min, max } => Some(format!("{min}-{max}")),
            _ => None,
        },
        _ => None,
    }
}

fn quoted(text: &str) -> String {
    format!("\"{}\"", text.replace('\\', "\\\\").replace('"', "\\\""))
}

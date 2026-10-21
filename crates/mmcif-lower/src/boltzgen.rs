use mmcif_ir::{
    Constraint, Element, Identity, Length, Scope, SecondaryStructure, SegmentId, Slot, Symbolic,
    SymbolicAddress,
};
use mmcif_spec::{Design, Diagnostic, SegmentKind, SegmentSource, SourceId, Span, one_letter};

use crate::capability::{self, Capabilities};
use crate::ranges;
use crate::text;

pub const CAPABILITIES: Capabilities = Capabilities {
    name: "boltzgen",
    schema: "boltzgen design specification, read 2026-10-21",
    designed_identity: true,
    bonds: true,
    distances: false,
    secondary_structure: true,
    ligands: true,
    inputs: None,
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

    let structures = structures(source, design, &mut problems);
    let mut out = String::from("entities:\n");
    let mut position = 0;
    while let Some(segment) = design.segments.get(position) {
        match segment.kind {
            SegmentKind::File => {
                let target = segment.target();
                let mut group = Vec::new();
                while let Some(next) = design.segments.get(position) {
                    if next.kind != SegmentKind::File || next.target() != target {
                        break;
                    }
                    group.push(next);
                    position += 1;
                }
                file(
                    &mut out,
                    source,
                    design,
                    target,
                    &group,
                    &structures,
                    &mut problems,
                );
            }
            SegmentKind::Protein => {
                protein(&mut out, source, design, segment, &mut problems);
                position += 1;
            }
            SegmentKind::Ligand => {
                ligand(&mut out, source, design, segment, &mut problems);
                position += 1;
            }
        }
    }
    constraints(&mut out, source, design, &mut problems);
    match problems.is_empty() {
        true => Ok(out),
        false => Err(problems),
    }
}

type Structures = Vec<(SegmentId, &'static str, Vec<u32>)>;

fn structures(source: SourceId, design: &Design, problems: &mut Vec<Diagnostic>) -> Structures {
    let mut found: Structures = Vec::new();
    for (index, constraint) in design.spec.constraints().iter().enumerate() {
        let Constraint::SecondaryStructure { scope, kind } = constraint else {
            continue;
        };
        let span = design.constraint_span(index).unwrap_or(Span::new(0, 0));
        let Some(key) = structure_key(*kind) else {
            problems.push(refusal(
                source,
                span,
                "this secondary structure has no keyword".to_string(),
            ));
            continue;
        };
        let addresses: &[SymbolicAddress] = match scope {
            Scope::Address(address) => std::slice::from_ref(address),
            Scope::Group(id) => match design.spec.group(id) {
                Ok(members) => members,
                Err(error) => {
                    problems.push(refusal(source, span, error.to_string()));
                    continue;
                }
            },
            _ => {
                problems.push(refusal(
                    source,
                    span,
                    "this constraint covers something that cannot be named".to_string(),
                ));
                continue;
            }
        };
        for address in addresses {
            let Some(label) = design.index_of(address) else {
                problems.push(refusal(
                    source,
                    design.span_of(address).unwrap_or(Span::new(0, 0)),
                    format!("{address} has no residue index"),
                ));
                continue;
            };
            match found
                .iter_mut()
                .find(|(segment, known, _)| segment == &address.segment && *known == key)
            {
                Some((_, _, labels)) => labels.push(label),
                None => found.push((address.segment.clone(), key, vec![label])),
            }
        }
    }
    found
}

fn file(
    out: &mut String,
    source: SourceId,
    design: &Design,
    target: Option<usize>,
    group: &[&SegmentSource],
    structures: &Structures,
    problems: &mut Vec<Diagnostic>,
) {
    let path = target.and_then(|index| design.targets.get(index));
    let Some(path) = path else {
        problems.push(refusal(
            source,
            group
                .first()
                .map(|segment| segment.span)
                .unwrap_or(Span::new(0, 0)),
            "a chain from a file has no file".to_string(),
        ));
        return;
    };
    out.push_str("  - file:\n");
    out.push_str(&format!("      path: {}\n", text::yaml(&path.declared)));
    out.push_str("      include:\n");
    for segment in group {
        chain(out, &segment.id, &labels(segment), "res_index");
    }
    let redesigned: Vec<(&SegmentSource, Vec<u32>)> = group
        .iter()
        .map(|segment| (*segment, crate::slots::redesigned(design, segment)))
        .filter(|(_, labels)| !labels.is_empty())
        .collect();
    if !redesigned.is_empty() {
        out.push_str("      design:\n");
        for (segment, labels) in &redesigned {
            chain(out, &segment.id, labels, "res_index");
        }
    }
    let mut lines = Vec::new();
    for segment in group {
        let mut own: Vec<&(SegmentId, &'static str, Vec<u32>)> = structures
            .iter()
            .filter(|(id, _, _)| id == &segment.id)
            .collect();
        own.sort_by_key(|(_, key, _)| rank(key));
        if own.is_empty() {
            continue;
        }
        let mut body = format!(
            "        - chain:\n            id: {}\n",
            text::yaml(segment.id.as_str())
        );
        for (_, key, labels) in own {
            body.push_str(&format!("            {key}: {}\n", ranges::write(labels)));
        }
        lines.push(body);
    }
    if !lines.is_empty() {
        out.push_str("      secondary_structure:\n");
        for line in lines {
            out.push_str(&line);
        }
    }
}

fn chain(out: &mut String, id: &SegmentId, labels: &[u32], key: &str) {
    out.push_str(&format!(
        "        - chain:\n            id: {}\n",
        text::yaml(id.as_str())
    ));
    if !labels.is_empty() {
        out.push_str(&format!("            {key}: {}\n", ranges::write(labels)));
    }
}

fn labels(segment: &SegmentSource) -> Vec<u32> {
    segment.labels().collect()
}

fn protein(
    out: &mut String,
    source: SourceId,
    design: &Design,
    segment: &SegmentSource,
    problems: &mut Vec<Diagnostic>,
) {
    let Ok(definition) = design.spec.segment(&segment.id) else {
        return;
    };
    let mut sequence = String::new();
    for (position, element) in definition.elements().iter().enumerate() {
        let span = segment
            .elements
            .get(position)
            .map(|info| info.span)
            .unwrap_or(segment.span);
        match piece(element) {
            Ok(text) => sequence.push_str(&text),
            Err(problem) => problems.push(refusal(source, span, problem)),
        }
    }
    out.push_str(&format!(
        "  - protein:\n      id: {}\n      sequence: {}\n",
        text::yaml(segment.id.as_str()),
        text::yaml(&sequence)
    ));
}

fn piece(element: &Element) -> Result<String, String> {
    let length = element.length();
    match &element.slot().identity {
        Identity::Any => Ok(match length {
            Length::Exact(count) => count.to_string(),
            Length::Range { min, max } => format!("{min}..{max}"),
            _ => return Err("this run has no length".to_string()),
        }),
        Identity::Exactly(comp) => {
            let Some(letter) = one_letter(comp) else {
                return Err(format!("{comp} has no one letter code"));
            };
            match length {
                Length::Exact(count) => Ok(letter.to_string().repeat(count as usize)),
                _ => Err(format!(
                    "a run of {comp} cannot have a variable length in a sequence"
                )),
            }
        }
        Identity::OneOf(_) => {
            Err("a set of allowed residues cannot be written as a sequence".to_string())
        }
        _ => Err("this residue cannot be written as a sequence".to_string()),
    }
}

fn ligand(
    out: &mut String,
    source: SourceId,
    design: &Design,
    segment: &SegmentSource,
    problems: &mut Vec<Diagnostic>,
) {
    let comp = design
        .spec
        .segment(&segment.id)
        .ok()
        .and_then(|definition| definition.elements().first().map(Element::slot))
        .and_then(component);
    match comp {
        Some(comp) => out.push_str(&format!(
            "  - ligand:\n      id: {}\n      ccd: {}\n",
            text::yaml(segment.id.as_str()),
            text::yaml(&comp)
        )),
        None => problems.push(refusal(
            source,
            segment.span,
            format!("ligand {} names no component", segment.id),
        )),
    }
}

fn component(slot: &Slot) -> Option<String> {
    match &slot.identity {
        Identity::Exactly(comp) => Some(comp.as_str().to_string()),
        _ => None,
    }
}

fn constraints(
    out: &mut String,
    source: SourceId,
    design: &Design,
    problems: &mut Vec<Diagnostic>,
) {
    let bonds: Vec<(usize, &Constraint<Symbolic>)> = design
        .spec
        .constraints()
        .iter()
        .enumerate()
        .filter(|(_, constraint)| matches!(constraint, Constraint::Bond { .. }))
        .collect();
    if bonds.is_empty() {
        return;
    }
    out.push_str("constraints:\n");
    for (index, constraint) in bonds {
        let Constraint::Bond { first, second } = constraint else {
            continue;
        };
        let span = design.constraint_span(index).unwrap_or(Span::new(0, 0));
        let (Some(one), Some(two)) = (
            atom(design, &first.slot, first.atom.as_str()),
            atom(design, &second.slot, second.atom.as_str()),
        ) else {
            problems.push(refusal(
                source,
                span,
                "a bond names a residue that has no index".to_string(),
            ));
            continue;
        };
        out.push_str(&format!(
            "  - bond:\n      atom1: {one}\n      atom2: {two}\n"
        ));
    }
}

fn atom(design: &Design, address: &SymbolicAddress, name: &str) -> Option<String> {
    let index = design.index_of(address)?;
    Some(format!(
        "[{}, {index}, {}]",
        text::yaml_flow(address.segment.as_str()),
        text::yaml_flow(name)
    ))
}

fn rank(key: &str) -> usize {
    ["loop", "helix", "sheet"]
        .iter()
        .position(|known| *known == key)
        .unwrap_or(usize::MAX)
}

fn structure_key(kind: SecondaryStructure) -> Option<&'static str> {
    match kind {
        SecondaryStructure::Helix => Some("helix"),
        SecondaryStructure::Strand => Some("sheet"),
        SecondaryStructure::Loop => Some("loop"),
        _ => None,
    }
}

fn refusal(source: SourceId, span: Span, message: String) -> Diagnostic {
    Diagnostic::error("cannot-lower", message).at(source, span, "this cannot be written out")
}

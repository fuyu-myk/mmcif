use std::path::Path;

use mmcif_ir::{
    AtomAddress, Constraint, Element, ElementIndex, Geometry, GroupId, Identity, IrError, Length,
    ResidueRef, Scope, Segment, SegmentId, Slot, Spec, Symbolic,
    SymbolicAddress,
};
use mmcif_model::AtomName;

use crate::diagnostic::{Diagnostic, Report};
use crate::document::{
    AtomRef, Document, Entity, EntityKind, FileEntity, Protein, Ranges, Selection, Spanned,
    StructureSelection,
};
use crate::resolve::{self, Residue, Target};
use crate::sequence::Token;
use crate::span::{SourceId, Span};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Addressing {
    Positional,
    Label,
}

#[derive(Debug, Clone)]
pub struct ElementInfo {
    pub span: Span,
    pub label: Option<u32>,
    pub target: Option<usize>,
}

#[derive(Debug, Clone)]
pub struct SegmentInfo {
    pub id: SegmentId,
    pub span: Span,
    pub addressing: Addressing,
    pub elements: Vec<ElementInfo>,
}

impl SegmentInfo {
    fn locate(&self, index: u32, spec: &Spec) -> Option<(ElementIndex, u32)> {
        match self.addressing {
            Addressing::Label => self
                .elements
                .iter()
                .position(|element| element.label == Some(index))
                .map(|element| (ElementIndex::new(element as u32), 0)),
            Addressing::Positional => {
                let segment = spec.segment(&self.id).ok()?;
                let mut seen = 0;
                for (position, element) in segment.elements().iter().enumerate() {
                    let count = element.length().min();
                    if index <= seen + count {
                        return Some((ElementIndex::new(position as u32), index - seen - 1));
                    }
                    seen += count;
                }
                None
            }
        }
    }

    fn minimum(&self, spec: &Spec) -> u32 {
        spec.segment(&self.id)
            .map(|segment| segment.length_bounds().0.min(u64::from(u32::MAX)) as u32)
            .unwrap_or(0)
    }
}

pub struct Lowered {
    pub spec: Spec,
    pub segments: Vec<SegmentInfo>,
    pub constraint_spans: Vec<Span>,
    pub targets: Vec<Target>,
}

impl Lowered {
    pub fn segment(&self, id: &str) -> Option<&SegmentInfo> {
        self.segments
            .iter()
            .find(|segment| segment.id.as_str() == id)
    }
}

pub fn lower(source: SourceId, base: &Path, document: &Document, report: &mut Report) -> Lowered {
    let mut lowered = Lowered {
        spec: Spec::new(),
        segments: Vec::new(),
        constraint_spans: Vec::new(),
        targets: Vec::new(),
    };
    let mut groups = 0;
    for entity in &document.entities {
        match &entity.kind {
            EntityKind::Protein(protein) => {
                designed(source, &mut lowered, entity, protein, report);
            }
            EntityKind::Ligand(ligand) => {
                let slot = Slot::new(Identity::Exactly(ligand.ccd.value.as_str().into()), Geometry::Free);
                push(
                    source,
                    &mut lowered,
                    &ligand.id,
                    entity.span,
                    Addressing::Positional,
                    vec![(Element::Slot(slot), ElementInfo {
                        span: ligand.ccd.span,
                        label: None,
                        target: None,
                    })],
                    report,
                );
            }
            EntityKind::File(file) => {
                from_file(source, base, &mut lowered, file, report, &mut groups);
            }
        }
    }
    for bond in &document.bonds {
        if let Some(constraint) = bond_constraint(source, &lowered, bond, report) {
            match lowered.spec.push_constraint(constraint) {
                Ok(()) => lowered.constraint_spans.push(bond.span),
                Err(error) => report.push(ir_error(source, bond.span, error)),
            }
        }
    }
    lowered
}

fn designed(
    source: SourceId,
    lowered: &mut Lowered,
    entity: &Entity,
    protein: &Protein,
    report: &mut Report,
) {
    let mut elements = Vec::new();
    for piece in &protein.sequence.value {
        let info = ElementInfo {
            span: piece.span,
            label: None,
            target: None,
        };
        let element = match &piece.token {
            Token::Residue(comp) => Element::Slot(Slot::new(
                Identity::Exactly(comp.clone()),
                Geometry::Free,
            )),
            Token::Designed { min, max } => match Length::range(*min, *max) {
                Ok(length) => Element::run(Slot::de_novo(), length),
                Err(error) => {
                    report.push(ir_error(source, piece.span, error));
                    continue;
                }
            },
        };
        elements.push((element, info));
    }
    push(
        source,
        lowered,
        &protein.id,
        entity.span,
        Addressing::Positional,
        elements,
        report,
    );
}

fn from_file(
    source: SourceId,
    base: &Path,
    lowered: &mut Lowered,
    file: &FileEntity,
    report: &mut Report,
    groups: &mut usize,
) {
    let Some(loaded) = resolve::load(source, base, &file.path, report) else {
        return;
    };
    let index = lowered.targets.len();
    let target = &loaded;

    let chains: Vec<&Selection> = match file.include.is_empty() {
        false => file.include.iter().collect(),
        true => Vec::new(),
    };
    let whole: Vec<Selection> = match chains.is_empty() {
        false => Vec::new(),
        true => every_chain(target),
    };
    let selections: Vec<&Selection> = match chains.is_empty() {
        false => chains,
        true => whole.iter().collect(),
    };

    for selection in selections {
        let Some(residues) = resolve::chain(source, target, &selection.chain, report) else {
            continue;
        };
        let residues = resolve::select(
            source,
            target,
            &selection.chain,
            residues,
            selection.residues.as_ref(),
            report,
        );
        let designed = design_ranges(file, &selection.chain.value);
        let elements: Vec<(Element, ElementInfo)> = residues
            .iter()
            .map(|residue| {
                let anchored = !designed.iter().any(|ranges| ranges.contains(residue.seq.get()));
                let slot = match anchored {
                    true => Slot::context(
                        residue.comp.clone(),
                        ResidueRef::Label(residue.key.clone()),
                    ),
                    false => Slot::de_novo(),
                };
                (
                    Element::Slot(slot),
                    ElementInfo {
                        span: selection.span,
                        label: Some(residue.seq.get()),
                        target: Some(index),
                    },
                )
            })
            .collect();
        push(
            source,
            lowered,
            &selection.chain,
            selection.span,
            Addressing::Label,
            elements,
            report,
        );
        check_design(source, file, &selection.chain.value, &residues, report);
    }
    lowered.targets.push(loaded);
    for selection in &file.structure {
        structure(source, lowered, selection, report, groups);
    }
}

fn every_chain(target: &Target) -> Vec<Selection> {
    let Some(model) = target.structure.first_model() else {
        return Vec::new();
    };
    model
        .chains()
        .filter(|chain| chain.operators().is_empty())
        .map(|chain| Selection {
            span: Span::new(0, 0),
            chain: Spanned::new(Span::new(0, 0), chain.label_asym().as_str().to_string()),
            residues: None,
        })
        .collect()
}

fn design_ranges<'a>(file: &'a FileEntity, chain: &str) -> Vec<&'a Ranges> {
    file.design
        .iter()
        .filter(|selection| selection.chain.value == chain)
        .filter_map(|selection| selection.residues.as_ref().map(|ranges| &ranges.value))
        .collect()
}

fn check_design(
    source: SourceId,
    file: &FileEntity,
    chain: &str,
    residues: &[Residue],
    report: &mut Report,
) {
    for selection in file
        .design
        .iter()
        .filter(|selection| selection.chain.value == chain)
    {
        let Some(ranges) = &selection.residues else {
            continue;
        };
        if !residues
            .iter()
            .any(|residue| ranges.value.contains(residue.seq.get()))
        {
            report.push(
                Diagnostic::error(
                    "design-outside-selection",
                    format!("no selected residue of chain {chain} is redesigned"),
                )
                .at(source, ranges.span, "none of these residues are included")
                .help("every redesigned residue must also be included"),
            );
        }
    }
}

fn structure(
    source: SourceId,
    lowered: &mut Lowered,
    selection: &StructureSelection,
    report: &mut Report,
    groups: &mut usize,
) {
    let Some(segment) = lowered.segment(&selection.chain.value).cloned() else {
        report.push(
            Diagnostic::error(
                "unknown-chain",
                format!("chain {} is not included", selection.chain.value),
            )
            .at(source, selection.chain.span, "not part of this entity")
            .help("secondary structure applies to included chains"),
        );
        return;
    };
    for (kind, ranges) in &selection.regions {
        let mut members = Vec::new();
        let mut fixed = Vec::new();
        for (position, element) in segment.elements.iter().enumerate() {
            let Some(label) = element.label else { continue };
            if !ranges.value.contains(label) {
                continue;
            }
            let index = ElementIndex::new(position as u32);
            let address = SymbolicAddress::at(segment.id.clone(), index, 0);
            if lowered
                .spec
                .segment(&segment.id)
                .ok()
                .and_then(|definition| definition.element(index).ok())
                .map(|element| element.slot().geometry.anchor().is_some())
                .unwrap_or(false)
            {
                fixed.push(label);
            }
            members.push(address);
        }
        if members.is_empty() {
            report.push(
                Diagnostic::error(
                    "unknown-residue",
                    format!("chain {} has no such residues", selection.chain.value),
                )
                .at(source, ranges.span, "selects no included residue"),
            );
            continue;
        }
        if !fixed.is_empty() {
            report.push(
                Diagnostic::warning(
                    "structure-on-fixed-residue",
                    format!("{} residues already hold their structure", fixed.len()),
                )
                .at(source, ranges.span, "these residues come from the file")
                .help("secondary structure guides designed residues"),
            );
        }
        *groups += 1;
        let id = GroupId::new(format!("ss{groups}:{}", segment.id));
        if let Err(error) = lowered.spec.define_group(id.clone(), members) {
            report.push(ir_error(source, ranges.span, error));
            continue;
        }
        let constraint = Constraint::<Symbolic>::SecondaryStructure {
            scope: Scope::Group(id),
            kind: *kind,
        };
        match lowered.spec.push_constraint(constraint) {
            Ok(()) => lowered.constraint_spans.push(ranges.span),
            Err(error) => report.push(ir_error(source, ranges.span, error)),
        }
    }
}

fn push(
    source: SourceId,
    lowered: &mut Lowered,
    id: &Spanned<String>,
    span: Span,
    addressing: Addressing,
    elements: Vec<(Element, ElementInfo)>,
    report: &mut Report,
) {
    if elements.is_empty() {
        return;
    }
    if let Some(known) = lowered.segment(&id.value) {
        report.push(
            Diagnostic::error("duplicate-chain", format!("chain {} is declared twice", id.value))
                .at(source, known.span, "first declared here")
                .at(source, id.span, "declared again here"),
        );
        return;
    }
    let segment = SegmentId::new(id.value.as_str());
    let (definitions, infos): (Vec<Element>, Vec<ElementInfo>) = elements.into_iter().unzip();
    match lowered
        .spec
        .push_segment(Segment::new(segment.clone(), definitions))
    {
        Ok(_) => lowered.segments.push(SegmentInfo {
            id: segment,
            span,
            addressing,
            elements: infos,
        }),
        Err(error) => report.push(ir_error(source, id.span, error)),
    }
}

fn bond_constraint(
    source: SourceId,
    lowered: &Lowered,
    bond: &crate::document::Bond,
    report: &mut Report,
) -> Option<Constraint<Symbolic>> {
    let first = address(source, lowered, &bond.first, report)?;
    let second = address(source, lowered, &bond.second, report)?;
    Some(Constraint::Bond {
        first: AtomAddress::new(first, AtomName::new(bond.first.atom.value.as_str())),
        second: AtomAddress::new(second, AtomName::new(bond.second.atom.value.as_str())),
    })
}

fn address(
    source: SourceId,
    lowered: &Lowered,
    atom: &AtomRef,
    report: &mut Report,
) -> Option<SymbolicAddress> {
    let Some(segment) = lowered.segment(&atom.chain.value) else {
        let mut names: Vec<&str> = lowered
            .segments
            .iter()
            .map(|segment| segment.id.as_str())
            .collect();
        names.sort_unstable();
        report.push(
            Diagnostic::error(
                "unknown-chain",
                format!("chain {} is not declared", atom.chain.value),
            )
            .at(source, atom.chain.span, "no entity declares this chain")
            .help(format!("declared chains: {}", names.join(", "))),
        );
        return None;
    };
    match segment.locate(atom.index.value, &lowered.spec) {
        Some((element, offset)) => {
            Some(SymbolicAddress::at(segment.id.clone(), element, offset))
        }
        None => {
            let help = match segment.addressing {
                Addressing::Positional => format!(
                    "chain {} holds {} residues at its shortest length",
                    segment.id,
                    segment.minimum(&lowered.spec)
                ),
                Addressing::Label => {
                    let labels: Vec<u32> = segment
                        .elements
                        .iter()
                        .filter_map(|element| element.label)
                        .collect();
                    match (labels.first(), labels.last()) {
                        (Some(first), Some(last)) => format!(
                            "chain {} covers residues {first} to {last}",
                            segment.id
                        ),
                        _ => format!("chain {} holds no residues", segment.id),
                    }
                }
            };
            report.push(
                Diagnostic::error(
                    "unknown-residue",
                    format!(
                        "chain {} has no residue {}",
                        atom.chain.value, atom.index.value
                    ),
                )
                .at(source, atom.index.span, "out of range")
                .help(help),
            );
            None
        }
    }
}

fn ir_error(source: SourceId, span: Span, error: IrError) -> Diagnostic {
    Diagnostic::error("invalid-spec", error.to_string()).at(source, span, "rejected by the model")
}

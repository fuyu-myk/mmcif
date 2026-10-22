use std::path::Path;

use mmcif_ir::{
    AtomAddress, Constraint, Element, ElementIndex, Geometry, GroupId, Identity, IrError, Length,
    ResidueRef, Scope, Segment, SegmentId, Slot, Spec, Symbolic, SymbolicAddress,
};
use mmcif_model::AtomName;

use crate::design::{ElementSource, SegmentKind, SegmentSource};
use crate::diagnostic::{Diagnostic, Report};
use crate::document::{
    AtomRef, Document, Entity, EntityKind, FileEntity, Protein, Ranges, Selection, Spanned,
    StructureSelection,
};
use crate::resolve::{self, Residue, Target};
use crate::sequence::Token;
use crate::span::{SourceId, Span};

/// How a residue index in a specification picks a residue of a chain.
///
/// # Examples
///
/// ```
/// use mmcif_ir::SegmentId;
/// use mmcif_spec::{Addressing, check_text};
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
/// ";
/// let checked = check_text("design.yaml", yaml, &base);
/// let design = checked.design.as_ref().unwrap();
///
/// let addressing = |id: &str| design.segment(&SegmentId::new(id)).unwrap().kind.addressing();
/// assert_eq!(addressing("A"), Addressing::Label);
/// assert_eq!(addressing("S"), Addressing::Positional);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Addressing {
    /// Indices count positions along the chain from 1, with every run of
    /// designed residues at its minimum length.
    Positional,
    /// Indices are the `label_seq_id` numbers of the structure file.
    Label,
}

fn locate(segment: &SegmentSource, index: u32, spec: &Spec) -> Option<(ElementIndex, u32)> {
    match segment.kind.addressing() {
        Addressing::Label => segment
            .elements
            .iter()
            .position(|element| element.label == Some(index))
            .map(|element| (ElementIndex::new(element as u32), 0)),
        Addressing::Positional => {
            let definition = spec.segment(&segment.id).ok()?;
            let index = u64::from(index);
            let mut seen = 0u64;
            for (position, element) in definition.elements().iter().enumerate() {
                let count = u64::from(element.length().min());
                if index <= seen + count {
                    return Some((
                        ElementIndex::new(position as u32),
                        (index - seen - 1) as u32,
                    ));
                }
                seen += count;
            }
            None
        }
    }
}

fn minimum(segment: &SegmentSource, spec: &Spec) -> u32 {
    spec.segment(&segment.id)
        .map(|definition| definition.length_bounds().0.min(u64::from(u32::MAX)) as u32)
        .unwrap_or(0)
}

pub struct Lowered {
    pub spec: Spec,
    pub segments: Vec<SegmentSource>,
    pub constraint_spans: Vec<Span>,
    pub targets: Vec<Target>,
}

impl Lowered {
    pub fn segment(&self, id: &str) -> Option<&SegmentSource> {
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
                let slot = Slot::new(
                    Identity::Exactly(ligand.ccd.value.as_str().into()),
                    Geometry::Free,
                );
                push(
                    source,
                    &mut lowered,
                    &ligand.id,
                    entity.span,
                    SegmentKind::Ligand,
                    vec![(
                        Element::Slot(slot),
                        ElementSource {
                            span: ligand.ccd.span,
                            label: None,
                            target: None,
                        },
                    )],
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
        let info = ElementSource {
            span: piece.span,
            label: None,
            target: None,
        };
        let element = match &piece.token {
            Token::Residue(comp) => {
                Element::Slot(Slot::new(Identity::Exactly(comp.clone()), Geometry::Free))
            }
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
        SegmentKind::Protein,
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
    let own = lowered.segments.len();
    let target = &loaded;

    let whole = match file.include.is_empty() {
        true => every_chain(target, file.path.span),
        false => Vec::new(),
    };
    let selections: Vec<&Selection> = match file.include.is_empty() {
        true => whole.iter().collect(),
        false => file.include.iter().collect(),
    };
    let included: Vec<&str> = selections
        .iter()
        .map(|selection| selection.chain.value.as_str())
        .collect();
    for selection in &file.design {
        if !included.contains(&selection.chain.value.as_str()) {
            report.push(
                Diagnostic::error(
                    "unknown-chain",
                    format!("chain {} is not included", selection.chain.value),
                )
                .at(source, selection.chain.span, "not part of this entity")
                .help("a redesigned chain must also be included"),
            );
        }
    }

    let explicit = !file.include.is_empty();
    for selection in selections {
        let Some(residues) = resolve::chain(source, target, &selection.chain, report) else {
            continue;
        };
        if residues.is_empty() {
            if explicit {
                report.push(resolve::empty_chain(source, target, &selection.chain));
            }
            continue;
        }
        let residues = resolve::select(
            source,
            &selection.chain,
            residues,
            selection.residues.as_ref(),
            report,
        );
        let designed = design_ranges(file, &selection.chain.value);
        let elements: Vec<(Element, ElementSource)> = residues
            .iter()
            .map(|residue| {
                let anchored = !designed
                    .iter()
                    .any(|ranges| ranges.contains(residue.seq.get()));
                let slot = match anchored {
                    true => {
                        Slot::context(residue.comp.clone(), ResidueRef::Label(residue.key.clone()))
                    }
                    false => Slot::de_novo(),
                };
                (
                    Element::Slot(slot),
                    ElementSource {
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
            SegmentKind::File,
            elements,
            report,
        );
        check_design(source, file, &selection.chain.value, &residues, report);
        if let Some(ranges) = &selection.residues {
            numbering(source, target, &selection.chain.value, ranges, report);
        }
    }
    for selection in &file.design {
        if let Some(ranges) = &selection.residues {
            numbering(source, target, &selection.chain.value, ranges, report);
        }
    }
    for selection in &file.structure {
        for (_, ranges) in &selection.regions {
            numbering(source, target, &selection.chain.value, ranges, report);
        }
    }
    lowered.targets.push(loaded);
    for selection in &file.structure {
        structure(source, lowered, selection, own, report, groups);
    }
}

fn written(ranges: &Ranges) -> Vec<u32> {
    match ranges {
        Ranges::All => Vec::new(),
        Ranges::List(list) => list
            .iter()
            .flat_map(|range| {
                std::iter::once(range.min).chain(range.max.filter(|max| *max != range.min))
            })
            .collect(),
    }
}

fn numbering(
    source: SourceId,
    target: &Target,
    chain: &str,
    ranges: &Spanned<Ranges>,
    report: &mut Report,
) {
    numbered(
        source,
        target,
        chain,
        ranges.span,
        &written(&ranges.value),
        report,
    );
}

fn numbered(
    source: SourceId,
    target: &Target,
    chain: &str,
    span: Span,
    indices: &[u32],
    report: &mut Report,
) {
    let differences = resolve::author_differences(target, chain, indices);
    if differences.is_empty() {
        return;
    }
    let reverse = differences.iter().find_map(|(index, _)| {
        resolve::label_for_author(target, chain, *index).map(|found| (*index, found))
    });
    let help = match reverse {
        Some((index, (author_chain, label))) => format!(
            "indices here count label numbers; author {author_chain} {index} is label {label}"
        ),
        None => "indices here count label numbers".to_string(),
    };
    let labels: Vec<String> = differences.into_iter().map(|(_, text)| text).collect();
    report.push(
        Diagnostic::warning(
            "author-numbering-differs",
            format!("chain {chain} is numbered differently by its authors"),
        )
        .at(source, span, labels.join(", "))
        .help(help),
    );
}

fn every_chain(target: &Target, span: Span) -> Vec<Selection> {
    let Some(model) = target.structure.first_model() else {
        return Vec::new();
    };
    model
        .chains()
        .filter(|chain| chain.operators().is_empty())
        .map(|chain| Selection {
            span,
            chain: Spanned::new(span, chain.label_asym().as_str().to_string()),
            residues: None,
        })
        .collect()
}

fn design_ranges(file: &FileEntity, chain: &str) -> Vec<Ranges> {
    file.design
        .iter()
        .filter(|selection| selection.chain.value == chain)
        .map(|selection| match &selection.residues {
            Some(ranges) => ranges.value.clone(),
            None => Ranges::All,
        })
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
        let (ranges, span) = match &selection.residues {
            Some(ranges) => (ranges.value.clone(), ranges.span),
            None => (Ranges::All, selection.chain.span),
        };
        if !residues
            .iter()
            .any(|residue| ranges.contains(residue.seq.get()))
        {
            report.push(
                Diagnostic::error(
                    "design-outside-selection",
                    format!("no selected residue of chain {chain} is redesigned"),
                )
                .at(source, span, "none of these residues are included")
                .help("every redesigned residue must also be included"),
            );
        }
    }
}

fn structure(
    source: SourceId,
    lowered: &mut Lowered,
    selection: &StructureSelection,
    own: usize,
    report: &mut Report,
    groups: &mut usize,
) {
    let Lowered {
        spec,
        segments,
        constraint_spans,
        ..
    } = lowered;
    let Some(segment) = segments
        .get(own..)
        .unwrap_or_default()
        .iter()
        .find(|segment| segment.id.as_str() == selection.chain.value)
    else {
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
            if spec
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
        if let Err(error) = spec.define_group(id.clone(), members) {
            report.push(ir_error(source, ranges.span, error));
            continue;
        }
        let constraint = Constraint::<Symbolic>::SecondaryStructure {
            scope: Scope::Group(id),
            kind: *kind,
        };
        match spec.push_constraint(constraint) {
            Ok(()) => constraint_spans.push(ranges.span),
            Err(error) => report.push(ir_error(source, ranges.span, error)),
        }
    }
}

fn push(
    source: SourceId,
    lowered: &mut Lowered,
    id: &Spanned<String>,
    span: Span,
    kind: SegmentKind,
    elements: Vec<(Element, ElementSource)>,
    report: &mut Report,
) {
    if elements.is_empty() {
        return;
    }
    if let Some(known) = lowered.segment(&id.value) {
        report.push(
            Diagnostic::error(
                "duplicate-chain",
                format!("chain {} is declared twice", id.value),
            )
            .at(source, known.span, "first declared here")
            .at(source, id.span, "declared again here"),
        );
        return;
    }
    let segment = SegmentId::new(id.value.as_str());
    let (definitions, infos): (Vec<Element>, Vec<ElementSource>) = elements.into_iter().unzip();
    match lowered
        .spec
        .push_segment(Segment::new(segment.clone(), definitions))
    {
        Ok(_) => lowered.segments.push(SegmentSource {
            id: segment,
            kind,
            span,
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
    match locate(segment, atom.index.value, &lowered.spec) {
        Some((element, offset)) => {
            if let Some(target) = segment
                .target()
                .and_then(|index| lowered.targets.get(index))
            {
                numbered(
                    source,
                    target,
                    segment.id.as_str(),
                    atom.index.span,
                    &[atom.index.value],
                    report,
                );
            }
            Some(SymbolicAddress::at(segment.id.clone(), element, offset))
        }
        None => {
            let help = match segment.kind.addressing() {
                Addressing::Positional => format!(
                    "chain {} holds {} residues at its shortest length",
                    segment.id,
                    minimum(segment, &lowered.spec)
                ),
                Addressing::Label => {
                    let labels: Vec<u32> = segment
                        .elements
                        .iter()
                        .filter_map(|element| element.label)
                        .collect();
                    match (labels.first(), labels.last()) {
                        (Some(first), Some(last)) => {
                            format!("chain {} covers residues {first} to {last}", segment.id)
                        }
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

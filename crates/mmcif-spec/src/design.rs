use std::path::PathBuf;

use mmcif_ir::{ElementIndex, SegmentId, Spec, SymbolicAddress};

use crate::diagnostic::Diagnostic;
use crate::lower::{Addressing, Lowered};
use crate::span::Span;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SegmentKind {
    Protein,
    Ligand,
    File,
}

impl SegmentKind {
    pub fn addressing(self) -> Addressing {
        match self {
            SegmentKind::File => Addressing::Label,
            SegmentKind::Protein | SegmentKind::Ligand => Addressing::Positional,
        }
    }
}

#[derive(Debug, Clone)]
pub struct ElementSource {
    pub span: Span,
    pub label: Option<u32>,
    pub target: Option<usize>,
}

#[derive(Debug, Clone)]
pub struct SegmentSource {
    pub id: SegmentId,
    pub kind: SegmentKind,
    pub span: Span,
    pub elements: Vec<ElementSource>,
}

impl SegmentSource {
    pub fn target(&self) -> Option<usize> {
        self.elements.iter().find_map(|element| element.target)
    }

    pub fn labels(&self) -> impl Iterator<Item = u32> + use<'_> {
        self.elements.iter().filter_map(|element| element.label)
    }

    pub fn element(&self, index: ElementIndex) -> Option<&ElementSource> {
        self.elements.get(index.as_usize())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TargetPath {
    pub declared: String,
    pub resolved: PathBuf,
}

#[derive(Debug, Clone)]
pub struct Design {
    pub spec: Spec,
    pub segments: Vec<SegmentSource>,
    pub targets: Vec<TargetPath>,
    constraint_spans: Vec<Span>,
    skipped: Vec<Diagnostic>,
}

impl Design {
    pub fn skipped(&self) -> &[Diagnostic] {
        &self.skipped
    }

    pub fn note_skipped(&mut self, diagnostics: impl IntoIterator<Item = Diagnostic>) {
        self.skipped.extend(diagnostics);
    }

    pub fn target_of(&self, segment: &SegmentSource) -> Option<&TargetPath> {
        self.targets.get(segment.target()?)
    }

    pub fn segment(&self, id: &SegmentId) -> Option<&SegmentSource> {
        self.segments.iter().find(|segment| &segment.id == id)
    }

    pub fn span_of(&self, address: &SymbolicAddress) -> Option<Span> {
        let segment = self.segment(&address.segment)?;
        segment
            .element(address.element)
            .map(|element| element.span)
            .or(Some(segment.span))
    }

    pub fn constraint_span(&self, index: usize) -> Option<Span> {
        self.constraint_spans.get(index).copied()
    }

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

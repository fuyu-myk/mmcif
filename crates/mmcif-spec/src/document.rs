use mmcif_ir::SecondaryStructure;

use crate::sequence::Piece;
use crate::span::Span;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Spanned<T> {
    pub span: Span,
    pub value: T,
}

impl<T> Spanned<T> {
    pub fn new(span: Span, value: T) -> Spanned<T> {
        Spanned { span, value }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Document {
    pub entities: Vec<Entity>,
    pub bonds: Vec<Bond>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Entity {
    pub span: Span,
    pub kind: EntityKind,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EntityKind {
    Protein(Protein),
    Ligand(Ligand),
    File(FileEntity),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Protein {
    pub id: Spanned<String>,
    pub sequence: Spanned<Vec<Piece>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ligand {
    pub id: Spanned<String>,
    pub ccd: Spanned<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileEntity {
    pub path: Spanned<String>,
    pub include: Vec<Selection>,
    pub design: Vec<Selection>,
    pub structure: Vec<StructureSelection>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Selection {
    pub span: Span,
    pub chain: Spanned<String>,
    pub residues: Option<Spanned<Ranges>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructureSelection {
    pub chain: Spanned<String>,
    pub regions: Vec<(SecondaryStructure, Spanned<Ranges>)>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Ranges {
    All,
    List(Vec<Range>),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Range {
    pub min: u32,
    pub max: Option<u32>,
}

impl Ranges {
    pub fn contains(&self, index: u32) -> bool {
        match self {
            Ranges::All => true,
            Ranges::List(ranges) => ranges.iter().any(|range| range.contains(index)),
        }
    }
}

impl Range {
    pub fn contains(self, index: u32) -> bool {
        index >= self.min && self.max.map(|max| index <= max).unwrap_or(true)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bond {
    pub span: Span,
    pub first: AtomRef,
    pub second: AtomRef,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomRef {
    pub span: Span,
    pub chain: Spanned<String>,
    pub index: Spanned<u32>,
    pub atom: Spanned<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ranges_cover_open_and_closed_ends() {
        let ranges = Ranges::List(vec![
            Range {
                min: 2,
                max: Some(50),
            },
            Range { min: 55, max: None },
        ]);
        assert!(ranges.contains(2) && ranges.contains(50) && ranges.contains(900));
        assert!(!ranges.contains(1) && !ranges.contains(52));
    }

    #[test]
    fn all_covers_everything_without_an_upper_bound() {
        assert!(Ranges::All.contains(1) && Ranges::All.contains(4096));
    }

    #[test]
    fn a_closed_list_stops_at_its_last_index() {
        let ranges = Ranges::List(vec![Range {
            min: 10,
            max: Some(13),
        }]);
        assert!(ranges.contains(13) && !ranges.contains(14));
    }
}

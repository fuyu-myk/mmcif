use crate::ids::{AtomName, CompId, EntityId, LabelAsymId, LabelSeqId};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EntityKind {
    Polymer,
    NonPolymer,
    Branched,
    Water,
    Macrolide,
    Other(String),
}

impl EntityKind {
    pub fn parse(text: &str) -> EntityKind {
        match text.to_ascii_lowercase().as_str() {
            "polymer" => EntityKind::Polymer,
            "non-polymer" => EntityKind::NonPolymer,
            "branched" => EntityKind::Branched,
            "water" => EntityKind::Water,
            "macrolide" => EntityKind::Macrolide,
            other => EntityKind::Other(other.to_owned()),
        }
    }

    pub fn as_mmcif_str(&self) -> &str {
        match self {
            EntityKind::Polymer => "polymer",
            EntityKind::NonPolymer => "non-polymer",
            EntityKind::Branched => "branched",
            EntityKind::Water => "water",
            EntityKind::Macrolide => "macrolide",
            EntityKind::Other(other) => other,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolymerType {
    PeptideL,
    PeptideD,
    Dna,
    Rna,
    DnaRnaHybrid,
    SaccharideD,
    SaccharideL,
    PeptideNucleicAcid,
    CyclicPseudoPeptide,
    Other(String),
}

const POLYMER_TYPE_NAMES: [(PolymerType, &str); 9] = [
    (PolymerType::PeptideL, "polypeptide(L)"),
    (PolymerType::PeptideD, "polypeptide(D)"),
    (PolymerType::Dna, "polydeoxyribonucleotide"),
    (PolymerType::Rna, "polyribonucleotide"),
    (
        PolymerType::DnaRnaHybrid,
        "polydeoxyribonucleotide/polyribonucleotide hybrid",
    ),
    (PolymerType::SaccharideD, "polysaccharide(D)"),
    (PolymerType::SaccharideL, "polysaccharide(L)"),
    (PolymerType::PeptideNucleicAcid, "peptide nucleic acid"),
    (PolymerType::CyclicPseudoPeptide, "cyclic-pseudo-peptide"),
];

impl PolymerType {
    pub fn parse(text: &str) -> PolymerType {
        POLYMER_TYPE_NAMES
            .iter()
            .find(|(_, name)| name.eq_ignore_ascii_case(text))
            .map_or_else(
                || PolymerType::Other(text.to_owned()),
                |(kind, _)| kind.clone(),
            )
    }

    pub fn as_mmcif_str(&self) -> &str {
        match self {
            PolymerType::Other(other) => other,
            known => POLYMER_TYPE_NAMES
                .iter()
                .find(|(kind, _)| kind == known)
                .map(|(_, name)| *name)
                .expect("every named polymer type is in the table"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequencePosition {
    pub(crate) seq: LabelSeqId,
    pub(crate) comps: Vec<CompId>,
}

impl SequencePosition {
    pub fn seq(&self) -> LabelSeqId {
        self.seq
    }

    pub fn comps(&self) -> &[CompId] {
        &self.comps
    }

    pub fn is_heterogeneous(&self) -> bool {
        self.comps.len() > 1
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polymer {
    pub(crate) kind: PolymerType,
    pub(crate) canonical_sequence: Option<String>,
    pub(crate) sequence: Vec<SequencePosition>,
}

impl Polymer {
    pub fn kind(&self) -> &PolymerType {
        &self.kind
    }

    pub fn canonical_sequence(&self) -> Option<&str> {
        self.canonical_sequence.as_deref()
    }

    pub fn sequence(&self) -> &[SequencePosition] {
        &self.sequence
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BranchMonomer {
    pub(crate) num: u32,
    pub(crate) comp: CompId,
}

impl BranchMonomer {
    pub fn num(&self) -> u32 {
        self.num
    }

    pub fn comp(&self) -> &CompId {
        &self.comp
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BranchAtom {
    pub(crate) num: u32,
    pub(crate) comp: CompId,
    pub(crate) atom: AtomName,
    pub(crate) leaving_atom: Option<AtomName>,
}

impl BranchAtom {
    pub fn num(&self) -> u32 {
        self.num
    }

    pub fn comp(&self) -> &CompId {
        &self.comp
    }

    pub fn atom(&self) -> &AtomName {
        &self.atom
    }

    pub fn leaving_atom(&self) -> Option<&AtomName> {
        self.leaving_atom.as_ref()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BranchLink {
    pub(crate) first: BranchAtom,
    pub(crate) second: BranchAtom,
    pub(crate) order: Option<String>,
}

impl BranchLink {
    pub fn first(&self) -> &BranchAtom {
        &self.first
    }

    pub fn second(&self) -> &BranchAtom {
        &self.second
    }

    pub fn order(&self) -> Option<&str> {
        self.order.as_deref()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Branch {
    pub(crate) kind: Option<String>,
    pub(crate) monomers: Vec<BranchMonomer>,
    pub(crate) links: Vec<BranchLink>,
}

impl Branch {
    pub fn kind(&self) -> Option<&str> {
        self.kind.as_deref()
    }

    pub fn monomers(&self) -> &[BranchMonomer] {
        &self.monomers
    }

    pub fn links(&self) -> &[BranchLink] {
        &self.links
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Entity {
    pub(crate) id: EntityId,
    pub(crate) kind: EntityKind,
    pub(crate) description: Option<String>,
    pub(crate) polymer: Option<Polymer>,
    pub(crate) branch: Option<Branch>,
    pub(crate) asym_ids: Vec<LabelAsymId>,
}

impl Entity {
    pub fn id(&self) -> &EntityId {
        &self.id
    }

    pub fn kind(&self) -> &EntityKind {
        &self.kind
    }

    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }

    pub fn polymer(&self) -> Option<&Polymer> {
        self.polymer.as_ref()
    }

    pub fn branch(&self) -> Option<&Branch> {
        self.branch.as_ref()
    }

    pub fn asym_ids(&self) -> &[LabelAsymId] {
        &self.asym_ids
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polymer_types_round_trip_through_mmcif_names() {
        for (kind, name) in POLYMER_TYPE_NAMES {
            assert_eq!(PolymerType::parse(name), kind);
            assert_eq!(kind.as_mmcif_str(), name);
        }
        assert_eq!(
            PolymerType::parse("something new"),
            PolymerType::Other("something new".into())
        );
    }

    #[test]
    fn entity_kinds_are_case_insensitive() {
        assert_eq!(EntityKind::parse("Non-polymer"), EntityKind::NonPolymer);
        assert_eq!(EntityKind::parse("WATER").as_mmcif_str(), "water");
    }
}

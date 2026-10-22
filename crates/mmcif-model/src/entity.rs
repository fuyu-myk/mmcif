use crate::ids::{AtomName, CompId, EntityId, LabelAsymId, LabelSeqId};

/// The type of an entity, from `_entity.type`.
///
/// # Examples
///
/// ```
/// use mmcif_model::EntityKind;
///
/// let kind = EntityKind::parse("non-polymer");
/// assert_eq!(kind, EntityKind::NonPolymer);
/// assert_eq!(kind.as_mmcif_str(), "non-polymer");
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EntityKind {
    /// A polymer, `polymer`.
    Polymer,
    /// A non-polymer, `non-polymer`.
    NonPolymer,
    /// A branched entity such as an oligosaccharide, `branched`.
    Branched,
    /// Water, `water`.
    Water,
    /// A macrolide, `macrolide`.
    Macrolide,
    /// Any other type, in lower case. Empty when the file gives no type.
    Other(String),
}

impl EntityKind {
    /// Parses an `_entity.type` value, ignoring case.
    ///
    /// A value that names none of the known types becomes
    /// [`EntityKind::Other`], in lower case.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::EntityKind;
    ///
    /// assert_eq!(EntityKind::parse("WATER"), EntityKind::Water);
    /// assert_eq!(EntityKind::parse("Glycan"), EntityKind::Other("glycan".to_owned()));
    /// ```
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

    /// Returns the `_entity.type` value for this type.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::EntityKind;
    ///
    /// assert_eq!(EntityKind::Branched.as_mmcif_str(), "branched");
    /// ```
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

/// The type of a polymer, from `_entity_poly.type`.
///
/// # Examples
///
/// ```
/// use mmcif_model::PolymerType;
///
/// let kind = PolymerType::parse("polypeptide(L)");
/// assert_eq!(kind, PolymerType::PeptideL);
/// assert_eq!(kind.as_mmcif_str(), "polypeptide(L)");
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolymerType {
    /// A polypeptide of L-amino acids, `polypeptide(L)`.
    PeptideL,
    /// A polypeptide of D-amino acids, `polypeptide(D)`.
    PeptideD,
    /// DNA, `polydeoxyribonucleotide`.
    Dna,
    /// RNA, `polyribonucleotide`.
    Rna,
    /// A DNA and RNA hybrid,
    /// `polydeoxyribonucleotide/polyribonucleotide hybrid`.
    DnaRnaHybrid,
    /// A polysaccharide of D-sugars, `polysaccharide(D)`.
    SaccharideD,
    /// A polysaccharide of L-sugars, `polysaccharide(L)`.
    SaccharideL,
    /// A peptide nucleic acid, `peptide nucleic acid`.
    PeptideNucleicAcid,
    /// A cyclic pseudo-peptide, `cyclic-pseudo-peptide`.
    CyclicPseudoPeptide,
    /// Any other type, as given in the file. Empty when the file gives no
    /// type.
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
    /// Parses an `_entity_poly.type` value, ignoring case.
    ///
    /// A value that names none of the known types becomes
    /// [`PolymerType::Other`] unchanged.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::PolymerType;
    ///
    /// assert_eq!(PolymerType::parse("POLYRIBONUCLEOTIDE"), PolymerType::Rna);
    /// assert_eq!(PolymerType::parse("Other"), PolymerType::Other("Other".to_owned()));
    /// ```
    pub fn parse(text: &str) -> PolymerType {
        POLYMER_TYPE_NAMES
            .iter()
            .find(|(_, name)| name.eq_ignore_ascii_case(text))
            .map_or_else(
                || PolymerType::Other(text.to_owned()),
                |(kind, _)| kind.clone(),
            )
    }

    /// Returns the `_entity_poly.type` value for this type.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::PolymerType;
    ///
    /// assert_eq!(PolymerType::Dna.as_mmcif_str(), "polydeoxyribonucleotide");
    /// ```
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

/// One position of a polymer sequence, from `_entity_poly_seq`.
///
/// A position has more than one component when the sequence is
/// microheterogeneous.
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let polymer = structure.entity("1").unwrap().polymer().unwrap();
/// let position = &polymer.sequence()[1];
///
/// assert_eq!(position.seq().get(), 2);
/// assert!(position.is_heterogeneous());
/// let comps: Vec<&str> = position.comps().iter().map(|comp| comp.as_str()).collect();
/// assert_eq!(comps, ["SER", "THR"]);
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequencePosition {
    pub(crate) seq: LabelSeqId,
    pub(crate) comps: Vec<CompId>,
}

impl SequencePosition {
    /// Returns the label residue number of the position.
    pub fn seq(&self) -> LabelSeqId {
        self.seq
    }

    /// Returns the components at the position, in file order.
    pub fn comps(&self) -> &[CompId] {
        &self.comps
    }

    /// Returns `true` if the position has more than one component.
    pub fn is_heterogeneous(&self) -> bool {
        self.comps.len() > 1
    }
}

/// The polymer description of an entity, from `_entity_poly` and
/// `_entity_poly_seq`.
///
/// # Examples
///
/// ```
/// use mmcif_model::PolymerType;
///
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let polymer = structure.entity("1").unwrap().polymer().unwrap();
///
/// assert_eq!(polymer.kind(), &PolymerType::PeptideL);
/// assert_eq!(polymer.canonical_sequence(), Some("GSAKVT"));
/// assert_eq!(polymer.sequence().len(), 6);
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polymer {
    pub(crate) kind: PolymerType,
    pub(crate) canonical_sequence: Option<String>,
    pub(crate) sequence: Vec<SequencePosition>,
}

impl Polymer {
    /// Returns the polymer type.
    pub fn kind(&self) -> &PolymerType {
        &self.kind
    }

    /// Returns the one-letter canonical sequence from
    /// `_entity_poly.pdbx_seq_one_letter_code_can` with whitespace removed,
    /// if the file gives one.
    pub fn canonical_sequence(&self) -> Option<&str> {
        self.canonical_sequence.as_deref()
    }

    /// Returns the sequence positions, in file order.
    pub fn sequence(&self) -> &[SequencePosition] {
        &self.sequence
    }
}

/// One monomer of a branched entity, from `_pdbx_entity_branch_list`.
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let branch = structure.entity("2").unwrap().branch().unwrap();
///
/// let monomers: Vec<(u32, &str)> = branch
///     .monomers()
///     .iter()
///     .map(|monomer| (monomer.num(), monomer.comp().as_str()))
///     .collect();
/// assert_eq!(monomers, [(1, "NAG"), (2, "BMA")]);
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BranchMonomer {
    pub(crate) num: u32,
    pub(crate) comp: CompId,
}

impl BranchMonomer {
    /// Returns the monomer number within the branched entity.
    pub fn num(&self) -> u32 {
        self.num
    }

    /// Returns the component id of the monomer.
    pub fn comp(&self) -> &CompId {
        &self.comp
    }
}

/// One end of a bond in a branched entity, from `_pdbx_entity_branch_link`.
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let link = &structure.entity("2").unwrap().branch().unwrap().links()[0];
/// let first = link.first();
///
/// assert_eq!(first.num(), 2);
/// assert_eq!(first.comp().as_str(), "BMA");
/// assert_eq!(first.atom().as_str(), "C1");
/// assert_eq!(first.leaving_atom().unwrap().as_str(), "O1");
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BranchAtom {
    pub(crate) num: u32,
    pub(crate) comp: CompId,
    pub(crate) atom: AtomName,
    pub(crate) leaving_atom: Option<AtomName>,
}

impl BranchAtom {
    /// Returns the number of the monomer that holds the atom.
    pub fn num(&self) -> u32 {
        self.num
    }

    /// Returns the component id of the monomer that holds the atom.
    pub fn comp(&self) -> &CompId {
        &self.comp
    }

    /// Returns the name of the bonded atom.
    pub fn atom(&self) -> &AtomName {
        &self.atom
    }

    /// Returns the name of the atom that leaves when the bond forms, if the
    /// file gives one.
    pub fn leaving_atom(&self) -> Option<&AtomName> {
        self.leaving_atom.as_ref()
    }
}

/// A bond between two monomers of a branched entity, from
/// `_pdbx_entity_branch_link`.
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let link = &structure.entity("2").unwrap().branch().unwrap().links()[0];
///
/// assert_eq!(link.first().comp().as_str(), "BMA");
/// assert_eq!(link.second().comp().as_str(), "NAG");
/// assert_eq!(link.order(), Some("sing"));
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BranchLink {
    pub(crate) first: BranchAtom,
    pub(crate) second: BranchAtom,
    pub(crate) order: Option<String>,
}

impl BranchLink {
    /// Returns the end given by the `_1` items.
    pub fn first(&self) -> &BranchAtom {
        &self.first
    }

    /// Returns the end given by the `_2` items.
    pub fn second(&self) -> &BranchAtom {
        &self.second
    }

    /// Returns the bond order from `value_order`, such as `sing`, if the file
    /// gives one.
    pub fn order(&self) -> Option<&str> {
        self.order.as_deref()
    }
}

/// The branched description of an entity, from `_pdbx_entity_branch`,
/// `_pdbx_entity_branch_list` and `_pdbx_entity_branch_link`.
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let branch = structure.entity("2").unwrap().branch().unwrap();
///
/// assert_eq!(branch.kind(), Some("oligosaccharide"));
/// assert_eq!(branch.monomers().len(), 2);
/// assert_eq!(branch.links().len(), 1);
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Branch {
    pub(crate) kind: Option<String>,
    pub(crate) monomers: Vec<BranchMonomer>,
    pub(crate) links: Vec<BranchLink>,
}

impl Branch {
    /// Returns the branch type, such as `oligosaccharide`, if the file gives
    /// one.
    pub fn kind(&self) -> Option<&str> {
        self.kind.as_deref()
    }

    /// Returns the monomers, in file order.
    pub fn monomers(&self) -> &[BranchMonomer] {
        &self.monomers
    }

    /// Returns the bonds between monomers, in file order.
    pub fn links(&self) -> &[BranchLink] {
        &self.links
    }
}

/// A chemically distinct part of a structure, from `_entity`.
///
/// A chain names its entity through [`ChainRef::entity`].
///
/// [`ChainRef::entity`]: crate::ChainRef::entity
///
/// # Examples
///
/// ```
/// use mmcif_model::EntityKind;
///
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let entity = structure.entity("1").unwrap();
///
/// assert_eq!(entity.kind(), &EntityKind::Polymer);
/// assert_eq!(entity.description(), Some("Heavy chain"));
/// assert!(entity.polymer().is_some());
/// assert!(entity.branch().is_none());
/// assert_eq!(entity.asym_ids()[0].as_str(), "AA");
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
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
    /// Returns the entity id.
    pub fn id(&self) -> &EntityId {
        &self.id
    }

    /// Returns the entity type.
    pub fn kind(&self) -> &EntityKind {
        &self.kind
    }

    /// Returns the description from `_entity.pdbx_description`, if the file
    /// gives one.
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }

    /// Returns the polymer description, if the entity has one.
    pub fn polymer(&self) -> Option<&Polymer> {
        self.polymer.as_ref()
    }

    /// Returns the branched description, if the entity has one.
    pub fn branch(&self) -> Option<&Branch> {
        self.branch.as_ref()
    }

    /// Returns the label ids of the chains that belong to the entity.
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

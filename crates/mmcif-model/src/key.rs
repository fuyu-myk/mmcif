use std::fmt;

use crate::ids::{AuthAsymId, AuthSeqId, CompId, InsCode, LabelAsymId, LabelSeqId, OperatorId};

/// A residue address in author numbering: author chain, residue number and
/// insertion code.
///
/// It displays as `chain:number` followed by the insertion code, if any.
///
/// # Examples
///
/// ```
/// use mmcif_model::{AuthResidueKey, AuthSeqId, InsCode};
///
/// let key = AuthResidueKey {
///     chain: "H".into(),
///     seq: AuthSeqId::new(100),
///     ins_code: InsCode::new('A'),
/// };
/// assert_eq!(key.to_string(), "H:100A");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AuthResidueKey {
    /// The author chain id.
    pub chain: AuthAsymId,
    /// The author residue number.
    pub seq: AuthSeqId,
    /// The insertion code.
    pub ins_code: InsCode,
}

impl fmt::Display for AuthResidueKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}{}", self.chain, self.seq, self.ins_code)
    }
}

/// A polymer residue address in label numbering: label chain and label
/// residue number.
///
/// It displays as `chain/number`.
///
/// # Examples
///
/// ```
/// use mmcif_model::{LabelResidueKey, LabelSeqId};
///
/// let key = LabelResidueKey {
///     asym: "B".into(),
///     seq: LabelSeqId::new(105).unwrap(),
/// };
/// assert_eq!(key.to_string(), "B/105");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct LabelResidueKey {
    /// The label chain id.
    pub asym: LabelAsymId,
    /// The label residue number.
    pub seq: LabelSeqId,
}

impl fmt::Display for LabelResidueKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/{}", self.asym, self.seq)
    }
}

/// The full address of a residue in a model: label chain, label and author
/// residue numbers, insertion code, component, and the assembly operators
/// that placed its chain.
///
/// It displays as `chain/label/author/component`, with `.` for a missing
/// label number and the insertion code after the author number. A residue of
/// an assembled chain adds `@` and its operators, separated by commas, as in
/// `AA/2/100/SER@2`. Keys come from [`ResidueRef::key`](crate::ResidueRef::key).
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let model = structure.first_model().unwrap();
///
/// let alanine = model.chain("AA").unwrap().residues().nth(2).unwrap();
/// let key = alanine.key();
/// assert_eq!(key.to_string(), "AA/3/100A/ALA");
/// assert_eq!(model.residue(&key), Some(alanine));
///
/// let sugar = model.chain("B").unwrap().residues().next().unwrap();
/// assert_eq!(sugar.key().to_string(), "B/./1/NAG");
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[non_exhaustive]
pub struct ResidueKey {
    /// The label chain id.
    pub asym: LabelAsymId,
    /// The label residue number, if the residue has one.
    pub label_seq: Option<LabelSeqId>,
    /// The author residue number.
    pub auth_seq: AuthSeqId,
    /// The insertion code.
    pub ins_code: InsCode,
    /// The chemical component id.
    pub comp: CompId,
    /// The assembly operators that placed the residue's chain, in the order
    /// [`ChainRef::operators`](crate::ChainRef::operators) gives them. Empty in
    /// a deposited model.
    pub operators: Box<[OperatorId]>,
}

impl fmt::Display for ResidueKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/", self.asym)?;
        match self.label_seq {
            Some(seq) => write!(f, "{seq}")?,
            None => f.write_str(".")?,
        }
        write!(f, "/{}{}/{}", self.auth_seq, self.ins_code, self.comp)?;
        for (index, operator) in self.operators.iter().enumerate() {
            f.write_str(if index == 0 { "@" } else { "," })?;
            write!(f, "{operator}")?;
        }
        Ok(())
    }
}

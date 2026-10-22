use std::collections::{HashMap, HashSet};

use crate::ids::{AuthAsymId, LabelAsymId};
use crate::key::{AuthResidueKey, LabelResidueKey};

/// An error from converting between label and author numbering with
/// [`Numbering`].
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum NumberingError {
    /// No chain has this label chain id.
    #[error("no chain with label_asym_id {0}")]
    UnknownLabelChain(LabelAsymId),
    /// No chain has this author chain id.
    #[error("no chain with auth_asym_id {0}")]
    UnknownAuthChain(AuthAsymId),
    /// The label chain exists but has no polymer position with this number.
    #[error("{0} is not a position in its polymer sequence")]
    NotInSequence(LabelResidueKey),
    /// The author address belongs to a residue outside any polymer.
    #[error("{0} is not a polymer residue and has no label_seq_id")]
    NotPolymer(AuthResidueKey),
    /// The author chain exists but has no residue at this address.
    #[error("no residue at {0}")]
    NotFound(AuthResidueKey),
    /// The author address belongs to a polymer residue in a file with no
    /// label numbering.
    #[error("{0} is a polymer residue but the file carries no label numbering")]
    Unnumbered(AuthResidueKey),
    /// The file gives no author numbering for this label position.
    #[error("{0} has no author numbering in this file")]
    NoAuthorNumbering(LabelResidueKey),
    /// The author address maps to more than one label position.
    #[error("{key} maps to {} label positions", candidates.len())]
    Ambiguous {
        /// The author address.
        key: AuthResidueKey,
        /// Every label position the address maps to, in the order read.
        candidates: Vec<LabelResidueKey>,
    },
}

/// The mapping between the label and author residue numbering of a
/// structure.
///
/// For mmCIF it is read from `_pdbx_poly_seq_scheme`, `_pdbx_branch_scheme`
/// and the residues of the first model. When a label position is given more
/// than once, the first author address read for it is kept. The residues of
/// a PDB file have author addresses only.
///
/// # Examples
///
/// ```
/// use mmcif_model::{AuthResidueKey, AuthSeqId, InsCode, LabelResidueKey, LabelSeqId};
///
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let numbering = structure.numbering();
///
/// let label = LabelResidueKey {
///     asym: "AA".into(),
///     seq: LabelSeqId::new(3).unwrap(),
/// };
/// let auth = AuthResidueKey {
///     chain: "H1".into(),
///     seq: AuthSeqId::new(100),
///     ins_code: InsCode::new('A'),
/// };
/// assert_eq!(numbering.to_auth(&label)?, auth);
/// assert_eq!(numbering.to_label(&auth)?, label);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Numbering {
    label_to_auth: HashMap<LabelResidueKey, AuthResidueKey>,
    auth_to_label: HashMap<AuthResidueKey, Vec<LabelResidueKey>>,
    label_chains: HashSet<LabelAsymId>,
    auth_chains: HashSet<AuthAsymId>,
    non_polymer: HashSet<AuthResidueKey>,
    unnumbered: HashSet<AuthResidueKey>,
    label_only: HashSet<LabelAsymId>,
    branch: HashMap<(LabelAsymId, u32), AuthResidueKey>,
}

impl Numbering {
    pub(crate) fn insert_polymer(&mut self, label: LabelResidueKey, auth: AuthResidueKey) {
        self.label_chains.insert(label.asym.clone());
        self.auth_chains.insert(auth.chain.clone());
        let candidates = self.auth_to_label.entry(auth.clone()).or_default();
        if !candidates.contains(&label) {
            candidates.push(label.clone());
        }
        self.label_to_auth.entry(label).or_insert(auth);
    }

    pub(crate) fn insert_non_polymer(&mut self, label_asym: LabelAsymId, auth: AuthResidueKey) {
        self.label_chains.insert(label_asym);
        self.auth_chains.insert(auth.chain.clone());
        self.non_polymer.insert(auth);
    }

    pub(crate) fn insert_unnumbered_polymer(
        &mut self,
        label_asym: LabelAsymId,
        auth: AuthResidueKey,
    ) {
        self.label_chains.insert(label_asym);
        self.auth_chains.insert(auth.chain.clone());
        self.unnumbered.insert(auth);
    }

    pub(crate) fn insert_label_only(&mut self, label_asym: LabelAsymId) {
        self.label_chains.insert(label_asym.clone());
        self.label_only.insert(label_asym);
    }

    pub(crate) fn insert_branch(&mut self, asym: LabelAsymId, num: u32, auth: AuthResidueKey) {
        self.branch.insert((asym, num), auth);
    }

    /// Returns the author address of monomer `num` of branched label chain
    /// `asym`, from `_pdbx_branch_scheme`, or `None` if the file gives none.
    ///
    /// The insertion code of the address is always [`InsCode::NONE`].
    ///
    /// [`InsCode::NONE`]: crate::InsCode::NONE
    ///
    /// # Examples
    ///
    /// ```
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let numbering = structure.numbering();
    ///
    /// let auth = numbering.branch_to_auth(&"B".into(), 2).unwrap();
    /// assert_eq!(auth.to_string(), "H1:2");
    /// assert!(numbering.branch_to_auth(&"B".into(), 3).is_none());
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn branch_to_auth(&self, asym: &LabelAsymId, num: u32) -> Option<&AuthResidueKey> {
        self.branch.get(&(asym.clone(), num))
    }

    /// Returns the author address of the polymer residue at label position
    /// `key`.
    ///
    /// # Errors
    ///
    /// Returns [`NumberingError::NoAuthorNumbering`] if the file gives no
    /// author numbering for `key`, [`NumberingError::NotInSequence`] if the
    /// label chain exists but has no polymer position at `key`, and
    /// [`NumberingError::UnknownLabelChain`] if no chain has the label chain
    /// id of `key`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::{LabelResidueKey, LabelSeqId, NumberingError};
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let numbering = structure.numbering();
    /// let label = |seq| LabelResidueKey {
    ///     asym: "AA".into(),
    ///     seq: LabelSeqId::new(seq).unwrap(),
    /// };
    ///
    /// assert_eq!(numbering.to_auth(&label(1))?.to_string(), "H1:99");
    /// assert_eq!(
    ///     numbering.to_auth(&label(7)),
    ///     Err(NumberingError::NotInSequence(label(7)))
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn to_auth(&self, key: &LabelResidueKey) -> Result<AuthResidueKey, NumberingError> {
        if let Some(auth) = self.label_to_auth.get(key) {
            return Ok(auth.clone());
        }
        if self.label_only.contains(&key.asym) {
            Err(NumberingError::NoAuthorNumbering(key.clone()))
        } else if self.label_chains.contains(&key.asym) {
            Err(NumberingError::NotInSequence(key.clone()))
        } else {
            Err(NumberingError::UnknownLabelChain(key.asym.clone()))
        }
    }

    /// Returns the label position of the polymer residue at author address
    /// `key`.
    ///
    /// # Errors
    ///
    /// Returns [`NumberingError::Ambiguous`] if `key` maps to more than one
    /// label position, [`NumberingError::NotPolymer`] if `key` belongs to a
    /// residue outside any polymer, [`NumberingError::Unnumbered`] if `key`
    /// belongs to a polymer residue in a file with no label numbering, such
    /// as a PDB file, [`NumberingError::NotFound`] if the author chain exists
    /// but has no residue at `key`, and [`NumberingError::UnknownAuthChain`]
    /// if no chain has the author chain id of `key`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::{AuthResidueKey, AuthSeqId, InsCode, NumberingError};
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let numbering = structure.numbering();
    /// let auth = |seq| AuthResidueKey {
    ///     chain: "H1".into(),
    ///     seq: AuthSeqId::new(seq),
    ///     ins_code: InsCode::NONE,
    /// };
    ///
    /// assert_eq!(numbering.to_label(&auth(101))?.to_string(), "AA/4");
    /// assert_eq!(
    ///     numbering.to_label(&auth(201)),
    ///     Err(NumberingError::NotPolymer(auth(201)))
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn to_label(&self, key: &AuthResidueKey) -> Result<LabelResidueKey, NumberingError> {
        match self.auth_to_label.get(key).map(Vec::as_slice) {
            Some([single]) => Ok(single.clone()),
            Some(candidates) => Err(NumberingError::Ambiguous {
                key: key.clone(),
                candidates: candidates.to_vec(),
            }),
            None if self.non_polymer.contains(key) => Err(NumberingError::NotPolymer(key.clone())),
            None if self.unnumbered.contains(key) => Err(NumberingError::Unnumbered(key.clone())),
            None if self.auth_chains.contains(&key.chain) => {
                Err(NumberingError::NotFound(key.clone()))
            }
            None => Err(NumberingError::UnknownAuthChain(key.chain.clone())),
        }
    }

    /// Returns every label position with an author address, paired with that
    /// address, in no particular order.
    ///
    /// # Examples
    ///
    /// ```
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let numbering = structure.numbering();
    ///
    /// assert_eq!(numbering.polymer_positions().count(), 6);
    /// for (label, auth) in numbering.polymer_positions() {
    ///     assert_eq!(numbering.to_auth(label).as_ref(), Ok(auth));
    /// }
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn polymer_positions(&self) -> impl Iterator<Item = (&LabelResidueKey, &AuthResidueKey)> {
        self.label_to_auth.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ids::{AuthSeqId, InsCode, LabelSeqId};

    fn label(asym: &str, seq: u32) -> LabelResidueKey {
        LabelResidueKey {
            asym: asym.into(),
            seq: LabelSeqId::new(seq).unwrap(),
        }
    }

    fn auth(chain: &str, seq: i32, ins: Option<char>) -> AuthResidueKey {
        AuthResidueKey {
            chain: chain.into(),
            seq: AuthSeqId::new(seq),
            ins_code: ins.map_or(InsCode::NONE, InsCode::new),
        }
    }

    fn numbering() -> Numbering {
        let mut numbering = Numbering::default();
        numbering.insert_polymer(label("A", 1), auth("H", 100, None));
        numbering.insert_polymer(label("A", 2), auth("H", 100, Some('A')));
        numbering.insert_polymer(label("A", 3), auth("H", 101, None));
        numbering.insert_non_polymer("B".into(), auth("H", 201, None));
        numbering
    }

    #[test]
    fn converts_both_ways_with_insertion_codes() {
        let numbering = numbering();
        assert_eq!(
            numbering.to_auth(&label("A", 2)).unwrap(),
            auth("H", 100, Some('A'))
        );
        assert_eq!(
            numbering.to_label(&auth("H", 100, Some('A'))).unwrap(),
            label("A", 2)
        );
        assert_eq!(
            numbering.to_label(&auth("H", 100, None)).unwrap(),
            label("A", 1)
        );
    }

    #[test]
    fn fails_cleanly_where_undefined() {
        let numbering = numbering();
        assert_eq!(
            numbering.to_auth(&label("A", 9)),
            Err(NumberingError::NotInSequence(label("A", 9)))
        );
        assert_eq!(
            numbering.to_auth(&label("Z", 1)),
            Err(NumberingError::UnknownLabelChain("Z".into()))
        );
        assert_eq!(
            numbering.to_label(&auth("H", 201, None)),
            Err(NumberingError::NotPolymer(auth("H", 201, None)))
        );
        assert_eq!(
            numbering.to_label(&auth("H", 999, None)),
            Err(NumberingError::NotFound(auth("H", 999, None)))
        );
        assert_eq!(
            numbering.to_label(&auth("Q", 1, None)),
            Err(NumberingError::UnknownAuthChain("Q".into()))
        );
    }

    #[test]
    fn duplicate_author_numbers_are_ambiguous() {
        let mut numbering = numbering();
        numbering.insert_polymer(label("A", 4), auth("H", 101, None));
        assert!(matches!(
            numbering.to_label(&auth("H", 101, None)),
            Err(NumberingError::Ambiguous { candidates, .. }) if candidates.len() == 2
        ));
    }
}

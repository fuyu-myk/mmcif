use std::collections::{HashMap, HashSet};

use crate::ids::{AuthAsymId, LabelAsymId};
use crate::key::{AuthResidueKey, LabelResidueKey};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum NumberingError {
    #[error("no chain with label_asym_id {0}")]
    UnknownLabelChain(LabelAsymId),
    #[error("no chain with auth_asym_id {0}")]
    UnknownAuthChain(AuthAsymId),
    #[error("{0} is not a position in its polymer sequence")]
    NotInSequence(LabelResidueKey),
    #[error("{0} is not a polymer residue and has no label_seq_id")]
    NotPolymer(AuthResidueKey),
    #[error("no residue at {0}")]
    NotFound(AuthResidueKey),
    #[error("{0} is a polymer residue but the file carries no label numbering")]
    Unnumbered(AuthResidueKey),
    #[error("{key} maps to {} label positions", candidates.len())]
    Ambiguous {
        key: AuthResidueKey,
        candidates: Vec<LabelResidueKey>,
    },
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Numbering {
    label_to_auth: HashMap<LabelResidueKey, AuthResidueKey>,
    auth_to_label: HashMap<AuthResidueKey, Vec<LabelResidueKey>>,
    label_chains: HashSet<LabelAsymId>,
    auth_chains: HashSet<AuthAsymId>,
    non_polymer: HashSet<AuthResidueKey>,
    unnumbered: HashSet<AuthResidueKey>,
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

    pub(crate) fn insert_branch(&mut self, asym: LabelAsymId, num: u32, auth: AuthResidueKey) {
        self.branch.insert((asym, num), auth);
    }

    pub fn branch_to_auth(&self, asym: &LabelAsymId, num: u32) -> Option<&AuthResidueKey> {
        self.branch.get(&(asym.clone(), num))
    }

    pub fn to_auth(&self, key: &LabelResidueKey) -> Result<AuthResidueKey, NumberingError> {
        if let Some(auth) = self.label_to_auth.get(key) {
            return Ok(auth.clone());
        }
        if self.label_chains.contains(&key.asym) {
            Err(NumberingError::NotInSequence(key.clone()))
        } else {
            Err(NumberingError::UnknownLabelChain(key.asym.clone()))
        }
    }

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

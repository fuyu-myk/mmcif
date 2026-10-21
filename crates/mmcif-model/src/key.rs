use std::fmt;

use crate::ids::{AuthAsymId, AuthSeqId, CompId, InsCode, LabelAsymId, LabelSeqId};

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AuthResidueKey {
    pub chain: AuthAsymId,
    pub seq: AuthSeqId,
    pub ins_code: InsCode,
}

impl fmt::Display for AuthResidueKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}{}", self.chain, self.seq, self.ins_code)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct LabelResidueKey {
    pub asym: LabelAsymId,
    pub seq: LabelSeqId,
}

impl fmt::Display for LabelResidueKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/{}", self.asym, self.seq)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ResidueKey {
    pub asym: LabelAsymId,
    pub label_seq: Option<LabelSeqId>,
    pub auth_seq: AuthSeqId,
    pub ins_code: InsCode,
    pub comp: CompId,
}

impl fmt::Display for ResidueKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/", self.asym)?;
        match self.label_seq {
            Some(seq) => write!(f, "{seq}")?,
            None => f.write_str(".")?,
        }
        write!(f, "/{}{}/{}", self.auth_seq, self.ins_code, self.comp)
    }
}

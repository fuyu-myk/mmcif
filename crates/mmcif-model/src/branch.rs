use crate::entity::{BranchAtom, BranchLink};
use crate::ids::{EntityId, LabelAsymId};
use crate::model::Model;
use crate::structure::Structure;
use crate::view::{AtomRef, ChainRef, ResidueRef};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum BranchError {
    #[error("no chain {0}")]
    UnknownChain(LabelAsymId),
    #[error("chain {0} has no branched entity")]
    NotBranched(LabelAsymId),
    #[error("chain {asym} has no author numbering for branch monomer {num}")]
    MissingBranchNumbering { asym: LabelAsymId, num: u32 },
    #[error("entity {0} is not defined")]
    UnknownEntity(EntityId),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinkedAtom<'a> {
    pub residue: ResidueRef<'a>,
    pub atom: AtomRef<'a>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BranchBond<'a> {
    pub link: &'a BranchLink,
    pub first: Option<LinkedAtom<'a>>,
    pub second: Option<LinkedAtom<'a>>,
}

impl Structure {
    pub fn branch_bonds<'a>(
        &'a self,
        model: &'a Model,
        asym: &str,
    ) -> Result<Vec<BranchBond<'a>>, BranchError> {
        let chain = model
            .chain(asym)
            .ok_or_else(|| BranchError::UnknownChain(LabelAsymId::new(asym)))?;
        let entity_id = chain
            .entity()
            .ok_or_else(|| BranchError::NotBranched(chain.label_asym().clone()))?;
        let entity = self
            .entity(entity_id.as_str())
            .ok_or_else(|| BranchError::UnknownEntity(entity_id.clone()))?;
        let branch = entity
            .branch()
            .ok_or_else(|| BranchError::NotBranched(chain.label_asym().clone()))?;
        branch
            .links()
            .iter()
            .map(|link| {
                Ok(BranchBond {
                    link,
                    first: self.locate(chain, &link.first)?,
                    second: self.locate(chain, &link.second)?,
                })
            })
            .collect()
    }

    fn locate<'a>(
        &self,
        chain: ChainRef<'a>,
        target: &BranchAtom,
    ) -> Result<Option<LinkedAtom<'a>>, BranchError> {
        let auth = self
            .numbering
            .branch_to_auth(chain.label_asym(), target.num)
            .ok_or_else(|| BranchError::MissingBranchNumbering {
                asym: chain.label_asym().clone(),
                num: target.num,
            })?;
        Ok(chain
            .residues_at_auth(auth.seq, auth.ins_code)
            .filter(|residue| residue.comp() == &target.comp)
            .find_map(|residue| {
                residue
                    .atom(target.atom.as_str())
                    .map(|atom| LinkedAtom { residue, atom })
            }))
    }
}

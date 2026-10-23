use crate::entity::{BranchAtom, BranchLink};
use crate::ids::{EntityId, LabelAsymId};
use crate::model::Model;
use crate::structure::Structure;
use crate::view::{AtomRef, ChainRef, ResidueRef};

/// An error from [`Structure::branch_bonds`].
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum BranchError {
    /// The model has no chain with this label id and no operators applied.
    #[error("no chain {0}")]
    UnknownChain(LabelAsymId),
    /// The chain has no entity, or its entity is not branched.
    #[error("chain {0} has no branched entity")]
    NotBranched(LabelAsymId),
    /// `_pdbx_branch_scheme` gives no author numbering for a monomer of the
    /// chain.
    #[error("chain {asym} has no author numbering for branch monomer {num}")]
    MissingBranchNumbering {
        /// The label chain id.
        asym: LabelAsymId,
        /// The monomer number.
        num: u32,
    },
    /// The chain names an entity that the structure does not define.
    #[error("entity {0} is not defined")]
    UnknownEntity(EntityId),
}

/// An atom of a branched chain, together with its residue.
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let model = structure.first_model().unwrap();
/// let bonds = structure.branch_bonds(model, "B")?;
/// let linked = bonds[0].second.unwrap();
///
/// assert_eq!(linked.residue.comp().as_str(), "NAG");
/// assert_eq!(linked.atom.name().as_str(), "O4");
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinkedAtom<'a> {
    /// The residue that holds the atom.
    pub residue: ResidueRef<'a>,
    /// The atom.
    pub atom: AtomRef<'a>,
}

/// A bond of a branched entity, resolved to atoms of a model.
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let model = structure.first_model().unwrap();
/// let bond = structure.branch_bonds(model, "B")?[0];
///
/// assert_eq!(bond.link.order(), Some("sing"));
/// assert_eq!(bond.first.unwrap().atom.name().as_str(), "C1");
/// assert_eq!(bond.second.unwrap().atom.name().as_str(), "O4");
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BranchBond<'a> {
    /// The bond as the entity defines it.
    pub link: &'a BranchLink,
    /// The atom at [`BranchLink::first`], or `None` if the model does not
    /// have it.
    pub first: Option<LinkedAtom<'a>>,
    /// The atom at [`BranchLink::second`], or `None` if the model does not
    /// have it.
    pub second: Option<LinkedAtom<'a>>,
}

impl Structure {
    /// Returns the bonds of the branched entity of chain `asym`, resolved to
    /// atoms of `model`.
    ///
    /// The bonds are in the order of [`Branch::links`]. Each end is found
    /// through the chain's `_pdbx_branch_scheme` numbering and must match
    /// the component and atom name the bond gives.
    ///
    /// [`Branch::links`]: crate::Branch::links
    ///
    /// # Errors
    ///
    /// Returns [`BranchError::UnknownChain`] if `model` has no chain `asym`
    /// with no operators applied, [`BranchError::NotBranched`] if the chain
    /// has no entity or its entity is not branched,
    /// [`BranchError::UnknownEntity`] if the structure does not define the
    /// chain's entity, and [`BranchError::MissingBranchNumbering`] if a
    /// monomer of a bond has no `_pdbx_branch_scheme` row for the chain.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::BranchError;
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let model = structure.first_model().unwrap();
    ///
    /// let bonds = structure.branch_bonds(model, "B")?;
    /// let first = bonds[0].first.unwrap();
    /// assert_eq!(first.residue.comp().as_str(), "BMA");
    /// assert_eq!(first.atom.name().as_str(), "C1");
    ///
    /// assert_eq!(
    ///     structure.branch_bonds(model, "AA"),
    ///     Err(BranchError::NotBranched("AA".into()))
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
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
            .branch_to_auth(chain.label_asym().as_str(), target.num)
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

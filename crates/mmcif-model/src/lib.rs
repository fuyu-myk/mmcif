//! A typed model of macromolecular structures read from mmCIF, BinaryCIF and
//! PDB files.
//!
//! A [`Structure`] holds the models, entities, assemblies and residue
//! numbering of one entry. Each [`Model`] holds chains, residues and atoms,
//! reached through the [`ChainRef`], [`ResidueRef`] and [`AtomRef`] views.
//! Author addresses ([`AuthResidueKey`]) and label addresses
//! ([`LabelResidueKey`]) are separate types, converted with [`Numbering`].
//!
//! # Examples
//!
//! ```
//! use mmcif_model::{AuthResidueKey, AuthSeqId, InsCode};
//!
//! # let text = include_str!("../tests/fixtures/pathologies.cif");
//! let structure = mmcif_model::parse_mmcif(text)?;
//! let model = structure.first_model().unwrap();
//!
//! let chain = model.chain("AA").unwrap();
//! assert_eq!(chain.auth_asym().as_str(), "H1");
//! assert_eq!(chain.residue_count(), 5);
//!
//! let key = AuthResidueKey {
//!     chain: "H1".into(),
//!     seq: AuthSeqId::new(100),
//!     ins_code: InsCode::new('A'),
//! };
//! assert_eq!(structure.numbering().to_label(&key)?.to_string(), "AA/3");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

#![warn(missing_docs)]

mod altloc;
mod assembly;
mod branch;
mod entity;
mod error;
mod expression;
mod handle;
mod ids;
mod key;
mod model;
mod numbering;
mod pdb;
mod read;
mod remark350;
mod structure;
mod view;

pub use altloc::AltlocPolicy;
pub use assembly::{Assembly, AssemblyError, Generator, Operator, Transform};
pub use branch::{BranchBond, BranchError, LinkedAtom};
pub use entity::{
    Branch, BranchAtom, BranchLink, BranchMonomer, Entity, EntityKind, Polymer, PolymerType,
    SequencePosition,
};
pub use error::ReadError;
pub use expression::{InvalidExpression, OperatorExpression};
pub use handle::{AtomHandle, ChainHandle, ResidueHandle};
pub use ids::{
    AltLoc, AssemblyId, AtomName, AuthAsymId, AuthSeqId, CompId, Element, EntityId, InsCode,
    LabelAsymId, LabelSeqId, ModelNumber, OperatorId,
};
pub use key::{AuthResidueKey, LabelResidueKey, ResidueKey};
pub use model::{AtomIndex, ChainIndex, Model, ResidueIndex};
pub use numbering::{Numbering, NumberingError};
pub use pdb::{parse_pdb, read_pdb};
pub use read::{
    Format, Loaded, ReadOptions, parse_bcif, parse_mmcif, read, read_bcif, read_mmcif, read_with,
};
pub use structure::Structure;
pub use view::{AtomRef, ChainRef, ResidueRef};

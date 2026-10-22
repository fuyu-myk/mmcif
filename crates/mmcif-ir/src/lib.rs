//! The shape of a design problem. Every type here may change until 1.0.
//!
//! A [`Spec`] holds segments of slots, named groups of slots and
//! constraints between them, all addressed symbolically. Instantiating it
//! with chosen [`Lengths`] numbers every slot, resolves each constraint to
//! slot indices and gives each slot a [`Determinacy`].
//!
//! # Examples
//!
//! ```
//! use mmcif_ir::{
//!     Corner, Element, Length, Lengths, ResidueRef, Segment, SegmentId, Slot, SlotIndex, Spec,
//! };
//! use mmcif_model::{CompId, LabelResidueKey, LabelSeqId};
//!
//! let anchor = ResidueRef::Label(LabelResidueKey {
//!     asym: "A".into(),
//!     seq: LabelSeqId::new(10).unwrap(),
//! });
//!
//! let mut spec = Spec::new();
//! spec.push_segment(Segment::new(
//!     SegmentId::new("target"),
//!     [Element::Slot(Slot::context(CompId::new("GLY"), anchor))],
//! ))?;
//! spec.push_segment(Segment::new(
//!     SegmentId::new("binder"),
//!     [Element::run(Slot::de_novo(), Length::range(8, 12)?)],
//! ))?;
//!
//! let design = spec.instantiate(&Lengths::new(&spec, [10])?)?;
//!
//! assert_eq!(design.slot_count(), 11);
//! assert_eq!(design.determinacy(SlotIndex::new(0))?.corner(), Some(Corner::Context));
//! assert_eq!(design.determinacy(SlotIndex::new(1))?.corner(), Some(Corner::DeNovo));
//! # Ok::<(), mmcif_ir::IrError>(())
//! ```

#![warn(missing_docs)]

mod constraint;
mod error;
mod ids;
mod instantiate;
mod segment;
mod slot;
mod spec;

pub use constraint::{
    Addressing, AtomAddress, Bounds, Constraint, Resolved, Scope, SecondaryStructure, Symbolic,
    SymbolicAddress,
};
pub use error::IrError;
pub use ids::{
    ElementIndex, GroupId, ResidueRef, SegmentId, SegmentIndex, SlotIndex, SlotSpan, TokenIndex,
    TokenLayout,
};
pub use instantiate::{Instantiated, Lengths, MAX_SLOTS, Origin};
pub use segment::{Element, Length, Segment};
pub use slot::{AlphabetSet, Corner, Determinacy, Geometry, Identity, Level, Slot};
pub use spec::Spec;

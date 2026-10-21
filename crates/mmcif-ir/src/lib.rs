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
    ElementIndex, GroupId, ResidueRef, SegmentId, SegmentIndex, SlotIndex, TokenIndex, TokenLayout,
};
pub use instantiate::{Instantiated, Lengths, Origin};
pub use segment::{Element, Length, Segment};
pub use slot::{AlphabetSet, Corner, Determinacy, Geometry, Identity, Level, Slot};
pub use spec::Spec;

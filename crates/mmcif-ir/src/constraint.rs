use std::fmt;

use mmcif_model::AtomName;

use crate::error::IrError;
use crate::ids::{ElementIndex, GroupId, SegmentId, SlotIndex};

pub trait Addressing {
    type One: Clone + fmt::Debug + PartialEq;
    type Many: Clone + fmt::Debug + PartialEq;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Symbolic;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Resolved;

impl Addressing for Symbolic {
    type One = SymbolicAddress;
    type Many = Scope;
}

impl Addressing for Resolved {
    type One = SlotIndex;
    type Many = Vec<SlotIndex>;
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SymbolicAddress {
    pub segment: SegmentId,
    pub element: ElementIndex,
    pub offset: Option<u32>,
}

impl SymbolicAddress {
    pub fn whole(segment: SegmentId, element: ElementIndex) -> SymbolicAddress {
        SymbolicAddress {
            segment,
            element,
            offset: None,
        }
    }

    pub fn at(segment: SegmentId, element: ElementIndex, offset: u32) -> SymbolicAddress {
        SymbolicAddress {
            segment,
            element,
            offset: Some(offset),
        }
    }
}

impl fmt::Display for SymbolicAddress {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}[{}]", self.segment, self.element)?;
        match self.offset {
            Some(offset) => write!(f, "+{offset}"),
            None => Ok(()),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Scope {
    Address(SymbolicAddress),
    Group(GroupId),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomAddress<T> {
    pub slot: T,
    pub atom: AtomName,
}

impl<T> AtomAddress<T> {
    pub fn new(slot: T, atom: AtomName) -> AtomAddress<T> {
        AtomAddress { slot, atom }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds {
    min: f64,
    max: f64,
}

impl Bounds {
    pub fn new(min: f64, max: f64) -> Result<Bounds, IrError> {
        if !min.is_finite() || !max.is_finite() || min < 0.0 || min > max {
            return Err(IrError::InvalidBounds { min, max });
        }
        Ok(Bounds { min, max })
    }

    pub fn at_most(max: f64) -> Result<Bounds, IrError> {
        Bounds::new(0.0, max)
    }

    pub fn min(self) -> f64 {
        self.min
    }

    pub fn max(self) -> f64 {
        self.max
    }

    pub fn contains(self, value: f64) -> bool {
        self.min <= value && value <= self.max
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum SecondaryStructure {
    Helix,
    Strand,
    Loop,
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Constraint<A: Addressing> {
    Distance {
        first: AtomAddress<A::One>,
        second: AtomAddress<A::One>,
        bounds: Bounds,
    },
    SecondaryStructure {
        scope: A::Many,
        kind: SecondaryStructure,
    },
}

impl Constraint<Resolved> {
    pub fn slots(&self) -> Vec<SlotIndex> {
        match self {
            Constraint::Distance { first, second, .. } => vec![first.slot, second.slot],
            Constraint::SecondaryStructure { scope, .. } => scope.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounds_reject_descending_negative_and_non_finite_pairs() {
        assert!(Bounds::new(2.0, 8.0).is_ok());
        assert!(Bounds::at_most(4.5).unwrap().contains(0.0));
        assert!(matches!(
            Bounds::new(8.0, 2.0),
            Err(IrError::InvalidBounds { .. })
        ));
        assert!(matches!(
            Bounds::new(-1.0, 2.0),
            Err(IrError::InvalidBounds { .. })
        ));
        assert!(matches!(
            Bounds::new(0.0, f64::NAN),
            Err(IrError::InvalidBounds { .. })
        ));
        assert!(matches!(
            Bounds::new(0.0, f64::INFINITY),
            Err(IrError::InvalidBounds { .. })
        ));
    }

    #[test]
    fn bounds_include_their_endpoints() {
        let bounds = Bounds::new(2.0, 8.0).unwrap();
        assert!(bounds.contains(2.0) && bounds.contains(8.0));
        assert!(!bounds.contains(1.999) && !bounds.contains(8.001));
        assert_eq!((bounds.min(), bounds.max()), (2.0, 8.0));
    }

    #[test]
    fn symbolic_addresses_print_their_position() {
        let whole = SymbolicAddress::whole(SegmentId::new("binder"), ElementIndex::new(2));
        assert_eq!(whole.to_string(), "binder[2]");
        assert_eq!(whole.offset, None);
        let offset = SymbolicAddress::at(SegmentId::new("binder"), ElementIndex::new(2), 5);
        assert_eq!(offset.to_string(), "binder[2]+5");
    }

    #[test]
    fn resolved_constraints_list_the_slots_they_touch() {
        let distance = Constraint::<Resolved>::Distance {
            first: AtomAddress::new(SlotIndex::new(3), AtomName::new("SG")),
            second: AtomAddress::new(SlotIndex::new(9), AtomName::new("SG")),
            bounds: Bounds::new(1.8, 2.2).unwrap(),
        };
        assert_eq!(distance.slots(), [SlotIndex::new(3), SlotIndex::new(9)]);
        let helix = Constraint::<Resolved>::SecondaryStructure {
            scope: vec![SlotIndex::new(1), SlotIndex::new(2)],
            kind: SecondaryStructure::Helix,
        };
        assert_eq!(helix.slots(), [SlotIndex::new(1), SlotIndex::new(2)]);
    }
}

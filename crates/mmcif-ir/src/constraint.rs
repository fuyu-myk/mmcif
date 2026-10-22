use std::fmt;

use mmcif_model::AtomName;

use crate::error::IrError;
use crate::ids::{ElementIndex, GroupId, SegmentId, SlotIndex};

/// The way a [`Constraint`] names its slots.
///
/// [`Symbolic`] names slots by their place in a [`Spec`](crate::Spec), and
/// [`Resolved`] names them by [`SlotIndex`] after instantiation.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{Addressing, Constraint, Resolved, SecondaryStructure, SlotIndex};
///
/// fn describe<A: Addressing>(constraint: &Constraint<A>) -> &'static str {
///     match constraint {
///         Constraint::Bond { .. } => "bond",
///         Constraint::Distance { .. } => "distance",
///         Constraint::SecondaryStructure { .. } => "secondary structure",
///         _ => "other",
///     }
/// }
///
/// let helix = Constraint::<Resolved>::SecondaryStructure {
///     scope: vec![SlotIndex::new(0), SlotIndex::new(1)],
///     kind: SecondaryStructure::Helix,
/// };
///
/// assert_eq!(describe(&helix), "secondary structure");
/// ```
pub trait Addressing {
    /// The type that names one slot.
    type One: Clone + fmt::Debug + PartialEq;
    /// The type that names a set of slots.
    type Many: Clone + fmt::Debug + PartialEq;
}

/// Addressing by [`SymbolicAddress`] and [`Scope`], used by the constraints
/// of a [`Spec`](crate::Spec).
///
/// # Examples
///
/// ```
/// use mmcif_ir::{Constraint, GroupId, Scope, SecondaryStructure, Symbolic};
///
/// let helix = Constraint::<Symbolic>::SecondaryStructure {
///     scope: Scope::Group(GroupId::new("core")),
///     kind: SecondaryStructure::Helix,
/// };
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Symbolic;

/// Addressing by [`SlotIndex`], used by the constraints of an
/// [`Instantiated`](crate::Instantiated) specification.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{AtomAddress, Constraint, Resolved, SlotIndex};
/// use mmcif_model::AtomName;
///
/// let bond = Constraint::<Resolved>::Bond {
///     first: AtomAddress::new(SlotIndex::new(4), AtomName::new("SG")),
///     second: AtomAddress::new(SlotIndex::new(9), AtomName::new("SG")),
/// };
///
/// assert_eq!(bond.slots(), [SlotIndex::new(4), SlotIndex::new(9)]);
/// ```
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

/// A place in a [`Spec`](crate::Spec): a whole element, or the slot at an
/// offset within an element.
///
/// An address displays as `segment[element]`, followed by `+offset` when it
/// has an offset.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{ElementIndex, SegmentId, SymbolicAddress};
///
/// let run = SymbolicAddress::whole(SegmentId::new("binder"), ElementIndex::new(2));
/// let slot = SymbolicAddress::at(SegmentId::new("binder"), ElementIndex::new(2), 3);
///
/// assert_eq!(run.to_string(), "binder[2]");
/// assert_eq!(slot.to_string(), "binder[2]+3");
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SymbolicAddress {
    /// The segment the address points into.
    pub segment: SegmentId,
    /// The element within the segment.
    pub element: ElementIndex,
    /// The offset of one slot within the element, or `None` for the whole
    /// element.
    pub offset: Option<u32>,
}

impl SymbolicAddress {
    /// Creates an address of a whole element.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{ElementIndex, SegmentId, SymbolicAddress};
    ///
    /// let run = SymbolicAddress::whole(SegmentId::new("binder"), ElementIndex::new(0));
    ///
    /// assert_eq!(run.offset, None);
    /// ```
    pub fn whole(segment: SegmentId, element: ElementIndex) -> SymbolicAddress {
        SymbolicAddress {
            segment,
            element,
            offset: None,
        }
    }

    /// Creates an address of the slot at `offset` within an element.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{ElementIndex, SegmentId, SymbolicAddress};
    ///
    /// let slot = SymbolicAddress::at(SegmentId::new("binder"), ElementIndex::new(0), 3);
    ///
    /// assert_eq!(slot.offset, Some(3));
    /// ```
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

/// The slots a symbolic secondary structure constraint covers.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{
///     Constraint, ElementIndex, Scope, SecondaryStructure, SegmentId, Symbolic, SymbolicAddress,
/// };
///
/// let run = SymbolicAddress::whole(SegmentId::new("binder"), ElementIndex::new(0));
/// let strand = Constraint::<Symbolic>::SecondaryStructure {
///     scope: Scope::Address(run),
///     kind: SecondaryStructure::Strand,
/// };
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Scope {
    /// The slots one address names.
    Address(SymbolicAddress),
    /// The slots of a group defined in the specification.
    Group(GroupId),
}

/// An atom of a slot, named by the slot and the atom name.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{AtomAddress, SlotIndex};
/// use mmcif_model::AtomName;
///
/// let sulfur = AtomAddress::new(SlotIndex::new(4), AtomName::new("SG"));
///
/// assert_eq!(sulfur.slot, SlotIndex::new(4));
/// assert_eq!(sulfur.atom, AtomName::new("SG"));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomAddress<T> {
    /// The slot that holds the atom.
    pub slot: T,
    /// The name of the atom.
    pub atom: AtomName,
}

impl<T> AtomAddress<T> {
    /// Creates an address of the atom named `atom` in `slot`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{AtomAddress, ElementIndex, SegmentId, SymbolicAddress};
    /// use mmcif_model::AtomName;
    ///
    /// let slot = SymbolicAddress::whole(SegmentId::new("loop"), ElementIndex::new(1));
    /// let sulfur = AtomAddress::new(slot.clone(), AtomName::new("SG"));
    ///
    /// assert_eq!(sulfur.slot, slot);
    /// ```
    pub fn new(slot: T, atom: AtomName) -> AtomAddress<T> {
        AtomAddress { slot, atom }
    }
}

/// An inclusive range of distances.
///
/// Both ends are finite, the minimum is not negative, and the minimum is not
/// greater than the maximum.
///
/// # Examples
///
/// ```
/// use mmcif_ir::Bounds;
///
/// let disulfide = Bounds::new(1.8, 2.2)?;
///
/// assert_eq!((disulfide.min(), disulfide.max()), (1.8, 2.2));
/// assert!(disulfide.contains(2.0));
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds {
    min: f64,
    max: f64,
}

impl Bounds {
    /// Creates bounds from `min` to `max` inclusive.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::InvalidBounds`] if `min` or `max` is NaN or
    /// infinite, if `min` is negative, or if `min` is greater than `max`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Bounds, IrError};
    ///
    /// assert_eq!(Bounds::new(2.0, 2.0)?.max(), 2.0);
    ///
    /// assert_eq!(
    ///     Bounds::new(8.0, 2.0),
    ///     Err(IrError::InvalidBounds { min: 8.0, max: 2.0 })
    /// );
    /// assert!(Bounds::new(-1.0, 2.0).is_err());
    /// assert!(Bounds::new(0.0, f64::INFINITY).is_err());
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn new(min: f64, max: f64) -> Result<Bounds, IrError> {
        if !min.is_finite() || !max.is_finite() || min < 0.0 || min > max {
            return Err(IrError::InvalidBounds { min, max });
        }
        Ok(Bounds { min, max })
    }

    /// Creates bounds from zero to `max` inclusive.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::InvalidBounds`] if `max` is negative, NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::Bounds;
    ///
    /// let contact = Bounds::at_most(8.0)?;
    ///
    /// assert_eq!((contact.min(), contact.max()), (0.0, 8.0));
    /// assert!(Bounds::at_most(-1.0).is_err());
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn at_most(max: f64) -> Result<Bounds, IrError> {
        Bounds::new(0.0, max)
    }

    /// Returns the minimum distance.
    pub fn min(self) -> f64 {
        self.min
    }

    /// Returns the maximum distance.
    pub fn max(self) -> f64 {
        self.max
    }

    /// Returns `true` if `value` lies between the minimum and the maximum,
    /// inclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::Bounds;
    ///
    /// let bounds = Bounds::new(2.0, 8.0)?;
    ///
    /// assert!(bounds.contains(2.0) && bounds.contains(8.0));
    /// assert!(!bounds.contains(8.5));
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn contains(self, value: f64) -> bool {
        self.min <= value && value <= self.max
    }
}

/// A kind of secondary structure.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{Constraint, Resolved, SecondaryStructure, SlotIndex};
///
/// let helix = Constraint::<Resolved>::SecondaryStructure {
///     scope: vec![SlotIndex::new(0), SlotIndex::new(1)],
///     kind: SecondaryStructure::Helix,
/// };
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum SecondaryStructure {
    /// A helix.
    Helix,
    /// A strand.
    Strand,
    /// A loop.
    Loop,
}

/// A constraint on the atoms or slots of a design, with slots named through
/// `A`.
///
/// A [`Spec`](crate::Spec) holds `Constraint<Symbolic>` values, and
/// [`Spec::instantiate`](crate::Spec::instantiate) resolves them to
/// `Constraint<Resolved>` values.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{
///     AtomAddress, Bounds, Constraint, Element, ElementIndex, Geometry, Identity, Length,
///     Lengths, Segment, SegmentId, Slot, SlotIndex, Spec, Symbolic, SymbolicAddress,
/// };
/// use mmcif_model::{AtomName, CompId};
///
/// let cysteine = Slot::new(Identity::Exactly(CompId::new("CYS")), Geometry::Free);
/// let mut spec = Spec::new();
/// spec.push_segment(Segment::new(
///     SegmentId::new("loop"),
///     [
///         Element::Slot(cysteine.clone()),
///         Element::run(Slot::de_novo(), Length::range(2, 4)?),
///         Element::Slot(cysteine),
///     ],
/// ))?;
/// let sulfur = |element| {
///     let slot = SymbolicAddress::whole(SegmentId::new("loop"), ElementIndex::new(element));
///     AtomAddress::new(slot, AtomName::new("SG"))
/// };
/// spec.push_constraint(Constraint::<Symbolic>::Distance {
///     first: sulfur(0),
///     second: sulfur(2),
///     bounds: Bounds::new(1.8, 2.2)?,
/// })?;
///
/// let design = spec.instantiate(&Lengths::new(&spec, [3])?)?;
///
/// assert_eq!(design.constraints()[0].slots(), [SlotIndex::new(0), SlotIndex::new(4)]);
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Constraint<A: Addressing> {
    /// A bond between two atoms.
    Bond {
        /// The first atom.
        first: AtomAddress<A::One>,
        /// The second atom.
        second: AtomAddress<A::One>,
    },
    /// A distance between two atoms that lies within bounds.
    Distance {
        /// The first atom.
        first: AtomAddress<A::One>,
        /// The second atom.
        second: AtomAddress<A::One>,
        /// The allowed distance between the two atoms.
        bounds: Bounds,
    },
    /// A secondary structure taken by a set of slots.
    SecondaryStructure {
        /// The slots the constraint covers.
        scope: A::Many,
        /// The kind of secondary structure.
        kind: SecondaryStructure,
    },
}

impl Constraint<Resolved> {
    /// Returns the slots the constraint names.
    ///
    /// For a bond or a distance these are the slot of the first atom and the
    /// slot of the second, in that order. For a secondary structure they are
    /// the slots of the scope.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{AtomAddress, Bounds, Constraint, Resolved, SlotIndex};
    /// use mmcif_model::AtomName;
    ///
    /// let distance = Constraint::<Resolved>::Distance {
    ///     first: AtomAddress::new(SlotIndex::new(9), AtomName::new("CA")),
    ///     second: AtomAddress::new(SlotIndex::new(3), AtomName::new("CA")),
    ///     bounds: Bounds::at_most(8.0)?,
    /// };
    ///
    /// assert_eq!(distance.slots(), [SlotIndex::new(9), SlotIndex::new(3)]);
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn slots(&self) -> Vec<SlotIndex> {
        match self {
            Constraint::Bond { first, second } => vec![first.slot, second.slot],
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

use std::collections::BTreeSet;

use crate::constraint::{Constraint, Scope, Symbolic, SymbolicAddress};
use crate::error::IrError;
use crate::ids::{ElementIndex, GroupId, SegmentId, SegmentIndex};
use crate::segment::{Length, Segment};
use crate::slot::Determinacy;

/// A design specification: segments of slots, named groups of addresses and
/// constraints between slots.
///
/// Segments, groups and constraints are kept in the order they were added,
/// and each address is checked against the segments when it is added.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{
///     Corner, Element, ElementIndex, Length, Segment, SegmentId, Slot, Spec, SymbolicAddress,
/// };
///
/// let mut spec = Spec::new();
/// spec.push_segment(Segment::new(
///     SegmentId::new("binder"),
///     [Element::run(Slot::de_novo(), Length::range(8, 12)?)],
/// ))?;
/// let binder = SymbolicAddress::whole(SegmentId::new("binder"), ElementIndex::new(0));
///
/// assert_eq!(spec.length_bounds(), (8, 12));
/// assert_eq!(spec.determinacy(&binder)?.corner(), Some(Corner::DeNovo));
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Spec {
    segments: Vec<Segment>,
    groups: Vec<(GroupId, Vec<SymbolicAddress>)>,
    constraints: Vec<Constraint<Symbolic>>,
}

impl Spec {
    /// Creates an empty specification.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::Spec;
    ///
    /// let spec = Spec::new();
    ///
    /// assert!(spec.segments().is_empty());
    /// assert_eq!(spec.length_bounds(), (0, 0));
    /// ```
    pub fn new() -> Spec {
        Spec::default()
    }

    /// Appends `segment` and returns its index.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::DuplicateSegment`] if the specification already
    /// holds a segment with the same id.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{IrError, Segment, SegmentId, SegmentIndex, Spec};
    ///
    /// let mut spec = Spec::new();
    ///
    /// assert_eq!(
    ///     spec.push_segment(Segment::new(SegmentId::new("binder"), []))?,
    ///     SegmentIndex::new(0)
    /// );
    /// assert_eq!(
    ///     spec.push_segment(Segment::new(SegmentId::new("binder"), [])),
    ///     Err(IrError::DuplicateSegment(SegmentId::new("binder")))
    /// );
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn push_segment(&mut self, segment: Segment) -> Result<SegmentIndex, IrError> {
        if self.segments.iter().any(|known| known.id() == segment.id()) {
            return Err(IrError::DuplicateSegment(segment.id().clone()));
        }
        self.segments.push(segment);
        Ok(SegmentIndex::new((self.segments.len() - 1) as u32))
    }

    /// Defines a group named `id` whose members are `members`, in the order
    /// given.
    ///
    /// The specification is unchanged when an error is returned.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::DuplicateGroup`] if a group named `id` is already
    /// defined, [`IrError::UnknownSegment`] if a member names a segment the
    /// specification does not hold, [`IrError::UnknownElement`] if a member
    /// names an element its segment does not have, and
    /// [`IrError::OffsetOutOfRange`] if the offset of a member is not below
    /// the maximum length of its element.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{
    ///     Element, ElementIndex, GroupId, IrError, Length, Segment, SegmentId, Slot, Spec,
    ///     SymbolicAddress,
    /// };
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("loop"),
    ///     [Element::run(Slot::de_novo(), Length::range(2, 4)?)],
    /// ))?;
    /// let first = SymbolicAddress::at(SegmentId::new("loop"), ElementIndex::new(0), 0);
    /// let past_end = SymbolicAddress::at(SegmentId::new("loop"), ElementIndex::new(0), 4);
    ///
    /// spec.define_group(GroupId::new("start"), [first])?;
    /// assert!(matches!(
    ///     spec.define_group(GroupId::new("end"), [past_end]),
    ///     Err(IrError::OffsetOutOfRange { offset: 4, .. })
    /// ));
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn define_group(
        &mut self,
        id: GroupId,
        members: impl IntoIterator<Item = SymbolicAddress>,
    ) -> Result<(), IrError> {
        if self.groups.iter().any(|(known, _)| known == &id) {
            return Err(IrError::DuplicateGroup(id));
        }
        let members: Vec<SymbolicAddress> = members.into_iter().collect();
        for member in &members {
            self.check_address(member)?;
        }
        self.groups.push((id, members));
        Ok(())
    }

    /// Appends `constraint` after checking the addresses it names.
    ///
    /// A secondary structure constraint may cover a whole element of any
    /// length. The specification is unchanged when an error is returned.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownSegment`], [`IrError::UnknownElement`] or
    /// [`IrError::OffsetOutOfRange`] if an address in the constraint fails
    /// the checks of [`define_group`](Spec::define_group),
    /// [`IrError::AddressNotSingular`] if a bond or distance atom names a
    /// whole element whose maximum length is not one, and
    /// [`IrError::UnknownGroup`] if a secondary structure constraint covers a
    /// group that is not defined.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{
    ///     AtomAddress, Bounds, Constraint, Element, ElementIndex, IrError, Length, Segment,
    ///     SegmentId, Slot, Spec, Symbolic, SymbolicAddress,
    /// };
    /// use mmcif_model::AtomName;
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("binder"),
    ///     [Element::run(Slot::de_novo(), Length::exact(10)?)],
    /// ))?;
    /// let binder = SegmentId::new("binder");
    /// let run = SymbolicAddress::whole(binder.clone(), ElementIndex::new(0));
    /// let third = SymbolicAddress::at(binder.clone(), ElementIndex::new(0), 2);
    /// let eighth = SymbolicAddress::at(binder, ElementIndex::new(0), 7);
    /// let alpha_carbon = |slot| AtomAddress::new(slot, AtomName::new("CA"));
    ///
    /// spec.push_constraint(Constraint::<Symbolic>::Distance {
    ///     first: alpha_carbon(third.clone()),
    ///     second: alpha_carbon(eighth),
    ///     bounds: Bounds::at_most(10.0)?,
    /// })?;
    /// assert_eq!(spec.constraints().len(), 1);
    ///
    /// let whole_run = Constraint::<Symbolic>::Distance {
    ///     first: alpha_carbon(run),
    ///     second: alpha_carbon(third),
    ///     bounds: Bounds::at_most(10.0)?,
    /// };
    /// assert!(matches!(
    ///     spec.push_constraint(whole_run),
    ///     Err(IrError::AddressNotSingular { slots: 10, .. })
    /// ));
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn push_constraint(&mut self, constraint: Constraint<Symbolic>) -> Result<(), IrError> {
        match &constraint {
            Constraint::Bond { first, second } | Constraint::Distance { first, second, .. } => {
                for address in [&first.slot, &second.slot] {
                    self.check_address(address)?;
                    self.check_singular(address)?;
                }
            }
            Constraint::SecondaryStructure { scope, .. } => self.check_scope(scope)?,
        }
        self.constraints.push(constraint);
        Ok(())
    }

    /// Returns the segments, in the order they were pushed.
    pub fn segments(&self) -> &[Segment] {
        &self.segments
    }

    /// Returns the index of the segment named `id`.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownSegment`] if no segment is named `id`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Segment, SegmentId, SegmentIndex, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(SegmentId::new("target"), []))?;
    /// spec.push_segment(Segment::new(SegmentId::new("binder"), []))?;
    ///
    /// assert_eq!(spec.segment_index(&SegmentId::new("binder"))?, SegmentIndex::new(1));
    /// assert!(spec.segment_index(&SegmentId::new("linker")).is_err());
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn segment_index(&self, id: &SegmentId) -> Result<SegmentIndex, IrError> {
        self.segments
            .iter()
            .position(|segment| segment.id() == id)
            .map(|index| SegmentIndex::new(index as u32))
            .ok_or_else(|| IrError::UnknownSegment(id.clone()))
    }

    /// Returns the segment named `id`.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownSegment`] if no segment is named `id`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, IrError, Segment, SegmentId, Slot, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(SegmentId::new("loop"), [Element::Slot(Slot::de_novo())]))?;
    ///
    /// assert_eq!(spec.segment(&SegmentId::new("loop"))?.elements().len(), 1);
    /// assert_eq!(
    ///     spec.segment(&SegmentId::new("linker")),
    ///     Err(IrError::UnknownSegment(SegmentId::new("linker")))
    /// );
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn segment(&self, id: &SegmentId) -> Result<&Segment, IrError> {
        let index = self.segment_index(id)?;
        Ok(&self.segments[index.as_usize()])
    }

    /// Returns an iterator over the groups and their members, in the order
    /// the groups were defined.
    pub fn groups(&self) -> impl ExactSizeIterator<Item = (&GroupId, &[SymbolicAddress])> {
        self.groups
            .iter()
            .map(|(id, members)| (id, members.as_slice()))
    }

    /// Returns the members of the group named `id`.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownGroup`] if no group is named `id`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{
    ///     Element, ElementIndex, GroupId, IrError, Segment, SegmentId, Slot, Spec,
    ///     SymbolicAddress,
    /// };
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(SegmentId::new("loop"), [Element::Slot(Slot::de_novo())]))?;
    /// let member = SymbolicAddress::whole(SegmentId::new("loop"), ElementIndex::new(0));
    /// spec.define_group(GroupId::new("core"), [member.clone()])?;
    ///
    /// assert_eq!(spec.group(&GroupId::new("core"))?, [member]);
    /// assert_eq!(
    ///     spec.group(&GroupId::new("rim")),
    ///     Err(IrError::UnknownGroup(GroupId::new("rim")))
    /// );
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn group(&self, id: &GroupId) -> Result<&[SymbolicAddress], IrError> {
        self.groups
            .iter()
            .find(|(known, _)| known == id)
            .map(|(_, members)| members.as_slice())
            .ok_or_else(|| IrError::UnknownGroup(id.clone()))
    }

    /// Returns the constraints, in the order they were pushed.
    pub fn constraints(&self) -> &[Constraint<Symbolic>] {
        &self.constraints
    }

    /// Returns the smallest and the largest number of slots the
    /// specification can produce, as `(min, max)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Length, Segment, SegmentId, Slot, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(SegmentId::new("target"), [Element::Slot(Slot::de_novo())]))?;
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("binder"),
    ///     [Element::run(Slot::de_novo(), Length::range(8, 12)?)],
    /// ))?;
    ///
    /// assert_eq!(spec.length_bounds(), (9, 13));
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn length_bounds(&self) -> (u64, u64) {
        self.segments.iter().fold((0, 0), |(min, max), segment| {
            let (segment_min, segment_max) = segment.length_bounds();
            (min + segment_min, max + segment_max)
        })
    }

    /// Returns `true` if any segment of the specification is symbolic.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Length, Segment, SegmentId, Slot, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(SegmentId::new("target"), [Element::Slot(Slot::de_novo())]))?;
    /// assert!(!spec.is_symbolic());
    ///
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("binder"),
    ///     [Element::run(Slot::de_novo(), Length::range(8, 12)?)],
    /// ))?;
    /// assert!(spec.is_symbolic());
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn is_symbolic(&self) -> bool {
        self.segments.iter().any(Segment::is_symbolic)
    }

    /// Returns an iterator over the symbolic elements, with the index of
    /// their segment, their own index and their length.
    ///
    /// Elements come in segment order and then element order. This is the
    /// order of the counts in [`Lengths`](crate::Lengths).
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, ElementIndex, Length, Segment, SegmentId, SegmentIndex, Slot, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("loop"),
    ///     [
    ///         Element::Slot(Slot::de_novo()),
    ///         Element::run(Slot::de_novo(), Length::range(3, 5)?),
    ///     ],
    /// ))?;
    /// let symbolic: Vec<_> = spec.symbolic_elements().collect();
    ///
    /// assert_eq!(
    ///     symbolic,
    ///     [(SegmentIndex::new(0), ElementIndex::new(1), Length::range(3, 5)?)]
    /// );
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn symbolic_elements(
        &self,
    ) -> impl Iterator<Item = (SegmentIndex, ElementIndex, Length)> + use<'_> {
        self.segments
            .iter()
            .enumerate()
            .flat_map(|(segment, definition)| {
                definition
                    .elements()
                    .iter()
                    .enumerate()
                    .filter(|(_, element)| element.is_symbolic())
                    .map(move |(element, definition)| {
                        (
                            SegmentIndex::new(segment as u32),
                            ElementIndex::new(element as u32),
                            definition.length(),
                        )
                    })
            })
    }

    /// Returns the determinacy of the slots `address` names.
    ///
    /// The geometry of a free slot is
    /// [`Level::Restricted`](crate::Level::Restricted) when a constraint
    /// names it: one of the atoms of the constraint, its scope address, or a
    /// member of its scope group addresses the whole element of the slot or
    /// the offset of the slot. An address of a whole element returns the
    /// determinacy all slots of the element share.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownSegment`] if no segment is named by
    /// `address`, [`IrError::UnknownElement`] if the segment has no element
    /// at its index, [`IrError::OffsetOutOfRange`] if its offset is not below
    /// the maximum length of the element, and [`IrError::MixedDeterminacy`]
    /// if `address` names a whole element and constraints name some, but not
    /// all, of the offsets below its maximum length.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{
    ///     AtomAddress, Bounds, Constraint, Element, ElementIndex, IrError, Length, Level,
    ///     Segment, SegmentId, Slot, Spec, Symbolic, SymbolicAddress,
    /// };
    /// use mmcif_model::AtomName;
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("binder"),
    ///     [Element::run(Slot::de_novo(), Length::exact(10)?)],
    /// ))?;
    /// let at = |offset| SymbolicAddress::at(SegmentId::new("binder"), ElementIndex::new(0), offset);
    /// let alpha_carbon = |offset| AtomAddress::new(at(offset), AtomName::new("CA"));
    /// spec.push_constraint(Constraint::<Symbolic>::Distance {
    ///     first: alpha_carbon(3),
    ///     second: alpha_carbon(4),
    ///     bounds: Bounds::at_most(8.0)?,
    /// })?;
    ///
    /// assert_eq!(spec.determinacy(&at(3))?.geometry, Level::Restricted);
    /// assert_eq!(spec.determinacy(&at(7))?.geometry, Level::Free);
    ///
    /// let run = SymbolicAddress::whole(SegmentId::new("binder"), ElementIndex::new(0));
    /// assert!(matches!(
    ///     spec.determinacy(&run),
    ///     Err(IrError::MixedDeterminacy { .. })
    /// ));
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn determinacy(&self, address: &SymbolicAddress) -> Result<Determinacy, IrError> {
        self.check_address(address)?;
        let element = self.segment(&address.segment)?.element(address.element)?;
        let constrained = self.constrained_at(address, element.length().max())?;
        Ok(element.slot().determinacy(constrained))
    }

    fn constrained_at(&self, address: &SymbolicAddress, slots: u32) -> Result<bool, IrError> {
        let (whole, offsets) = self.constraint_extent(&address.segment, address.element);
        match address.offset {
            Some(offset) => Ok(whole || offsets.contains(&offset)),
            None if whole || offsets.len() as u32 == slots => Ok(true),
            None if offsets.is_empty() => Ok(false),
            None => Err(IrError::MixedDeterminacy {
                segment: address.segment.clone(),
                element: address.element,
            }),
        }
    }

    fn constraint_extent(
        &self,
        segment: &SegmentId,
        element: ElementIndex,
    ) -> (bool, BTreeSet<u32>) {
        let mut whole = false;
        let mut offsets = BTreeSet::new();
        for address in self.constraint_addresses() {
            if &address.segment != segment || address.element != element {
                continue;
            }
            match address.offset {
                None => whole = true,
                Some(offset) => {
                    offsets.insert(offset);
                }
            }
        }
        (whole, offsets)
    }

    fn constraint_addresses(&self) -> Vec<&SymbolicAddress> {
        let mut found = Vec::new();
        for constraint in &self.constraints {
            match constraint {
                Constraint::Bond { first, second } | Constraint::Distance { first, second, .. } => {
                    found.push(&first.slot);
                    found.push(&second.slot);
                }
                Constraint::SecondaryStructure { scope, .. } => match scope {
                    Scope::Address(address) => found.push(address),
                    Scope::Group(id) => {
                        if let Ok(members) = self.group(id) {
                            found.extend(members.iter());
                        }
                    }
                },
            }
        }
        found
    }

    fn check_address(&self, address: &SymbolicAddress) -> Result<(), IrError> {
        let element = self.segment(&address.segment)?.element(address.element)?;
        match address.offset {
            Some(offset) if offset >= element.length().max() => Err(IrError::OffsetOutOfRange {
                segment: address.segment.clone(),
                element: address.element,
                offset,
            }),
            _ => Ok(()),
        }
    }

    fn check_singular(&self, address: &SymbolicAddress) -> Result<(), IrError> {
        let element = self.segment(&address.segment)?.element(address.element)?;
        let slots = element.length().max();
        if address.offset.is_none() && slots != 1 {
            return Err(IrError::AddressNotSingular {
                segment: address.segment.clone(),
                element: address.element,
                slots,
            });
        }
        Ok(())
    }

    fn check_scope(&self, scope: &Scope) -> Result<(), IrError> {
        match scope {
            Scope::Address(address) => self.check_address(address),
            Scope::Group(id) => self.group(id).map(|_| ()),
        }
    }
}

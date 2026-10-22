use crate::constraint::{AtomAddress, Constraint, Resolved, Scope, SymbolicAddress};
use crate::error::IrError;
use crate::ids::{ElementIndex, GroupId, SegmentId, SegmentIndex, SlotIndex, SlotSpan};
use crate::segment::Length;
use crate::slot::{Determinacy, Slot};
use crate::spec::Spec;

/// The largest number of slots an instantiation may hold.
pub const MAX_SLOTS: u32 = 1 << 24;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Choice {
    segment: SegmentIndex,
    element: ElementIndex,
    count: u32,
}

/// A chosen count for each symbolic element of a [`Spec`].
///
/// The counts follow the order of [`Spec::symbolic_elements`].
///
/// # Examples
///
/// ```
/// use mmcif_ir::{Element, Length, Lengths, Segment, SegmentId, Slot, Spec};
///
/// let mut spec = Spec::new();
/// spec.push_segment(Segment::new(
///     SegmentId::new("loop"),
///     [
///         Element::run(Slot::de_novo(), Length::range(3, 5)?),
///         Element::Slot(Slot::de_novo()),
///         Element::run(Slot::de_novo(), Length::range(1, 4)?),
///     ],
/// ))?;
///
/// assert_eq!(Lengths::new(&spec, [4, 2])?.counts(), [4, 2]);
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Lengths {
    choices: Vec<Choice>,
}

impl Lengths {
    /// Chooses the minimum count for each symbolic element of `spec`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Length, Lengths, Segment, SegmentId, Slot, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("loop"),
    ///     [
    ///         Element::run(Slot::de_novo(), Length::range(3, 5)?),
    ///         Element::run(Slot::de_novo(), Length::range(0, 4)?),
    ///     ],
    /// ))?;
    ///
    /// assert_eq!(Lengths::minimum(&spec).counts(), [3, 0]);
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn minimum(spec: &Spec) -> Lengths {
        Lengths::chosen(spec, Length::min)
    }

    /// Chooses the maximum count for each symbolic element of `spec`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Length, Lengths, Segment, SegmentId, Slot, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("loop"),
    ///     [
    ///         Element::run(Slot::de_novo(), Length::range(3, 5)?),
    ///         Element::run(Slot::de_novo(), Length::range(0, 4)?),
    ///     ],
    /// ))?;
    ///
    /// assert_eq!(Lengths::maximum(&spec).counts(), [5, 4]);
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn maximum(spec: &Spec) -> Lengths {
        Lengths::chosen(spec, Length::max)
    }

    fn chosen(spec: &Spec, pick: impl Fn(Length) -> u32) -> Lengths {
        Lengths {
            choices: spec
                .symbolic_elements()
                .map(|(segment, element, length)| Choice {
                    segment,
                    element,
                    count: pick(length),
                })
                .collect(),
        }
    }

    /// Chooses `counts` for the symbolic elements of `spec`, in the order of
    /// [`Spec::symbolic_elements`].
    ///
    /// # Errors
    ///
    /// Returns [`IrError::LengthCountMismatch`] if the number of counts
    /// differs from the number of symbolic elements, and
    /// [`IrError::LengthOutOfRange`] if a count lies outside the range of its
    /// element.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, IrError, Length, Lengths, Segment, SegmentId, Slot, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("loop"),
    ///     [Element::run(Slot::de_novo(), Length::range(3, 5)?)],
    /// ))?;
    ///
    /// assert_eq!(Lengths::new(&spec, [4])?.counts(), [4]);
    /// assert_eq!(
    ///     Lengths::new(&spec, [4, 4]),
    ///     Err(IrError::LengthCountMismatch { expected: 1, actual: 2 })
    /// );
    /// assert_eq!(
    ///     Lengths::new(&spec, [9]),
    ///     Err(IrError::LengthOutOfRange { chosen: 9, min: 3, max: 5 })
    /// );
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn new(spec: &Spec, counts: impl IntoIterator<Item = u32>) -> Result<Lengths, IrError> {
        let counts: Vec<u32> = counts.into_iter().collect();
        let expected = spec.symbolic_elements().count();
        if counts.len() != expected {
            return Err(IrError::LengthCountMismatch {
                expected,
                actual: counts.len(),
            });
        }
        let lengths = Lengths {
            choices: spec
                .symbolic_elements()
                .zip(counts)
                .map(|((segment, element, _), count)| Choice {
                    segment,
                    element,
                    count,
                })
                .collect(),
        };
        lengths.check(spec)?;
        Ok(lengths)
    }

    /// Returns the chosen counts, in the order of
    /// [`Spec::symbolic_elements`].
    pub fn counts(&self) -> Vec<u32> {
        self.choices.iter().map(|choice| choice.count).collect()
    }

    fn check(&self, spec: &Spec) -> Result<(), IrError> {
        let expected = spec.symbolic_elements().count();
        if self.choices.len() != expected {
            return Err(IrError::LengthCountMismatch {
                expected,
                actual: self.choices.len(),
            });
        }
        for (choice, (segment, element, length)) in
            self.choices.iter().zip(spec.symbolic_elements())
        {
            if choice.segment != segment || choice.element != element {
                return Err(IrError::LengthsFromAnotherSpec {
                    expected: address_of(spec, segment, element),
                    actual: address_of(spec, choice.segment, choice.element),
                });
            }
            if !length.contains(choice.count) {
                return Err(IrError::LengthOutOfRange {
                    chosen: choice.count,
                    min: length.min(),
                    max: length.max(),
                });
            }
        }
        Ok(())
    }

    fn total_slots(&self, spec: &Spec) -> u64 {
        let fixed: u64 = spec
            .segments()
            .iter()
            .flat_map(|segment| segment.elements())
            .filter(|element| !element.is_symbolic())
            .map(|element| u64::from(element.length().min()))
            .sum();
        fixed
            + self
                .choices
                .iter()
                .map(|choice| u64::from(choice.count))
                .sum::<u64>()
    }
}

fn address_of(spec: &Spec, segment: SegmentIndex, element: ElementIndex) -> SymbolicAddress {
    let id = spec
        .segments()
        .get(segment.as_usize())
        .map(|segment| segment.id().clone())
        .unwrap_or_else(|| SegmentId::new(segment.to_string()));
    SymbolicAddress::whole(id, element)
}

/// The element and the offset an instantiated slot was produced from.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{
///     Element, ElementIndex, Length, Lengths, Origin, Segment, SegmentId, SegmentIndex, Slot,
///     SlotIndex, Spec,
/// };
///
/// let mut spec = Spec::new();
/// spec.push_segment(Segment::new(
///     SegmentId::new("binder"),
///     [
///         Element::Slot(Slot::de_novo()),
///         Element::run(Slot::de_novo(), Length::range(2, 4)?),
///     ],
/// ))?;
/// let design = spec.instantiate(&Lengths::new(&spec, [3])?)?;
///
/// assert_eq!(
///     design.origin(SlotIndex::new(2))?,
///     Origin {
///         segment: SegmentIndex::new(0),
///         element: ElementIndex::new(1),
///         offset: 1,
///     }
/// );
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Origin {
    /// The index of the segment the slot belongs to.
    pub segment: SegmentIndex,
    /// The index of the element that produced the slot.
    pub element: ElementIndex,
    /// The position of the slot among the slots its element produced.
    pub offset: u32,
}

/// A specification with every length chosen, its slots numbered and its
/// groups and constraints resolved to slot indices.
///
/// Slots are numbered from zero in segment order, then element order, then
/// offset order. Instantiations are created by
/// [`Spec::instantiate`].
///
/// # Examples
///
/// ```
/// use mmcif_ir::{Corner, Element, Length, Lengths, ResidueRef, Segment, SegmentId, Slot, Spec};
/// use mmcif_model::{CompId, LabelResidueKey, LabelSeqId};
///
/// let anchor = ResidueRef::Label(LabelResidueKey {
///     asym: "A".into(),
///     seq: LabelSeqId::new(10).unwrap(),
/// });
/// let mut spec = Spec::new();
/// spec.push_segment(Segment::new(
///     SegmentId::new("target"),
///     [Element::Slot(Slot::context(CompId::new("GLY"), anchor))],
/// ))?;
/// spec.push_segment(Segment::new(
///     SegmentId::new("binder"),
///     [Element::run(Slot::de_novo(), Length::range(8, 12)?)],
/// ))?;
///
/// let design = spec.instantiate(&Lengths::new(&spec, [10])?)?;
/// let binder = design.segment_slots(&SegmentId::new("binder"))?;
///
/// assert_eq!(design.slot_count(), 11);
/// assert_eq!(binder.len(), 10);
/// assert_eq!(design.determinacy(binder.start())?.corner(), Some(Corner::DeNovo));
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Instantiated {
    slots: Vec<Slot>,
    origins: Vec<Origin>,
    segment_ids: Vec<SegmentId>,
    segment_slots: Vec<SlotSpan>,
    element_slots: Vec<Vec<SlotSpan>>,
    groups: Vec<(GroupId, Vec<SlotIndex>)>,
    constraints: Vec<Constraint<Resolved>>,
    constrained: Vec<bool>,
}

impl Spec {
    /// Instantiates the specification with the counts in `lengths`.
    ///
    /// Each element produces as many copies of its slot as its exact length
    /// or its chosen count, and a count of zero produces no slot. Each group
    /// and each secondary structure scope resolves to the slots its addresses
    /// name, in ascending order without repeats. Each bond or distance atom
    /// resolves to one slot.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::LengthCountMismatch`] if `lengths` holds a
    /// different number of counts than the specification has symbolic
    /// elements, [`IrError::LengthsFromAnotherSpec`] if `lengths` was chosen
    /// for a specification whose symbolic elements sit at other addresses,
    /// [`IrError::LengthOutOfRange`] if a count lies outside the range of its
    /// element, [`IrError::TooManySlots`] if more than [`MAX_SLOTS`] slots
    /// would be produced, [`IrError::OffsetOutOfRange`] if an offset in a
    /// group or a constraint is not below the number of slots its element
    /// produced, and [`IrError::AddressNotSingular`] if a bond or distance
    /// atom names a whole element that did not produce exactly one slot.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{
    ///     AtomAddress, Bounds, Constraint, Element, ElementIndex, IrError, Length, Lengths,
    ///     Segment, SegmentId, Slot, SlotIndex, Spec, Symbolic, SymbolicAddress,
    /// };
    /// use mmcif_model::AtomName;
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("loop"),
    ///     [Element::run(Slot::de_novo(), Length::range(1, 6)?)],
    /// ))?;
    /// let alpha_carbon = |offset| {
    ///     let slot = SymbolicAddress::at(SegmentId::new("loop"), ElementIndex::new(0), offset);
    ///     AtomAddress::new(slot, AtomName::new("CA"))
    /// };
    /// spec.push_constraint(Constraint::<Symbolic>::Distance {
    ///     first: alpha_carbon(0),
    ///     second: alpha_carbon(4),
    ///     bounds: Bounds::at_most(10.0)?,
    /// })?;
    ///
    /// let design = spec.instantiate(&Lengths::new(&spec, [5])?)?;
    /// assert_eq!(design.slot_count(), 5);
    /// assert_eq!(design.constraints()[0].slots(), [SlotIndex::new(0), SlotIndex::new(4)]);
    ///
    /// assert!(matches!(
    ///     spec.instantiate(&Lengths::new(&spec, [3])?),
    ///     Err(IrError::OffsetOutOfRange { offset: 4, .. })
    /// ));
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn instantiate(&self, lengths: &Lengths) -> Result<Instantiated, IrError> {
        lengths.check(self)?;
        let total = lengths.total_slots(self);
        if total > u64::from(MAX_SLOTS) {
            return Err(IrError::TooManySlots(MAX_SLOTS));
        }
        let mut slots = Vec::with_capacity(total as usize);
        let mut origins = Vec::with_capacity(total as usize);
        let mut element_slots = Vec::with_capacity(self.segments().len());
        let mut segment_slots = Vec::with_capacity(self.segments().len());
        let mut choices = lengths.choices.iter();
        for (segment_index, segment) in self.segments().iter().enumerate() {
            let segment_start = slots.len() as u32;
            let mut spans = Vec::with_capacity(segment.elements().len());
            for (element_index, element) in segment.elements().iter().enumerate() {
                let count = if element.is_symbolic() {
                    choices
                        .next()
                        .expect("lengths cover every symbolic element")
                        .count
                } else {
                    element.length().min()
                };
                let element_start = slots.len() as u32;
                for offset in 0..count {
                    slots.push(element.slot().clone());
                    origins.push(Origin {
                        segment: SegmentIndex::new(segment_index as u32),
                        element: ElementIndex::new(element_index as u32),
                        offset,
                    });
                }
                spans.push(SlotSpan::new(element_start, slots.len() as u32));
            }
            element_slots.push(spans);
            segment_slots.push(SlotSpan::new(segment_start, slots.len() as u32));
        }

        let groups = self
            .groups()
            .map(|(id, members)| {
                let mut resolved = Vec::new();
                for member in members {
                    resolved.extend(resolve(self, &element_slots, member)?);
                }
                Ok((id.clone(), sorted(resolved)))
            })
            .collect::<Result<Vec<_>, IrError>>()?;

        let constraints = self
            .constraints()
            .iter()
            .map(|constraint| resolve_constraint(self, &element_slots, &groups, constraint))
            .collect::<Result<Vec<_>, IrError>>()?;

        let mut constrained = vec![false; slots.len()];
        for slot in constraints.iter().flat_map(Constraint::slots) {
            constrained[slot.as_usize()] = true;
        }

        Ok(Instantiated {
            slots,
            origins,
            segment_ids: self
                .segments()
                .iter()
                .map(|segment| segment.id().clone())
                .collect(),
            segment_slots,
            element_slots,
            groups,
            constraints,
            constrained,
        })
    }
}

fn sorted(mut slots: Vec<SlotIndex>) -> Vec<SlotIndex> {
    slots.sort_unstable();
    slots.dedup();
    slots
}

fn resolve_within(
    elements: &[SlotSpan],
    address: &SymbolicAddress,
) -> Result<Vec<SlotIndex>, IrError> {
    let span = elements
        .get(address.element.as_usize())
        .copied()
        .ok_or_else(|| IrError::UnknownElement {
            segment: address.segment.clone(),
            element: address.element,
        })?;
    match address.offset {
        None => Ok(span.iter().collect()),
        Some(offset) if (offset as usize) < span.len() => {
            Ok(vec![SlotIndex::new(span.start().get() + offset)])
        }
        Some(offset) => Err(IrError::OffsetOutOfRange {
            segment: address.segment.clone(),
            element: address.element,
            offset,
        }),
    }
}

fn resolve(
    spec: &Spec,
    element_slots: &[Vec<SlotSpan>],
    address: &SymbolicAddress,
) -> Result<Vec<SlotIndex>, IrError> {
    let segment = spec.segment_index(&address.segment)?;
    resolve_within(&element_slots[segment.as_usize()], address)
}

fn resolve_one(
    spec: &Spec,
    element_slots: &[Vec<SlotSpan>],
    address: &SymbolicAddress,
) -> Result<SlotIndex, IrError> {
    let resolved = resolve(spec, element_slots, address)?;
    match resolved.as_slice() {
        [slot] => Ok(*slot),
        slots => Err(IrError::AddressNotSingular {
            segment: address.segment.clone(),
            element: address.element,
            slots: slots.len() as u32,
        }),
    }
}

fn resolve_constraint(
    spec: &Spec,
    element_slots: &[Vec<SlotSpan>],
    groups: &[(GroupId, Vec<SlotIndex>)],
    constraint: &Constraint<crate::constraint::Symbolic>,
) -> Result<Constraint<Resolved>, IrError> {
    match constraint {
        Constraint::Bond { first, second } => Ok(Constraint::Bond {
            first: AtomAddress::new(
                resolve_one(spec, element_slots, &first.slot)?,
                first.atom.clone(),
            ),
            second: AtomAddress::new(
                resolve_one(spec, element_slots, &second.slot)?,
                second.atom.clone(),
            ),
        }),
        Constraint::Distance {
            first,
            second,
            bounds,
        } => Ok(Constraint::Distance {
            first: AtomAddress::new(
                resolve_one(spec, element_slots, &first.slot)?,
                first.atom.clone(),
            ),
            second: AtomAddress::new(
                resolve_one(spec, element_slots, &second.slot)?,
                second.atom.clone(),
            ),
            bounds: *bounds,
        }),
        Constraint::SecondaryStructure { scope, kind } => {
            let slots = match scope {
                Scope::Address(address) => sorted(resolve(spec, element_slots, address)?),
                Scope::Group(id) => groups
                    .iter()
                    .find(|(known, _)| known == id)
                    .map(|(_, slots)| slots.clone())
                    .ok_or_else(|| IrError::UnknownGroup(id.clone()))?,
            };
            Ok(Constraint::SecondaryStructure {
                scope: slots,
                kind: *kind,
            })
        }
    }
}

impl Instantiated {
    /// Returns the number of slots.
    pub fn slot_count(&self) -> usize {
        self.slots.len()
    }

    /// Returns `true` if the instantiation holds no slots.
    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }

    /// Returns the slots, in index order.
    pub fn slots(&self) -> &[Slot] {
        &self.slots
    }

    /// Returns the slot at `index`.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownSlot`] if `index` is not below the slot
    /// count.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, IrError, Length, Lengths, Segment, SegmentId, Slot, SlotIndex, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("binder"),
    ///     [Element::run(Slot::de_novo(), Length::range(2, 4)?)],
    /// ))?;
    /// let design = spec.instantiate(&Lengths::new(&spec, [3])?)?;
    ///
    /// assert_eq!(design.slot(SlotIndex::new(2))?, &Slot::de_novo());
    /// assert_eq!(
    ///     design.slot(SlotIndex::new(3)),
    ///     Err(IrError::UnknownSlot(SlotIndex::new(3)))
    /// );
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn slot(&self, index: SlotIndex) -> Result<&Slot, IrError> {
        self.slots
            .get(index.as_usize())
            .ok_or(IrError::UnknownSlot(index))
    }

    /// Returns the element and the offset the slot at `index` was produced
    /// from.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownSlot`] if `index` is not below the slot
    /// count.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, ElementIndex, Length, Lengths, Segment, SegmentId, Slot, SlotIndex, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("binder"),
    ///     [
    ///         Element::Slot(Slot::de_novo()),
    ///         Element::run(Slot::de_novo(), Length::range(2, 4)?),
    ///     ],
    /// ))?;
    /// let design = spec.instantiate(&Lengths::new(&spec, [3])?)?;
    /// let origin = design.origin(SlotIndex::new(3))?;
    ///
    /// assert_eq!((origin.element, origin.offset), (ElementIndex::new(1), 2));
    /// assert!(design.origin(SlotIndex::new(4)).is_err());
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn origin(&self, index: SlotIndex) -> Result<Origin, IrError> {
        self.origins
            .get(index.as_usize())
            .copied()
            .ok_or(IrError::UnknownSlot(index))
    }

    /// Returns the address of the slot at `index` in the specification, with
    /// its offset.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownSlot`] if `index` is not below the slot
    /// count.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{
    ///     Element, ElementIndex, Length, Lengths, Segment, SegmentId, Slot, SlotIndex, Spec,
    ///     SymbolicAddress,
    /// };
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("binder"),
    ///     [
    ///         Element::Slot(Slot::de_novo()),
    ///         Element::run(Slot::de_novo(), Length::range(2, 4)?),
    ///     ],
    /// ))?;
    /// let design = spec.instantiate(&Lengths::new(&spec, [3])?)?;
    /// let source = design.source(SlotIndex::new(3))?;
    ///
    /// assert_eq!(
    ///     source,
    ///     SymbolicAddress::at(SegmentId::new("binder"), ElementIndex::new(1), 2)
    /// );
    /// assert_eq!(design.resolve(&source)?, [SlotIndex::new(3)]);
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn source(&self, index: SlotIndex) -> Result<SymbolicAddress, IrError> {
        let origin = self.origin(index)?;
        Ok(SymbolicAddress::at(
            self.segment_ids[origin.segment.as_usize()].clone(),
            origin.element,
            origin.offset,
        ))
    }

    /// Returns the ids of the segments, in the order of the specification.
    pub fn segment_ids(&self) -> &[SegmentId] {
        &self.segment_ids
    }

    /// Returns the span of the slots the segment named `id` produced.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownSegment`] if no segment is named `id`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Length, Lengths, Segment, SegmentId, Slot, SlotIndex, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(SegmentId::new("target"), [Element::Slot(Slot::de_novo())]))?;
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("binder"),
    ///     [Element::run(Slot::de_novo(), Length::range(8, 12)?)],
    /// ))?;
    /// let design = spec.instantiate(&Lengths::minimum(&spec))?;
    /// let binder = design.segment_slots(&SegmentId::new("binder"))?;
    ///
    /// assert_eq!((binder.start(), binder.len()), (SlotIndex::new(1), 8));
    /// assert!(design.segment_slots(&SegmentId::new("linker")).is_err());
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn segment_slots(&self, id: &SegmentId) -> Result<SlotSpan, IrError> {
        Ok(self.segment_slots[self.segment_position(id)?])
    }

    /// Returns an iterator over the groups and their slots, in the order the
    /// groups were defined.
    ///
    /// The slots of each group are in ascending order without repeats.
    pub fn groups(&self) -> impl ExactSizeIterator<Item = (&GroupId, &[SlotIndex])> {
        self.groups.iter().map(|(id, slots)| (id, slots.as_slice()))
    }

    /// Returns the slots of the group named `id`, in ascending order without
    /// repeats.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownGroup`] if no group is named `id`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{
    ///     Element, ElementIndex, GroupId, Length, Lengths, Segment, SegmentId, Slot, SlotIndex,
    ///     Spec, SymbolicAddress,
    /// };
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("loop"),
    ///     [Element::run(Slot::de_novo(), Length::exact(4)?)],
    /// ))?;
    /// let loop_id = SegmentId::new("loop");
    /// spec.define_group(
    ///     GroupId::new("ends"),
    ///     [
    ///         SymbolicAddress::at(loop_id.clone(), ElementIndex::new(0), 3),
    ///         SymbolicAddress::at(loop_id.clone(), ElementIndex::new(0), 0),
    ///         SymbolicAddress::at(loop_id, ElementIndex::new(0), 3),
    ///     ],
    /// )?;
    /// let design = spec.instantiate(&Lengths::minimum(&spec))?;
    ///
    /// assert_eq!(
    ///     design.group(&GroupId::new("ends"))?,
    ///     [SlotIndex::new(0), SlotIndex::new(3)]
    /// );
    /// assert!(design.group(&GroupId::new("core")).is_err());
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn group(&self, id: &GroupId) -> Result<&[SlotIndex], IrError> {
        self.groups
            .iter()
            .find(|(known, _)| known == id)
            .map(|(_, slots)| slots.as_slice())
            .ok_or_else(|| IrError::UnknownGroup(id.clone()))
    }

    /// Returns the resolved constraints, in the order they were pushed to the
    /// specification.
    pub fn constraints(&self) -> &[Constraint<Resolved>] {
        &self.constraints
    }

    /// Returns `true` if a resolved constraint names the slot at `index`.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownSlot`] if `index` is not below the slot
    /// count.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{
    ///     AtomAddress, Constraint, Element, ElementIndex, Length, Lengths, Segment, SegmentId,
    ///     Slot, SlotIndex, Spec, Symbolic, SymbolicAddress,
    /// };
    /// use mmcif_model::AtomName;
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("loop"),
    ///     [Element::run(Slot::de_novo(), Length::exact(6)?)],
    /// ))?;
    /// let sulfur = |offset| {
    ///     let slot = SymbolicAddress::at(SegmentId::new("loop"), ElementIndex::new(0), offset);
    ///     AtomAddress::new(slot, AtomName::new("SG"))
    /// };
    /// spec.push_constraint(Constraint::<Symbolic>::Bond {
    ///     first: sulfur(1),
    ///     second: sulfur(4),
    /// })?;
    /// let design = spec.instantiate(&Lengths::minimum(&spec))?;
    ///
    /// assert!(design.is_constrained(SlotIndex::new(4))?);
    /// assert!(!design.is_constrained(SlotIndex::new(2))?);
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn is_constrained(&self, index: SlotIndex) -> Result<bool, IrError> {
        self.constrained
            .get(index.as_usize())
            .copied()
            .ok_or(IrError::UnknownSlot(index))
    }

    /// Returns the determinacy of the slot at `index`.
    ///
    /// The geometry of a free slot is
    /// [`Level::Restricted`](crate::Level::Restricted) when the slot
    /// [is constrained](Instantiated::is_constrained).
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownSlot`] if `index` is not below the slot
    /// count.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{
    ///     AtomAddress, Constraint, Corner, Determinacy, Element, ElementIndex, Geometry,
    ///     Identity, Length, Lengths, Level, Segment, SegmentId, Slot, SlotIndex, Spec, Symbolic,
    ///     SymbolicAddress,
    /// };
    /// use mmcif_model::{AtomName, CompId};
    ///
    /// let cysteine = Slot::new(Identity::Exactly(CompId::new("CYS")), Geometry::Free);
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("loop"),
    ///     [
    ///         Element::Slot(cysteine.clone()),
    ///         Element::run(Slot::de_novo(), Length::exact(3)?),
    ///         Element::Slot(cysteine),
    ///     ],
    /// ))?;
    /// let sulfur = |element| {
    ///     let slot = SymbolicAddress::whole(SegmentId::new("loop"), ElementIndex::new(element));
    ///     AtomAddress::new(slot, AtomName::new("SG"))
    /// };
    /// spec.push_constraint(Constraint::<Symbolic>::Bond {
    ///     first: sulfur(0),
    ///     second: sulfur(2),
    /// })?;
    /// let design = spec.instantiate(&Lengths::minimum(&spec))?;
    ///
    /// assert_eq!(
    ///     design.determinacy(SlotIndex::new(0))?,
    ///     Determinacy::new(Level::Fixed, Level::Restricted)
    /// );
    /// assert_eq!(design.determinacy(SlotIndex::new(1))?.corner(), Some(Corner::DeNovo));
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn determinacy(&self, index: SlotIndex) -> Result<Determinacy, IrError> {
        Ok(self.slot(index)?.determinacy(self.is_constrained(index)?))
    }

    /// Returns the slots `address` names, in ascending order.
    ///
    /// An address of a whole element returns every slot the element
    /// produced, which is none for a count of zero.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownSegment`] if no segment is named by
    /// `address`, [`IrError::UnknownElement`] if the segment has no element
    /// at its index, and [`IrError::OffsetOutOfRange`] if its offset is not
    /// below the number of slots the element produced.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{
    ///     Element, ElementIndex, IrError, Length, Lengths, Segment, SegmentId, Slot, SlotIndex,
    ///     Spec, SymbolicAddress,
    /// };
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("binder"),
    ///     [
    ///         Element::Slot(Slot::de_novo()),
    ///         Element::run(Slot::de_novo(), Length::range(2, 4)?),
    ///     ],
    /// ))?;
    /// let design = spec.instantiate(&Lengths::new(&spec, [2])?)?;
    /// let binder = SegmentId::new("binder");
    ///
    /// let run = SymbolicAddress::whole(binder.clone(), ElementIndex::new(1));
    /// assert_eq!(design.resolve(&run)?, [SlotIndex::new(1), SlotIndex::new(2)]);
    ///
    /// let past_end = SymbolicAddress::at(binder, ElementIndex::new(1), 2);
    /// assert!(matches!(
    ///     design.resolve(&past_end),
    ///     Err(IrError::OffsetOutOfRange { offset: 2, .. })
    /// ));
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn resolve(&self, address: &SymbolicAddress) -> Result<Vec<SlotIndex>, IrError> {
        let segment = self.segment_position(&address.segment)?;
        resolve_within(&self.element_slots[segment], address)
    }

    fn segment_position(&self, id: &SegmentId) -> Result<usize, IrError> {
        self.segment_ids
            .iter()
            .position(|known| known == id)
            .ok_or_else(|| IrError::UnknownSegment(id.clone()))
    }
}

use crate::constraint::{AtomAddress, Constraint, Resolved, Scope, SymbolicAddress};
use crate::error::IrError;
use crate::ids::{ElementIndex, GroupId, SegmentId, SegmentIndex, SlotIndex, SlotSpan};
use crate::segment::Length;
use crate::slot::{Determinacy, Slot};
use crate::spec::Spec;

pub const MAX_SLOTS: u32 = 1 << 24;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Choice {
    segment: SegmentIndex,
    element: ElementIndex,
    count: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Lengths {
    choices: Vec<Choice>,
}

impl Lengths {
    pub fn minimum(spec: &Spec) -> Lengths {
        Lengths::chosen(spec, Length::min)
    }

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Origin {
    pub segment: SegmentIndex,
    pub element: ElementIndex,
    pub offset: u32,
}

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
    pub fn slot_count(&self) -> usize {
        self.slots.len()
    }

    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }

    pub fn slots(&self) -> &[Slot] {
        &self.slots
    }

    pub fn slot(&self, index: SlotIndex) -> Result<&Slot, IrError> {
        self.slots
            .get(index.as_usize())
            .ok_or(IrError::UnknownSlot(index))
    }

    pub fn origin(&self, index: SlotIndex) -> Result<Origin, IrError> {
        self.origins
            .get(index.as_usize())
            .copied()
            .ok_or(IrError::UnknownSlot(index))
    }

    pub fn source(&self, index: SlotIndex) -> Result<SymbolicAddress, IrError> {
        let origin = self.origin(index)?;
        Ok(SymbolicAddress::at(
            self.segment_ids[origin.segment.as_usize()].clone(),
            origin.element,
            origin.offset,
        ))
    }

    pub fn segment_ids(&self) -> &[SegmentId] {
        &self.segment_ids
    }

    pub fn segment_slots(&self, id: &SegmentId) -> Result<SlotSpan, IrError> {
        Ok(self.segment_slots[self.segment_position(id)?])
    }

    pub fn groups(&self) -> impl ExactSizeIterator<Item = (&GroupId, &[SlotIndex])> {
        self.groups.iter().map(|(id, slots)| (id, slots.as_slice()))
    }

    pub fn group(&self, id: &GroupId) -> Result<&[SlotIndex], IrError> {
        self.groups
            .iter()
            .find(|(known, _)| known == id)
            .map(|(_, slots)| slots.as_slice())
            .ok_or_else(|| IrError::UnknownGroup(id.clone()))
    }

    pub fn constraints(&self) -> &[Constraint<Resolved>] {
        &self.constraints
    }

    pub fn is_constrained(&self, index: SlotIndex) -> Result<bool, IrError> {
        self.constrained
            .get(index.as_usize())
            .copied()
            .ok_or(IrError::UnknownSlot(index))
    }

    pub fn determinacy(&self, index: SlotIndex) -> Result<Determinacy, IrError> {
        Ok(self.slot(index)?.determinacy(self.is_constrained(index)?))
    }

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

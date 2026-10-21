use std::ops::Range;

use crate::constraint::{AtomAddress, Constraint, Resolved, Scope, SymbolicAddress};
use crate::error::IrError;
use crate::ids::{ElementIndex, GroupId, SegmentId, SegmentIndex, SlotIndex};
use crate::slot::{Determinacy, Slot};
use crate::spec::Spec;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Lengths {
    counts: Vec<u32>,
}

impl Lengths {
    pub fn minimum(spec: &Spec) -> Lengths {
        Lengths {
            counts: spec
                .symbolic_elements()
                .map(|(_, _, length)| length.min())
                .collect(),
        }
    }

    pub fn maximum(spec: &Spec) -> Lengths {
        Lengths {
            counts: spec
                .symbolic_elements()
                .map(|(_, _, length)| length.max())
                .collect(),
        }
    }

    pub fn new(spec: &Spec, counts: impl IntoIterator<Item = u32>) -> Result<Lengths, IrError> {
        let lengths = Lengths {
            counts: counts.into_iter().collect(),
        };
        lengths.check(spec)?;
        Ok(lengths)
    }

    pub fn counts(&self) -> &[u32] {
        &self.counts
    }

    fn check(&self, spec: &Spec) -> Result<(), IrError> {
        let expected = spec.symbolic_elements().count();
        if self.counts.len() != expected {
            return Err(IrError::LengthCountMismatch {
                expected,
                actual: self.counts.len(),
            });
        }
        for (&chosen, (_, _, length)) in self.counts.iter().zip(spec.symbolic_elements()) {
            if !length.contains(chosen) {
                return Err(IrError::LengthOutOfRange {
                    chosen,
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
                .counts
                .iter()
                .map(|&count| u64::from(count))
                .sum::<u64>()
    }
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
    segment_slots: Vec<Range<u32>>,
    element_slots: Vec<Vec<Range<u32>>>,
    groups: Vec<(GroupId, Vec<SlotIndex>)>,
    constraints: Vec<Constraint<Resolved>>,
    constrained: Vec<bool>,
}

impl Spec {
    pub fn instantiate(&self, lengths: &Lengths) -> Result<Instantiated, IrError> {
        lengths.check(self)?;
        let total = lengths.total_slots(self);
        if total > u64::from(u32::MAX) {
            return Err(IrError::TooManySlots(u32::MAX));
        }
        let mut slots = Vec::with_capacity(total as usize);
        let mut origins = Vec::with_capacity(total as usize);
        let mut element_slots = Vec::with_capacity(self.segments().len());
        let mut segment_slots = Vec::with_capacity(self.segments().len());
        let mut choices = lengths.counts.iter();
        for (segment_index, segment) in self.segments().iter().enumerate() {
            let segment_start = slots.len() as u32;
            let mut ranges = Vec::with_capacity(segment.elements().len());
            for (element_index, element) in segment.elements().iter().enumerate() {
                let count = if element.is_symbolic() {
                    *choices
                        .next()
                        .expect("lengths cover every symbolic element")
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
                ranges.push(element_start..slots.len() as u32);
            }
            element_slots.push(ranges);
            segment_slots.push(segment_start..slots.len() as u32);
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
    elements: &[Range<u32>],
    address: &SymbolicAddress,
) -> Result<Vec<SlotIndex>, IrError> {
    let range =
        elements
            .get(address.element.as_usize())
            .ok_or_else(|| IrError::UnknownElement {
                segment: address.segment.clone(),
                element: address.element,
            })?;
    match address.offset {
        None => Ok((range.start..range.end).map(SlotIndex::new).collect()),
        Some(offset) if offset < range.end - range.start => {
            Ok(vec![SlotIndex::new(range.start + offset)])
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
    element_slots: &[Vec<Range<u32>>],
    address: &SymbolicAddress,
) -> Result<Vec<SlotIndex>, IrError> {
    let segment = spec.segment_index(&address.segment)?;
    resolve_within(&element_slots[segment.as_usize()], address)
}

fn resolve_one(
    spec: &Spec,
    element_slots: &[Vec<Range<u32>>],
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
    element_slots: &[Vec<Range<u32>>],
    groups: &[(GroupId, Vec<SlotIndex>)],
    constraint: &Constraint<crate::constraint::Symbolic>,
) -> Result<Constraint<Resolved>, IrError> {
    match constraint {
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

    pub fn segment_slots(&self, id: &SegmentId) -> Result<Range<usize>, IrError> {
        let index = self
            .segment_ids
            .iter()
            .position(|known| known == id)
            .ok_or_else(|| IrError::UnknownSegment(id.clone()))?;
        let range = &self.segment_slots[index];
        Ok(range.start as usize..range.end as usize)
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
        let index = self
            .segment_ids
            .iter()
            .position(|known| known == &address.segment)
            .ok_or_else(|| IrError::UnknownSegment(address.segment.clone()))?;
        resolve_within(&self.element_slots[index], address)
    }
}

use std::collections::BTreeSet;

use crate::constraint::{Constraint, Scope, Symbolic, SymbolicAddress};
use crate::error::IrError;
use crate::ids::{ElementIndex, GroupId, SegmentId, SegmentIndex};
use crate::segment::{Length, Segment};
use crate::slot::Determinacy;

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Spec {
    segments: Vec<Segment>,
    groups: Vec<(GroupId, Vec<SymbolicAddress>)>,
    constraints: Vec<Constraint<Symbolic>>,
}

impl Spec {
    pub fn new() -> Spec {
        Spec::default()
    }

    pub fn push_segment(&mut self, segment: Segment) -> Result<SegmentIndex, IrError> {
        if self.segments.iter().any(|known| known.id() == segment.id()) {
            return Err(IrError::DuplicateSegment(segment.id().clone()));
        }
        self.segments.push(segment);
        Ok(SegmentIndex::new((self.segments.len() - 1) as u32))
    }

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

    pub fn segments(&self) -> &[Segment] {
        &self.segments
    }

    pub fn segment_index(&self, id: &SegmentId) -> Result<SegmentIndex, IrError> {
        self.segments
            .iter()
            .position(|segment| segment.id() == id)
            .map(|index| SegmentIndex::new(index as u32))
            .ok_or_else(|| IrError::UnknownSegment(id.clone()))
    }

    pub fn segment(&self, id: &SegmentId) -> Result<&Segment, IrError> {
        let index = self.segment_index(id)?;
        Ok(&self.segments[index.as_usize()])
    }

    pub fn groups(&self) -> impl ExactSizeIterator<Item = (&GroupId, &[SymbolicAddress])> {
        self.groups
            .iter()
            .map(|(id, members)| (id, members.as_slice()))
    }

    pub fn group(&self, id: &GroupId) -> Result<&[SymbolicAddress], IrError> {
        self.groups
            .iter()
            .find(|(known, _)| known == id)
            .map(|(_, members)| members.as_slice())
            .ok_or_else(|| IrError::UnknownGroup(id.clone()))
    }

    pub fn constraints(&self) -> &[Constraint<Symbolic>] {
        &self.constraints
    }

    pub fn length_bounds(&self) -> (u64, u64) {
        self.segments.iter().fold((0, 0), |(min, max), segment| {
            let (segment_min, segment_max) = segment.length_bounds();
            (min + segment_min, max + segment_max)
        })
    }

    pub fn is_symbolic(&self) -> bool {
        self.segments.iter().any(Segment::is_symbolic)
    }

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
                Constraint::Bond { first, second }
                | Constraint::Distance { first, second, .. } => {
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

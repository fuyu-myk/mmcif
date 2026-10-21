use crate::error::IrError;
use crate::ids::{ElementIndex, SegmentId};
use crate::slot::Slot;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Length {
    Exact(u32),
    Range { min: u32, max: u32 },
}

impl Length {
    pub fn exact(count: u32) -> Result<Length, IrError> {
        if count == 0 {
            return Err(IrError::EmptyRun);
        }
        Ok(Length::Exact(count))
    }

    pub fn range(min: u32, max: u32) -> Result<Length, IrError> {
        if min > max {
            return Err(IrError::DescendingRange { min, max });
        }
        if max == 0 {
            return Err(IrError::EmptyRun);
        }
        if min == max {
            return Ok(Length::Exact(min));
        }
        Ok(Length::Range { min, max })
    }

    pub fn min(self) -> u32 {
        match self {
            Length::Exact(count) => count,
            Length::Range { min, .. } => min,
        }
    }

    pub fn max(self) -> u32 {
        match self {
            Length::Exact(count) => count,
            Length::Range { max, .. } => max,
        }
    }

    pub fn contains(self, count: u32) -> bool {
        self.min() <= count && count <= self.max()
    }

    pub fn is_symbolic(self) -> bool {
        self.min() != self.max()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Element {
    Slot(Slot),
    Run { slot: Slot, length: Length },
}

impl Element {
    pub fn run(slot: Slot, length: Length) -> Element {
        match length {
            Length::Exact(1) => Element::Slot(slot),
            _ => Element::Run { slot, length },
        }
    }

    pub fn slot(&self) -> &Slot {
        match self {
            Element::Slot(slot) => slot,
            Element::Run { slot, .. } => slot,
        }
    }

    pub fn length(&self) -> Length {
        match self {
            Element::Slot(_) => Length::Exact(1),
            Element::Run { length, .. } => *length,
        }
    }

    pub fn is_symbolic(&self) -> bool {
        self.length().is_symbolic()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    id: SegmentId,
    elements: Vec<Element>,
}

impl Segment {
    pub fn new(id: SegmentId, elements: impl IntoIterator<Item = Element>) -> Segment {
        Segment {
            id,
            elements: elements.into_iter().collect(),
        }
    }

    pub fn id(&self) -> &SegmentId {
        &self.id
    }

    pub fn elements(&self) -> &[Element] {
        &self.elements
    }

    pub fn element(&self, index: ElementIndex) -> Result<&Element, IrError> {
        self.elements
            .get(index.as_usize())
            .ok_or_else(|| IrError::UnknownElement {
                segment: self.id.clone(),
                element: index,
            })
    }

    pub fn length_bounds(&self) -> (u64, u64) {
        self.elements.iter().fold((0, 0), |(min, max), element| {
            (
                min + u64::from(element.length().min()),
                max + u64::from(element.length().max()),
            )
        })
    }

    pub fn is_symbolic(&self) -> bool {
        self.elements.iter().any(Element::is_symbolic)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::slot::Identity;
    use mmcif_model::CompId;

    fn cysteine() -> Slot {
        Slot::new(
            Identity::Exactly(CompId::new("CYS")),
            crate::slot::Geometry::Free,
        )
    }

    #[test]
    fn degenerate_ranges_collapse_to_exact_lengths() {
        assert_eq!(Length::range(4, 4).unwrap(), Length::Exact(4));
        assert!(!Length::range(4, 4).unwrap().is_symbolic());
        assert!(Length::range(2, 5).unwrap().is_symbolic());
        assert_eq!(
            Length::range(5, 2),
            Err(IrError::DescendingRange { min: 5, max: 2 })
        );
        assert_eq!(Length::range(0, 0), Err(IrError::EmptyRun));
        assert_eq!(Length::exact(0), Err(IrError::EmptyRun));
        assert_eq!(Length::range(0, 3).unwrap().min(), 0);
    }

    #[test]
    fn a_run_of_one_is_a_single_slot() {
        assert_eq!(
            Element::run(Slot::de_novo(), Length::Exact(1)),
            Element::Slot(Slot::de_novo())
        );
        assert_eq!(Element::Slot(Slot::de_novo()).length(), Length::Exact(1));
    }

    #[test]
    fn a_disulfide_pattern_bounds_its_own_length() {
        let segment = Segment::new(
            SegmentId::new("binder"),
            [
                Element::run(Slot::de_novo(), Length::range(3, 5).unwrap()),
                Element::Slot(cysteine()),
                Element::run(Slot::de_novo(), Length::exact(2).unwrap()),
                Element::Slot(cysteine()),
                Element::run(Slot::de_novo(), Length::range(1, 4).unwrap()),
            ],
        );
        assert_eq!(
            segment.length_bounds(),
            (3 + 1 + 2 + 1 + 1, 5 + 1 + 2 + 1 + 4)
        );
        assert!(segment.is_symbolic());
        assert_eq!(segment.elements().len(), 5);
        assert_eq!(
            segment.element(ElementIndex::new(1)).unwrap().slot(),
            &cysteine()
        );
    }

    #[test]
    fn an_unknown_element_names_its_segment() {
        let segment = Segment::new(SegmentId::new("a"), [Element::Slot(Slot::de_novo())]);
        assert_eq!(
            segment.element(ElementIndex::new(3)),
            Err(IrError::UnknownElement {
                segment: SegmentId::new("a"),
                element: ElementIndex::new(3)
            })
        );
    }

    #[test]
    fn a_fully_concrete_segment_is_not_symbolic() {
        let segment = Segment::new(
            SegmentId::new("target"),
            [Element::run(Slot::de_novo(), Length::exact(9).unwrap())],
        );
        assert!(!segment.is_symbolic());
        assert_eq!(segment.length_bounds(), (9, 9));
    }
}

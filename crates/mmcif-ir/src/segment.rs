use crate::error::IrError;
use crate::ids::{ElementIndex, SegmentId};
use crate::slot::Slot;

/// The number of slots an [`Element`] produces, either exact or an inclusive
/// range.
///
/// A range is a symbolic length: the count is chosen when the specification
/// is instantiated with [`Lengths`](crate::Lengths).
///
/// # Examples
///
/// ```
/// use mmcif_ir::Length;
///
/// let length = Length::range(3, 5)?;
///
/// assert_eq!((length.min(), length.max()), (3, 5));
/// assert!(length.contains(4));
/// assert!(length.is_symbolic());
/// assert!(!Length::exact(4)?.is_symbolic());
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Length {
    min: u32,
    max: u32,
}

impl Length {
    const ONE: Length = Length { min: 1, max: 1 };

    /// Creates a length of exactly `count` slots.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::EmptyRun`] if `count` is zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{IrError, Length};
    ///
    /// let length = Length::exact(4)?;
    /// assert_eq!((length.min(), length.max()), (4, 4));
    /// assert_eq!(Length::exact(0), Err(IrError::EmptyRun));
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn exact(count: u32) -> Result<Length, IrError> {
        Length::range(count, count)
    }

    /// Creates a length of `min` to `max` slots inclusive.
    ///
    /// The minimum may be zero.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::DescendingRange`] if `min` is greater than `max`,
    /// and [`IrError::EmptyRun`] if `max` is zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{IrError, Length};
    ///
    /// let length = Length::range(0, 4)?;
    /// assert_eq!((length.min(), length.max()), (0, 4));
    /// assert_eq!(Length::range(4, 4)?, Length::exact(4)?);
    /// assert_eq!(
    ///     Length::range(5, 2),
    ///     Err(IrError::DescendingRange { min: 5, max: 2 })
    /// );
    /// assert_eq!(Length::range(0, 0), Err(IrError::EmptyRun));
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn range(min: u32, max: u32) -> Result<Length, IrError> {
        if min > max {
            return Err(IrError::DescendingRange { min, max });
        }
        if max == 0 {
            return Err(IrError::EmptyRun);
        }
        Ok(Length { min, max })
    }

    /// Returns the smallest number of slots the length allows.
    pub fn min(self) -> u32 {
        self.min
    }

    /// Returns the largest number of slots the length allows.
    pub fn max(self) -> u32 {
        self.max
    }

    /// Returns `true` if `count` lies between the minimum and the maximum,
    /// inclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::Length;
    ///
    /// let length = Length::range(3, 5)?;
    ///
    /// assert!(length.contains(3) && length.contains(5));
    /// assert!(!length.contains(6));
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn contains(self, count: u32) -> bool {
        self.min() <= count && count <= self.max()
    }

    /// Returns `true` if the minimum and the maximum differ.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::Length;
    ///
    /// assert!(Length::range(3, 5)?.is_symbolic());
    /// assert!(!Length::range(4, 4)?.is_symbolic());
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn is_symbolic(self) -> bool {
        self.min() != self.max()
    }
}

/// An entry of a [`Segment`]: a single slot, or a run of copies of one slot.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{Element, Length, Slot};
///
/// let run = Element::run(Slot::de_novo(), Length::range(3, 5)?);
///
/// assert_eq!(run.slot(), &Slot::de_novo());
/// assert_eq!(run.length(), Length::range(3, 5)?);
/// assert_eq!(Element::Slot(Slot::de_novo()).length(), Length::exact(1)?);
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Element {
    /// A single slot.
    Slot(Slot),
    /// Copies of one slot, made by [`Element::run`].
    #[non_exhaustive]
    Run {
        /// The slot each copy repeats.
        slot: Slot,
        /// The number of copies.
        length: Length,
    },
}

impl Element {
    /// Creates an element that repeats `slot` `length` times.
    ///
    /// Returns [`Element::Slot`] when `length` is exactly one, and
    /// [`Element::Run`] otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Length, Slot};
    ///
    /// let single = Element::run(Slot::de_novo(), Length::exact(1)?);
    /// assert_eq!(single, Element::Slot(Slot::de_novo()));
    ///
    /// let run = Element::run(Slot::de_novo(), Length::exact(4)?);
    /// assert!(matches!(run, Element::Run { .. }));
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn run(slot: Slot, length: Length) -> Element {
        if length == Length::ONE {
            Element::Slot(slot)
        } else {
            Element::Run { slot, length }
        }
    }

    /// Returns the slot the element holds or repeats.
    pub fn slot(&self) -> &Slot {
        match self {
            Element::Slot(slot) => slot,
            Element::Run { slot, .. } => slot,
        }
    }

    /// Returns the length of the element, which is exactly one for
    /// [`Element::Slot`].
    pub fn length(&self) -> Length {
        match self {
            Element::Slot(_) => Length::ONE,
            Element::Run { length, .. } => *length,
        }
    }

    /// Returns `true` if the length of the element is symbolic.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Length, Slot};
    ///
    /// assert!(Element::run(Slot::de_novo(), Length::range(3, 5)?).is_symbolic());
    /// assert!(!Element::Slot(Slot::de_novo()).is_symbolic());
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn is_symbolic(&self) -> bool {
        self.length().is_symbolic()
    }
}

/// A named, ordered list of elements.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{Element, Length, Segment, SegmentId, Slot};
///
/// let segment = Segment::new(
///     SegmentId::new("binder"),
///     [
///         Element::run(Slot::de_novo(), Length::range(3, 5)?),
///         Element::Slot(Slot::de_novo()),
///     ],
/// );
///
/// assert_eq!(segment.id().as_str(), "binder");
/// assert_eq!(segment.elements().len(), 2);
/// assert_eq!(segment.length_bounds(), (4, 6));
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    id: SegmentId,
    elements: Vec<Element>,
}

impl Segment {
    /// Creates a segment from its id and its elements, in order.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Segment, SegmentId, Slot};
    ///
    /// let segment = Segment::new(SegmentId::new("loop"), [Element::Slot(Slot::de_novo())]);
    ///
    /// assert_eq!(segment.elements(), [Element::Slot(Slot::de_novo())]);
    /// ```
    pub fn new(id: SegmentId, elements: impl IntoIterator<Item = Element>) -> Segment {
        Segment {
            id,
            elements: elements.into_iter().collect(),
        }
    }

    /// Returns the id of the segment.
    pub fn id(&self) -> &SegmentId {
        &self.id
    }

    /// Returns the elements of the segment, in order.
    pub fn elements(&self) -> &[Element] {
        &self.elements
    }

    /// Returns the element at `index`.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownElement`] if the segment has no element at
    /// `index`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, ElementIndex, IrError, Segment, SegmentId, Slot};
    ///
    /// let segment = Segment::new(SegmentId::new("loop"), [Element::Slot(Slot::de_novo())]);
    ///
    /// assert_eq!(segment.element(ElementIndex::new(0))?.slot(), &Slot::de_novo());
    /// assert_eq!(
    ///     segment.element(ElementIndex::new(1)),
    ///     Err(IrError::UnknownElement {
    ///         segment: SegmentId::new("loop"),
    ///         element: ElementIndex::new(1),
    ///     })
    /// );
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn element(&self, index: ElementIndex) -> Result<&Element, IrError> {
        self.elements
            .get(index.as_usize())
            .ok_or_else(|| IrError::UnknownElement {
                segment: self.id.clone(),
                element: index,
            })
    }

    /// Returns the smallest and the largest number of slots the segment can
    /// produce, as `(min, max)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Length, Segment, SegmentId, Slot};
    ///
    /// let segment = Segment::new(
    ///     SegmentId::new("loop"),
    ///     [
    ///         Element::run(Slot::de_novo(), Length::range(0, 4)?),
    ///         Element::run(Slot::de_novo(), Length::exact(2)?),
    ///     ],
    /// );
    ///
    /// assert_eq!(segment.length_bounds(), (2, 6));
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn length_bounds(&self) -> (u64, u64) {
        self.elements.iter().fold((0, 0), |(min, max), element| {
            (
                min + u64::from(element.length().min()),
                max + u64::from(element.length().max()),
            )
        })
    }

    /// Returns `true` if any element of the segment is symbolic.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Length, Segment, SegmentId, Slot};
    ///
    /// let fixed = Segment::new(
    ///     SegmentId::new("target"),
    ///     [Element::run(Slot::de_novo(), Length::exact(9)?)],
    /// );
    /// let open = Segment::new(
    ///     SegmentId::new("binder"),
    ///     [Element::run(Slot::de_novo(), Length::range(8, 12)?)],
    /// );
    ///
    /// assert!(!fixed.is_symbolic());
    /// assert!(open.is_symbolic());
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn is_symbolic(&self) -> bool {
        self.elements.iter().any(Element::is_symbolic)
    }
}

/// ```compile_fail
/// use mmcif_ir::Length;
///
/// let empty = Length::Exact(0);
/// ```
///
/// ```compile_fail
/// use mmcif_ir::Length;
///
/// let descending = Length::Range { min: 5, max: 2 };
/// ```
///
/// ```compile_fail,E0639
/// use mmcif_ir::{Element, Length, Slot};
///
/// let single = Element::Run {
///     slot: Slot::de_novo(),
///     length: Length::exact(1).unwrap(),
/// };
/// ```
#[cfg(doctest)]
struct LengthsAndRunsAreValidated;

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
        assert_eq!(Length::range(4, 4).unwrap(), Length::exact(4).unwrap());
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
            Element::run(Slot::de_novo(), Length::exact(1).unwrap()),
            Element::Slot(Slot::de_novo())
        );
        assert_eq!(
            Element::Slot(Slot::de_novo()).length(),
            Length::exact(1).unwrap()
        );
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

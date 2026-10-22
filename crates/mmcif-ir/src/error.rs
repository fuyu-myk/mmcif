use crate::constraint::SymbolicAddress;
use crate::ids::{ElementIndex, GroupId, SegmentId, SlotIndex, TokenIndex};

/// An error from building, querying or instantiating a [`Spec`](crate::Spec),
/// or from building or querying a [`TokenLayout`](crate::TokenLayout).
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum IrError {
    /// A set of allowed components names no component.
    #[error("an allowed identity set must name at least one component")]
    EmptyAlphabet,
    /// A length allows no slot: an exact length of zero or a range whose
    /// maximum is zero.
    #[error("a run must allow at least one slot")]
    EmptyRun,
    /// A length range has a minimum greater than its maximum.
    #[error("length range {min}..={max} is descending")]
    DescendingRange {
        /// The minimum of the range.
        min: u32,
        /// The maximum of the range.
        max: u32,
    },
    /// Distance bounds have a NaN or infinite end, a negative minimum, or a
    /// minimum greater than the maximum.
    #[error("distance bounds {min}..={max} are not a finite ascending range")]
    InvalidBounds {
        /// The minimum of the bounds.
        min: f64,
        /// The maximum of the bounds.
        max: f64,
    },
    /// A segment with this id is already in the specification.
    #[error("segment {0} is already defined")]
    DuplicateSegment(SegmentId),
    /// A group with this id is already in the specification.
    #[error("group {0} is already defined")]
    DuplicateGroup(GroupId),
    /// No segment has this id.
    #[error("no segment {0}")]
    UnknownSegment(SegmentId),
    /// No group has this id.
    #[error("no group {0}")]
    UnknownGroup(GroupId),
    /// A segment has no element at this index.
    #[error("segment {segment} has no element {element}")]
    UnknownElement {
        /// The segment that was searched.
        segment: SegmentId,
        /// The index that has no element.
        element: ElementIndex,
    },
    /// An offset lies at or past the number of slots its element holds.
    #[error("segment {segment} element {element} has no slot at offset {offset}")]
    OffsetOutOfRange {
        /// The segment of the address.
        segment: SegmentId,
        /// The element of the address.
        element: ElementIndex,
        /// The offset that has no slot.
        offset: u32,
    },
    /// An address that must name one slot names a whole element that holds
    /// a different number of slots.
    #[error(
        "segment {segment} element {element} holds {slots} slots where one address is required"
    )]
    AddressNotSingular {
        /// The segment of the address.
        segment: SegmentId,
        /// The element of the address.
        element: ElementIndex,
        /// The maximum length of the element in a specification, or the
        /// number of slots it produced in an instantiation.
        slots: u32,
    },
    /// The number of chosen lengths differs from the number of symbolic
    /// elements in the specification.
    #[error("expected {expected} chosen lengths, got {actual}")]
    LengthCountMismatch {
        /// The number of symbolic elements in the specification.
        expected: usize,
        /// The number of chosen lengths.
        actual: usize,
    },
    /// A chosen length lies outside the range of its element.
    #[error("length {chosen} is outside {min}..={max}")]
    LengthOutOfRange {
        /// The chosen length.
        chosen: u32,
        /// The minimum length of the element.
        min: u32,
        /// The maximum length of the element.
        max: u32,
    },
    /// No slot has this index.
    #[error("no slot at index {0}")]
    UnknownSlot(SlotIndex),
    /// No token has this index.
    #[error("no token at index {0}")]
    UnknownToken(TokenIndex),
    /// A slot exists but has no token in the layout.
    #[error("slot {0} has no token")]
    SlotNotTokenized(SlotIndex),
    /// A slot appears more than once in the order of a token layout.
    #[error("slot {0} appears twice in the token layout")]
    DuplicateToken(SlotIndex),
    /// The number of slots exceeds the given maximum.
    #[error("the specification holds more than {0} slots")]
    TooManySlots(u32),
    /// An address names a whole element whose slots differ in determinacy.
    #[error("segment {segment} element {element} holds slots of differing determinacy")]
    MixedDeterminacy {
        /// The segment of the address.
        segment: SegmentId,
        /// The element of the address.
        element: ElementIndex,
    },
    /// Chosen lengths belong to a specification whose symbolic elements sit
    /// at different addresses.
    #[error("chosen lengths describe {actual} but this specification expects {expected}")]
    LengthsFromAnotherSpec {
        /// The symbolic element this specification has at the first
        /// mismatched position.
        expected: SymbolicAddress,
        /// The element the lengths were chosen for at that position.
        actual: SymbolicAddress,
    },
}

use crate::constraint::SymbolicAddress;
use crate::ids::{ElementIndex, GroupId, SegmentId, SlotIndex, TokenIndex};

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum IrError {
    #[error("an allowed identity set must name at least one component")]
    EmptyAlphabet,
    #[error("a run must allow at least one slot")]
    EmptyRun,
    #[error("length range {min}..={max} is descending")]
    DescendingRange { min: u32, max: u32 },
    #[error("distance bounds {min}..={max} are not a finite ascending range")]
    InvalidBounds { min: f64, max: f64 },
    #[error("segment {0} is already defined")]
    DuplicateSegment(SegmentId),
    #[error("group {0} is already defined")]
    DuplicateGroup(GroupId),
    #[error("no segment {0}")]
    UnknownSegment(SegmentId),
    #[error("no group {0}")]
    UnknownGroup(GroupId),
    #[error("segment {segment} has no element {element}")]
    UnknownElement {
        segment: SegmentId,
        element: ElementIndex,
    },
    #[error("segment {segment} element {element} has no slot at offset {offset}")]
    OffsetOutOfRange {
        segment: SegmentId,
        element: ElementIndex,
        offset: u32,
    },
    #[error(
        "segment {segment} element {element} holds {slots} slots where one address is required"
    )]
    AddressNotSingular {
        segment: SegmentId,
        element: ElementIndex,
        slots: u32,
    },
    #[error("expected {expected} chosen lengths, got {actual}")]
    LengthCountMismatch { expected: usize, actual: usize },
    #[error("length {chosen} is outside {min}..={max}")]
    LengthOutOfRange { chosen: u32, min: u32, max: u32 },
    #[error("no slot at index {0}")]
    UnknownSlot(SlotIndex),
    #[error("no token at index {0}")]
    UnknownToken(TokenIndex),
    #[error("slot {0} has no token")]
    SlotNotTokenized(SlotIndex),
    #[error("slot {0} appears twice in the token layout")]
    DuplicateToken(SlotIndex),
    #[error("the specification holds more than {0} slots")]
    TooManySlots(u32),
    #[error("segment {segment} element {element} holds slots of differing determinacy")]
    MixedDeterminacy {
        segment: SegmentId,
        element: ElementIndex,
    },
    #[error("chosen lengths describe {actual} but this specification expects {expected}")]
    LengthsFromAnotherSpec {
        expected: SymbolicAddress,
        actual: SymbolicAddress,
    },
}

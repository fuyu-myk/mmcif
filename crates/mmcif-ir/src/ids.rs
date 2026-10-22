use std::fmt;

use mmcif_model::{AuthResidueKey, LabelResidueKey};

use crate::error::IrError;

macro_rules! index_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
        pub struct $name(u32);

        impl $name {
            /// Creates an index with the given value.
            ///
            /// # Examples
            ///
            /// ```
            #[doc = concat!("use mmcif_ir::", stringify!($name), ";")]
            ///
            #[doc = concat!("let index = ", stringify!($name), "::new(3);")]
            ///
            /// assert_eq!(index.get(), 3);
            /// assert_eq!(index.to_string(), "3");
            /// ```
            pub fn new(value: u32) -> Self {
                $name(value)
            }

            /// Returns the value of the index.
            pub fn get(self) -> u32 {
                self.0
            }

            /// Returns the value of the index as a `usize`.
            pub fn as_usize(self) -> usize {
                self.0 as usize
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                self.0.fmt(f)
            }
        }
    };
}

index_id! {
    /// The position of a segment in a [`Spec`](crate::Spec), counted from
    /// zero in the order the segments were pushed.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Segment, SegmentId, SegmentIndex, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(SegmentId::new("target"), []))?;
    /// let binder = spec.push_segment(Segment::new(SegmentId::new("binder"), []))?;
    ///
    /// assert_eq!(binder, SegmentIndex::new(1));
    /// assert_eq!(spec.segments()[binder.as_usize()].id().as_str(), "binder");
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    SegmentIndex
}

index_id! {
    /// The position of an element in a [`Segment`](crate::Segment), counted
    /// from zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, ElementIndex, Segment, SegmentId, Slot};
    ///
    /// let segment = Segment::new(SegmentId::new("loop"), [Element::Slot(Slot::de_novo())]);
    ///
    /// assert!(segment.element(ElementIndex::new(0)).is_ok());
    /// assert!(segment.element(ElementIndex::new(1)).is_err());
    /// ```
    ElementIndex
}

index_id! {
    /// The position of a slot in an [`Instantiated`](crate::Instantiated)
    /// specification, counted from zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Length, Lengths, Segment, SegmentId, Slot, SlotIndex, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(
    ///     SegmentId::new("binder"),
    ///     [Element::run(Slot::de_novo(), Length::exact(3)?)],
    /// ))?;
    /// let design = spec.instantiate(&Lengths::minimum(&spec))?;
    ///
    /// assert!(design.slot(SlotIndex::new(2)).is_ok());
    /// assert!(design.slot(SlotIndex::new(3)).is_err());
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    SlotIndex
}

index_id! {
    /// The position of a token in a [`TokenLayout`], counted from zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{SlotIndex, TokenIndex, TokenLayout};
    ///
    /// let layout = TokenLayout::new(3, [SlotIndex::new(2), SlotIndex::new(0)])?;
    ///
    /// assert_eq!(layout.slot_for(TokenIndex::new(0))?, SlotIndex::new(2));
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    TokenIndex
}

macro_rules! name_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
        pub struct $name(Box<str>);

        impl $name {
            /// Creates an id from a name.
            ///
            /// # Examples
            ///
            /// ```
            #[doc = concat!("use mmcif_ir::", stringify!($name), ";")]
            ///
            #[doc = concat!("let id = ", stringify!($name), "::new(\"core\");")]
            ///
            /// assert_eq!(id.as_str(), "core");
            /// assert_eq!(id.to_string(), "core");
            /// ```
            pub fn new(value: impl Into<Box<str>>) -> Self {
                $name(value.into())
            }

            /// Returns the name as a string slice.
            pub fn as_str(&self) -> &str {
                &self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.0)
            }
        }

        impl From<&str> for $name {
            fn from(value: &str) -> Self {
                $name::new(value)
            }
        }
    };
}

name_id! {
    /// The name of a [`Segment`](crate::Segment), unique within a
    /// [`Spec`](crate::Spec).
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Segment, SegmentId, Spec};
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(SegmentId::new("binder"), []))?;
    ///
    /// assert!(spec.segment(&SegmentId::from("binder")).is_ok());
    /// assert!(spec.push_segment(Segment::new(SegmentId::new("binder"), [])).is_err());
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    SegmentId
}

name_id! {
    /// The name of a group of addresses, unique within a
    /// [`Spec`](crate::Spec).
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{
    ///     Element, ElementIndex, GroupId, Segment, SegmentId, Slot, Spec, SymbolicAddress,
    /// };
    ///
    /// let mut spec = Spec::new();
    /// spec.push_segment(Segment::new(SegmentId::new("loop"), [Element::Slot(Slot::de_novo())]))?;
    /// let member = SymbolicAddress::whole(SegmentId::new("loop"), ElementIndex::new(0));
    /// spec.define_group(GroupId::new("core"), [member])?;
    ///
    /// assert_eq!(spec.group(&GroupId::from("core"))?.len(), 1);
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    GroupId
}

/// A residue in an existing structure, named by its author numbering or its
/// label numbering.
///
/// A reference displays as `auth` or `label` followed by its key.
///
/// # Examples
///
/// ```
/// use mmcif_ir::ResidueRef;
/// use mmcif_model::{LabelResidueKey, LabelSeqId};
///
/// let residue = ResidueRef::Label(LabelResidueKey {
///     asym: "A".into(),
///     seq: LabelSeqId::new(3).unwrap(),
/// });
///
/// assert!(residue.as_label().is_some());
/// assert!(residue.as_auth().is_none());
/// assert_eq!(residue.to_string(), "label A/3");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ResidueRef {
    /// A residue named by its author chain, sequence number and insertion
    /// code.
    Auth(AuthResidueKey),
    /// A residue named by its label asym id and label sequence number.
    Label(LabelResidueKey),
}

impl ResidueRef {
    /// Returns the author key if the residue is named by author numbering.
    pub fn as_auth(&self) -> Option<&AuthResidueKey> {
        match self {
            ResidueRef::Auth(key) => Some(key),
            ResidueRef::Label(_) => None,
        }
    }

    /// Returns the label key if the residue is named by label numbering.
    pub fn as_label(&self) -> Option<&LabelResidueKey> {
        match self {
            ResidueRef::Label(key) => Some(key),
            ResidueRef::Auth(_) => None,
        }
    }
}

impl fmt::Display for ResidueRef {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ResidueRef::Auth(key) => write!(f, "auth {key}"),
            ResidueRef::Label(key) => write!(f, "label {key}"),
        }
    }
}

/// A half-open range of slot indices, from [`start`](SlotSpan::start)
/// inclusive to [`end`](SlotSpan::end) exclusive.
///
/// Spans are returned by
/// [`Instantiated::segment_slots`](crate::Instantiated::segment_slots).
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
///     [Element::run(Slot::de_novo(), Length::exact(3)?)],
/// ))?;
/// let design = spec.instantiate(&Lengths::minimum(&spec))?;
///
/// let binder = design.segment_slots(&SegmentId::new("binder"))?;
/// assert_eq!((binder.start(), binder.end()), (SlotIndex::new(1), SlotIndex::new(4)));
/// assert_eq!(binder.len(), 3);
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SlotSpan {
    start: u32,
    end: u32,
}

impl SlotSpan {
    pub(crate) fn new(start: u32, end: u32) -> SlotSpan {
        SlotSpan { start, end }
    }

    /// Returns the inclusive start of the span.
    pub fn start(self) -> SlotIndex {
        SlotIndex::new(self.start)
    }

    /// Returns the exclusive end of the span.
    pub fn end(self) -> SlotIndex {
        SlotIndex::new(self.end)
    }

    /// Returns the number of slots in the span.
    pub fn len(self) -> usize {
        (self.end - self.start) as usize
    }

    /// Returns `true` if the span holds no slots.
    pub fn is_empty(self) -> bool {
        self.start == self.end
    }

    /// Returns `true` if `slot` lies in the span.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Length, Lengths, Segment, SegmentId, Slot, SlotIndex, Spec};
    ///
    /// # let mut spec = Spec::new();
    /// # spec.push_segment(Segment::new(SegmentId::new("target"), [Element::Slot(Slot::de_novo())]))?;
    /// # spec.push_segment(Segment::new(
    /// #     SegmentId::new("binder"),
    /// #     [Element::run(Slot::de_novo(), Length::exact(3)?)],
    /// # ))?;
    /// # let design = spec.instantiate(&Lengths::minimum(&spec))?;
    /// let binder = design.segment_slots(&SegmentId::new("binder"))?;
    ///
    /// assert!(binder.contains(SlotIndex::new(1)));
    /// assert!(binder.contains(SlotIndex::new(3)));
    /// assert!(!binder.contains(SlotIndex::new(4)));
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn contains(self, slot: SlotIndex) -> bool {
        self.start <= slot.get() && slot.get() < self.end
    }

    /// Returns an iterator over the slots in the span, in ascending order.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Element, Length, Lengths, Segment, SegmentId, Slot, SlotIndex, Spec};
    ///
    /// # let mut spec = Spec::new();
    /// # spec.push_segment(Segment::new(SegmentId::new("target"), [Element::Slot(Slot::de_novo())]))?;
    /// # spec.push_segment(Segment::new(
    /// #     SegmentId::new("binder"),
    /// #     [Element::run(Slot::de_novo(), Length::exact(3)?)],
    /// # ))?;
    /// # let design = spec.instantiate(&Lengths::minimum(&spec))?;
    /// let binder = design.segment_slots(&SegmentId::new("binder"))?;
    /// let slots: Vec<u32> = binder.iter().map(SlotIndex::get).collect();
    ///
    /// assert_eq!(slots, [1, 2, 3]);
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn iter(self) -> impl ExactSizeIterator<Item = SlotIndex> {
        (self.start..self.end).map(SlotIndex::new)
    }
}

/// A numbering of a subset of slots as tokens.
///
/// Tokens are numbered from zero in the order their slots were given. Each
/// slot has at most one token, and a slot left out of the order has none.
/// The default layout covers no slots.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{SlotIndex, TokenIndex, TokenLayout};
///
/// let order = [SlotIndex::new(2), SlotIndex::new(0), SlotIndex::new(3)];
/// let layout = TokenLayout::new(4, order)?;
///
/// assert_eq!(layout.token_count(), 3);
/// assert_eq!(layout.token_for(SlotIndex::new(0))?, TokenIndex::new(1));
/// assert_eq!(layout.slot_for(TokenIndex::new(0))?, SlotIndex::new(2));
/// assert!(layout.token_for(SlotIndex::new(1)).is_err());
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TokenLayout {
    slots: Vec<Option<TokenIndex>>,
    tokens: Vec<SlotIndex>,
}

impl TokenLayout {
    /// Creates a layout over `slot_count` slots that gives a token to each
    /// slot in `order`, numbering the tokens from zero in that order.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::TooManySlots`] with `u32::MAX` if `slot_count`
    /// exceeds `u32::MAX`, [`IrError::UnknownSlot`] if a slot in `order` is
    /// not below `slot_count`, and [`IrError::DuplicateToken`] if a slot
    /// appears in `order` more than once.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{IrError, SlotIndex, TokenLayout};
    ///
    /// let layout = TokenLayout::new(3, [SlotIndex::new(1), SlotIndex::new(2)])?;
    /// assert_eq!(layout.token_count(), 2);
    ///
    /// assert_eq!(
    ///     TokenLayout::new(3, [SlotIndex::new(3)]),
    ///     Err(IrError::UnknownSlot(SlotIndex::new(3)))
    /// );
    /// assert_eq!(
    ///     TokenLayout::new(3, [SlotIndex::new(1), SlotIndex::new(1)]),
    ///     Err(IrError::DuplicateToken(SlotIndex::new(1)))
    /// );
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn new(
        slot_count: usize,
        order: impl IntoIterator<Item = SlotIndex>,
    ) -> Result<TokenLayout, IrError> {
        if u32::try_from(slot_count).is_err() {
            return Err(IrError::TooManySlots(u32::MAX));
        }
        let mut slots = vec![None; slot_count];
        let mut tokens = Vec::new();
        for slot in order {
            let entry = slots
                .get_mut(slot.as_usize())
                .ok_or(IrError::UnknownSlot(slot))?;
            if entry.is_some() {
                return Err(IrError::DuplicateToken(slot));
            }
            *entry = Some(TokenIndex::new(tokens.len() as u32));
            tokens.push(slot);
        }
        Ok(TokenLayout { slots, tokens })
    }

    /// Returns the number of tokens.
    pub fn token_count(&self) -> usize {
        self.tokens.len()
    }

    /// Returns the token of `slot`.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownSlot`] if `slot` is not below the slot count
    /// of the layout, and [`IrError::SlotNotTokenized`] if `slot` has no
    /// token.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{IrError, SlotIndex, TokenIndex, TokenLayout};
    ///
    /// let layout = TokenLayout::new(3, [SlotIndex::new(2)])?;
    ///
    /// assert_eq!(layout.token_for(SlotIndex::new(2))?, TokenIndex::new(0));
    /// assert_eq!(
    ///     layout.token_for(SlotIndex::new(0)),
    ///     Err(IrError::SlotNotTokenized(SlotIndex::new(0)))
    /// );
    /// assert_eq!(
    ///     layout.token_for(SlotIndex::new(3)),
    ///     Err(IrError::UnknownSlot(SlotIndex::new(3)))
    /// );
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn token_for(&self, slot: SlotIndex) -> Result<TokenIndex, IrError> {
        match self.slots.get(slot.as_usize()) {
            None => Err(IrError::UnknownSlot(slot)),
            Some(None) => Err(IrError::SlotNotTokenized(slot)),
            Some(Some(token)) => Ok(*token),
        }
    }

    /// Returns the slot that `token` stands for.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::UnknownToken`] if `token` is not below the token
    /// count.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{IrError, SlotIndex, TokenIndex, TokenLayout};
    ///
    /// let layout = TokenLayout::new(3, [SlotIndex::new(2)])?;
    ///
    /// assert_eq!(layout.slot_for(TokenIndex::new(0))?, SlotIndex::new(2));
    /// assert_eq!(
    ///     layout.slot_for(TokenIndex::new(1)),
    ///     Err(IrError::UnknownToken(TokenIndex::new(1)))
    /// );
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn slot_for(&self, token: TokenIndex) -> Result<SlotIndex, IrError> {
        self.tokens
            .get(token.as_usize())
            .copied()
            .ok_or(IrError::UnknownToken(token))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use mmcif_model::{AuthSeqId, InsCode, LabelSeqId};

    fn slots(values: [u32; 3]) -> Vec<SlotIndex> {
        values.into_iter().map(SlotIndex::new).collect()
    }

    #[test]
    fn token_layout_converts_both_ways() {
        let layout = TokenLayout::new(4, slots([2, 0, 3])).unwrap();
        assert_eq!(layout.token_count(), 3);
        assert_eq!(
            layout.token_for(SlotIndex::new(0)).unwrap(),
            TokenIndex::new(1)
        );
        assert_eq!(
            layout.slot_for(TokenIndex::new(0)).unwrap(),
            SlotIndex::new(2)
        );
        for token in 0..layout.token_count() as u32 {
            let token = TokenIndex::new(token);
            let slot = layout.slot_for(token).unwrap();
            assert_eq!(layout.token_for(slot).unwrap(), token);
        }
    }

    #[test]
    fn token_layout_rejects_an_unrepresentable_slot_count() {
        assert_eq!(
            TokenLayout::new(usize::MAX, []),
            Err(IrError::TooManySlots(u32::MAX))
        );
    }

    #[test]
    fn slot_spans_iterate_and_contain_their_slots() {
        let span = SlotSpan::new(2, 5);
        assert_eq!(span.len(), 3);
        assert!(!span.is_empty());
        assert_eq!(span.start(), SlotIndex::new(2));
        assert_eq!(span.end(), SlotIndex::new(5));
        assert!(span.contains(SlotIndex::new(4)) && !span.contains(SlotIndex::new(5)));
        assert_eq!(span.iter().count(), 3);
        assert!(SlotSpan::new(2, 2).is_empty());
    }

    #[test]
    fn token_layout_rejects_unknown_and_repeated_slots() {
        assert_eq!(
            TokenLayout::new(2, slots([0, 1, 2])),
            Err(IrError::UnknownSlot(SlotIndex::new(2)))
        );
        assert_eq!(
            TokenLayout::new(3, slots([0, 1, 1])),
            Err(IrError::DuplicateToken(SlotIndex::new(1)))
        );
    }

    #[test]
    fn untokenized_and_missing_slots_fail_separately() {
        let layout = TokenLayout::new(3, [SlotIndex::new(0)]).unwrap();
        assert_eq!(
            layout.token_for(SlotIndex::new(1)),
            Err(IrError::SlotNotTokenized(SlotIndex::new(1)))
        );
        assert_eq!(
            layout.token_for(SlotIndex::new(9)),
            Err(IrError::UnknownSlot(SlotIndex::new(9)))
        );
        assert_eq!(
            layout.slot_for(TokenIndex::new(1)),
            Err(IrError::UnknownToken(TokenIndex::new(1)))
        );
    }

    #[test]
    fn residue_refs_do_not_coerce_between_numbering_schemes() {
        let auth = ResidueRef::Auth(AuthResidueKey {
            chain: "H".into(),
            seq: AuthSeqId::new(100),
            ins_code: InsCode::new('A'),
        });
        let label = ResidueRef::Label(LabelResidueKey {
            asym: "A".into(),
            seq: LabelSeqId::new(3).unwrap(),
        });
        assert!(auth.as_auth().is_some() && auth.as_label().is_none());
        assert!(label.as_label().is_some() && label.as_auth().is_none());
        assert_eq!(auth.to_string(), "auth H:100A");
        assert_eq!(label.to_string(), "label A/3");
    }
}

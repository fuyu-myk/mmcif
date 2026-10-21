use std::fmt;

use mmcif_model::{AuthResidueKey, LabelResidueKey};

use crate::error::IrError;

macro_rules! index_id {
    ($name:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
        pub struct $name(u32);

        impl $name {
            pub fn new(value: u32) -> Self {
                $name(value)
            }

            pub fn get(self) -> u32 {
                self.0
            }

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

index_id!(SegmentIndex);
index_id!(ElementIndex);
index_id!(SlotIndex);
index_id!(TokenIndex);

macro_rules! name_id {
    ($name:ident) => {
        #[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
        pub struct $name(Box<str>);

        impl $name {
            pub fn new(value: impl Into<Box<str>>) -> Self {
                $name(value.into())
            }

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

name_id!(SegmentId);
name_id!(GroupId);

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ResidueRef {
    Auth(AuthResidueKey),
    Label(LabelResidueKey),
}

impl ResidueRef {
    pub fn as_auth(&self) -> Option<&AuthResidueKey> {
        match self {
            ResidueRef::Auth(key) => Some(key),
            ResidueRef::Label(_) => None,
        }
    }

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SlotSpan {
    start: u32,
    end: u32,
}

impl SlotSpan {
    pub(crate) fn new(start: u32, end: u32) -> SlotSpan {
        SlotSpan { start, end }
    }

    pub fn start(self) -> SlotIndex {
        SlotIndex::new(self.start)
    }

    pub fn end(self) -> SlotIndex {
        SlotIndex::new(self.end)
    }

    pub fn len(self) -> usize {
        (self.end - self.start) as usize
    }

    pub fn is_empty(self) -> bool {
        self.start == self.end
    }

    pub fn contains(self, slot: SlotIndex) -> bool {
        self.start <= slot.get() && slot.get() < self.end
    }

    pub fn iter(self) -> impl ExactSizeIterator<Item = SlotIndex> {
        (self.start..self.end).map(SlotIndex::new)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TokenLayout {
    slots: Vec<Option<TokenIndex>>,
    tokens: Vec<SlotIndex>,
}

impl TokenLayout {
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

    pub fn token_count(&self) -> usize {
        self.tokens.len()
    }

    pub fn token_for(&self, slot: SlotIndex) -> Result<TokenIndex, IrError> {
        match self.slots.get(slot.as_usize()) {
            None => Err(IrError::UnknownSlot(slot)),
            Some(None) => Err(IrError::SlotNotTokenized(slot)),
            Some(Some(token)) => Ok(*token),
        }
    }

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

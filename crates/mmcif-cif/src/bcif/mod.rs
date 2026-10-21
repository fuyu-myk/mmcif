mod array;
mod blocks;
mod encoding;
mod error;
mod msgpack;
mod numeric;
mod record;
mod strings;
#[cfg(test)]
mod test_files;
#[cfg(test)]
mod test_pack;
mod typed;

use std::path::Path;

use crate::document::Document;
use crate::error::Error;

pub use error::BcifError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BcifLimits {
    pub max_elements: usize,
}

impl Default for BcifLimits {
    fn default() -> Self {
        BcifLimits {
            max_elements: 1 << 28,
        }
    }
}

pub(crate) struct Budget {
    remaining: usize,
    limit: usize,
}

impl Budget {
    pub fn new(limits: BcifLimits) -> Self {
        Budget {
            remaining: limits.max_elements,
            limit: limits.max_elements,
        }
    }

    pub fn take(&mut self, count: usize) -> Result<(), BcifError> {
        self.remaining = self
            .remaining
            .checked_sub(count)
            .ok_or(BcifError::LimitExceeded { limit: self.limit })?;
        Ok(())
    }
}

impl Document {
    pub fn from_bcif(bytes: &[u8]) -> Result<Document, BcifError> {
        Document::from_bcif_with_limits(bytes, BcifLimits::default())
    }

    pub fn from_bcif_with_limits(bytes: &[u8], limits: BcifLimits) -> Result<Document, BcifError> {
        blocks::decode_document(bytes, limits)
    }

    pub fn open_bcif(path: impl AsRef<Path>) -> Result<Document, Error> {
        let bytes = std::fs::read(path)?;
        Ok(Document::from_bcif(&bytes)?)
    }
}

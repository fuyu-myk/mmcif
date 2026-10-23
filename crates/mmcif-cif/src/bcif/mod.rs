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

/// Limits on how much data BinaryCIF decoding may create.
///
/// The default allows `2^28` elements and `2^31` bytes of text.
///
/// # Examples
///
/// ```
/// use mmcif_cif::{BcifLimits, Document};
///
/// # let bytes = b"\x81\xaadataBlocks\x91\x82\xa6header\xa41ABC\xaacategories\x91\x83\
/// #     \xa4name\xa5_atom\xa8rowCount\x02\xa7columns\x91\x82\xa4name\xa2id\xa4data\
/// #     \x82\xa4data\xc4\x02\x01\x02\xa8encoding\x91\x82\xa4kind\xa9ByteArray\xa4type\x01";
/// let limits = BcifLimits {
///     max_text_bytes: 1 << 20,
///     ..BcifLimits::default()
/// };
/// let document = Document::from_bcif_with_limits(bytes, limits)?;
///
/// assert_eq!(document.len(), 1);
/// # Ok::<(), mmcif_cif::BcifError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BcifLimits {
    /// The maximum number of elements decoding may create, counting the
    /// entries of MessagePack arrays and maps and the elements of decoded
    /// column data.
    pub max_elements: usize,
    /// The maximum number of bytes of text the decoded document may store for
    /// block names, tags and values held as text.
    pub max_text_bytes: usize,
}

impl Default for BcifLimits {
    fn default() -> Self {
        BcifLimits {
            max_elements: 1 << 28,
            max_text_bytes: 1 << 31,
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
    /// Decodes a BinaryCIF document from `bytes`, within the default
    /// [`BcifLimits`].
    ///
    /// Each data block becomes a block named by its header. A category with
    /// more than one row becomes a loop, a category with one row becomes
    /// pairs, and a category with no rows or no columns is left out. Masked
    /// values become `.` or `?`, and strings take the kind chosen by
    /// [`ValueKind::for_text`](crate::ValueKind::for_text). Floating-point
    /// numbers are written in decimal notation with the fewest digits that
    /// read back to the same value, and fixed-point numbers with a factor of
    /// `10^n` are written with `n` decimal places.
    ///
    /// # Errors
    ///
    /// Returns:
    ///
    /// - [`Truncated`], [`UnsupportedMarker`], [`TooDeep`] or
    ///   [`TrailingBytes`] if `bytes` is not one well-formed MessagePack value,
    ///   and [`InvalidUtf8`] if a string in it is not UTF-8.
    /// - [`MissingField`] or [`UnexpectedType`] if a field the format requires
    ///   is missing or does not have the type it requires.
    /// - [`UnknownEncoding`], [`UnsupportedDataType`], [`UnexpectedArray`],
    ///   [`InvalidParameter`], [`MisalignedBytes`], [`Overflow`],
    ///   [`OddRunLength`], [`UnterminatedPacking`], [`InvalidStringOffsets`]
    ///   or [`InvalidStringIndex`] if the encoded data of a column cannot be
    ///   decoded.
    /// - [`LengthMismatch`] if decoded data has a different number of elements
    ///   than its encoding or category declares, and [`InvalidMask`] if a mask
    ///   value is not 0, 1 or 2.
    /// - [`LimitExceeded`] or [`TextLimitExceeded`] if decoding goes beyond the
    ///   limits.
    /// - [`UnrepresentableValue`] if a string contains a line break followed by
    ///   `;`.
    /// - [`InvalidBlockName`] if a block header is empty or contains
    ///   whitespace, [`InvalidTag`] if a category or column name makes a tag
    ///   with whitespace, [`DuplicateTag`] if a block has a tag twice in any
    ///   ASCII case, and [`TooLarge`] if the decoded document is larger than a
    ///   [`Document`] can hold.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{BcifError, Document};
    ///
    /// # let bytes = b"\x81\xaadataBlocks\x91\x82\xa6header\xa41ABC\xaacategories\x91\x83\
    /// #     \xa4name\xa5_atom\xa8rowCount\x02\xa7columns\x91\x82\xa4name\xa2id\xa4data\
    /// #     \x82\xa4data\xc4\x02\x01\x02\xa8encoding\x91\x82\xa4kind\xa9ByteArray\xa4type\x01";
    /// let document = Document::from_bcif(bytes)?;
    /// let block = document.block(0).unwrap();
    ///
    /// assert_eq!(block.name(), "1ABC");
    /// assert_eq!(block.category("atom").unwrap().row_count(), 2);
    /// assert_eq!(document.to_cif_string(), "data_1ABC\n#\nloop_\n_atom.id\n1\n2\n#\n");
    ///
    /// assert_eq!(
    ///     Document::from_bcif(&[0x80]).unwrap_err(),
    ///     BcifError::MissingField("dataBlocks")
    /// );
    /// # Ok::<(), BcifError>(())
    /// ```
    ///
    /// [`Truncated`]: BcifError::Truncated
    /// [`UnsupportedMarker`]: BcifError::UnsupportedMarker
    /// [`TooDeep`]: BcifError::TooDeep
    /// [`TrailingBytes`]: BcifError::TrailingBytes
    /// [`InvalidUtf8`]: BcifError::InvalidUtf8
    /// [`MissingField`]: BcifError::MissingField
    /// [`UnexpectedType`]: BcifError::UnexpectedType
    /// [`UnknownEncoding`]: BcifError::UnknownEncoding
    /// [`UnsupportedDataType`]: BcifError::UnsupportedDataType
    /// [`UnexpectedArray`]: BcifError::UnexpectedArray
    /// [`InvalidParameter`]: BcifError::InvalidParameter
    /// [`MisalignedBytes`]: BcifError::MisalignedBytes
    /// [`Overflow`]: BcifError::Overflow
    /// [`OddRunLength`]: BcifError::OddRunLength
    /// [`UnterminatedPacking`]: BcifError::UnterminatedPacking
    /// [`InvalidStringOffsets`]: BcifError::InvalidStringOffsets
    /// [`InvalidStringIndex`]: BcifError::InvalidStringIndex
    /// [`LengthMismatch`]: BcifError::LengthMismatch
    /// [`InvalidMask`]: BcifError::InvalidMask
    /// [`LimitExceeded`]: BcifError::LimitExceeded
    /// [`TextLimitExceeded`]: BcifError::TextLimitExceeded
    /// [`UnrepresentableValue`]: BcifError::UnrepresentableValue
    /// [`InvalidBlockName`]: BcifError::InvalidBlockName
    /// [`InvalidTag`]: BcifError::InvalidTag
    /// [`DuplicateTag`]: BcifError::DuplicateTag
    /// [`TooLarge`]: BcifError::TooLarge
    pub fn from_bcif(bytes: &[u8]) -> Result<Document, BcifError> {
        Document::from_bcif_with_limits(bytes, BcifLimits::default())
    }

    /// Decodes a BinaryCIF document from `bytes`, within `limits`.
    ///
    /// The document is built as by [`Document::from_bcif`].
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Document::from_bcif`], with
    /// [`BcifError::LimitExceeded`] and [`BcifError::TextLimitExceeded`]
    /// checked against `limits`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{BcifError, BcifLimits, Document};
    ///
    /// # let bytes = b"\x81\xaadataBlocks\x91\x82\xa6header\xa41ABC\xaacategories\x91\x83\
    /// #     \xa4name\xa5_atom\xa8rowCount\x02\xa7columns\x91\x82\xa4name\xa2id\xa4data\
    /// #     \x82\xa4data\xc4\x02\x01\x02\xa8encoding\x91\x82\xa4kind\xa9ByteArray\xa4type\x01";
    /// let limits = BcifLimits {
    ///     max_elements: 8,
    ///     ..BcifLimits::default()
    /// };
    ///
    /// assert_eq!(
    ///     Document::from_bcif_with_limits(bytes, limits).unwrap_err(),
    ///     BcifError::LimitExceeded { limit: 8 }
    /// );
    /// ```
    pub fn from_bcif_with_limits(bytes: &[u8], limits: BcifLimits) -> Result<Document, BcifError> {
        blocks::decode_document(bytes, limits)
    }

    /// Reads the BinaryCIF file at `path` and decodes it within the default
    /// [`BcifLimits`].
    ///
    /// The file is read with [`read_file`](crate::read_file), so it may be
    /// gzip-compressed when the `gzip` feature is on.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] if the file cannot be read, [`Error::Compressed`]
    /// if it is compressed and the `gzip` feature is off, and [`Error::Bcif`]
    /// with the errors of [`Document::from_bcif`] if its contents cannot be
    /// decoded.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use mmcif_cif::Document;
    ///
    /// let document = Document::open_bcif("1abc.bcif")?;
    /// let atoms = document.block(0).unwrap().category("atom_site").unwrap();
    /// println!("{} atoms", atoms.row_count());
    /// # Ok::<(), mmcif_cif::Error>(())
    /// ```
    pub fn open_bcif(path: impl AsRef<Path>) -> Result<Document, Error> {
        Ok(Document::from_bcif(&crate::read_file(path)?)?)
    }
}

use std::fs::File;
use std::path::Path;

use memmap2::Mmap;

use crate::error::{Error, ParseError, ParseErrorKind};

pub(crate) struct MappedText(Mmap);

impl MappedText {
    #[allow(unsafe_code)]
    pub unsafe fn open(path: &Path) -> Result<Source, Error> {
        let file = File::open(path)?;
        if !file.metadata()?.is_file() || file.metadata()?.len() == 0 {
            return Ok(Source::Owned(read_text(std::fs::read(path)?)?));
        }
        // SAFETY: the caller guarantees the file is neither modified nor truncated while mapped.
        let map = unsafe { Mmap::map(&file)? };
        if let Err(error) = std::str::from_utf8(&map) {
            return Err(
                ParseError::at(&map, error.valid_up_to(), ParseErrorKind::InvalidUtf8).into(),
            );
        }
        Ok(Source::Mapped(MappedText(map)))
    }

    fn text(&self) -> &str {
        #[allow(unsafe_code)]
        // SAFETY: `open` validated the mapping as UTF-8 and the caller guarantees it does not change.
        unsafe {
            std::str::from_utf8_unchecked(&self.0)
        }
    }
}

pub(crate) enum Source {
    Owned(String),
    Mapped(MappedText),
}

impl Source {
    pub fn text(&self) -> &str {
        match self {
            Source::Owned(text) => text,
            Source::Mapped(mapped) => mapped.text(),
        }
    }
}

pub(crate) fn read_text(bytes: Vec<u8>) -> Result<String, ParseError> {
    String::from_utf8(bytes).map_err(|error| {
        let offset = error.utf8_error().valid_up_to();
        ParseError::at(error.as_bytes(), offset, ParseErrorKind::InvalidUtf8)
    })
}

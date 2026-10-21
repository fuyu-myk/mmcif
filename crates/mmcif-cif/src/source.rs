use std::fs::File;
use std::io;
use std::path::Path;

use memmap2::Mmap;

pub(crate) enum Source {
    Owned(String),
    Mapped(Mmap),
}

impl Source {
    pub fn map(path: &Path) -> io::Result<Source> {
        let file = File::open(path)?;
        if file.metadata()?.len() == 0 {
            return Ok(Source::Owned(String::new()));
        }
        #[allow(unsafe_code)]
        // SAFETY: the mapping is read-only and documents assume the file is not modified while open.
        let map = unsafe { Mmap::map(&file)? };
        Ok(Source::Mapped(map))
    }

    pub fn bytes(&self) -> &[u8] {
        match self {
            Source::Owned(text) => text.as_bytes(),
            Source::Mapped(map) => map,
        }
    }

    pub fn text(&self) -> &str {
        match self {
            Source::Owned(text) => text,
            #[allow(unsafe_code)]
            // SAFETY: mapped sources are only exposed after `Document::from_source` validated them as UTF-8.
            Source::Mapped(map) => unsafe { std::str::from_utf8_unchecked(map) },
        }
    }
}

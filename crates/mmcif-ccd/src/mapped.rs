use std::fs::File;
use std::io;
use std::path::Path;

use memmap2::Mmap;

pub(crate) enum FileBytes {
    Owned(Vec<u8>),
    Mapped(Mmap),
}

impl FileBytes {
    pub fn map(path: &Path) -> io::Result<FileBytes> {
        let file = File::open(path)?;
        if file.metadata()?.len() == 0 {
            return Ok(FileBytes::Owned(Vec::new()));
        }
        #[allow(unsafe_code)]
        // SAFETY: the mapping is read-only and readers assume the file is not modified while mapped.
        let map = unsafe { Mmap::map(&file)? };
        Ok(FileBytes::Mapped(map))
    }

    pub fn as_slice(&self) -> &[u8] {
        match self {
            FileBytes::Owned(bytes) => bytes,
            FileBytes::Mapped(map) => map,
        }
    }
}

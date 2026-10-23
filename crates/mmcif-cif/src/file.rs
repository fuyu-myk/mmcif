use std::path::Path;

use crate::error::Error;

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Reads the file at `path`, decompressing it if it is gzip-compressed.
///
/// A file that begins with the gzip magic bytes `1f 8b` is compressed; no
/// CIF, BinaryCIF or PDB file begins with them. Decompressing needs the
/// `gzip` feature.
///
/// # Errors
///
/// Returns [`Error::Io`] if the file cannot be read or its compressed data is
/// corrupt, and [`Error::Compressed`] if it is compressed and the `gzip`
/// feature is off.
///
/// # Examples
///
/// ```no_run
/// let bytes = mmcif_cif::read_file("1abc.cif.gz")?;
/// let document = mmcif_cif::Document::from_bytes(bytes)?;
/// println!("{} blocks", document.len());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn read_file(path: impl AsRef<Path>) -> Result<Vec<u8>, Error> {
    let bytes = std::fs::read(path)?;
    if bytes.starts_with(&GZIP_MAGIC) {
        decompress(&bytes)
    } else {
        Ok(bytes)
    }
}

#[cfg(feature = "gzip")]
fn decompress(bytes: &[u8]) -> Result<Vec<u8>, Error> {
    use std::io::Read;

    let mut decompressed = Vec::new();
    flate2::read::MultiGzDecoder::new(bytes).read_to_end(&mut decompressed)?;
    Ok(decompressed)
}

#[cfg(not(feature = "gzip"))]
fn decompress(_bytes: &[u8]) -> Result<Vec<u8>, Error> {
    Err(Error::Compressed)
}

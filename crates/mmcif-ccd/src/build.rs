use std::io::{self, Write};
use std::path::Path;

use mmcif_cif::stream_blocks;

use crate::codec::{put_str, put_varint, varint_len};
use crate::component::Component;
use crate::error::CcdError;
use crate::format::{DIRECTORY_ENTRY_LEN, FORMAT_VERSION, IdKey, MAGIC, id_key, key_id};
use crate::mapped::FileBytes;
use crate::record::encode_record;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BuildSummary {
    components: usize,
    atoms: usize,
    bonds: usize,
    bytes: u64,
}

impl BuildSummary {
    pub fn components(&self) -> usize {
        self.components
    }

    pub fn atoms(&self) -> usize {
        self.atoms
    }

    pub fn bonds(&self) -> usize {
        self.bonds
    }

    pub fn bytes(&self) -> u64 {
        self.bytes
    }
}

struct PendingRecord {
    key: IdKey,
    start: usize,
    end: usize,
}

impl PendingRecord {
    fn len(&self) -> u64 {
        (self.end - self.start) as u64
    }
}

pub struct IndexBuilder {
    release: String,
    pending: Vec<PendingRecord>,
    records: Vec<u8>,
    atoms: usize,
    bonds: usize,
}

impl IndexBuilder {
    pub fn new(release: &str) -> Self {
        IndexBuilder {
            release: release.to_owned(),
            pending: Vec::new(),
            records: Vec::new(),
            atoms: 0,
            bonds: 0,
        }
    }

    pub fn add(&mut self, component: &Component) -> Result<(), CcdError> {
        let key = id_key(component.id())?;
        let start = self.records.len();
        encode_record(component, &mut self.records);
        self.pending.push(PendingRecord {
            key,
            start,
            end: self.records.len(),
        });
        self.atoms += component.atoms().len();
        self.bonds += component.bonds().len();
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.pending.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pending.is_empty()
    }

    pub fn finish(mut self, mut out: impl Write) -> Result<BuildSummary, CcdError> {
        self.pending.sort_unstable_by_key(|record| record.key);
        if let Some(pair) = self
            .pending
            .windows(2)
            .find(|pair| pair[0].key == pair[1].key)
        {
            let id = key_id(&pair[0].key).unwrap_or_default().to_owned();
            return Err(CcdError::DuplicateId { id });
        }
        let header = self.header();
        out.write_all(&header)?;
        let mut prefix = Vec::with_capacity(varint_len(u64::MAX));
        let mut records_len = 0;
        for record in &self.pending {
            prefix.clear();
            put_varint(&mut prefix, record.len());
            out.write_all(&prefix)?;
            out.write_all(&self.records[record.start..record.end])?;
            records_len += (prefix.len() + record.end - record.start) as u64;
        }
        out.flush()?;
        Ok(BuildSummary {
            components: self.pending.len(),
            atoms: self.atoms,
            bonds: self.bonds,
            bytes: header.len() as u64 + records_len,
        })
    }

    fn header(&self) -> Vec<u8> {
        let mut header = Vec::with_capacity(self.pending.len() * DIRECTORY_ENTRY_LEN);
        header.extend_from_slice(MAGIC);
        header.extend_from_slice(&FORMAT_VERSION.to_le_bytes());
        put_str(&mut header, &self.release);
        put_varint(&mut header, self.pending.len() as u64);
        let mut offset = 0u64;
        for record in &self.pending {
            header.extend_from_slice(&record.key);
            header.extend_from_slice(&offset.to_le_bytes());
            offset += varint_len(record.len()) as u64 + record.len();
        }
        header
    }
}

pub fn build_index(
    source_text: &str,
    release: &str,
    out: impl Write,
) -> Result<BuildSummary, CcdError> {
    let mut builder = IndexBuilder::new(release);
    for streamed in stream_blocks(source_text) {
        builder.add(&Component::from_block(streamed?.block())?)?;
    }
    builder.finish(out)
}

pub fn build_index_from_file(
    source: impl AsRef<Path>,
    release: &str,
    out: impl Write,
) -> Result<BuildSummary, CcdError> {
    let bytes = FileBytes::map(source.as_ref())?;
    let text = std::str::from_utf8(bytes.as_slice())
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    build_index(text, release, out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format::MAX_ID_LEN;

    fn block(id: &str, atoms: &[&str]) -> String {
        let mut text = format!("data_{id}\n_chem_comp.id {id}\n");
        if !atoms.is_empty() {
            text.push_str("loop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n");
            for atom in atoms {
                text.push_str(&format!("{atom} C\n"));
            }
        }
        text
    }

    fn build(text: &str) -> Result<(BuildSummary, Vec<u8>), CcdError> {
        let mut out = Vec::new();
        let summary = build_index(text, "2026-10", &mut out)?;
        Ok((summary, out))
    }

    #[test]
    fn writes_header_and_sorted_directory() {
        let text = [
            block("ZN", &["ZN"]),
            block("ALA", &["N", "CA"]),
            block("A1AAA", &[]),
        ]
        .concat();
        let (summary, bytes) = build(&text).unwrap();
        assert_eq!(
            (summary.components(), summary.atoms(), summary.bonds()),
            (3, 3, 0)
        );
        assert_eq!(summary.bytes(), bytes.len() as u64);
        assert_eq!(&bytes[..8], MAGIC);
        assert_eq!(&bytes[8..12], &FORMAT_VERSION.to_le_bytes());
        assert_eq!(&bytes[12..20], b"\x072026-10");
        assert_eq!(bytes[20], 3);
        let directory: Vec<&str> = bytes[21..21 + 3 * DIRECTORY_ENTRY_LEN]
            .chunks(DIRECTORY_ENTRY_LEN)
            .map(|entry| key_id(&entry[..MAX_ID_LEN]).unwrap())
            .collect();
        assert_eq!(directory, ["A1AAA", "ALA", "ZN"]);
    }

    #[test]
    fn empty_sources_build_empty_indexes() {
        let (summary, bytes) = build("").unwrap();
        assert_eq!(summary.components(), 0);
        assert_eq!(bytes.len(), 8 + 4 + 8 + 1);
    }

    #[test]
    fn rejects_duplicate_ids() {
        let text = [block("ALA", &[]), block("GLY", &[]), block("ALA", &["N"])].concat();
        assert!(matches!(build(&text), Err(CcdError::DuplicateId { ref id }) if id == "ALA"));
    }

    #[test]
    fn rejects_ids_that_do_not_fit() {
        let long = block("ABCDEFGHI", &[]);
        assert!(matches!(build(&long), Err(CcdError::IdTooLong { ref id }) if id == "ABCDEFGHI"));
        let spaced = "data_X\n_chem_comp.id 'A B'\n";
        assert!(matches!(build(spaced), Err(CcdError::InvalidId { .. })));
        let eight = block("ABCDEFGH", &[]);
        assert_eq!(build(&eight).unwrap().0.components(), 1);
    }

    #[test]
    fn reports_parse_and_component_errors() {
        assert!(matches!(
            build("data_X\n_chem_comp.id 'open\n"),
            Err(CcdError::Parse(_))
        ));
        assert!(matches!(
            build("data_X\n_x.y 1\n"),
            Err(CcdError::MissingItem { .. })
        ));
    }

    #[test]
    fn reports_missing_source_files() {
        let missing = Path::new(env!("CARGO_MANIFEST_DIR")).join("no-such-file.cif");
        let result = build_index_from_file(missing, "r", Vec::new());
        assert!(matches!(result, Err(CcdError::Io(_))));
    }
}

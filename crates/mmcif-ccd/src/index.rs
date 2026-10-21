use std::cmp::Ordering;
use std::fmt;
use std::path::Path;

use crate::codec::Reader;
use crate::component::Component;
use crate::error::{CcdError, corrupt};
use crate::format::{DIRECTORY_ENTRY_LEN, FORMAT_VERSION, MAGIC, MAX_ID_LEN, id_key, key_id};
use crate::record::decode_record;

const OUT_OF_BOUNDS: &str = "record offset is out of bounds";

struct Layout {
    release: String,
    count: usize,
    directory_start: usize,
    records_start: usize,
}

impl Layout {
    fn read(data: &[u8]) -> Result<Layout, CcdError> {
        let mut reader = Reader::new(data);
        if reader.take(MAGIC.len()).ok() != Some(MAGIC.as_slice()) {
            return Err(CcdError::BadMagic);
        }
        let found = reader.u32()?;
        if found != FORMAT_VERSION {
            return Err(CcdError::UnsupportedVersion { found });
        }
        let release = reader.str()?.to_owned();
        let count = reader.count(DIRECTORY_ENTRY_LEN)?;
        let directory_start = data.len() - reader.remaining();
        validate_directory(reader.take(count * DIRECTORY_ENTRY_LEN)?)?;
        Ok(Layout {
            release,
            count,
            directory_start,
            records_start: data.len() - reader.remaining(),
        })
    }
}

fn validate_directory(directory: &[u8]) -> Result<(), CcdError> {
    let mut previous: Option<&[u8]> = None;
    let (entries, _) = directory.as_chunks::<DIRECTORY_ENTRY_LEN>();
    for entry in entries {
        let key = &entry[..MAX_ID_LEN];
        if key_id(key).is_none() {
            return Err(corrupt("invalid id in directory"));
        }
        if previous.is_some_and(|previous| previous >= key) {
            return Err(corrupt("directory ids are not sorted and unique"));
        }
        previous = Some(key);
    }
    Ok(())
}

pub struct CcdIndex {
    bytes: Vec<u8>,
    layout: Layout,
}

impl CcdIndex {
    pub fn open(path: impl AsRef<Path>) -> Result<CcdIndex, CcdError> {
        CcdIndex::from_bytes(std::fs::read(path)?)
    }

    pub fn from_bytes(bytes: Vec<u8>) -> Result<CcdIndex, CcdError> {
        let layout = Layout::read(&bytes)?;
        Ok(CcdIndex { bytes, layout })
    }

    pub fn release(&self) -> &str {
        &self.layout.release
    }

    pub fn len(&self) -> usize {
        self.layout.count
    }

    pub fn is_empty(&self) -> bool {
        self.layout.count == 0
    }

    pub fn ids(&self) -> impl Iterator<Item = &str> + '_ {
        (0..self.len()).filter_map(|index| self.entry(index).get(..MAX_ID_LEN).and_then(key_id))
    }

    pub fn get(&self, id: &str) -> Result<Option<Component>, CcdError> {
        let Ok(key) = id_key(id) else {
            return Ok(None);
        };
        let Some(index) = self.find(&key) else {
            return Ok(None);
        };
        let mut entry = Reader::new(self.entry(index));
        entry.take(MAX_ID_LEN)?;
        let offset = usize::try_from(entry.u64()?).map_err(|_| corrupt(OUT_OF_BOUNDS))?;
        let records = self
            .data()
            .get(self.layout.records_start..)
            .unwrap_or_default();
        let mut record = Reader::new(records.get(offset..).ok_or(corrupt(OUT_OF_BOUNDS))?);
        let len = record.length()?;
        let component = decode_record(record.take(len)?)?;
        if component.id() != id {
            return Err(corrupt("record id does not match directory"));
        }
        Ok(Some(component))
    }

    fn data(&self) -> &[u8] {
        &self.bytes
    }

    fn entry(&self, index: usize) -> &[u8] {
        let start = self.layout.directory_start + index * DIRECTORY_ENTRY_LEN;
        self.data()
            .get(start..start + DIRECTORY_ENTRY_LEN)
            .unwrap_or_default()
    }

    fn find(&self, key: &[u8]) -> Option<usize> {
        let (mut low, mut high) = (0, self.len());
        while low < high {
            let middle = low + (high - low) / 2;
            match self.entry(middle).get(..MAX_ID_LEN)?.cmp(key) {
                Ordering::Less => low = middle + 1,
                Ordering::Greater => high = middle,
                Ordering::Equal => return Some(middle),
            }
        }
        None
    }
}

impl fmt::Debug for CcdIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CcdIndex")
            .field("release", &self.release())
            .field("len", &self.len())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::build::IndexBuilder;
    use crate::record::tests::sample_component;

    const IDS: [&str; 4] = ["NAG", "ALA", "A1AAA", "ZN"];

    fn index_bytes() -> Vec<u8> {
        let mut builder = IndexBuilder::new("2026-10-20");
        for id in IDS {
            builder.add(&sample_component(id)).unwrap();
        }
        let mut out = Vec::new();
        builder.finish(&mut out).unwrap();
        out
    }

    fn is_index_error(error: &CcdError) -> bool {
        matches!(
            error,
            CcdError::BadMagic
                | CcdError::UnsupportedVersion { .. }
                | CcdError::CorruptIndex { .. }
        )
    }

    #[test]
    fn looks_components_up_by_id() {
        let index = CcdIndex::from_bytes(index_bytes()).unwrap();
        assert_eq!(index.release(), "2026-10-20");
        assert_eq!(index.len(), 4);
        assert_eq!(
            index.ids().collect::<Vec<_>>(),
            ["A1AAA", "ALA", "NAG", "ZN"]
        );
        for id in IDS {
            assert_eq!(index.get(id).unwrap(), Some(sample_component(id)));
        }
        for missing in ["", "AL", "ALAA", "A", "ZZZ", "TOO-LONG-ID", "ala"] {
            assert_eq!(index.get(missing).unwrap(), None, "{missing}");
        }
    }

    #[test]
    fn rejects_bad_magic() {
        let mut bytes = index_bytes();
        bytes[0] = b'X';
        assert!(matches!(
            CcdIndex::from_bytes(bytes),
            Err(CcdError::BadMagic)
        ));
        assert!(matches!(
            CcdIndex::from_bytes(Vec::new()),
            Err(CcdError::BadMagic)
        ));
    }

    #[test]
    fn rejects_other_format_versions() {
        let mut bytes = index_bytes();
        bytes[8..12].copy_from_slice(&2u32.to_le_bytes());
        let error = CcdIndex::from_bytes(bytes).unwrap_err();
        assert!(matches!(error, CcdError::UnsupportedVersion { found: 2 }));
    }

    #[test]
    fn rejects_counts_beyond_the_file() {
        let mut bytes = MAGIC.to_vec();
        bytes.extend_from_slice(&FORMAT_VERSION.to_le_bytes());
        bytes.push(0);
        bytes.extend_from_slice(&[0xFF; 9]);
        bytes.push(0x01);
        let error = CcdIndex::from_bytes(bytes).unwrap_err();
        assert!(matches!(error, CcdError::CorruptIndex { .. }));
    }

    #[test]
    fn rejects_unsorted_directories() {
        let mut bytes = index_bytes();
        let directory = 8 + 4 + 11 + 1;
        let (first, second) = bytes[directory..].split_at_mut(DIRECTORY_ENTRY_LEN);
        first.swap_with_slice(&mut second[..DIRECTORY_ENTRY_LEN]);
        let error = CcdIndex::from_bytes(bytes).unwrap_err();
        assert!(matches!(error, CcdError::CorruptIndex { .. }));
    }

    #[test]
    fn detects_records_that_do_not_match_their_entry() {
        let mut bytes = index_bytes();
        let directory = 8 + 4 + 11 + 1;
        let offsets = directory + MAX_ID_LEN;
        let second = offsets + DIRECTORY_ENTRY_LEN;
        let first_offset: [u8; 8] = bytes[offsets..offsets + 8].try_into().unwrap();
        bytes.copy_within(second..second + 8, offsets);
        bytes[second..second + 8].copy_from_slice(&first_offset);
        let index = CcdIndex::from_bytes(bytes).unwrap();
        let error = index.get("A1AAA").unwrap_err();
        assert!(matches!(error, CcdError::CorruptIndex { .. }));
    }

    #[test]
    fn truncated_indexes_fail_cleanly() {
        let bytes = index_bytes();
        for len in 0..bytes.len() {
            let Ok(index) = CcdIndex::from_bytes(bytes[..len].to_vec()) else {
                continue;
            };
            assert_eq!(index.len(), IDS.len(), "length {len}");
            let mut failures = 0;
            for id in IDS {
                match index.get(id) {
                    Ok(found) => assert_eq!(found, Some(sample_component(id)), "length {len}"),
                    Err(error) => {
                        assert!(is_index_error(&error), "length {len}: {error}");
                        failures += 1;
                    }
                }
            }
            assert!(failures > 0, "length {len} decoded every record");
        }
    }

    #[test]
    fn corrupted_bytes_never_panic() {
        let bytes = index_bytes();
        for offset in 0..bytes.len() {
            for mask in [0x01, 0x80, 0xFF] {
                let mut corrupted = bytes.clone();
                corrupted[offset] ^= mask;
                let index = match CcdIndex::from_bytes(corrupted) {
                    Ok(index) => index,
                    Err(error) => {
                        assert!(is_index_error(&error), "offset {offset}: {error}");
                        continue;
                    }
                };
                let ids: Vec<String> = index.ids().map(str::to_owned).collect();
                for id in ids.iter().map(String::as_str).chain(IDS) {
                    if let Err(error) = index.get(id) {
                        assert!(is_index_error(&error), "offset {offset}: {error}");
                    }
                }
            }
        }
    }
}

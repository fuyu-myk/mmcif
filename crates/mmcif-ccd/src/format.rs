use crate::error::CcdError;

pub(crate) const MAGIC: &[u8; 8] = b"MMCCDIDX";
pub(crate) const FORMAT_VERSION: u32 = 1;
pub(crate) const MAX_ID_LEN: usize = 8;
pub(crate) const DIRECTORY_ENTRY_LEN: usize = MAX_ID_LEN + 8;

pub(crate) type IdKey = [u8; MAX_ID_LEN];

pub(crate) fn id_key(id: &str) -> Result<IdKey, CcdError> {
    if id.len() > MAX_ID_LEN {
        return Err(CcdError::IdTooLong { id: id.to_owned() });
    }
    if id.is_empty() || !id.bytes().all(|byte| byte.is_ascii_graphic()) {
        return Err(CcdError::InvalidId { id: id.to_owned() });
    }
    let mut key = [0; MAX_ID_LEN];
    key[..id.len()].copy_from_slice(id.as_bytes());
    Ok(key)
}

pub(crate) fn key_id(key: &[u8]) -> Option<&str> {
    let len = key.iter().position(|&byte| byte == 0).unwrap_or(key.len());
    let (id, padding) = key.split_at(len);
    if id.is_empty() || !id.iter().all(u8::is_ascii_graphic) || padding.iter().any(|&b| b != 0) {
        return None;
    }
    std::str::from_utf8(id).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keys_pad_ids_with_zeros() {
        assert_eq!(&id_key("ALA").unwrap(), b"ALA\0\0\0\0\0");
        assert_eq!(&id_key("ABCDEFGH").unwrap(), b"ABCDEFGH");
        assert_eq!(key_id(b"ALA\0\0\0\0\0"), Some("ALA"));
        assert_eq!(key_id(b"ABCDEFGH"), Some("ABCDEFGH"));
    }

    #[test]
    fn rejects_unrepresentable_ids() {
        assert!(matches!(
            id_key("ABCDEFGHI"),
            Err(CcdError::IdTooLong { .. })
        ));
        assert!(matches!(id_key(""), Err(CcdError::InvalidId { .. })));
        assert!(matches!(id_key("A B"), Err(CcdError::InvalidId { .. })));
        assert!(matches!(id_key("A\0"), Err(CcdError::InvalidId { .. })));
    }

    #[test]
    fn rejects_malformed_keys() {
        assert_eq!(key_id(&[0; 8]), None);
        assert_eq!(key_id(b"A\0B\0\0\0\0\0"), None);
        assert_eq!(key_id(b"A \0\0\0\0\0\0"), None);
        assert_eq!(key_id(&[0xC3, 0xA9, 0, 0, 0, 0, 0, 0]), None);
    }

    #[test]
    fn key_order_matches_id_order() {
        let mut ids = ["ZN", "A1AAA", "ALA", "AL", "HOH", "0G"];
        let mut keys: Vec<IdKey> = ids.iter().map(|id| id_key(id).unwrap()).collect();
        ids.sort_unstable();
        keys.sort_unstable();
        let decoded: Vec<&str> = keys.iter().filter_map(|key| key_id(key)).collect();
        assert_eq!(decoded, ids);
    }
}

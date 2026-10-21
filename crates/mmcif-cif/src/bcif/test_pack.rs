use super::encoding;
use super::msgpack::{self, MsgValue};
use super::{BcifError, BcifLimits, Budget};

#[derive(Debug, Clone)]
pub(crate) enum Pack {
    Nil,
    Bool(bool),
    Int(i64),
    UInt(u64),
    Float(f64),
    Str(String),
    Bin(Vec<u8>),
    Array(Vec<Pack>),
    Map(Vec<(String, Pack)>),
}

pub(crate) fn obj<'k>(entries: impl IntoIterator<Item = (&'k str, Pack)>) -> Pack {
    Pack::Map(
        entries
            .into_iter()
            .map(|(key, value)| (key.to_owned(), value))
            .collect(),
    )
}

pub(crate) fn text(value: &str) -> Pack {
    Pack::Str(value.to_owned())
}

struct Markers {
    fix: Option<(u8, usize)>,
    len8: Option<u8>,
    len16: u8,
    len32: u8,
}

const STR: Markers = Markers {
    fix: Some((0xa0, 32)),
    len8: Some(0xd9),
    len16: 0xda,
    len32: 0xdb,
};
const BIN: Markers = Markers {
    fix: None,
    len8: Some(0xc4),
    len16: 0xc5,
    len32: 0xc6,
};
const ARRAY: Markers = Markers {
    fix: Some((0x90, 16)),
    len8: None,
    len16: 0xdc,
    len32: 0xdd,
};
const MAP: Markers = Markers {
    fix: Some((0x80, 16)),
    len8: None,
    len16: 0xde,
    len32: 0xdf,
};

fn header(out: &mut Vec<u8>, len: usize, markers: Markers) {
    match markers {
        Markers {
            fix: Some((base, limit)),
            ..
        } if len < limit => out.push(base | len as u8),
        Markers {
            len8: Some(marker), ..
        } if len <= 0xff => out.extend([marker, len as u8]),
        Markers { len16, .. } if len <= 0xffff => {
            out.push(len16);
            out.extend((len as u16).to_be_bytes());
        }
        Markers { len32, .. } => {
            out.push(len32);
            out.extend(u32::try_from(len).unwrap().to_be_bytes());
        }
    }
}

impl Pack {
    pub fn encode(&self) -> Vec<u8> {
        let mut out = Vec::new();
        self.write(&mut out);
        out
    }

    fn write(&self, out: &mut Vec<u8>) {
        match self {
            Pack::Nil => out.push(0xc0),
            Pack::Bool(value) => out.push(if *value { 0xc3 } else { 0xc2 }),
            Pack::Int(value) => write_int(out, *value),
            Pack::UInt(value) => {
                out.push(0xcf);
                out.extend(value.to_be_bytes());
            }
            Pack::Float(value) => {
                out.push(0xcb);
                out.extend(value.to_be_bytes());
            }
            Pack::Str(value) => {
                header(out, value.len(), STR);
                out.extend(value.as_bytes());
            }
            Pack::Bin(value) => {
                header(out, value.len(), BIN);
                out.extend(value);
            }
            Pack::Array(items) => {
                header(out, items.len(), ARRAY);
                items.iter().for_each(|item| item.write(out));
            }
            Pack::Map(entries) => {
                header(out, entries.len(), MAP);
                for (key, value) in entries {
                    Pack::Str(key.clone()).write(out);
                    value.write(out);
                }
            }
        }
    }
}

fn write_int(out: &mut Vec<u8>, value: i64) {
    if let Ok(small) = u8::try_from(value)
        && small < 0x80
    {
        out.push(small);
    } else if (-32..0).contains(&value) {
        out.extend(i8::try_from(value).unwrap().to_be_bytes());
    } else if let Ok(byte) = i8::try_from(value) {
        out.push(0xd0);
        out.extend(byte.to_be_bytes());
    } else if let Ok(short) = i16::try_from(value) {
        out.push(0xd1);
        out.extend(short.to_be_bytes());
    } else if let Ok(word) = i32::try_from(value) {
        out.push(0xd2);
        out.extend(word.to_be_bytes());
    } else {
        out.push(0xd3);
        out.extend(value.to_be_bytes());
    }
}

pub(crate) const INT8: i64 = 1;
pub(crate) const INT16: i64 = 2;
pub(crate) const INT32: i64 = 3;
pub(crate) const UINT8: i64 = 4;
pub(crate) const FLOAT32: i64 = 32;
pub(crate) const FLOAT64: i64 = 33;

pub(crate) fn byte_array(code: i64) -> Pack {
    obj([("kind", text("ByteArray")), ("type", Pack::Int(code))])
}

pub(crate) fn fixed_point(factor: Pack, src_type: i64) -> Pack {
    obj([
        ("kind", text("FixedPoint")),
        ("factor", factor),
        ("srcType", Pack::Int(src_type)),
    ])
}

pub(crate) fn run_length(src_size: i64) -> Pack {
    obj([
        ("kind", text("RunLength")),
        ("srcType", Pack::Int(INT32)),
        ("srcSize", Pack::Int(src_size)),
    ])
}

pub(crate) fn delta(origin: i64) -> Pack {
    obj([
        ("kind", text("Delta")),
        ("origin", Pack::Int(origin)),
        ("srcType", Pack::Int(INT32)),
    ])
}

pub(crate) fn integer_packing(byte_count: i64, unsigned: bool, src_size: i64) -> Pack {
    obj([
        ("kind", text("IntegerPacking")),
        ("byteCount", Pack::Int(byte_count)),
        ("isUnsigned", Pack::Bool(unsigned)),
        ("srcSize", Pack::Int(src_size)),
    ])
}

pub(crate) fn string_array(strings: &str, offsets: &[i32], data_encoding: Vec<Pack>) -> Pack {
    obj([
        ("kind", text("StringArray")),
        ("dataEncoding", Pack::Array(data_encoding)),
        ("stringData", text(strings)),
        ("offsetEncoding", Pack::Array(vec![byte_array(INT32)])),
        ("offsets", Pack::Bin(i32_bytes(offsets))),
    ])
}

pub(crate) fn i8_bytes(values: &[i8]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}

pub(crate) fn i16_bytes(values: &[i16]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}

pub(crate) fn i32_bytes(values: &[i32]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}

pub(crate) fn decode_texts_with(
    data: &[u8],
    encodings: Vec<Pack>,
    limits: BcifLimits,
) -> Result<Vec<String>, BcifError> {
    let bytes = Pack::Array(encodings).encode();
    let mut budget = Budget::new(limits);
    let parsed = msgpack::parse(&bytes, &mut budget)?;
    let MsgValue::Array(list) = &parsed else {
        unreachable!("the test encoder wrote an array");
    };
    let array = encoding::decode(data, list, &mut budget)?;
    let mut buffer = String::new();
    Ok((0..array.len())
        .map(|row| array.text(row, &mut buffer).unwrap().to_owned())
        .collect())
}

pub(crate) fn decode_texts(data: &[u8], encodings: Vec<Pack>) -> Result<Vec<String>, BcifError> {
    decode_texts_with(data, encodings, BcifLimits::default())
}

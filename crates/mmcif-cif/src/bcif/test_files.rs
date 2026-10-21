use super::test_pack::{
    FLOAT32, FLOAT64, INT8, INT32, Pack, UINT8, byte_array, delta, fixed_point, i8_bytes,
    i32_bytes, obj, run_length, string_array, text,
};

pub(crate) fn encoded(data: Vec<u8>, encodings: Vec<Pack>) -> Pack {
    obj([
        ("data", Pack::Bin(data)),
        ("encoding", Pack::Array(encodings)),
    ])
}

pub(crate) fn column(name: &str, data: Pack, mask: Option<&[u8]>) -> Pack {
    let mask = mask.map_or(Pack::Nil, |mask| {
        encoded(mask.to_vec(), vec![byte_array(UINT8)])
    });
    obj([("name", text(name)), ("data", data), ("mask", mask)])
}

pub(crate) fn strings(values: &[&str]) -> Pack {
    let mut unique: Vec<&str> = Vec::new();
    let indices: Vec<i32> = values
        .iter()
        .map(|value| match unique.iter().position(|seen| seen == value) {
            Some(index) => index as i32,
            None => {
                unique.push(value);
                unique.len() as i32 - 1
            }
        })
        .collect();
    let offsets: Vec<i32> = std::iter::once(0)
        .chain(unique.iter().scan(0, |end, value| {
            *end += value.len() as i32;
            Some(*end)
        }))
        .collect();
    let encoding = string_array(&unique.concat(), &offsets, vec![byte_array(INT32)]);
    encoded(i32_bytes(&indices), vec![encoding])
}

pub(crate) fn category(name: &str, rows: i64, columns: Vec<Pack>) -> Pack {
    obj([
        ("name", text(name)),
        ("rowCount", Pack::Int(rows)),
        ("columns", Pack::Array(columns)),
    ])
}

pub(crate) fn file(header: &str, categories: Vec<Pack>) -> Vec<u8> {
    let block = obj([
        ("header", text(header)),
        ("categories", Pack::Array(categories)),
    ]);
    obj([
        ("version", text("0.3.0")),
        ("encoder", text("test")),
        ("dataBlocks", Pack::Array(vec![block])),
    ])
    .encode()
}

fn floats32(values: &[f32]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}

pub(crate) fn sample() -> Vec<u8> {
    let entry = category("_entry", 1, vec![column("id", strings(&["1ABC"]), None)]);
    let atoms = vec![
        column(
            "id",
            encoded(
                i32_bytes(&[0, 1, 1, 3]),
                vec![delta(1), run_length(4), byte_array(INT32)],
            ),
            None,
        ),
        column(
            "label_atom_id",
            strings(&["N", "CA", "O5'", "two words"]),
            None,
        ),
        column(
            "label_alt_id",
            strings(&["", "", "A", ""]),
            Some(&[1, 1, 0, 2]),
        ),
        column(
            "Cartn_x",
            encoded(
                i32_bytes(&[12345, -500, 0, 7]),
                vec![fixed_point(Pack::Int(1000), FLOAT32), byte_array(INT32)],
            ),
            None,
        ),
        column(
            "occupancy",
            encoded(floats32(&[1.0, 0.5, 0.25, 0.1]), vec![byte_array(FLOAT32)]),
            None,
        ),
        column(
            "B_iso_or_equiv",
            encoded(
                [10.5f64, 0.1, 3.0, 1e-7]
                    .iter()
                    .flat_map(|v| v.to_le_bytes())
                    .collect(),
                vec![byte_array(FLOAT64)],
            ),
            None,
        ),
        column(
            "pdbx_formal_charge",
            encoded(i8_bytes(&[0, 0, -1, 1]), vec![byte_array(INT8)]),
            Some(&[0, 2, 0, 0]),
        ),
    ];
    file(
        "TEST",
        vec![
            entry,
            category("_empty", 0, Vec::new()),
            category("_atom_site", 4, atoms),
        ],
    )
}

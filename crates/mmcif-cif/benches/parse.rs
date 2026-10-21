use std::hint::black_box;
use std::path::PathBuf;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use mmcif_cif::Document;

const SIZED_IDS: [&str; 3] = ["1crn", "4hhb", "1aon"];
const WRITE_ID: &str = "4hhb";

fn corpus_dir() -> PathBuf {
    let dir = std::env::var_os("MMCIF_CORPUS_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/../../corpus/data"))
        });
    assert!(
        dir.is_dir(),
        "corpus directory {} is missing: run `harness/.venv/bin/python corpus/fetch.py` or set MMCIF_CORPUS_DIR",
        dir.display()
    );
    dir
}

fn corpus_bytes(id: &str, extension: &str) -> Vec<u8> {
    let path = corpus_dir().join(format!("{id}.{extension}"));
    std::fs::read(&path).unwrap_or_else(|error| panic!("corpus entry {}: {error}", path.display()))
}

fn corpus_text(id: &str) -> String {
    String::from_utf8(corpus_bytes(id, "cif"))
        .unwrap_or_else(|error| panic!("corpus entry {id}.cif: {error}"))
}

fn parse_text(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse");
    for id in SIZED_IDS {
        let text = corpus_text(id);
        group.throughput(Throughput::Bytes(text.len() as u64));
        group.bench_with_input(BenchmarkId::new("text", id), text.as_str(), |b, text| {
            b.iter(|| Document::parse(black_box(text)).expect("corpus entry parses"))
        });
    }
    group.finish();
}

fn decode_bcif(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse");
    for id in SIZED_IDS {
        let bytes = corpus_bytes(id, "bcif");
        group.throughput(Throughput::Bytes(bytes.len() as u64));
        group.bench_with_input(
            BenchmarkId::new("bcif", id),
            bytes.as_slice(),
            |b, bytes| {
                b.iter(|| Document::from_bcif(black_box(bytes)).expect("corpus entry decodes"))
            },
        );
    }
    group.finish();
}

fn write_text(c: &mut Criterion) {
    let document = Document::parse(corpus_text(WRITE_ID)).expect("corpus entry parses");
    let written_len = document.to_cif_string().len();
    let mut group = c.benchmark_group("write");
    group.throughput(Throughput::Bytes(written_len as u64));
    group.bench_function(BenchmarkId::new("text", WRITE_ID), |b| {
        b.iter(|| black_box(&document).to_cif_string())
    });
    group.finish();
}

criterion_group!(benches, parse_text, decode_bcif, write_text);
criterion_main!(benches);

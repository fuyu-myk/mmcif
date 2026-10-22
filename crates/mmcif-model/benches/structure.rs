use std::hint::black_box;
use std::path::{Path, PathBuf};

use criterion::{
    BenchmarkId, Criterion, SamplingMode, Throughput, criterion_group, criterion_main,
};
use mmcif_model::{AltlocPolicy, Structure, read_bcif, read_mmcif};

const SIZED_IDS: [&str; 3] = ["1crn", "4hhb", "1aon"];
const ALTLOC_ID: &str = "2vb1";
const ASSEMBLY_IDS: [&str; 2] = ["1stp", "2buk"];
const ASSEMBLY: &str = "1";
const SWEEP_SAMPLES: usize = 10;

type Reader = fn(&Path) -> Structure;

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

fn corpus_path(id: &str) -> PathBuf {
    let path = corpus_dir().join(format!("{id}.cif"));
    assert!(path.is_file(), "corpus entry {} is missing", path.display());
    path
}

fn corpus_files(extension: &str) -> Vec<PathBuf> {
    let dir = corpus_dir();
    let mut files: Vec<PathBuf> = std::fs::read_dir(&dir)
        .unwrap_or_else(|error| panic!("corpus directory {}: {error}", dir.display()))
        .map(|entry| entry.expect("corpus directory entry is readable").path())
        .filter(|path| path.extension().is_some_and(|e| e == extension))
        .collect();
    files.sort();
    assert!(
        !files.is_empty(),
        "no .{extension} files in {}",
        dir.display()
    );
    files
}

fn file_len(path: &Path) -> u64 {
    std::fs::metadata(path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()))
        .len()
}

fn read_text(path: &Path) -> Structure {
    read_mmcif(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

fn read_binary(path: &Path) -> Structure {
    read_bcif(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

fn read_structures(c: &mut Criterion) {
    let mut group = c.benchmark_group("read");
    for id in SIZED_IDS {
        let path = corpus_path(id);
        group.throughput(Throughput::Bytes(file_len(&path)));
        group.bench_with_input(BenchmarkId::new("text", id), path.as_path(), |b, path| {
            b.iter(|| read_text(black_box(path)))
        });
        let binary = path.with_extension("bcif");
        group.throughput(Throughput::Bytes(file_len(&binary)));
        group.bench_with_input(BenchmarkId::new("bcif", id), binary.as_path(), |b, path| {
            b.iter(|| read_binary(black_box(path)))
        });
    }
    group.finish();
}

fn pick_single_conformer(c: &mut Criterion) {
    let structure = read_text(&corpus_path(ALTLOC_ID));
    let model = structure.first_model().expect("altloc entry has a model");
    let mut group = c.benchmark_group("altloc");
    group.bench_function(BenchmarkId::new("highest_occupancy", ALTLOC_ID), |b| {
        b.iter(|| black_box(model).single_conformer(AltlocPolicy::HighestOccupancy))
    });
    group.finish();
}

fn expand_assemblies(c: &mut Criterion) {
    let mut group = c.benchmark_group("assemble");
    for id in ASSEMBLY_IDS {
        let structure = read_text(&corpus_path(id));
        let model = structure.first_model().expect("assembly entry has a model");
        let assembly = structure
            .assemblies()
            .iter()
            .find(|assembly| assembly.id() == ASSEMBLY)
            .unwrap_or_else(|| panic!("{id} has no assembly {ASSEMBLY}"));
        group.bench_function(BenchmarkId::from_parameter(id), |b| {
            b.iter(|| {
                structure
                    .assemble(black_box(assembly.id()), black_box(model))
                    .expect("assembly expands")
            })
        });
    }
    group.finish();
}

fn sweep_corpus(c: &mut Criterion) {
    let readers: [(&str, Reader); 2] = [("cif", read_text), ("bcif", read_binary)];
    let mut group = c.benchmark_group("sweep");
    group
        .sample_size(SWEEP_SAMPLES)
        .sampling_mode(SamplingMode::Flat);
    for (extension, read) in readers {
        let files = corpus_files(extension);
        group.throughput(Throughput::Bytes(
            files.iter().map(|path| file_len(path)).sum(),
        ));
        group.bench_with_input(
            BenchmarkId::from_parameter(extension),
            files.as_slice(),
            |b, files| {
                b.iter(|| {
                    for path in files {
                        black_box(read(black_box(path)));
                    }
                })
            },
        );
    }
    group.finish();
}

criterion_group!(
    benches,
    read_structures,
    pick_single_conformer,
    expand_assemblies,
    sweep_corpus
);
criterion_main!(benches);

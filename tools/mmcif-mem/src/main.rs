mod allocator;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use mmcif_model::{Structure, read_bcif, read_mmcif};
use serde_json::{Value, json};

use allocator::CountingAllocator;

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

const USAGE: &str = "usage: mmcif-mem [--corpus] [FILE.cif | FILE.bcif]...

Reads each structure, recording peak heap bytes above the level before the read and wall
time, and prints JSON: {\"files\": {NAME: {\"peak_bytes\", \"seconds\"}}, \"sweep_peak_bytes\"}.

  --corpus    add every .cif and .bcif in $MMCIF_CORPUS_DIR (default <workspace>/corpus/data)";

type Reader = fn(&Path) -> Result<Structure, String>;

struct Input {
    path: PathBuf,
    read: Reader,
}

struct Measurement {
    peak_bytes: usize,
    seconds: f64,
    absolute_peak_bytes: usize,
}

enum Command {
    Help,
    Measure(BTreeMap<String, Input>),
}

fn reader_for(path: &Path) -> Option<Reader> {
    match path.extension().and_then(|extension| extension.to_str()) {
        Some("cif") => Some(|path| read_mmcif(path).map_err(|error| error.to_string())),
        Some("bcif") => Some(|path| read_bcif(path).map_err(|error| error.to_string())),
        _ => None,
    }
}

fn corpus_dir() -> PathBuf {
    std::env::var_os("MMCIF_CORPUS_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/../../corpus/data")))
}

fn corpus_files() -> Result<Vec<PathBuf>, String> {
    let dir = corpus_dir();
    let entries = std::fs::read_dir(&dir).map_err(|error| {
        format!(
            "corpus directory {}: {error}; run `harness/.venv/bin/python corpus/fetch.py` or set MMCIF_CORPUS_DIR",
            dir.display()
        )
    })?;
    let mut files = Vec::new();
    for entry in entries {
        let path = entry
            .map_err(|error| format!("{}: {error}", dir.display()))?
            .path();
        if reader_for(&path).is_some() {
            files.push(path);
        }
    }
    if files.is_empty() {
        return Err(format!("no .cif or .bcif files in {}", dir.display()));
    }
    Ok(files)
}

fn add_input(inputs: &mut BTreeMap<String, Input>, path: PathBuf) -> Result<(), String> {
    let read = reader_for(&path)
        .ok_or_else(|| format!("{}: expected a .cif or .bcif file", path.display()))?;
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| format!("{}: file name is not UTF-8", path.display()))?
        .to_owned();
    if inputs.contains_key(&name) {
        return Err(format!(
            "{}: another input is also named {name}",
            path.display()
        ));
    }
    inputs.insert(name, Input { path, read });
    Ok(())
}

fn parse_args(args: impl Iterator<Item = String>) -> Result<Command, String> {
    let mut paths = Vec::new();
    for arg in args {
        match arg.as_str() {
            "-h" | "--help" => return Ok(Command::Help),
            "--corpus" => paths.extend(corpus_files()?),
            option if option.starts_with('-') => return Err(format!("unknown option {option}")),
            _ => paths.push(PathBuf::from(arg)),
        }
    }
    if paths.is_empty() {
        return Err("no input files".to_owned());
    }
    let mut inputs = BTreeMap::new();
    for path in paths {
        add_input(&mut inputs, path)?;
    }
    Ok(Command::Measure(inputs))
}

fn measure(input: &Input) -> Result<Measurement, String> {
    let baseline_bytes = allocator::reset_peak();
    let started = Instant::now();
    let structure = (input.read)(&input.path)?;
    let seconds = started.elapsed().as_secs_f64();
    drop(structure);
    let absolute_peak_bytes = allocator::peak_bytes();
    Ok(Measurement {
        peak_bytes: absolute_peak_bytes - baseline_bytes,
        seconds,
        absolute_peak_bytes,
    })
}

fn report(inputs: &BTreeMap<String, Input>) -> Result<Value, String> {
    let mut files = BTreeMap::new();
    let mut sweep_peak_bytes = 0;
    for (name, input) in inputs {
        let measurement =
            measure(input).map_err(|error| format!("{}: {error}", input.path.display()))?;
        sweep_peak_bytes = sweep_peak_bytes.max(measurement.absolute_peak_bytes);
        files.insert(
            name.as_str(),
            json!({ "peak_bytes": measurement.peak_bytes, "seconds": measurement.seconds }),
        );
    }
    Ok(json!({ "files": files, "sweep_peak_bytes": sweep_peak_bytes }))
}

fn main() -> ExitCode {
    let inputs = match parse_args(std::env::args().skip(1)) {
        Ok(Command::Measure(inputs)) => inputs,
        Ok(Command::Help) => {
            println!("{USAGE}");
            return ExitCode::SUCCESS;
        }
        Err(error) => {
            eprintln!("mmcif-mem: {error}\n\n{USAGE}");
            return ExitCode::from(2);
        }
    };
    match report(&inputs) {
        Ok(value) => {
            println!("{value:#}");
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("mmcif-mem: {error}");
            ExitCode::FAILURE
        }
    }
}

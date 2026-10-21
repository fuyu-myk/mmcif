use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::process::ExitCode;

use mmcif_ccd::{BuildSummary, CcdError, build_index_from_file};

const USAGE: &str = "usage: mmcif-ccd-build <components.cif> <out.idx> --release <tag>";

struct Args {
    source: String,
    output: String,
    release: String,
}

fn parse_args(mut args: impl Iterator<Item = String>) -> Option<Args> {
    let mut positional = Vec::new();
    let mut release = None;
    while let Some(arg) = args.next() {
        if arg == "--release" {
            release = Some(args.next()?);
        } else if arg.starts_with("--") {
            return None;
        } else {
            positional.push(arg);
        }
    }
    let [source, output] = <[String; 2]>::try_from(positional).ok()?;
    Some(Args {
        source,
        output,
        release: release?,
    })
}

fn run(args: &Args) -> Result<BuildSummary, CcdError> {
    let out = BufWriter::new(File::create(&args.output)?);
    build_index_from_file(&args.source, &args.release, out)
}

fn main() -> ExitCode {
    let Some(args) = parse_args(env::args().skip(1)) else {
        eprintln!("{USAGE}");
        return ExitCode::from(2);
    };
    match run(&args) {
        Ok(summary) => {
            println!(
                "indexed {} components ({} atoms, {} bonds) for release {} into {} ({} bytes)",
                summary.components(),
                summary.atoms(),
                summary.bonds(),
                args.release,
                args.output,
                summary.bytes(),
            );
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("mmcif-ccd-build: {}: {error}", args.source);
            ExitCode::FAILURE
        }
    }
}

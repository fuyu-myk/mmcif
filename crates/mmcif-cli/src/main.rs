use std::process::ExitCode;

use mmcif_ir::{Corner, Instantiated, SegmentId};
use mmcif_lower::{Backend, lower};
use mmcif_spec::{Checked, Report, check_path};

const COMMANDS: &str = "usage: mmcif <command> <specification>

commands:
  check <spec.yaml>              report problems in a design specification
  show <spec.yaml>               print what the specification resolves to
  lower <spec.yaml> --to <name>  write the specification for a backend

backends:";

const OPTIONS: &str = "options:
  -h, --help                     print this message
  -V, --version                  print the version";

fn help() -> String {
    let backends: String = Backend::all()
        .iter()
        .map(|backend| format!("  {:30} {}\n", backend.name(), backend.schema()))
        .collect();
    format!("{COMMANDS}\n{backends}\n{OPTIONS}")
}

fn main() -> ExitCode {
    let arguments: Vec<String> = std::env::args().skip(1).collect();
    match arguments.first().map(String::as_str) {
        None => usage("no command given"),
        Some("-h" | "--help") => {
            println!("{}", help());
            ExitCode::SUCCESS
        }
        Some("-V" | "--version") => {
            println!("mmcif {}", env!("CARGO_PKG_VERSION"));
            ExitCode::SUCCESS
        }
        Some(command @ ("check" | "show")) => match arguments.get(1) {
            None => usage(&format!("{command} needs a specification")),
            Some(_) if arguments.len() > 2 => usage(&format!(
                "{command} takes one specification, not {}",
                arguments.len() - 1
            )),
            Some(path) => run(command, path),
        },
        Some("lower") => match arguments.get(1..) {
            Some([path, flag, name]) if flag == "--to" => match Backend::named(name) {
                Some(backend) => emit(path, backend),
                None => usage(&format!(
                    "unknown backend {name}, expected {}",
                    Backend::names().join(" or ")
                )),
            },
            Some([_]) | Some([_, _]) => usage("lower needs a specification and --to <backend>"),
            _ => usage("lower takes one specification and --to <backend>"),
        },
        Some(command) => usage(&format!("unknown command {command}")),
    }
}

fn usage(problem: &str) -> ExitCode {
    eprintln!("mmcif: {problem}\n\n{}", help());
    ExitCode::from(2)
}

fn run(command: &str, path: &str) -> ExitCode {
    let checked = match check_path(path) {
        Ok(checked) => checked,
        Err(error) => {
            eprintln!("mmcif: {error}");
            return ExitCode::from(2);
        }
    };
    match command {
        "show" => show(&checked),
        _ => report(&checked),
    }
}

fn emit(path: &str, backend: Backend) -> ExitCode {
    let checked = match check_path(path) {
        Ok(checked) => checked,
        Err(error) => {
            eprintln!("mmcif: {error}");
            return ExitCode::from(2);
        }
    };
    let rendered = checked.render();
    if !rendered.is_empty() {
        eprint!("{rendered}");
    }
    let Some(design) = checked.design.as_ref().filter(|_| !checked.has_errors()) else {
        eprintln!("mmcif: the specification has problems, nothing was written");
        return ExitCode::FAILURE;
    };
    match lower(backend, checked.source, design) {
        Ok(text) => {
            print!("{text}");
            ExitCode::SUCCESS
        }
        Err(problems) => {
            let mut report = Report::new();
            report.extend(problems);
            eprint!("{}", checked.sources.render_all(&report));
            eprintln!(
                "mmcif: {} cannot express this specification",
                backend.name()
            );
            ExitCode::FAILURE
        }
    }
}

fn report(checked: &Checked) -> ExitCode {
    let rendered = checked.render();
    if !rendered.is_empty() {
        eprint!("{rendered}");
    }
    let errors = checked.report.errors().count();
    let warnings = checked.report.len() - errors;
    match (errors, warnings) {
        (0, 0) => {
            println!("no problems found");
            ExitCode::SUCCESS
        }
        (0, warnings) => {
            println!("{} passed", counted(warnings, "warning"));
            ExitCode::SUCCESS
        }
        (errors, 0) => {
            println!("{} found", counted(errors, "error"));
            ExitCode::FAILURE
        }
        (errors, warnings) => {
            println!(
                "{} and {} found",
                counted(errors, "error"),
                counted(warnings, "warning")
            );
            ExitCode::FAILURE
        }
    }
}

fn show(checked: &Checked) -> ExitCode {
    let code = report(checked);
    let Some(shortest) = checked.shortest.as_ref().filter(|_| !checked.has_errors()) else {
        return code;
    };
    println!();
    for id in shortest.segment_ids() {
        println!("{}", segment_line(shortest, id));
    }
    let constraints = shortest.constraints().len();
    let groups = shortest.groups().len();
    println!(
        "{} over {} slots, {}",
        counted(constraints, "constraint"),
        shortest.slot_count(),
        counted(groups, "group")
    );
    code
}

fn segment_line(shortest: &Instantiated, id: &SegmentId) -> String {
    let Ok(span) = shortest.segment_slots(id) else {
        return format!("{id}: unresolved");
    };
    let mut counts = [0usize; 5];
    for slot in span.iter() {
        let corner = shortest
            .determinacy(slot)
            .ok()
            .and_then(|slot| slot.corner());
        let position = match corner {
            Some(Corner::DeNovo) => 0,
            Some(Corner::FixedBackbone) => 1,
            Some(Corner::Flexible) => 2,
            Some(Corner::Context) => 3,
            _ => 4,
        };
        counts[position] += 1;
    }
    let names = [
        "de novo",
        "fixed backbone",
        "flexible",
        "context",
        "partial",
    ];
    let mix: Vec<String> = counts
        .iter()
        .zip(names)
        .filter(|(count, _)| **count > 0)
        .map(|(count, name)| format!("{count} {name}"))
        .collect();
    format!("{id}: {} ({})", counted(span.len(), "slot"), mix.join(", "))
}

fn counted(count: usize, name: &str) -> String {
    match count {
        1 => format!("1 {name}"),
        _ => format!("{count} {name}s"),
    }
}

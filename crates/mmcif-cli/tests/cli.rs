use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const TARGET: &str = "data_target
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
ATOM 1 N N  CYS A 1 1 0.000  0.000 0.000 1 A
ATOM 2 C CA CYS A 1 1 1.500  0.000 0.000 1 A
ATOM 3 S SG CYS A 1 1 2.000  1.000 0.000 1 A
ATOM 4 N N  GLY A 1 2 3.800  0.000 0.000 2 A
ATOM 5 C CA GLY A 1 2 5.000  0.000 0.000 2 A
ATOM 6 N N  CYS A 1 3 7.600  0.000 0.000 3 A
ATOM 7 C CA CYS A 1 3 8.000  0.000 0.000 3 A
ATOM 8 S SG CYS A 1 3 20.000 1.000 0.000 3 A
";

fn workspace(name: &str) -> PathBuf {
    let directory = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    std::fs::create_dir_all(&directory).expect("a writable test directory");
    std::fs::write(directory.join("target.cif"), TARGET).expect("the target is written");
    directory
}

fn spec(name: &str, yaml: &str) -> PathBuf {
    let path = workspace(name).join("design.yaml");
    std::fs::write(&path, yaml).expect("the specification is written");
    path
}

fn mmcif(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_mmcif"))
        .args(arguments)
        .output()
        .expect("the binary runs")
}

fn run(arguments: &[&str]) -> (i32, String, String) {
    let output = mmcif(arguments);
    (
        output.status.code().expect("the process exits normally"),
        String::from_utf8(output.stdout).expect("stdout is text"),
        String::from_utf8(output.stderr).expect("stderr is text"),
    )
}

const CLEAN: &str = "entities:
  - file:
      path: target.cif
      include:
        - chain:
            id: A
  - protein:
      id: S
      sequence: 10..14C6C3
constraints:
  - bond:
      atom1: [S, 11, SG]
      atom2: [S, 18, SG]
";

#[test]
fn a_clean_specification_passes_quietly() {
    let path = spec("clean", CLEAN);
    let (code, out, err) = run(&["check", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert_eq!(out, "no problems found\n");
    assert!(err.is_empty(), "{err}");
}

#[test]
fn an_infeasible_bond_fails_with_its_diagnostic() {
    let path = spec(
        "infeasible",
        concat!(
            "entities:\n",
            "  - file:\n",
            "      path: target.cif\n",
            "      include:\n",
            "        - chain:\n",
            "            id: A\n",
            "constraints:\n",
            "  - bond:\n",
            "      atom1: [A, 1, SG]\n",
            "      atom2: [A, 3, SG]\n",
        ),
    );
    let (code, out, err) = run(&["check", path.to_str().unwrap()]);
    assert_eq!(code, 1);
    assert_eq!(out, "1 error found\n");
    assert!(err.contains("error[infeasible-geometry]"), "{err}");
    assert!(err.contains("--> design.yaml:8:5"), "{err}");
}

#[test]
fn warnings_alone_still_pass() {
    let path = spec(
        "warnings",
        concat!(
            "entities:\n",
            "  - file:\n",
            "      path: target.cif\n",
            "      structure_groups:\n",
            "        - group:\n",
            "            visibility: 1\n",
            "            id: A\n",
        ),
    );
    let (code, out, err) = run(&["check", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert_eq!(out, "1 warning passed\n");
    assert!(err.contains("warning[unsupported-key]"), "{err}");
}

#[test]
fn show_prints_what_each_chain_resolves_to() {
    let path = spec("show", CLEAN);
    let (code, out, _) = run(&["show", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert_eq!(
        out,
        concat!(
            "no problems found\n",
            "\n",
            "A: 3 slots (3 context)\n",
            "S: 21 slots (19 de novo, 2 partial)\n",
            "1 constraint over 24 slots, 0 groups\n",
        )
    );
}

#[test]
fn show_on_a_broken_specification_reports_and_stops() {
    let path = spec(
        "broken",
        "entities:\n  - protein:\n      id: S\n      sequence: 0\n",
    );
    let (code, out, err) = run(&["show", path.to_str().unwrap()]);
    assert_eq!(code, 1);
    assert_eq!(out, "1 error found\n");
    assert!(err.contains("error[bad-sequence]"), "{err}");
}

#[test]
fn usage_problems_are_told_apart_from_specification_problems() {
    let path = spec("usage", CLEAN);
    let one = path.to_str().unwrap();
    assert_eq!(run(&[]).0, 2);
    assert_eq!(run(&["check"]).0, 2);
    assert_eq!(run(&["lower", one]).0, 2);
    assert_eq!(run(&["check", one, one]).0, 2);
    assert_eq!(run(&["check", "nowhere.yaml"]).0, 2);
    assert!(run(&["nonsense"]).2.contains("unknown command nonsense"));
    assert!(run(&["check", "nowhere.yaml"]).2.contains("nowhere.yaml"));
}

#[test]
fn help_and_version_are_available() {
    let (code, out, _) = run(&["--help"]);
    assert_eq!(code, 0);
    assert!(out.contains("usage: mmcif"));
    assert!(out.contains("check <spec.yaml>"));
    let (code, out, _) = run(&["-V"]);
    assert_eq!(code, 0);
    assert!(out.starts_with("mmcif "));
}

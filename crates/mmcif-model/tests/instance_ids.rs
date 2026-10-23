use mmcif_model::{AltlocPolicy, Model, Structure, parse_mmcif};

const FIXTURE: &str = include_str!("fixtures/pathologies.cif");
const IDENTITY: &str = "1 'identity operation' 1 0 0 0 0 1 0 0 0 0 1 0\n";

fn structure(text: &str) -> Structure {
    parse_mmcif(text).unwrap()
}

fn instance_ids(model: &Model) -> Vec<&str> {
    model.chains().map(|chain| chain.instance_id()).collect()
}

fn assembled_ids(text: &str, assembly: &str) -> Vec<String> {
    let structure = structure(text);
    let model = structure
        .assemble(assembly, structure.first_model().unwrap())
        .unwrap();
    instance_ids(&model)
        .into_iter()
        .map(str::to_owned)
        .collect()
}

#[test]
fn a_deposited_chain_is_named_by_its_author_id() {
    let structure = structure(FIXTURE);
    for model in structure.models() {
        for chain in model.chains() {
            assert_eq!(chain.instance_id(), chain.auth_asym().as_str());
        }
    }
}

#[test]
fn copies_append_the_ids_of_the_operators_that_placed_them() {
    assert_eq!(assembled_ids(FIXTURE, "2"), ["H1", "H1", "H1-2", "H1-2"]);
    assert_eq!(assembled_ids(FIXTURE, "3"), ["H1-2-3"]);
    assert_eq!(assembled_ids(FIXTURE, "1"), ["H1", "H1", "H1", "H1"]);
}

#[test]
fn an_identity_within_the_tolerance_keeps_the_author_id() {
    let nearly = FIXTURE.replace(
        IDENTITY,
        "1 'identity operation' 0.9999999999999999 0 0 0 0 1 0 0 0 0 1 0.0005\n",
    );
    assert_ne!(nearly, FIXTURE);
    assert_eq!(assembled_ids(&nearly, "2"), ["H1", "H1", "H1-2", "H1-2"]);

    let moved = FIXTURE.replace(
        IDENTITY,
        "1 'identity operation' 1 0 0 0 0 1 0 0 0 0 1 0.01\n",
    );
    assert_eq!(assembled_ids(&moved, "2"), ["H1-1", "H1-1", "H1-2", "H1-2"]);
}

#[test]
fn choosing_one_conformer_keeps_the_instance_ids() {
    let structure = structure(FIXTURE);
    let dimer = structure
        .assemble("2", structure.first_model().unwrap())
        .unwrap();
    let single = dimer.single_conformer(AltlocPolicy::HighestOccupancy);
    assert_eq!(instance_ids(&single), instance_ids(&dimer));
}

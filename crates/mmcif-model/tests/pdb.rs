use mmcif_model::{
    AltLoc, AuthResidueKey, AuthSeqId, EntityKind, InsCode, NumberingError, ReadError, parse_pdb,
    read_pdb,
};

const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/legacy.pdb");

fn auth(chain: &str, seq: i32) -> AuthResidueKey {
    AuthResidueKey {
        chain: chain.into(),
        seq: AuthSeqId::new(seq),
        ins_code: InsCode::NONE,
    }
}

#[test]
fn reads_models_and_header() {
    let structure = read_pdb(FIXTURE).unwrap();
    assert_eq!(structure.entry_id(), Some("1ABC"));
    let counts: Vec<usize> = structure.models().iter().map(|m| m.atom_count()).collect();
    assert_eq!(counts, [15, 1]);
    assert_eq!(structure.models()[1].number().get(), 2);
}

#[test]
fn splits_chains_into_polymer_ligand_and_water_label_chains() {
    let structure = read_pdb(FIXTURE).unwrap();
    let model = structure.first_model().unwrap();
    let chains: Vec<String> = model
        .chains()
        .iter()
        .map(|c| {
            format!(
                "{}:{}:{}",
                c.label_asym(),
                c.auth_asym(),
                c.residues().len()
            )
        })
        .collect();
    assert_eq!(
        chains,
        ["A:A:5", "B:B:1", "C:A:1", "D:A:1", "E:A:2", "F:B:1"]
    );
    let polymer = model.chain("A").unwrap();
    let comps: Vec<&str> = polymer
        .residues()
        .iter()
        .map(|r| r.comp().as_str())
        .collect();
    assert_eq!(comps, ["GLY", "SER", "ALA", "MSE", "LYS"]);
    assert!(polymer.residues()[3].is_hetero());
    assert_eq!(polymer.residues()[2].ins_code(), InsCode::new('A'));
    assert_eq!(
        polymer.residues()[1].altlocs(),
        [AltLoc::new('A'), AltLoc::new('B')]
    );
}

#[test]
fn synthesises_entities() {
    let structure = read_pdb(FIXTURE).unwrap();
    let summary: Vec<String> = structure
        .entities()
        .iter()
        .map(|e| {
            let asyms: Vec<&str> = e.asym_ids().iter().map(|a| a.as_str()).collect();
            format!("{}:{}:{}", e.id(), e.kind().as_mmcif_str(), asyms.join(","))
        })
        .collect();
    assert_eq!(
        summary,
        [
            "1:polymer:A",
            "2:polymer:B",
            "3:non-polymer:C",
            "4:non-polymer:D",
            "5:water:E,F"
        ]
    );
    let second_model_chain = &structure.models()[1].chains()[0];
    assert_eq!(second_model_chain.entity().unwrap().as_str(), "1");
    assert_eq!(structure.entity("5").unwrap().kind(), &EntityKind::Water);
}

#[test]
fn reads_atom_columns() {
    let structure = read_pdb(FIXTURE).unwrap();
    let model = structure.first_model().unwrap();
    let zinc = &model.chain("C").unwrap().residues()[0].atoms()[0];
    assert_eq!(zinc.element().symbol(), "Zn");
    assert_eq!(zinc.charge(), Some(2));
    assert_eq!(zinc.position(), [0.0, 0.0, 5.0]);
    let inferred = &model.chain("B").unwrap().residues()[0].atoms()[1];
    assert_eq!(inferred.element().symbol(), "C");
    let selenium = model.chain("A").unwrap().residues()[3].atom("SE").unwrap();
    assert_eq!(selenium.element().symbol(), "Se");
    assert_eq!(
        model.chain("A").unwrap().residues()[1].atoms()[1].occupancy(),
        0.4
    );
}

#[test]
fn legacy_files_have_no_label_numbering() {
    let structure = read_pdb(FIXTURE).unwrap();
    let numbering = structure.numbering();
    assert_eq!(
        numbering.to_label(&auth("A", 1)),
        Err(NumberingError::Unnumbered(auth("A", 1)))
    );
    assert_eq!(
        numbering.to_label(&auth("A", 101)),
        Err(NumberingError::NotPolymer(auth("A", 101)))
    );
    assert_eq!(
        numbering.to_label(&auth("Z", 1)),
        Err(NumberingError::UnknownAuthChain("Z".into()))
    );
    assert!(model_residue_has_no_label_seq(&structure));
}

fn model_residue_has_no_label_seq(structure: &mmcif_model::Structure) -> bool {
    structure
        .first_model()
        .unwrap()
        .chains()
        .iter()
        .flat_map(|c| c.residues())
        .all(|r| r.label_seq().is_none())
}

#[test]
fn reads_hybrid36_serials_and_residue_numbers() {
    let line = "ATOM  A0000  CA  GLY AA000       1.000   2.000   3.000  1.00 10.00           C  ";
    let structure = parse_pdb(line).unwrap();
    let chain = &structure.first_model().unwrap().chains()[0];
    assert_eq!(chain.residues()[0].auth_seq().get(), 10_000);
    assert_eq!(chain.residues()[0].atoms()[0].serial(), Some(100_000));
}

#[test]
fn reports_bad_lines() {
    let bad_x = "ATOM      1  CA  GLY A   1       x.000   2.000   3.000  1.00 10.00           C  ";
    assert!(matches!(
        parse_pdb(bad_x).unwrap_err(),
        ReadError::Pdb {
            line: 1,
            problem: "invalid x coordinate"
        }
    ));
    let bad_occupancy =
        "ATOM      1  CA  GLY A   1       1.000   2.000   3.000  abcd 10.00           C  ";
    assert!(matches!(
        parse_pdb(bad_occupancy).unwrap_err(),
        ReadError::Pdb {
            problem: "invalid occupancy",
            ..
        }
    ));
}

#[test]
fn ignores_files_without_atoms() {
    let structure = parse_pdb("HEADER    EMPTY\nEND\n").unwrap();
    assert!(structure.models().is_empty());
}

fn atom_line(record: &str, serial: u32, name: &str, comp: &str, chain: char, seq: i32) -> String {
    format!(
        "{record:<6}{serial:>5}  {name:<3} {comp:>3} {chain}{seq:>4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}  ",
        seq as f64,
        0.0,
        0.0,
        1.0,
        0.0,
        &name[..1]
    )
}

#[test]
fn atom_records_after_ter_start_a_new_polymer_chain() {
    let lines = [
        atom_line("ATOM", 1, "N", "ALA", ' ', 1),
        atom_line("ATOM", 2, "N", "GLY", ' ', 2),
        "TER".to_owned(),
        atom_line("ATOM", 3, "N", "SER", ' ', 1),
        atom_line("ATOM", 4, "N", "THR", ' ', 2),
        atom_line("ATOM", 5, "N", "VAL", ' ', 3),
        "TER".to_owned(),
        atom_line("HETATM", 6, "C1", "LIG", ' ', 4),
    ];
    let structure = parse_pdb(&lines.join("\n")).unwrap();
    let chains: Vec<String> = structure
        .first_model()
        .unwrap()
        .chains()
        .iter()
        .map(|c| format!("{}:{}", c.label_asym(), c.residues().len()))
        .collect();
    assert_eq!(chains, ["A:2", "B:3", "C:1"]);
    let kinds: Vec<&str> = structure
        .entities()
        .iter()
        .map(|e| e.kind().as_mmcif_str())
        .collect();
    assert_eq!(kinds, ["polymer", "polymer", "non-polymer"]);
}

#[test]
fn rejects_malformed_serials_and_charges() {
    let mut bad_charge = atom_line("HETATM", 1, "ZN", "ZN", 'A', 1);
    bad_charge.replace_range(78..80, "+2");
    assert!(matches!(
        parse_pdb(&bad_charge).unwrap_err(),
        ReadError::Pdb {
            problem: "invalid formal charge",
            ..
        }
    ));
    let mut bad_serial = atom_line("ATOM", 1, "N", "ALA", 'A', 1);
    bad_serial.replace_range(6..11, "1x2y3");
    assert!(matches!(
        parse_pdb(&bad_serial).unwrap_err(),
        ReadError::Pdb {
            problem: "invalid atom serial",
            ..
        }
    ));
}

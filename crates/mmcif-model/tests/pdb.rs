use mmcif_model::{
    AltLoc, AuthResidueKey, AuthSeqId, EntityKind, InsCode, NumberingError, ReadError, ReadOptions,
    Transform, parse_pdb, read_pdb,
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
    let comps: Vec<&str> = polymer.residues().map(|r| r.comp().as_str()).collect();
    assert_eq!(comps, ["GLY", "SER", "ALA", "MSE", "LYS"]);
    assert!(polymer.residues().nth(3).unwrap().is_hetero());
    assert_eq!(
        polymer.residues().nth(2).unwrap().ins_code(),
        InsCode::new('A')
    );
    assert_eq!(
        polymer.residues().nth(1).unwrap().altlocs(),
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
    let second_model_chain = &structure.models()[1].chains().next().unwrap();
    assert_eq!(second_model_chain.entity().unwrap().as_str(), "1");
    assert_eq!(structure.entity("5").unwrap().kind(), &EntityKind::Water);
}

#[test]
fn reads_atom_columns() {
    let structure = read_pdb(FIXTURE).unwrap();
    let model = structure.first_model().unwrap();
    let zinc = &model
        .chain("C")
        .unwrap()
        .residues()
        .next()
        .unwrap()
        .atoms()
        .next()
        .unwrap();
    assert_eq!(zinc.element().symbol(), "Zn");
    assert_eq!(zinc.charge(), Some(2));
    assert_eq!(zinc.position(), [0.0, 0.0, 5.0]);
    let inferred = &model
        .chain("B")
        .unwrap()
        .residues()
        .next()
        .unwrap()
        .atoms()
        .nth(1)
        .unwrap();
    assert_eq!(inferred.element().symbol(), "C");
    let selenium = model
        .chain("A")
        .unwrap()
        .residues()
        .nth(3)
        .unwrap()
        .atom("SE")
        .unwrap();
    assert_eq!(selenium.element().symbol(), "Se");
    assert_eq!(
        model
            .chain("A")
            .unwrap()
            .residues()
            .nth(1)
            .unwrap()
            .atoms()
            .nth(1)
            .unwrap()
            .occupancy(),
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
        .flat_map(|c| c.residues())
        .all(|r| r.label_seq().is_none())
}

#[test]
fn reads_hybrid36_serials_and_residue_numbers() {
    let line = "ATOM  A0000  CA  GLY AA000       1.000   2.000   3.000  1.00 10.00           C  ";
    let structure = parse_pdb(line).unwrap();
    let chain = &structure.first_model().unwrap().chains().next().unwrap();
    assert_eq!(chain.residues().next().unwrap().auth_seq().get(), 10_000);
    assert_eq!(
        chain
            .residues()
            .next()
            .unwrap()
            .atoms()
            .next()
            .unwrap()
            .serial(),
        Some(100_000)
    );
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
fn reads_a_formal_charge_with_its_sign_first() {
    let charges: Vec<Option<i8>> = ["+2", "-1"]
        .iter()
        .map(|field| {
            let mut line = atom_line("HETATM", 1, "ZN", "ZN", 'A', 1);
            line.replace_range(78..80, field);
            parse_pdb(&line).unwrap().first_model().unwrap().charges()[0]
        })
        .collect();
    assert_eq!(charges, [Some(2), Some(-1)]);
}

#[test]
fn reads_a_formal_charge_written_as_a_bare_digit() {
    let charges: Vec<Option<i8>> = [" 0", "0 ", " 2", "2+", "1-", "  "]
        .iter()
        .map(|field| {
            let mut line = atom_line("HETATM", 1, "ZN", "ZN", 'A', 1);
            line.replace_range(78..80, field);
            let structure = parse_pdb(&line).unwrap();
            structure.first_model().unwrap().charges()[0]
        })
        .collect();
    assert_eq!(
        charges,
        [Some(0), Some(0), Some(2), Some(2), Some(-1), None]
    );
}

#[test]
fn rejects_malformed_serials_and_charges() {
    let mut bad_charge = atom_line("HETATM", 1, "ZN", "ZN", 'A', 1);
    bad_charge.replace_range(78..80, "2x");
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

const REMARK_350: &str = "\
REMARK 350 BIOMOLECULE: 1
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC
REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: DIMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A,
REMARK 350                    AND CHAINS: B
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000       10.00000
REMARK 350   BIOMT2   2  0.000000 -1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
REMARK 350
REMARK 350 BIOMOLECULE: 2
REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: MONOMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350 APPLY THE FOLLOWING TO CHAINS: B
REMARK 350   BIOMT1   2  1.000000  0.000000  0.000000       20.00000
REMARK 350   BIOMT2   2  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
";

fn with_remark_350(remarks: &str) -> String {
    let legacy = std::fs::read_to_string(FIXTURE).unwrap();
    let (header, rest) = legacy.split_once('\n').unwrap();
    format!("{header}\n{remarks}{rest}")
}

#[test]
fn reads_assemblies_from_remark_350() {
    let structure = parse_pdb(&with_remark_350(REMARK_350)).unwrap();

    let operators: Vec<(&str, Transform)> = structure
        .operators()
        .iter()
        .map(|operator| (operator.id().as_str(), *operator.transform()))
        .collect();
    assert_eq!(operators.len(), 3);
    assert_eq!(operators[0], ("1", Transform::IDENTITY));
    assert_eq!(
        operators[1].1,
        Transform {
            rotation: [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [10.0, 0.0, 0.0],
        }
    );
    assert_eq!(operators[2].1.translation, [20.0, 0.0, 0.0]);

    let assemblies = structure.assemblies();
    let ids: Vec<&str> = assemblies.iter().map(|a| a.id().as_str()).collect();
    assert_eq!(ids, ["1", "2"]);
    let details: Vec<Option<&str>> = assemblies.iter().map(|a| a.details()).collect();
    assert_eq!(
        details,
        [
            Some("author_and_software_defined_assembly"),
            Some("software_defined_assembly")
        ]
    );
    let oligomers: Vec<Option<&str>> = assemblies.iter().map(|a| a.oligomeric_details()).collect();
    assert_eq!(oligomers, [Some("dimeric"), Some("monomeric")]);
    assert!(assemblies.iter().all(|a| a.oligomeric_count().is_none()));
    let generators = |index: usize| -> Vec<(String, Vec<String>)> {
        assemblies[index]
            .generators()
            .iter()
            .map(|g| {
                (
                    g.expression().to_owned(),
                    g.asym_ids().iter().map(|a| a.to_string()).collect(),
                )
            })
            .collect()
    };
    let strings = |ids: &[&str]| ids.iter().map(|id| id.to_string()).collect::<Vec<_>>();
    assert_eq!(
        generators(0),
        [("1,2".to_owned(), strings(&["A", "C", "D", "E", "B", "F"]))]
    );
    assert_eq!(
        generators(1),
        [
            ("1".to_owned(), strings(&["A", "C", "D", "E"])),
            ("3".to_owned(), strings(&["B", "F"])),
        ]
    );

    let model = structure.first_model().unwrap();
    let dimer = structure.assemble("1", model).unwrap();
    assert_eq!(dimer.atom_count(), 2 * model.atom_count());
    assert_eq!(structure.assembly_is_deposited("1", 1e-6), Ok(false));
    assert_eq!(structure.assembly_is_deposited("2", 1e-6), Ok(false));
}

#[test]
fn files_without_remark_350_have_no_assemblies() {
    let structure = read_pdb(FIXTURE).unwrap();
    assert!(structure.assemblies().is_empty() && structure.operators().is_empty());
}

#[test]
fn malformed_remark_350_fails_a_strict_read_and_is_skipped_by_a_lenient_one() {
    let bad_number = REMARK_350.replace(
        "BIOMT2   1  0.000000  1.000000",
        "BIOMT2   1  0.000000  x.000000",
    );
    let early = REMARK_350.replace("REMARK 350 APPLY THE FOLLOWING TO CHAINS: A,\n", "");
    let short = REMARK_350.replace(
        "REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000\n\
         REMARK 350\n",
        "REMARK 350\n",
    );
    for remarks in [bad_number, early, short] {
        assert_ne!(remarks, REMARK_350);
        let text = with_remark_350(&remarks);
        assert!(
            matches!(parse_pdb(&text), Err(ReadError::Pdb { .. })),
            "{remarks}"
        );

        let path =
            std::env::temp_dir().join(format!("mmcif-remark-350-{}.pdb", std::process::id()));
        std::fs::write(&path, &text).unwrap();
        let loaded = mmcif_model::read_with(&path, &ReadOptions::new().lenient(true)).unwrap();
        assert!(
            loaded.structure.assemblies().is_empty() && loaded.structure.operators().is_empty()
        );
        assert!(matches!(loaded.skipped.as_slice(), [ReadError::Pdb { .. }]));
        assert_eq!(loaded.structure.first_model().unwrap().atom_count(), 15);
    }
}

#[test]
fn a_blank_b_factor_is_none_and_a_non_finite_one_is_invalid() {
    let blank = atom_line("ATOM", 1, "N", "GLY", 'A', 1);
    let mut columns: Vec<char> = blank.chars().collect();
    columns.resize(80, ' ');
    let with_b = |text: &str| -> String {
        let mut line = columns.clone();
        for (offset, character) in format!("{text:>6}").chars().enumerate() {
            line[60 + offset] = character;
        }
        line.into_iter().collect::<String>() + "\n"
    };
    let structure = parse_pdb(&with_b("")).unwrap();
    let atom = structure.first_model().unwrap().atoms().next().unwrap();
    assert_eq!(atom.b_factor(), None);
    assert_eq!(
        parse_pdb(&with_b("12.50"))
            .unwrap()
            .first_model()
            .unwrap()
            .atoms()
            .next()
            .unwrap()
            .b_factor(),
        Some(12.5)
    );
    for value in ["nan", "inf"] {
        assert!(
            matches!(
                parse_pdb(&with_b(value)),
                Err(ReadError::Pdb { line: 1, .. })
            ),
            "{value}"
        );
    }
}

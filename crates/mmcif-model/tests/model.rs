use mmcif_model::{
    AltLoc, AltlocPolicy, AssemblyError, AssemblyId, AuthResidueKey, AuthSeqId, EntityKind,
    InsCode, LabelResidueKey, LabelSeqId, ModelNumber, NumberingError, PolymerType, ReadError,
    Structure, parse_mmcif, read_mmcif,
};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/pathologies.cif"
);

fn structure() -> Structure {
    read_mmcif(FIXTURE).unwrap()
}

fn auth(seq: i32, ins: Option<char>) -> AuthResidueKey {
    AuthResidueKey {
        chain: "H1".into(),
        seq: AuthSeqId::new(seq),
        ins_code: ins.map_or(InsCode::NONE, InsCode::new),
    }
}

fn label(seq: u32) -> LabelResidueKey {
    LabelResidueKey {
        asym: "AA".into(),
        seq: LabelSeqId::new(seq).unwrap(),
    }
}

#[test]
fn keeps_models_separate() {
    let structure = structure();
    assert_eq!(structure.entry_id(), Some("TEST"));
    let numbers: Vec<i32> = structure
        .models()
        .iter()
        .map(|m| m.number().get())
        .collect();
    assert_eq!(numbers, [1, 2]);
    assert_eq!(structure.first_model().unwrap().atom_count(), 16);
    assert_eq!(
        structure.model(ModelNumber::new(2)).unwrap().atom_count(),
        1
    );
}

#[test]
fn groups_atoms_into_label_chains_with_author_ids() {
    let structure = structure();
    let model = structure.first_model().unwrap();
    let labels: Vec<&str> = model
        .chains()
        .iter()
        .map(|c| c.label_asym().as_str())
        .collect();
    assert_eq!(labels, ["AA", "B", "C", "D"]);
    assert_eq!(model.auth_chain("H1").count(), 4);
    let waters = model.chain("D").unwrap();
    assert_eq!(waters.residues().len(), 2);
    assert_eq!(waters.entity().unwrap().as_str(), "4");
}

#[test]
fn insertion_codes_give_distinct_residues() {
    let structure = structure();
    let model = structure.first_model().unwrap();
    let plain: Vec<&str> = model
        .residues_at(&auth(100, None))
        .iter()
        .map(|(_, r)| r.comp().as_str())
        .collect();
    assert_eq!(plain, ["SER", "THR"]);
    let inserted = model.residues_at(&auth(100, Some('A')));
    assert_eq!(inserted.len(), 1);
    assert_eq!(inserted[0].1.comp().as_str(), "ALA");
    assert_eq!(inserted[0].1.label_seq(), LabelSeqId::new(3));
}

#[test]
fn microheterogeneity_keeps_both_components() {
    let structure = structure();
    let chain = structure.first_model().unwrap().chain("AA").unwrap();
    let at_two: Vec<&str> = chain
        .residues_at_label(LabelSeqId::new(2).unwrap())
        .map(|r| r.comp().as_str())
        .collect();
    assert_eq!(at_two, ["SER", "THR"]);
    let polymer = structure.entity("1").unwrap().polymer().unwrap();
    let position = &polymer.sequence()[1];
    assert!(position.is_heterogeneous());
    let comps: Vec<&str> = position.comps().iter().map(|c| c.as_str()).collect();
    assert_eq!(comps, ["SER", "THR"]);
}

#[test]
fn residue_keys_are_unique_and_resolvable() {
    let structure = structure();
    let model = structure.first_model().unwrap();
    let mut keys = Vec::new();
    for chain in model.chains() {
        for residue in chain.residues() {
            let key = chain.residue_key(residue);
            assert_eq!(model.residue(&key), Some(residue));
            keys.push(key);
        }
    }
    let total = keys.len();
    keys.sort();
    keys.dedup();
    assert_eq!(keys.len(), total);
}

#[test]
fn altloc_conformers_and_selection() {
    let structure = structure();
    let model = structure.first_model().unwrap();
    let chain = model.chain("AA").unwrap();
    let alanine = chain
        .residues_at_auth(AuthSeqId::new(100), InsCode::new('A'))
        .next()
        .unwrap();
    assert_eq!(alanine.altlocs(), [AltLoc::new('A'), AltLoc::new('B')]);
    let conformer_b: Vec<[f64; 3]> = alanine
        .conformer(AltLoc::new('B'))
        .map(|a| a.position())
        .collect();
    assert_eq!(conformer_b.len(), 3);
    assert_eq!(conformer_b[2], [4.5, 1.0, 3.0]);

    let comps = |policy| -> Vec<String> {
        model
            .single_conformer(policy)
            .chain("AA")
            .unwrap()
            .residues()
            .iter()
            .map(|r| format!("{}:{}", r.comp(), r.atoms().len()))
            .collect()
    };
    assert_eq!(
        comps(AltlocPolicy::First),
        ["SER:2", "ALA:3", "LYS:1", "VAL:1"]
    );
    assert_eq!(
        comps(AltlocPolicy::HighestOccupancy),
        ["SER:2", "ALA:3", "LYS:1", "VAL:1"]
    );
    assert_eq!(
        comps(AltlocPolicy::Prefer(AltLoc::new('B'))),
        ["THR:2", "ALA:3", "LYS:1", "VAL:1"]
    );
    assert_eq!(model.single_conformer(AltlocPolicy::First).atom_count(), 13);
}

#[test]
fn numbering_converts_between_label_and_author() {
    let structure = structure();
    let numbering = structure.numbering();
    assert_eq!(numbering.to_auth(&label(1)).unwrap(), auth(99, None));
    assert_eq!(numbering.to_auth(&label(2)).unwrap(), auth(100, None));
    assert_eq!(numbering.to_auth(&label(3)).unwrap(), auth(100, Some('A')));
    assert_eq!(numbering.to_auth(&label(6)).unwrap(), auth(103, None));
    assert_eq!(numbering.to_label(&auth(100, Some('A'))).unwrap(), label(3));
    assert_eq!(numbering.to_label(&auth(101, None)).unwrap(), label(4));
    assert_eq!(
        numbering.to_label(&auth(201, None)),
        Err(NumberingError::NotPolymer(auth(201, None)))
    );
    assert_eq!(
        numbering.to_auth(&label(7)),
        Err(NumberingError::NotInSequence(label(7)))
    );
}

#[test]
fn reads_entities() {
    let structure = structure();
    let kinds: Vec<&EntityKind> = structure.entities().iter().map(|e| e.kind()).collect();
    assert_eq!(
        kinds,
        [
            &EntityKind::Polymer,
            &EntityKind::Branched,
            &EntityKind::NonPolymer,
            &EntityKind::Water
        ]
    );
    let polymer = structure.entity("1").unwrap();
    assert_eq!(polymer.description(), Some("Heavy chain"));
    assert_eq!(polymer.polymer().unwrap().kind(), &PolymerType::PeptideL);
    assert_eq!(
        polymer.polymer().unwrap().canonical_sequence(),
        Some("GSAKVT")
    );
    assert_eq!(polymer.polymer().unwrap().sequence().len(), 6);
    let asyms: Vec<&str> = polymer.asym_ids().iter().map(|a| a.as_str()).collect();
    assert_eq!(asyms, ["AA"]);
}

#[test]
fn reads_atom_properties() {
    let structure = structure();
    let model = structure.first_model().unwrap();
    let zinc = model.chain("C").unwrap().residues()[0].atom("ZN").unwrap();
    assert_eq!(zinc.element().symbol(), "Zn");
    assert_eq!(zinc.charge(), Some(2));
    assert_eq!(zinc.serial(), Some(14));
    assert_eq!(zinc.b_factor(), 30.0);
    assert!(model.chain("C").unwrap().residues()[0].is_hetero());
}

#[test]
fn resolves_glycan_links() {
    let structure = structure();
    let model = structure.first_model().unwrap();
    let bonds = structure.branch_bonds(model, "B").unwrap();
    assert_eq!(bonds.len(), 1);
    let first = bonds[0].first.unwrap();
    let second = bonds[0].second.unwrap();
    assert_eq!(
        (first.residue.comp().as_str(), first.atom.name().as_str()),
        ("BMA", "C1")
    );
    assert_eq!(
        (second.residue.comp().as_str(), second.atom.name().as_str()),
        ("NAG", "O4")
    );
    assert_eq!(bonds[0].link.order(), Some("sing"));
    assert!(structure.branch_bonds(model, "AA").is_err());
}

#[test]
fn expands_assemblies() {
    let structure = structure();
    let model = structure.first_model().unwrap();

    let identity = structure.assemble(&AssemblyId::new("1"), model).unwrap();
    assert_eq!(identity.atom_count(), model.atom_count());
    let positions =
        |m: &mmcif_model::Model| m.atoms().map(|(_, _, a)| a.position()).collect::<Vec<_>>();
    assert_eq!(positions(&identity), positions(model));

    let dimer = structure.assemble(&AssemblyId::new("2"), model).unwrap();
    let copies: Vec<String> = dimer
        .chains()
        .iter()
        .map(|c| {
            format!(
                "{}x{}",
                c.label_asym(),
                c.operators()
                    .iter()
                    .map(|o| o.as_str())
                    .collect::<Vec<_>>()
                    .join("")
            )
        })
        .collect();
    assert_eq!(copies, ["AAx1", "Cx1", "AAx2", "Cx2"]);
    let rotated_zinc = dimer.chains()[3].residues()[0].atoms()[0].position();
    assert_eq!(rotated_zinc, [-0.0, -0.0, 5.0]);

    let product = structure.assemble(&AssemblyId::new("3"), model).unwrap();
    assert_eq!(product.chains().len(), 1);
    assert_eq!(
        product.chains()[0].residues()[0].atoms()[0].position(),
        [-11.0, -2.0, 3.0]
    );

    assert_eq!(
        structure.assemble(&AssemblyId::new("9"), model),
        Err(AssemblyError::UnknownAssembly(AssemblyId::new("9")))
    );
    assert_eq!(structure.assemblies()[1].operator_count(), 2);
}

fn with_atom_site(rows: &str) -> String {
    format!(
        "data_x\nloop_\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_comp_id\n\
         _atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.type_symbol\n_atom_site.label_seq_id\n\
         _atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n{rows}\n"
    )
}

#[test]
fn reports_malformed_values_with_location() {
    let error = parse_mmcif(&with_atom_site("A A GLY N . N 1 1 abc 0 0")).unwrap_err();
    assert!(matches!(
        error,
        ReadError::InvalidValue {
            category: "atom_site",
            item: "Cartn_x",
            row: 0,
            ..
        }
    ));
    let error = parse_mmcif(&with_atom_site("A A GLY N AB N 1 1 0 0 0")).unwrap_err();
    assert!(matches!(
        error,
        ReadError::InvalidValue {
            item: "label_alt_id",
            ..
        }
    ));
    let error = parse_mmcif(&with_atom_site("A A GLY N . N 0 1 0 0 0")).unwrap_err();
    assert!(matches!(
        error,
        ReadError::InvalidValue {
            item: "label_seq_id",
            ..
        }
    ));
    let error = parse_mmcif(&with_atom_site("A A GLY N . Q1 1 1 0 0 0")).unwrap_err();
    assert!(matches!(
        error,
        ReadError::InvalidValue {
            item: "type_symbol",
            ..
        }
    ));
}

#[test]
fn rejects_label_chains_split_across_author_chains() {
    let error = parse_mmcif(&with_atom_site(
        "A A GLY N . N 1 1 0 0 0\nA B GLY CA . C 1 1 0 0 0",
    ))
    .unwrap_err();
    assert!(matches!(error, ReadError::InconsistentAuthChain { .. }));
}

#[test]
fn requires_atom_site() {
    assert!(matches!(
        parse_mmcif("data_x\n_entry.id x\n").unwrap_err(),
        ReadError::MissingCategory("atom_site")
    ));
    assert!(matches!(
        parse_mmcif("").unwrap_err(),
        ReadError::NoDataBlock
    ));
}

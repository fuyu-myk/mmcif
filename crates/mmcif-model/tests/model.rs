use mmcif_model::{
    AltLoc, AltlocPolicy, AssemblyError, AssemblyId, AuthResidueKey, AuthSeqId, EntityKind,
    InsCode, LabelResidueKey, LabelSeqId, ModelNumber, NumberingError, OperatorId, PolymerType,
    ReadError, Structure, parse_mmcif, read_mmcif,
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
    let labels: Vec<&str> = model.chains().map(|c| c.label_asym().as_str()).collect();
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
        .map(|r| r.comp().as_str())
        .collect();
    assert_eq!(plain, ["SER", "THR"]);
    let inserted = model.residues_at(&auth(100, Some('A')));
    assert_eq!(inserted.len(), 1);
    assert_eq!(inserted[0].comp().as_str(), "ALA");
    assert_eq!(inserted[0].label_seq(), LabelSeqId::new(3));
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
            let key = residue.key();
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
    let zinc = model
        .chain("C")
        .unwrap()
        .residues()
        .next()
        .unwrap()
        .atom("ZN")
        .unwrap();
    assert_eq!(zinc.element().symbol(), "Zn");
    assert_eq!(zinc.charge(), Some(2));
    assert_eq!(zinc.serial(), Some(14));
    assert_eq!(zinc.b_factor(), Some(30.0));
    assert!(
        model
            .chain("C")
            .unwrap()
            .residues()
            .next()
            .unwrap()
            .is_hetero()
    );
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
    let positions = |m: &mmcif_model::Model| m.atoms().map(|a| a.position()).collect::<Vec<_>>();
    assert_eq!(positions(&identity), positions(model));

    let dimer = structure.assemble(&AssemblyId::new("2"), model).unwrap();
    let copies: Vec<String> = dimer
        .chains()
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
    let rotated_zinc = dimer
        .chains()
        .nth(3)
        .unwrap()
        .residues()
        .next()
        .unwrap()
        .atoms()
        .next()
        .unwrap()
        .position();
    assert_eq!(rotated_zinc, [-0.0, -0.0, 5.0]);

    let product = structure.assemble(&AssemblyId::new("3"), model).unwrap();
    assert_eq!(product.chains().len(), 1);
    assert_eq!(
        product
            .chains()
            .next()
            .unwrap()
            .residues()
            .next()
            .unwrap()
            .atoms()
            .next()
            .unwrap()
            .position(),
        [-11.0, -2.0, 3.0]
    );

    assert_eq!(
        structure.assemble(&AssemblyId::new("9"), model),
        Err(AssemblyError::UnknownAssembly(AssemblyId::new("9")))
    );
    assert_eq!(structure.assemblies()[1].operator_count(), 2);
}

fn replaced(text: &str, from: &str, to: &str) -> String {
    assert!(text.contains(from), "the fixture lacks {from:?}");
    text.replace(from, to)
}

fn with_oligomers(rows: [&str; 3]) -> String {
    let text = std::fs::read_to_string(FIXTURE).unwrap();
    let text = replaced(
        &text,
        "_pdbx_struct_assembly.details\n",
        "_pdbx_struct_assembly.details\n_pdbx_struct_assembly.oligomeric_details\n\
         _pdbx_struct_assembly.oligomeric_count\n",
    );
    let text = replaced(
        &text,
        "1 author_defined_assembly\n",
        &format!("1 author_defined_assembly {}\n", rows[0]),
    );
    let text = replaced(
        &text,
        "2 software_defined_assembly\n",
        &format!("2 software_defined_assembly {}\n", rows[1]),
    );
    replaced(
        &text,
        "3 software_defined_assembly\n",
        &format!("3 software_defined_assembly {}\n", rows[2]),
    )
}

#[test]
fn reads_oligomeric_details_and_counts() {
    let annotated = parse_mmcif(&with_oligomers(["monomeric 1", "dimeric 2", "? ?"])).unwrap();
    let assemblies = annotated.assemblies();
    assert_eq!(assemblies[0].oligomeric_details(), Some("monomeric"));
    assert_eq!(assemblies[0].oligomeric_count(), Some(1));
    assert_eq!(assemblies[1].oligomeric_count(), Some(2));
    assert_eq!(
        (
            assemblies[2].oligomeric_details(),
            assemblies[2].oligomeric_count()
        ),
        (None, None)
    );
    let plain = structure();
    let unannotated = &plain.assemblies()[0];
    assert_eq!(
        (
            unannotated.oligomeric_details(),
            unannotated.oligomeric_count()
        ),
        (None, None)
    );

    let error = parse_mmcif(&with_oligomers(["monomeric one", "dimeric 2", "? ?"])).unwrap_err();
    assert!(matches!(
        error,
        ReadError::InvalidValue {
            category: "pdbx_struct_assembly",
            item: "oligomeric_count",
            row: 0,
            ..
        }
    ));
}

#[test]
fn tells_whether_an_assembly_shows_the_deposited_arrangement() {
    let structure = structure();
    let deposited = |id: &str| structure.assembly_is_deposited(&AssemblyId::new(id), 1e-3);
    assert_eq!(deposited("1"), Ok(true));
    assert_eq!(deposited("2"), Ok(false));
    assert_eq!(deposited("3"), Ok(true));
    assert_eq!(
        deposited("9"),
        Err(AssemblyError::UnknownAssembly(AssemblyId::new("9")))
    );
}

#[test]
fn a_deposited_assembly_may_span_several_generators() {
    let text = std::fs::read_to_string(FIXTURE).unwrap();
    let deposited = |text: &str| {
        parse_mmcif(text)
            .unwrap()
            .assembly_is_deposited(&AssemblyId::new("1"), 1e-3)
    };
    let split = replaced(&text, "1 1 AA,B,C,D\n", "1 1 AA,B\n1 1 C,D\n");
    assert_eq!(deposited(&split), Ok(true));
    let rearranged = replaced(&text, "1 1 AA,B,C,D\n", "1 1 AA,B\n1 2 C,D\n");
    assert_eq!(deposited(&rearranged), Ok(false));
    let copied = replaced(&text, "1 1 AA,B,C,D\n", "1 1 AA,B,C,D\n1 1 AA\n");
    assert_eq!(deposited(&copied), Ok(false));
    let undefined = replaced(&text, "1 1 AA,B,C,D\n", "1 7 AA,B,C,D\n");
    assert_eq!(
        deposited(&undefined),
        Err(AssemblyError::UnknownOperator {
            assembly: AssemblyId::new("1"),
            operator: OperatorId::new("7"),
        })
    );
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

#[test]
fn highest_occupancy_compares_mean_occupancy_per_altloc() {
    let rows = "A A TRP N A N 1 1 0 0 0 0.40\nA A TRP CA A C 1 1 1 0 0 0.40\nA A TRP CB A C 1 1 2 0 0 0.40\n\
                A A GLY N B N 1 1 0 1 0 0.60";
    let text = format!(
        "data_x\nloop_\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_comp_id\n\
         _atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.type_symbol\n_atom_site.label_seq_id\n\
         _atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n{rows}\n"
    );
    let structure = parse_mmcif(&text).unwrap();
    let selected = structure
        .first_model()
        .unwrap()
        .single_conformer(AltlocPolicy::HighestOccupancy);
    let comps: Vec<&str> = selected
        .chains()
        .next()
        .unwrap()
        .residues()
        .map(|r| r.comp().as_str())
        .collect();
    assert_eq!(comps, ["GLY"]);
}

#[test]
fn files_without_author_columns_do_not_invent_author_numbering() {
    let text = "data_x\nloop_\n_atom_site.label_asym_id\n_atom_site.label_comp_id\n_atom_site.label_atom_id\n\
                _atom_site.type_symbol\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n\
                _atom_site.Cartn_z\nA GLY N N 1 0 0 0\n";
    let structure = parse_mmcif(text).unwrap();
    let chain = &structure.first_model().unwrap().chains().next().unwrap();
    assert_eq!(chain.auth_asym().as_str(), "A");
    assert_eq!(
        structure.numbering().to_auth(&LabelResidueKey {
            asym: "A".into(),
            seq: LabelSeqId::new(1).unwrap()
        }),
        Err(NumberingError::NoAuthorNumbering(LabelResidueKey {
            asym: "A".into(),
            seq: LabelSeqId::new(1).unwrap()
        }))
    );
}

#[test]
fn null_author_values_are_errors_when_the_column_exists() {
    let error = parse_mmcif(&with_atom_site("A A GLY N . N 1 ? 0 0 0")).unwrap_err();
    assert!(matches!(
        error,
        ReadError::InvalidValue {
            item: "auth_seq_id",
            ..
        }
    ));
    let error = parse_mmcif(&with_atom_site("A ? GLY N . N 1 1 0 0 0")).unwrap_err();
    assert!(matches!(
        error,
        ReadError::InvalidValue {
            item: "auth_asym_id",
            ..
        }
    ));
}

#[test]
fn interleaved_rows_are_grouped_into_contiguous_residues() {
    let rows = "A A GLY N . N 1 1 1 0 0\nB B ALA N . N 1 1 2 0 0\nA A GLY CA . C 1 1 3 0 0\n\
                A A SER N . N 2 2 4 0 0\nB B ALA CA . C 1 1 5 0 0";
    let structure = parse_mmcif(&with_atom_site(rows)).unwrap();
    let model = structure.first_model().unwrap();
    let layout: Vec<String> = model
        .residues()
        .map(|r| {
            let atoms: Vec<String> = r
                .atoms()
                .map(|a| format!("{}@{}", a.name(), a.position()[0]))
                .collect();
            format!(
                "{}:{}:{}",
                r.chain().label_asym(),
                r.comp(),
                atoms.join(",")
            )
        })
        .collect();
    assert_eq!(layout, ["A:GLY:N@1,CA@3", "A:SER:N@4", "B:ALA:N@2,CA@5"]);
    assert_eq!(model.positions().len(), 5);
    assert_eq!(
        model.chain("B").unwrap().positions(),
        [[2.0, 0.0, 0.0], [5.0, 0.0, 0.0]]
    );
    let atom = model.atoms().nth(4).unwrap();
    assert_eq!(atom.residue().chain().label_asym().as_str(), "B");
}

fn assembly_document(rows: usize, expression: &str) -> String {
    let atoms: String = (0..rows).map(|_| "A A GLY N . N 1 1 0 0 0\n").collect();
    format!(
        "{}_pdbx_struct_assembly_gen.assembly_id 1\n_pdbx_struct_assembly_gen.oper_expression '{expression}'\n\
         _pdbx_struct_assembly_gen.asym_id_list A\nloop_\n_pdbx_struct_oper_list.id\n\
         _pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n\
         _pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n\
         _pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n\
         _pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.vector[3]\n\
         1 1 0 0 0 1 0 0 0 1 0 0 0\n",
        with_atom_site(atoms.trim_end())
    )
}

#[test]
fn undefined_operators_are_reported_before_expansion() {
    let structure = parse_mmcif(&assembly_document(10, "(1-1000000)")).unwrap();
    let error = structure
        .assemble(&AssemblyId::new("1"), structure.first_model().unwrap())
        .unwrap_err();
    assert!(
        matches!(error, AssemblyError::UnknownOperator { ref operator, .. } if operator.as_str() == "2")
    );
}

#[test]
fn oversized_assemblies_fail_before_allocating() {
    let ones = vec!["1"; 1024].join(",");
    let structure = parse_mmcif(&assembly_document(300, &format!("({ones})({ones})"))).unwrap();
    let error = structure
        .assemble(&AssemblyId::new("1"), structure.first_model().unwrap())
        .unwrap_err();
    assert!(matches!(error, AssemblyError::TooLarge { atoms, .. } if atoms == 300 * 1024 * 1024));
}

#[test]
fn derived_models_equal_fresh_parses_of_the_same_atoms() {
    let structure = structure();
    let model = structure.first_model().unwrap();
    let selected = model.single_conformer(AltlocPolicy::Prefer(AltLoc::new('B')));
    let rows = "AA H1 THR N B N 2 100 1.100 2.100 3.100\nAA H1 THR CA B C 2 100 2.100 2.100 3.100";
    let fresh = parse_mmcif(&with_atom_site(rows)).unwrap();
    let selected_threonine: Vec<(String, [f64; 3])> = selected
        .chain("AA")
        .unwrap()
        .residues()
        .next()
        .unwrap()
        .atoms()
        .map(|a| (a.name().to_string(), a.position()))
        .collect();
    let fresh_threonine: Vec<(String, [f64; 3])> = fresh
        .first_model()
        .unwrap()
        .atoms()
        .map(|a| (a.name().to_string(), a.position()))
        .collect();
    assert_eq!(selected_threonine, fresh_threonine);
    let identity = selected.single_conformer(AltlocPolicy::First);
    assert_eq!(identity, selected);
}

fn with_b_factors(values: [&str; 2]) -> String {
    format!(
        "data_x\nloop_\n_atom_site.label_asym_id\n_atom_site.label_comp_id\n_atom_site.label_atom_id\n\
         _atom_site.type_symbol\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n\
         _atom_site.Cartn_z\n_atom_site.B_iso_or_equiv\nA GLY N N 1 0 0 0 {}\nA GLY CA C 1 1.5 0 0 {}\n",
        values[0], values[1]
    )
}

#[test]
fn a_missing_b_factor_is_none_and_nan() {
    for missing in [".", "?"] {
        let structure = parse_mmcif(&with_b_factors([missing, "12.5"])).unwrap();
        let model = structure.first_model().unwrap();
        let b_factors: Vec<Option<f64>> = model.atoms().map(|atom| atom.b_factor()).collect();
        assert_eq!(b_factors, [None, Some(12.5)]);
        assert!(model.b_factors()[0].is_nan());
        assert_eq!(model.b_factors()[1], 12.5);
    }
}

#[test]
fn a_b_factor_that_is_not_finite_is_invalid() {
    for value in ["nan", "inf", "-inf"] {
        let error = parse_mmcif(&with_b_factors(["1.0", value])).unwrap_err();
        assert!(
            matches!(
                error,
                ReadError::InvalidValue {
                    category: "atom_site",
                    item: "B_iso_or_equiv",
                    row: 1,
                    ..
                }
            ),
            "{value}: {error}"
        );
    }
}

#[test]
fn structures_with_missing_b_factors_equal_themselves() {
    let text = with_b_factors([".", "?"]);
    let structure = parse_mmcif(&text).unwrap();
    assert_eq!(structure, structure.clone());
    assert_eq!(structure, parse_mmcif(&text).unwrap());
    assert_ne!(
        structure,
        parse_mmcif(&with_b_factors([".", "1.0"])).unwrap()
    );
}

use std::fmt::Write;

use mmcif_model::{AltlocPolicy, LabelResidueKey, LabelSeqId, OperatorId, parse_mmcif};
use proptest::prelude::*;

const HEADER: &str = "loop_\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_comp_id\n\
    _atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.type_symbol\n_atom_site.label_seq_id\n\
    _atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n\
    _atom_site.Cartn_z\n_atom_site.occupancy\n";

const IDENTITY_ASSEMBLY: &str = "_pdbx_struct_assembly_gen.assembly_id 1\n\
    _pdbx_struct_assembly_gen.oper_expression 1\n_pdbx_struct_assembly_gen.asym_id_list A,B\n\
    loop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n\
    _pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n\
    _pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n\
    _pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.vector[2]\n\
    _pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 1 0 0 0 1 0 0 0\n";

#[derive(Debug, Clone)]
struct GeneratedResidue {
    auth_seq: i32,
    ins_code: Option<char>,
    atoms: Vec<(Option<char>, [f64; 3], f64)>,
}

fn residues() -> impl Strategy<Value = Vec<GeneratedResidue>> {
    let atom = (
        prop::option::weighted(0.3, prop::sample::select(vec!['A', 'B', 'C'])),
        prop::array::uniform3(-500.0f64..500.0),
        0.0f64..1.0,
    );
    let residue = (
        1i32..4,
        prop::option::weighted(0.3, prop::sample::select(vec!['A', 'B'])),
        prop::collection::vec(atom, 1..5),
    );
    (-50i32..50, prop::collection::vec(residue, 1..30)).prop_map(|(start, specs)| {
        let mut seq = start;
        specs
            .into_iter()
            .map(|(gap, ins_code, atoms)| {
                if ins_code.is_none() {
                    seq += gap;
                }
                GeneratedResidue {
                    auth_seq: seq,
                    ins_code,
                    atoms,
                }
            })
            .collect::<Vec<_>>()
    })
}

fn dedup_positions(mut residues: Vec<GeneratedResidue>) -> Vec<GeneratedResidue> {
    let mut seen = std::collections::HashSet::new();
    residues.retain(|r| seen.insert((r.auth_seq, r.ins_code)));
    residues
}

fn document(residues: &[GeneratedResidue], extra: &str) -> String {
    let mut text = format!("data_gen\n{HEADER}");
    for (index, residue) in residues.iter().enumerate() {
        for (atom_index, (altloc, position, occupancy)) in residue.atoms.iter().enumerate() {
            writeln!(
                text,
                "A H GLY X{atom_index} {} C {} {} {} {} {} {} {}",
                altloc.map_or(".".to_owned(), String::from),
                index + 1,
                residue.auth_seq,
                residue.ins_code.map_or("?".to_owned(), String::from),
                position[0],
                position[1],
                position[2],
                occupancy
            )
            .unwrap();
        }
    }
    writeln!(text, "B H HOH O . O . 900 ? 0 0 0 1").unwrap();
    text.push_str(extra);
    text
}

proptest! {
    #[test]
    fn numbering_round_trips_for_every_polymer_residue(generated in residues().prop_map(dedup_positions)) {
        let structure = parse_mmcif(&document(&generated, "")).unwrap();
        let model = structure.first_model().unwrap();
        let chain = model.chain("A").unwrap();
        prop_assert_eq!(chain.residues().len(), generated.len());
        for (index, residue) in chain.residues().enumerate() {
            let label = LabelResidueKey { asym: "A".into(), seq: LabelSeqId::new(index as u32 + 1).unwrap() };
            let auth = structure.numbering().to_auth(&label).unwrap();
            prop_assert_eq!(&auth, &residue.auth_key());
            prop_assert_eq!(structure.numbering().to_label(&auth).unwrap(), label);
        }
    }

    #[test]
    fn identity_assembly_is_a_no_op(generated in residues().prop_map(dedup_positions)) {
        let structure = parse_mmcif(&document(&generated, IDENTITY_ASSEMBLY)).unwrap();
        let model = structure.first_model().unwrap();
        let assembled = structure.assemble("1", model).unwrap();
        prop_assert_eq!(assembled.chains().len(), model.chains().len());
        for (copy, original) in assembled.chains().zip(model.chains()) {
            let addresses = |chain: mmcif_model::ChainRef<'_>| {
                chain
                    .residues()
                    .map(|r| {
                        let key = r.key();
                        (key.asym, key.label_seq, key.auth_seq, key.ins_code, key.comp)
                    })
                    .collect::<Vec<_>>()
            };
            prop_assert_eq!(addresses(copy), addresses(original));
            prop_assert!(copy.residues().all(|r| r.key().operators.as_ref() == [OperatorId::new("1")]));
            prop_assert_eq!(copy.positions(), original.positions());
        }
    }

    #[test]
    fn single_conformer_keeps_one_altloc_per_position(generated in residues().prop_map(dedup_positions)) {
        let structure = parse_mmcif(&document(&generated, "")).unwrap();
        let model = structure.first_model().unwrap();
        for policy in [AltlocPolicy::First, AltlocPolicy::HighestOccupancy] {
            let selected = model.single_conformer(policy);
            prop_assert!(selected.atom_count() <= model.atom_count());
            for residue in selected.chains().flat_map(|c| c.residues()) {
                prop_assert!(residue.altlocs().len() <= 1);
            }
            let untouched = |chain: mmcif_model::ChainRef<'_>| chain.residues().filter(|r| r.altlocs().is_empty()).count();
            prop_assert_eq!(
                untouched(selected.chain("A").unwrap()),
                untouched(model.chain("A").unwrap())
            );
        }
    }
}

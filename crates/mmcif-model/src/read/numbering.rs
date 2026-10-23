use mmcif_cif::Block;

use crate::error::ReadError;
use crate::ids::{AuthAsymId, AuthSeqId, InsCode, LabelAsymId, LabelSeqId};
use crate::key::{AuthResidueKey, LabelResidueKey};
use crate::model::Model;
use crate::numbering::Numbering;
use crate::read::Skipping;
use crate::read::table::Table;

pub(crate) fn read_numbering(
    block: &Block<'_>,
    models: &[Model],
    has_author_numbering: bool,
    skipping: &mut Skipping,
) -> Result<Numbering, ReadError> {
    let mut numbering = Numbering::default();
    for (label, auth) in skipping.or_default(read_polymer_scheme(block))? {
        numbering.insert_polymer(label, auth);
    }
    for (asym, num, auth) in skipping.or_default(read_branch_scheme(block))? {
        numbering.insert_branch(asym, num, auth);
    }
    if let Some(model) = models.first() {
        for chain in model.chains() {
            if !has_author_numbering {
                numbering.insert_label_only(chain.label_asym().clone());
                continue;
            }
            for residue in chain.residues() {
                let auth = residue.auth_key();
                match residue.label_seq() {
                    Some(seq) => numbering.insert_polymer(
                        LabelResidueKey {
                            asym: chain.label_asym().clone(),
                            seq,
                        },
                        auth,
                    ),
                    None => numbering.insert_non_polymer(chain.label_asym().clone(), auth),
                }
            }
        }
    }
    Ok(numbering)
}

fn read_polymer_scheme(
    block: &Block<'_>,
) -> Result<Vec<(LabelResidueKey, AuthResidueKey)>, ReadError> {
    let mut pairs = Vec::new();
    if let Some(table) = Table::optional(block, "pdbx_poly_seq_scheme") {
        let asyms = table.required_column("asym_id")?;
        let seqs = table.required_column("seq_id")?;
        let auth_seqs = table.required_column("pdb_seq_num")?;
        let ins_codes = table.column("pdb_ins_code");
        let strands = table.required_column("pdb_strand_id")?;
        for row in 0..table.rows() {
            let seq: u32 = seqs.required_parse(row)?;
            let label = LabelResidueKey {
                asym: LabelAsymId::new(asyms.required_text(row)?),
                seq: LabelSeqId::new(seq)
                    .ok_or_else(|| seqs.invalid(row, "0", "seq_id starts at 1"))?,
            };
            let auth = AuthResidueKey {
                chain: AuthAsymId::new(strands.required_text(row)?),
                seq: AuthSeqId::new(auth_seqs.required_parse(row)?),
                ins_code: ins_codes
                    .single_char(row)?
                    .map_or(InsCode::NONE, InsCode::new),
            };
            pairs.push((label, auth));
        }
    }
    Ok(pairs)
}

fn read_branch_scheme(
    block: &Block<'_>,
) -> Result<Vec<(LabelAsymId, u32, AuthResidueKey)>, ReadError> {
    let mut entries = Vec::new();
    if let Some(table) = Table::optional(block, "pdbx_branch_scheme") {
        let asyms = table.required_column("asym_id")?;
        let nums = table.required_column("num")?;
        let auth_seqs = table.required_column("pdb_seq_num")?;
        let strands = table.required_column("pdb_asym_id")?;
        for row in 0..table.rows() {
            let auth = AuthResidueKey {
                chain: AuthAsymId::new(strands.required_text(row)?),
                seq: AuthSeqId::new(auth_seqs.required_parse(row)?),
                ins_code: InsCode::NONE,
            };
            entries.push((
                LabelAsymId::new(asyms.required_text(row)?),
                nums.required_parse(row)?,
                auth,
            ));
        }
    }
    Ok(entries)
}

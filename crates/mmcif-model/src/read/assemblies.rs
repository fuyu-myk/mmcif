use mmcif_cif::Block;

use crate::assembly::{Assembly, Generator, Operator, Transform};
use crate::error::ReadError;
use crate::ids::{AssemblyId, LabelAsymId, OperatorId};
use crate::read::table::Table;

const MATRIX_ITEMS: [[&str; 3]; 3] = [
    ["matrix[1][1]", "matrix[1][2]", "matrix[1][3]"],
    ["matrix[2][1]", "matrix[2][2]", "matrix[2][3]"],
    ["matrix[3][1]", "matrix[3][2]", "matrix[3][3]"],
];
const VECTOR_ITEMS: [&str; 3] = ["vector[1]", "vector[2]", "vector[3]"];

pub(crate) fn read_operators(block: &Block<'_>) -> Result<Vec<Operator>, ReadError> {
    let Some(table) = Table::optional(block, "pdbx_struct_oper_list") else {
        return Ok(Vec::new());
    };
    let ids = table.required_column("id")?;
    let kinds = table.column("type");
    let matrix = MATRIX_ITEMS
        .iter()
        .map(|row| row.iter().map(|item| table.required_column(item)).collect())
        .collect::<Result<Vec<Vec<_>>, _>>()?;
    let vector = VECTOR_ITEMS
        .iter()
        .map(|item| table.required_column(item))
        .collect::<Result<Vec<_>, _>>()?;
    (0..table.rows())
        .map(|row| {
            let mut transform = Transform::IDENTITY;
            for (i, fields) in matrix.iter().enumerate() {
                for (j, field) in fields.iter().enumerate() {
                    transform.rotation[i][j] = field.required_parse(row)?;
                }
                transform.translation[i] = vector[i].required_parse(row)?;
            }
            Ok(Operator {
                id: OperatorId::new(ids.required_text(row)?),
                kind: kinds.text(row).map(str::to_owned),
                transform,
            })
        })
        .collect()
}

pub(crate) fn read_assemblies(block: &Block<'_>) -> Result<Vec<Assembly>, ReadError> {
    let mut assemblies: Vec<Assembly> = Vec::new();
    if let Some(table) = Table::optional(block, "pdbx_struct_assembly") {
        let ids = table.required_column("id")?;
        let details = table.column("details");
        for row in 0..table.rows() {
            assemblies.push(Assembly {
                id: AssemblyId::new(ids.required_text(row)?),
                details: details.text(row).map(str::to_owned),
                generators: Vec::new(),
            });
        }
    }
    let Some(table) = Table::optional(block, "pdbx_struct_assembly_gen") else {
        return Ok(assemblies);
    };
    let ids = table.required_column("assembly_id")?;
    let expressions = table.required_column("oper_expression")?;
    let asym_lists = table.required_column("asym_id_list")?;
    for row in 0..table.rows() {
        let id = ids.required_text(row)?;
        let asym_ids = asym_lists
            .required_text(row)?
            .split(',')
            .map(str::trim)
            .filter(|asym| !asym.is_empty())
            .map(LabelAsymId::new)
            .collect();
        let generator = Generator::new(expressions.required_text(row)?, asym_ids)?;
        let index = match assemblies
            .iter()
            .position(|assembly| assembly.id.as_str() == id)
        {
            Some(index) => index,
            None => {
                assemblies.push(Assembly::new(AssemblyId::new(id), Vec::new()));
                assemblies.len() - 1
            }
        };
        assemblies[index].generators.push(generator);
    }
    Ok(assemblies)
}

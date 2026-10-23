use crate::assembly::{Assembly, Generator, Operator, Transform};
use crate::error::ReadError;
use crate::ids::{AssemblyId, LabelAsymId, OperatorId};
use crate::model::Model;
use crate::pdb::Line;

const ROWS: usize = 3;

pub(crate) struct Biomolecule {
    id: String,
    author: Option<String>,
    software: Option<String>,
    groups: Vec<Group>,
}

struct Group {
    chains: Vec<String>,
    transforms: Vec<Rows>,
}

struct Rows {
    serial: String,
    transform: Transform,
    filled: usize,
}

pub(crate) fn parse(lines: &[Line<'_>]) -> Result<Vec<Biomolecule>, ReadError> {
    let mut biomolecules: Vec<Biomolecule> = Vec::new();
    for line in lines {
        let content = line.columns(10, 80)?.trim();
        if let Some(id) = content.strip_prefix("BIOMOLECULE:") {
            close_group(&biomolecules, line)?;
            biomolecules.push(Biomolecule {
                id: id.trim().to_owned(),
                author: None,
                software: None,
                groups: Vec::new(),
            });
        } else if let Some(unit) = content.strip_prefix("AUTHOR DETERMINED BIOLOGICAL UNIT:") {
            current(&mut biomolecules, line)?.author = Some(unit.trim().to_ascii_lowercase());
        } else if let Some(unit) = content.strip_prefix("SOFTWARE DETERMINED QUATERNARY STRUCTURE:")
        {
            current(&mut biomolecules, line)?.software = Some(unit.trim().to_ascii_lowercase());
        } else if let Some(chains) = content.strip_prefix("APPLY THE FOLLOWING TO CHAINS:") {
            close_group(&biomolecules, line)?;
            current(&mut biomolecules, line)?.groups.push(Group {
                chains: chain_list(chains),
                transforms: Vec::new(),
            });
        } else if let Some(chains) = content.strip_prefix("AND CHAINS:") {
            group(&mut biomolecules, line)?
                .chains
                .extend(chain_list(chains));
        } else if content.starts_with("BIOMT") {
            read_biomt(group(&mut biomolecules, line)?, line, content)?;
        }
    }
    if let Some(last) = lines.last() {
        close_group(&biomolecules, last)?;
    }
    Ok(biomolecules)
}

fn current<'b>(
    biomolecules: &'b mut [Biomolecule],
    line: &Line<'_>,
) -> Result<&'b mut Biomolecule, ReadError> {
    biomolecules
        .last_mut()
        .ok_or_else(|| line.error("REMARK 350 record before any BIOMOLECULE"))
}

fn group<'b>(
    biomolecules: &'b mut [Biomolecule],
    line: &Line<'_>,
) -> Result<&'b mut Group, ReadError> {
    current(biomolecules, line)?
        .groups
        .last_mut()
        .ok_or_else(|| line.error("BIOMT or AND CHAINS before APPLY THE FOLLOWING TO CHAINS"))
}

fn close_group(biomolecules: &[Biomolecule], line: &Line<'_>) -> Result<(), ReadError> {
    let Some(group) = biomolecules
        .last()
        .and_then(|biomolecule| biomolecule.groups.last())
    else {
        return Ok(());
    };
    match group.transforms.last() {
        None => Err(line.error("APPLY THE FOLLOWING TO CHAINS has no BIOMT transformation")),
        Some(rows) if rows.filled < ROWS => Err(line.error("BIOMT transformation lacks a row")),
        Some(_) => Ok(()),
    }
}

fn chain_list(text: &str) -> Vec<String> {
    text.split(',')
        .map(str::trim)
        .filter(|chain| !chain.is_empty())
        .map(str::to_owned)
        .collect()
}

fn read_biomt(group: &mut Group, line: &Line<'_>, content: &str) -> Result<(), ReadError> {
    let invalid = || line.error("invalid BIOMT record");
    let fields: Vec<&str> = content.split_whitespace().collect();
    let [label, serial, values @ ..] = fields.as_slice() else {
        return Err(invalid());
    };
    let row = match label.strip_prefix("BIOMT") {
        Some("1") => 0,
        Some("2") => 1,
        Some("3") => 2,
        _ => return Err(invalid()),
    };
    let numbers = values
        .iter()
        .map(|value| value.parse::<f64>().map_err(|_| invalid()))
        .collect::<Result<Vec<f64>, ReadError>>()?;
    let [x, y, z, shift] = numbers.as_slice() else {
        return Err(invalid());
    };
    if row == 0 {
        if group
            .transforms
            .last()
            .is_some_and(|rows| rows.filled < ROWS)
        {
            return Err(line.error("BIOMT transformation lacks a row"));
        }
        group.transforms.push(Rows {
            serial: (*serial).to_owned(),
            transform: Transform::IDENTITY,
            filled: 0,
        });
    }
    let rows = group
        .transforms
        .last_mut()
        .filter(|rows| rows.filled == row && rows.serial == *serial)
        .ok_or_else(|| line.error("BIOMT rows out of order"))?;
    rows.transform.rotation[row] = [*x, *y, *z];
    rows.transform.translation[row] = *shift;
    rows.filled += 1;
    Ok(())
}

pub(crate) fn build(
    biomolecules: Vec<Biomolecule>,
    model: Option<&Model>,
) -> (Vec<Operator>, Vec<Assembly>) {
    let mut operators: Vec<Operator> = Vec::new();
    let assemblies = biomolecules
        .into_iter()
        .map(|biomolecule| {
            let generators = biomolecule
                .groups
                .iter()
                .map(|group| {
                    let ids: Vec<String> = group
                        .transforms
                        .iter()
                        .map(|rows| operator_id(&mut operators, &rows.transform))
                        .collect();
                    Generator::new(&ids.join(","), label_chains(model, &group.chains))
                        .expect("operator ids joined by commas form a valid expression")
                })
                .collect();
            Assembly {
                details: details(&biomolecule),
                id: AssemblyId::new(biomolecule.id),
                oligomeric_details: biomolecule.author.or(biomolecule.software),
                oligomeric_count: None,
                generators,
            }
        })
        .collect();
    (operators, assemblies)
}

fn operator_id(operators: &mut Vec<Operator>, transform: &Transform) -> String {
    if let Some(existing) = operators
        .iter()
        .find(|operator| operator.transform == *transform)
    {
        return existing.id.to_string();
    }
    let id = (operators.len() + 1).to_string();
    operators.push(Operator {
        id: OperatorId::new(id.as_str()),
        kind: None,
        transform: *transform,
    });
    id
}

fn label_chains(model: Option<&Model>, authors: &[String]) -> Vec<LabelAsymId> {
    authors
        .iter()
        .flat_map(|author| {
            model
                .into_iter()
                .flat_map(|model| model.chains())
                .filter(move |chain| chain.auth_asym().as_str() == author)
                .map(|chain| chain.label_asym().clone())
        })
        .collect()
}

fn details(biomolecule: &Biomolecule) -> Option<String> {
    let details = match (&biomolecule.author, &biomolecule.software) {
        (Some(_), Some(_)) => "author_and_software_defined_assembly",
        (Some(_), None) => "author_defined_assembly",
        (None, Some(_)) => "software_defined_assembly",
        (None, None) => return None,
    };
    Some(details.to_owned())
}

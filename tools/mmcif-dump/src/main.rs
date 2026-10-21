use std::process::ExitCode;

use mmcif_model::{Model, Structure};
use serde_json::{Value, json};

fn round3(value: f64) -> f64 {
    (value * 1000.0).round() / 1000.0
}

fn code(value: Option<char>) -> String {
    value.map(String::from).unwrap_or_default()
}

fn dump_model(model: &Model) -> Value {
    let atoms: Vec<Value> = model
        .atoms()
        .map(|atom| {
            let residue = atom.residue();
            let chain = residue.chain();
            let [x, y, z] = atom.position();
            json!([
                chain.auth_asym().as_str(),
                chain.label_asym().as_str(),
                residue.auth_seq().get(),
                code(residue.ins_code().get()),
                residue.comp().as_str(),
                atom.name().as_str(),
                code(atom.altloc().get()),
                atom.element().symbol().to_ascii_uppercase(),
                round3(x),
                round3(y),
                round3(z),
                round3(atom.occupancy()),
                round3(atom.b_factor()),
            ])
        })
        .collect();
    let residues: Vec<Value> = model
        .chains()
        .flat_map(|chain| {
            chain.residues().map(move |residue| {
                json!([
                    chain.auth_asym().as_str(),
                    chain.label_asym().as_str(),
                    residue.auth_seq().get(),
                    code(residue.ins_code().get()),
                    residue.comp().as_str(),
                    residue.label_seq().map(|seq| seq.get()),
                ])
            })
        })
        .collect();
    json!({ "number": model.number().get(), "atoms": atoms, "residues": residues })
}

fn dump_assemblies(structure: &Structure) -> Result<Vec<Value>, String> {
    let Some(model) = structure.first_model() else {
        return Ok(Vec::new());
    };
    structure
        .assemblies()
        .iter()
        .map(|assembly| {
            let assembled = structure
                .assemble(assembly.id(), model)
                .map_err(|error| error.to_string())?;
            let copies: Vec<Value> = assembled
                .chains()
                .filter(|chain| chain.atom_count() > 0)
                .map(|chain| {
                    let positions = chain.positions();
                    let count = positions.len() as f64;
                    let centroid = [0, 1, 2]
                        .map(|axis| round3(positions.iter().map(|p| p[axis]).sum::<f64>() / count));
                    json!({ "atoms": positions.len(), "centroid": centroid })
                })
                .collect();
            Ok(json!({
                "id": assembly.id().as_str(),
                "atom_count": assembled.atom_count(),
                "copies": copies,
            }))
        })
        .collect()
}

fn dump(structure: &Structure) -> Result<Value, String> {
    let entities: Vec<Value> = structure
        .entities()
        .iter()
        .map(|entity| {
            let mut subchains: Vec<&str> =
                entity.asym_ids().iter().map(|asym| asym.as_str()).collect();
            subchains.sort_unstable();
            json!({
                "id": entity.id().as_str(),
                "type": entity.kind().as_mmcif_str(),
                "polymer_type": entity.polymer().map(|polymer| polymer.kind().as_mmcif_str()),
                "subchains": subchains,
            })
        })
        .collect();
    Ok(json!({
        "entry": structure.entry_id(),
        "models": structure.models().iter().map(|model| dump_model(model)).collect::<Vec<_>>(),
        "entities": entities,
        "assemblies": dump_assemblies(structure)?,
    }))
}

fn main() -> ExitCode {
    let Some(path) = std::env::args().nth(1) else {
        eprintln!("usage: mmcif-dump <structure file>");
        return ExitCode::from(2);
    };
    let result = mmcif_model::read(&path)
        .map_err(|error| error.to_string())
        .and_then(|structure| dump(&structure));
    match result {
        Ok(value) => {
            println!("{value}");
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("{path}: {error}");
            ExitCode::FAILURE
        }
    }
}

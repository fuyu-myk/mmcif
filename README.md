# mmcif

Correctness-first mmCIF, BinaryCIF and PDB structure IO with typed numbering, for Rust and Python.

Every residue in an mmCIF file has two addresses: the author numbering people read off papers and viewers, and the label numbering the file defines. They often disagree, and tools that keep only one put constraints and selections on the wrong residue without saying so. This project keeps both as separate types and converts between them only when you ask.

## Example

Heavy chain H of antibody 12e8 numbers its CDR residues with insertion codes. Author residue H 100A is label residue B 105: the chain, the number and the insertion code all differ.

```bash
curl -O https://files.rcsb.org/download/12e8.cif
```

```python
import mmcif_io

structure = mmcif_io.read("12e8.cif")
print(structure.to_label("H", 100, "A"))  # ('B', 105)
```

```rust
use mmcif_model::{AuthResidueKey, AuthSeqId, InsCode};

let structure = mmcif_model::read("12e8.cif")?;
let key = AuthResidueKey { chain: "H".into(), seq: AuthSeqId::new(100), ins_code: InsCode::new('A') };
assert_eq!(structure.numbering().to_label(&key)?.to_string(), "B/105");
```

## Install

```bash
pip install mmcif-io
```

```bash
cargo add mmcif-model
```

```bash
cargo install mmcif-cli
```

`mmcif-io` is the Python package, imported as `mmcif_io`. `mmcif-model` is the Rust entry point. `mmcif-cli` installs the `mmcif` command. The Rust crates and the Python package share one version number and release together.

## How it differs

| | mmcif | gemmi 0.7.5 | pdbtbx 0.12.0 |
| --- | --- | --- | --- |
| Reads BinaryCIF | Yes | No | No |
| Label and author numbering | Separate types, explicit conversion | Both stored on each residue | Author only, label used when author is missing |
| Biological assemblies | Expanded, including composite operator expressions | Expanded, but skips assemblies with unrecognised details and expands only the first group of a composite expression | Not expanded |
| Languages | Rust, Python | C++, Python | Rust |

Checked on 21 October 2026 against gemmi 0.7.5 and pdbtbx 0.12.0. gemmi is the broader toolkit. This project competes on numbering, BinaryCIF and assembly correctness, not on breadth.

## Correctness

Every push to main and every pull request runs the same entries through this project and gemmi, comparing atom counts, coordinates, occupancies, B-factors, residue keys, label numbering, entities and expanded assemblies. Any difference without a written explanation fails the build.

The corpus covers insertion codes, alternate conformations, microheterogeneity, numbering divergence, assembly operators, branched entities, multi-character chain ids and multiple models, plus a baseline set: 44 entries today. All 44 agree with gemmi except in assembly expansion on two entries, where 4,872 differences are each explained by a rule in [`harness/explained.toml`](https://github.com/fuyu-myk/mmcif/blob/main/harness/explained.toml). To reproduce:

```bash
pip install -r harness/requirements.txt
python corpus/fetch.py
python harness/run.py
```

For every corpus entry, writing a parsed document and parsing it again gives the same document. The CIF tokenizer, the BinaryCIF decoder, the PDB reader and structure building are fuzzed in CI.

## Design specifications

`mmcif check` reads a [BoltzGen](https://github.com/HannesStark/boltzgen) design specification and resolves every residue reference against the structure it names. BoltzGen addresses chains by label id, so writing the author id is an easy mistake. In 12e8, label chain H is a water chain:

```
error[no-polymer-residues]: chain H in 12e8.cif has no polymer residues
 --> binder.yaml:6:17
  |
6 |             id: H
  |                 ^ it holds HOH
  = help: author chain H is label chain B; polymer chains in the file: A, B, C, D
```

Where an index in a specification selects a different residue under author numbering, it warns and names the author residue the index selects. It also checks that every sampled length is satisfiable and that bonds and distances are geometrically possible against the coordinates in the file.

`mmcif lower spec.yaml --to boltzgen` and `mmcif lower spec.yaml --to rfd3` write a checked specification out for BoltzGen or [RFdiffusion3](https://github.com/RosettaCommons/foundry). Where a backend cannot express part of a specification, such as a disulfide bond for RFdiffusion3, lowering refuses and points at the line instead of dropping it. These types and formats may change before 1.0.

## Crates

| Crate | |
| --- | --- |
| [`mmcif-model`](https://crates.io/crates/mmcif-model) | Typed structure model: hierarchy, numbering, altlocs, entities, assemblies. Start here. |
| [`mmcif-cif`](https://crates.io/crates/mmcif-cif) | Lossless CIF tokenizer, table layer and BinaryCIF decoder |
| [`mmcif-ccd`](https://crates.io/crates/mmcif-ccd) | Compiled, versioned index of the Chemical Component Dictionary |
| [`mmcif-geom`](https://crates.io/crates/mmcif-geom) | Neighbour grid queries over atomic coordinates |
| [`mmcif-io`](https://pypi.org/project/mmcif-io/) | Python bindings for all of the above |
| [`mmcif-cli`](https://crates.io/crates/mmcif-cli) | The `mmcif` command |
| [`mmcif-ir`](https://crates.io/crates/mmcif-ir) | Design specification model: slots, determinacy, constraints |
| [`mmcif-spec`](https://crates.io/crates/mmcif-spec) | Design specification parser and checker |
| [`mmcif-lower`](https://crates.io/crates/mmcif-lower) | Lowering to BoltzGen and RFdiffusion3 |

## License

MIT or Apache-2.0, at your option.

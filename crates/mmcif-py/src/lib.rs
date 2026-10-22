mod ccd;
mod document;
mod errors;
mod geom;
mod hierarchy;
mod structure;

use pyo3::prelude::*;

/// Read mmCIF, BinaryCIF and PDB files into a typed structure model.
///
/// ``read`` returns a ``Structure`` holding models, chains, residues and
/// atoms, with coordinates as NumPy arrays. ``Document`` gives the CIF
/// tables of a file, ``NeighbourGrid`` answers distance queries over
/// coordinates, and ``CcdIndex`` reads a compiled Chemical Component
/// Dictionary.
///
/// The examples assume ``import mmcif_io``.
///
/// Examples
/// --------
/// >>> structure = mmcif_io.read("1crn.bcif")
/// >>> model = structure.models[0]
/// >>> model.atom_count
/// 327
/// >>> model.coordinates.shape
/// (327, 3)
#[pymodule]
fn mmcif_io(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    module.add("ReadError", py.get_type::<errors::ReadError>())?;
    module.add("NumberingError", py.get_type::<errors::NumberingError>())?;
    module.add("AssemblyError", py.get_type::<errors::AssemblyError>())?;
    module.add_class::<document::PyDocument>()?;
    module.add_class::<document::PyBlock>()?;
    module.add_class::<document::PyCategory>()?;
    module.add_class::<structure::PyStructure>()?;
    module.add_class::<structure::PyEntity>()?;
    module.add_class::<structure::PyAssembly>()?;
    module.add_class::<hierarchy::PyModel>()?;
    module.add_class::<hierarchy::PyChain>()?;
    module.add_class::<hierarchy::PyResidue>()?;
    module.add_class::<hierarchy::PyAtom>()?;
    module.add_class::<geom::PyNeighbourGrid>()?;
    module.add_class::<ccd::PyCcdIndex>()?;
    module.add_class::<ccd::PyComponent>()?;
    module.add_function(wrap_pyfunction!(structure::read, module)?)?;
    module.add_function(wrap_pyfunction!(structure::read_mmcif, module)?)?;
    module.add_function(wrap_pyfunction!(structure::read_bcif, module)?)?;
    module.add_function(wrap_pyfunction!(structure::read_pdb, module)?)?;
    module.add_function(wrap_pyfunction!(structure::parse_mmcif, module)?)?;
    module.add_function(wrap_pyfunction!(structure::parse_bcif, module)?)?;
    module.add_function(wrap_pyfunction!(structure::parse_pdb, module)?)?;
    Ok(())
}

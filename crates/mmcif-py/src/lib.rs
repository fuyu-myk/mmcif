mod ccd;
mod document;
mod errors;
mod hierarchy;
mod structure;

use pyo3::prelude::*;

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

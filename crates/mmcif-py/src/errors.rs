use std::fmt::Display;

use pyo3::PyErr;
use pyo3::exceptions::PyValueError;

pyo3::create_exception!(mmcif_io, ReadError, PyValueError);
pyo3::create_exception!(mmcif_io, NumberingError, PyValueError);
pyo3::create_exception!(mmcif_io, AssemblyError, PyValueError);

pub fn read_error(error: impl Display) -> PyErr {
    ReadError::new_err(error.to_string())
}

pub fn numbering_error(error: impl Display) -> PyErr {
    NumberingError::new_err(error.to_string())
}

pub fn assembly_error(error: impl Display) -> PyErr {
    AssemblyError::new_err(error.to_string())
}

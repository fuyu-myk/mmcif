use std::fmt::Display;
use std::io;
use std::path::Path;

use pyo3::PyErr;
use pyo3::exceptions::{
    PyFileNotFoundError, PyIsADirectoryError, PyOSError, PyPermissionError, PyValueError,
};

pyo3::create_exception!(
    mmcif_io,
    ReadError,
    PyValueError,
    "Raised when a file, text or bytes cannot be read as mmCIF, BinaryCIF, PDB\nor a component index.\n\nA subclass of ``ValueError``."
);
pyo3::create_exception!(
    mmcif_io,
    NumberingError,
    PyValueError,
    "Raised when a residue address has no counterpart in the other numbering\nscheme.\n\nA subclass of ``ValueError``."
);
pyo3::create_exception!(
    mmcif_io,
    AssemblyError,
    PyValueError,
    "Raised when a biological assembly cannot be built.\n\nA subclass of ``ValueError``."
);

pub fn read_error(error: impl Display) -> PyErr {
    ReadError::new_err(error.to_string())
}

pub fn numbering_error(error: impl Display) -> PyErr {
    NumberingError::new_err(error.to_string())
}

pub fn assembly_error(error: impl Display) -> PyErr {
    AssemblyError::new_err(error.to_string())
}

pub fn os_error(error: &io::Error, path: &Path) -> PyErr {
    let arguments = (
        error.raw_os_error().unwrap_or(0),
        error.to_string(),
        path.display().to_string(),
    );
    match error.kind() {
        io::ErrorKind::NotFound => PyFileNotFoundError::new_err(arguments),
        io::ErrorKind::PermissionDenied => PyPermissionError::new_err(arguments),
        io::ErrorKind::IsADirectory => PyIsADirectoryError::new_err(arguments),
        _ => PyOSError::new_err(arguments),
    }
}

pub fn cif_error(error: mmcif_cif::Error, path: &Path) -> PyErr {
    match error {
        mmcif_cif::Error::Io(io) => os_error(&io, path),
        other => read_error(other),
    }
}

pub fn model_error(error: mmcif_model::ReadError, path: &Path) -> PyErr {
    match error {
        mmcif_model::ReadError::Cif(cif) => cif_error(cif, path),
        other => read_error(other),
    }
}

pub fn ccd_error(error: mmcif_ccd::CcdError, path: &Path) -> PyErr {
    match error {
        mmcif_ccd::CcdError::Io(io) => os_error(&io, path),
        other => read_error(other),
    }
}

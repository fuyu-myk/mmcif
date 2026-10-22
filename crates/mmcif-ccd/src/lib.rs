//! A compiled, versioned index of the Chemical Component Dictionary.
//!
//! [`build_index`] and [`IndexBuilder`] compile dictionary entries into an
//! index for one dictionary release, and [`CcdIndex`] reads that index and
//! returns each [`Component`] with its atoms and bonds by id.
//!
//! # Examples
//!
//! ```
//! use mmcif_ccd::{BondOrder, CcdIndex, build_index};
//!
//! let dictionary = include_str!("../tests/fixtures/components-sample.cif");
//! let mut bytes = Vec::new();
//! let summary = build_index(dictionary, "2026-10-20", &mut bytes)?;
//! assert_eq!(summary.components(), 8);
//!
//! let index = CcdIndex::from_bytes(bytes)?;
//! assert_eq!(index.release(), "2026-10-20");
//!
//! let alanine = index.get("ALA")?.unwrap();
//! assert_eq!(alanine.name(), "ALANINE");
//! assert_eq!(alanine.atoms().len(), 13);
//! assert_eq!(alanine.bonds()[6].order(), BondOrder::Double);
//! # Ok::<(), mmcif_ccd::CcdError>(())
//! ```

#![warn(missing_docs)]

mod block;
mod build;
mod codec;
mod component;
mod error;
mod format;
mod index;
mod record;

pub use build::{BuildSummary, IndexBuilder, build_index, build_index_from_file};
pub use component::{AtomStereo, BondOrder, BondStereo, Component, ComponentAtom, ComponentBond};
pub use error::CcdError;
pub use index::CcdIndex;

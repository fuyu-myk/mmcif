//! Distance queries over points in three dimensions.
//!
//! A [`NeighbourGrid`] holds a set of points and answers radius,
//! nearest-point, pair and contact queries against them. Points are
//! identified by their position in the slice the grid was built from.
//!
//! # Examples
//!
//! ```
//! use mmcif_geom::NeighbourGrid;
//!
//! let points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
//! let grid = NeighbourGrid::new(&points, 4.0)?;
//!
//! assert_eq!(grid.within([0.0, 0.0, 0.0], 1.5)?, [0, 1]);
//! assert_eq!(grid.pairs_within(1.5)?, [(0, 1)]);
//! assert_eq!(grid.nearest([4.0, 0.0, 0.0], 10.0)?, Some((2, 1.0)));
//! # Ok::<(), mmcif_geom::GeomError>(())
//! ```

#![warn(missing_docs)]

mod cell;
mod check;
mod error;
mod grid;
mod index;

pub use error::GeomError;
pub use grid::NeighbourGrid;

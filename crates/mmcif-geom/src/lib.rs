//! Distance queries over points in three dimensions.
//!
//! A [`NeighbourGrid`] holds a set of points and answers radius,
//! nearest-point, pair and contact queries against them, each within a
//! [`Radius`]. Points are identified by their position in the order the grid
//! was built from.
//!
//! # Examples
//!
//! ```
//! use mmcif_geom::{NeighbourGrid, Radius};
//!
//! let points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
//! let grid = NeighbourGrid::new(&points, 4.0)?;
//! let radius = Radius::new(1.5)?;
//!
//! assert_eq!(grid.within([0.0, 0.0, 0.0], radius), [0, 1]);
//! assert_eq!(grid.pairs_within(radius), [(0, 1)]);
//! assert_eq!(grid.nearest([4.0, 0.0, 0.0], Radius::INFINITE), Some((2, 1.0)));
//! # Ok::<(), mmcif_geom::GeomError>(())
//! ```

#![warn(missing_docs)]

mod cell;
mod check;
mod error;
mod grid;
mod index;
mod radius;

pub use error::GeomError;
pub use grid::NeighbourGrid;
pub use radius::Radius;

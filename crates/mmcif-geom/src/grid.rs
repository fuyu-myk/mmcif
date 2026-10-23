mod query;

use std::borrow::Borrow;

use crate::cell::cell_key;
use crate::check;
use crate::error::GeomError;
use crate::index::CellIndex;

/// A set of points indexed by cubic cells for distance queries.
///
/// Points are identified by their position in the order given to
/// [`NeighbourGrid::new`]. Distances are Euclidean and every radius is
/// inclusive.
///
/// # Examples
///
/// ```
/// use mmcif_geom::{NeighbourGrid, Radius};
///
/// let grid = NeighbourGrid::new([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], 4.0)?;
///
/// assert_eq!(grid.len(), 2);
/// assert_eq!(grid.within([0.0, 0.0, 0.0], Radius::new(2.0)?), [0, 1]);
/// # Ok::<(), mmcif_geom::GeomError>(())
/// ```
#[derive(Debug, Clone)]
pub struct NeighbourGrid {
    cell_size: f64,
    sorted_points: Vec<[f64; 3]>,
    sorted_indices: Vec<u32>,
    index: CellIndex,
}

impl NeighbourGrid {
    /// Builds a grid over `points` with cubic cells of edge length `cell_size`.
    ///
    /// `points` is any sequence of points or references to them, such as a
    /// slice or an iterator. The grid keeps one copy of the points, in its own
    /// order.
    ///
    /// # Errors
    ///
    /// Returns [`GeomError::InvalidCellSize`] if `cell_size` is not positive
    /// and finite, [`GeomError::NonFiniteCoordinate`] if a point has a NaN or
    /// infinite coordinate, and [`GeomError::TooManyPoints`] if there are more
    /// than `u32::MAX` points.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_geom::{GeomError, NeighbourGrid};
    ///
    /// let grid = NeighbourGrid::new([[0.0, 0.0, 0.0]], 4.0)?;
    /// assert_eq!(grid.cell_size(), 4.0);
    ///
    /// assert_eq!(
    ///     NeighbourGrid::new([[0.0, 0.0, 0.0]], 0.0).unwrap_err(),
    ///     GeomError::InvalidCellSize(0.0)
    /// );
    /// # Ok::<(), GeomError>(())
    /// ```
    pub fn new<P>(points: impl IntoIterator<Item = P>, cell_size: f64) -> Result<Self, GeomError>
    where
        P: Borrow<[f64; 3]>,
    {
        let cell_size = check::cell_size(cell_size)?;
        let points = points
            .into_iter()
            .enumerate()
            .map(|(index, point)| {
                let point = *point.borrow();
                if check::is_finite(point) {
                    Ok(point)
                } else {
                    Err(GeomError::NonFiniteCoordinate { index })
                }
            })
            .collect::<Result<Vec<[f64; 3]>, GeomError>>()?;
        check::point_count(points.len())?;
        let mut keyed: Vec<([i64; 3], u32)> = points
            .iter()
            .zip(0u32..)
            .map(|(&point, index)| (cell_key(point, cell_size), index))
            .collect();
        keyed.sort_unstable();
        let sorted_indices: Vec<u32> = keyed.iter().map(|&(_, index)| index).collect();
        let sorted_points = sorted_indices
            .iter()
            .map(|&index| points[index as usize])
            .collect();
        Ok(Self {
            cell_size,
            sorted_points,
            sorted_indices,
            index: CellIndex::from_sorted(&keyed),
        })
    }

    /// Returns the number of points in the grid.
    pub fn len(&self) -> usize {
        self.sorted_points.len()
    }

    /// Returns `true` if the grid holds no points.
    pub fn is_empty(&self) -> bool {
        self.sorted_points.is_empty()
    }

    /// Returns the edge length of a grid cell.
    pub fn cell_size(&self) -> f64 {
        self.cell_size
    }
}

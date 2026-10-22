mod query;

use crate::cell::cell_key;
use crate::check;
use crate::error::GeomError;
use crate::index::CellIndex;

/// A set of points indexed by cubic cells for distance queries.
///
/// Points are identified by their position in the slice passed to
/// [`NeighbourGrid::new`]. Distances are Euclidean and every radius is
/// inclusive.
///
/// # Examples
///
/// ```
/// use mmcif_geom::NeighbourGrid;
///
/// let grid = NeighbourGrid::new(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], 4.0)?;
///
/// assert_eq!(grid.len(), 2);
/// assert_eq!(grid.within([0.0, 0.0, 0.0], 2.0)?, [0, 1]);
/// # Ok::<(), mmcif_geom::GeomError>(())
/// ```
#[derive(Debug, Clone)]
pub struct NeighbourGrid {
    cell_size: f64,
    points: Vec<[f64; 3]>,
    sorted_points: Vec<[f64; 3]>,
    sorted_indices: Vec<u32>,
    index: CellIndex,
}

impl NeighbourGrid {
    /// Builds a grid over `points` with cubic cells of edge length `cell_size`.
    ///
    /// The grid keeps its own copy of `points`.
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
    /// let grid = NeighbourGrid::new(&[[0.0, 0.0, 0.0]], 4.0)?;
    /// assert_eq!(grid.cell_size(), 4.0);
    ///
    /// assert_eq!(
    ///     NeighbourGrid::new(&[[0.0, 0.0, 0.0]], 0.0).unwrap_err(),
    ///     GeomError::InvalidCellSize(0.0)
    /// );
    /// # Ok::<(), GeomError>(())
    /// ```
    pub fn new(points: &[[f64; 3]], cell_size: f64) -> Result<Self, GeomError> {
        let cell_size = check::cell_size(cell_size)?;
        check::point_count(points.len())?;
        check::finite_points(points)?;
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
            points: points.to_vec(),
            sorted_points,
            sorted_indices,
            index: CellIndex::from_sorted(&keyed),
        })
    }

    /// Returns the number of points in the grid.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Returns `true` if the grid holds no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Returns the edge length of a grid cell.
    pub fn cell_size(&self) -> f64 {
        self.cell_size
    }

    /// Returns the points in the order they were given.
    pub fn points(&self) -> &[[f64; 3]] {
        &self.points
    }
}

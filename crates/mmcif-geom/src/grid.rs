mod query;

use crate::cell::cell_key;
use crate::check;
use crate::error::GeomError;
use crate::index::CellIndex;

#[derive(Debug, Clone)]
pub struct NeighbourGrid {
    cell_size: f64,
    points: Vec<[f64; 3]>,
    sorted_points: Vec<[f64; 3]>,
    sorted_indices: Vec<u32>,
    index: CellIndex,
}

impl NeighbourGrid {
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

    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    pub fn cell_size(&self) -> f64 {
        self.cell_size
    }

    pub fn points(&self) -> &[[f64; 3]] {
        &self.points
    }
}

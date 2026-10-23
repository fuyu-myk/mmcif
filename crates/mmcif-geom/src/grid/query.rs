use std::borrow::Borrow;
use std::ops::Range;

use super::NeighbourGrid;
use crate::cell::{CellRange, distance_squared};
use crate::check;
use crate::error::GeomError;
use crate::radius::Radius;

struct Probe {
    point: [f64; 3],
    radius_squared: f64,
    cells: CellRange,
}

impl Probe {
    fn new(point: [f64; 3], radius: f64, cell_size: f64) -> Self {
        let radius_squared = radius * radius;
        let cells = if radius_squared.is_finite() {
            CellRange::around(point, radius, cell_size)
        } else {
            CellRange::UNBOUNDED
        };
        Self {
            point,
            radius_squared,
            cells,
        }
    }
}

impl NeighbourGrid {
    /// Returns the indices of the points within `radius` of `point`, in
    /// ascending order.
    ///
    /// A `point` with a NaN or infinite coordinate finds nothing.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_geom::{NeighbourGrid, Radius};
    ///
    /// let points = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [9.0, 0.0, 0.0]];
    /// let grid = NeighbourGrid::new(&points, 4.0)?;
    /// let radius = Radius::new(2.0)?;
    ///
    /// assert_eq!(grid.within([1.0, 0.0, 0.0], radius), [0, 1]);
    /// assert!(grid.within([20.0, 0.0, 0.0], radius).is_empty());
    /// # Ok::<(), mmcif_geom::GeomError>(())
    /// ```
    pub fn within(&self, point: [f64; 3], radius: Radius) -> Vec<u32> {
        let mut found = Vec::new();
        if let Some(probe) = self.probe(point, radius.get()) {
            self.visit_matches(&probe, |index, _| found.push(index));
        }
        found.sort_unstable();
        found
    }

    /// Calls `visit` with the index and distance of each point within `radius`
    /// of `point`, in no particular order.
    ///
    /// A `point` with a NaN or infinite coordinate finds nothing, so `visit` is
    /// not called.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_geom::{NeighbourGrid, Radius};
    ///
    /// let grid = NeighbourGrid::new([[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]], 4.0)?;
    ///
    /// let mut found = Vec::new();
    /// grid.for_each_within([0.0, 0.0, 0.0], Radius::new(5.0)?, |index, distance| {
    ///     found.push((index, distance));
    /// });
    /// found.sort_by_key(|&(index, _)| index);
    ///
    /// assert_eq!(found, [(0, 0.0), (1, 5.0)]);
    /// # Ok::<(), mmcif_geom::GeomError>(())
    /// ```
    pub fn for_each_within(
        &self,
        point: [f64; 3],
        radius: Radius,
        mut visit: impl FnMut(u32, f64),
    ) {
        if let Some(probe) = self.probe(point, radius.get()) {
            self.visit_matches(&probe, |index, squared| visit(index, squared.sqrt()));
        }
    }

    /// Returns the index and distance of the point closest to `point`, if one
    /// lies within `max_radius`.
    ///
    /// When several points are equally close, the one with the lowest index
    /// is returned. A `point` with a NaN or infinite coordinate finds nothing.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_geom::{NeighbourGrid, Radius};
    ///
    /// let points = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
    /// let grid = NeighbourGrid::new(&points, 4.0)?;
    ///
    /// assert_eq!(grid.nearest([8.0, 0.0, 0.0], Radius::new(5.0)?), Some((1, 2.0)));
    /// assert_eq!(grid.nearest([5.0, 0.0, 0.0], Radius::new(5.0)?), Some((0, 5.0)));
    /// assert_eq!(grid.nearest([5.0, 0.0, 0.0], Radius::new(4.0)?), None);
    /// assert_eq!(grid.nearest([50.0, 0.0, 0.0], Radius::INFINITE), Some((1, 40.0)));
    /// # Ok::<(), mmcif_geom::GeomError>(())
    /// ```
    pub fn nearest(&self, point: [f64; 3], max_radius: Radius) -> Option<(u32, f64)> {
        let max_radius = max_radius.get();
        let mut radius = self.cell_size.min(max_radius);
        loop {
            let probe = self.probe(point, radius)?;
            if radius < max_radius && probe.cells.contains(self.index.occupied()) {
                radius = max_radius;
                continue;
            }
            let closest = self.closest(&probe);
            if closest.is_some() || radius >= max_radius {
                return closest.map(|(index, squared)| (index, squared.sqrt()));
            }
            radius = (radius * 2.0).min(max_radius);
        }
    }

    /// Returns every pair of points within `radius` of each other.
    ///
    /// Each pair is `(lower, higher)` by index, and the pairs are in ascending
    /// order.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_geom::{NeighbourGrid, Radius};
    ///
    /// let points = [[0.0, 0.0, 0.0], [9.0, 0.0, 0.0], [1.0, 0.0, 0.0], [9.5, 0.0, 0.0]];
    /// let grid = NeighbourGrid::new(&points, 4.0)?;
    ///
    /// assert_eq!(grid.pairs_within(Radius::new(1.0)?), [(0, 2), (1, 3)]);
    /// # Ok::<(), mmcif_geom::GeomError>(())
    /// ```
    pub fn pairs_within(&self, radius: Radius) -> Vec<(u32, u32)> {
        let radius = radius.get();
        let radius_squared = radius * radius;
        let mut pairs = Vec::new();
        let mut runs = Vec::new();
        for cell in self.index.cells() {
            let neighbourhood = cell
                .clone()
                .map(|position| {
                    Probe::new(self.sorted_points[position], radius, self.cell_size).cells
                })
                .fold(CellRange::EMPTY, CellRange::union);
            runs.clear();
            self.index
                .for_each_run(&neighbourhood, |run| runs.push(run));
            for position in cell {
                self.collect_later_pairs(position, radius_squared, &runs, &mut pairs);
            }
        }
        pairs.sort_unstable();
        pairs
    }

    /// Returns an `(other, point)` index pair for each point in `others` that
    /// lies within `radius` of a point in the grid.
    ///
    /// `others` is any sequence of points or references to them. Pairs are
    /// ordered by `other` and then by `point`. A point in `others` with a NaN
    /// or infinite coordinate is in no pair.
    ///
    /// # Errors
    ///
    /// Returns [`GeomError::TooManyPoints`] if `others` holds more than
    /// `u32::MAX` points.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_geom::{NeighbourGrid, Radius};
    ///
    /// let grid = NeighbourGrid::new([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]], 4.0)?;
    /// let others = [[10.5, 0.0, 0.0], [50.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
    ///
    /// assert_eq!(grid.contacts(&others, Radius::new(1.0)?)?, [(0, 1), (2, 0)]);
    /// # Ok::<(), mmcif_geom::GeomError>(())
    /// ```
    pub fn contacts<P>(
        &self,
        others: impl IntoIterator<Item = P>,
        radius: Radius,
    ) -> Result<Vec<(u32, u32)>, GeomError>
    where
        P: Borrow<[f64; 3]>,
    {
        let mut contacts = Vec::new();
        let mut found = Vec::new();
        let mut others = others.into_iter().enumerate();
        while let Some((index, point)) = others.next() {
            let Ok(other) = u32::try_from(index) else {
                return Err(GeomError::TooManyPoints(index + 1 + others.count()));
            };
            let Some(probe) = self.probe(*point.borrow(), radius.get()) else {
                continue;
            };
            found.clear();
            self.visit_matches(&probe, |index, _| found.push(index));
            found.sort_unstable();
            contacts.extend(found.iter().map(|&index| (other, index)));
        }
        Ok(contacts)
    }

    fn probe(&self, point: [f64; 3], radius: f64) -> Option<Probe> {
        check::is_finite(point).then(|| Probe::new(point, radius, self.cell_size))
    }

    fn closest(&self, probe: &Probe) -> Option<(u32, f64)> {
        let mut closest: Option<(u32, f64)> = None;
        self.visit_matches(probe, |index, squared| {
            if closest.is_none_or(|(best, best_squared)| (squared, index) < (best_squared, best)) {
                closest = Some((index, squared));
            }
        });
        closest
    }

    fn collect_later_pairs(
        &self,
        position: usize,
        radius_squared: f64,
        runs: &[Range<usize>],
        pairs: &mut Vec<(u32, u32)>,
    ) {
        let point = self.sorted_points[position];
        let first = self.sorted_indices[position];
        for run in runs {
            let later = run.start.max(position + 1)..run.end;
            self.visit_run(point, radius_squared, later, |second, _| {
                pairs.push((first.min(second), first.max(second)));
            });
        }
    }

    fn visit_matches(&self, probe: &Probe, mut visit: impl FnMut(u32, f64)) {
        self.index.for_each_run(&probe.cells, |run| {
            self.visit_run(probe.point, probe.radius_squared, run, &mut visit);
        });
    }

    fn visit_run(
        &self,
        point: [f64; 3],
        radius_squared: f64,
        run: Range<usize>,
        mut visit: impl FnMut(u32, f64),
    ) {
        if run.is_empty() {
            return;
        }
        let candidates = self.sorted_points[run.clone()]
            .iter()
            .zip(&self.sorted_indices[run]);
        for (&candidate, &index) in candidates {
            let squared = distance_squared(point, candidate);
            if squared <= radius_squared {
                visit(index, squared);
            }
        }
    }
}

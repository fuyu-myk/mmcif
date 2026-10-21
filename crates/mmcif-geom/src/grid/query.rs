use std::ops::Range;

use super::NeighbourGrid;
use crate::cell::{CellRange, distance_squared};
use crate::check;
use crate::error::GeomError;

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

    fn checked(point: [f64; 3], radius: f64, cell_size: f64) -> Result<Self, GeomError> {
        let radius = check::radius(radius)?;
        let point = check::query_point(point)?;
        Ok(Self::new(point, radius, cell_size))
    }
}

impl NeighbourGrid {
    pub fn within(&self, point: [f64; 3], radius: f64) -> Result<Vec<u32>, GeomError> {
        let probe = Probe::checked(point, radius, self.cell_size)?;
        let mut found = Vec::new();
        self.visit_matches(&probe, |index, _| found.push(index));
        found.sort_unstable();
        Ok(found)
    }

    pub fn for_each_within(
        &self,
        point: [f64; 3],
        radius: f64,
        mut visit: impl FnMut(u32, f64),
    ) -> Result<(), GeomError> {
        let probe = Probe::checked(point, radius, self.cell_size)?;
        self.visit_matches(&probe, |index, squared| visit(index, squared.sqrt()));
        Ok(())
    }

    pub fn nearest(
        &self,
        point: [f64; 3],
        max_radius: f64,
    ) -> Result<Option<(u32, f64)>, GeomError> {
        let max_radius = check::radius(max_radius)?;
        let point = check::query_point(point)?;
        let mut radius = self.cell_size.min(max_radius);
        loop {
            let probe = Probe::new(point, radius, self.cell_size);
            if radius < max_radius && probe.cells.contains(self.index.occupied()) {
                radius = max_radius;
                continue;
            }
            let closest = self.closest(&probe);
            if closest.is_some() || radius >= max_radius {
                return Ok(closest.map(|(index, squared)| (index, squared.sqrt())));
            }
            radius = (radius * 2.0).min(max_radius);
        }
    }

    pub fn pairs_within(&self, radius: f64) -> Result<Vec<(u32, u32)>, GeomError> {
        let radius = check::radius(radius)?;
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
        Ok(pairs)
    }

    pub fn contacts(&self, others: &[[f64; 3]], radius: f64) -> Result<Vec<(u32, u32)>, GeomError> {
        let radius = check::radius(radius)?;
        check::point_count(others.len())?;
        check::finite_points(others)?;
        let mut contacts = Vec::new();
        let mut found = Vec::new();
        for (&point, other) in others.iter().zip(0u32..) {
            let probe = Probe::new(point, radius, self.cell_size);
            found.clear();
            self.visit_matches(&probe, |index, _| found.push(index));
            found.sort_unstable();
            contacts.extend(found.iter().map(|&index| (other, index)));
        }
        Ok(contacts)
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

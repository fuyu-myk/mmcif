use std::ops::Range;

use crate::cell::CellRange;

#[derive(Debug, Clone)]
pub(crate) struct CellIndex {
    columns: Vec<[i64; 2]>,
    column_cell_starts: Vec<u32>,
    cell_layers: Vec<i64>,
    cell_point_starts: Vec<u32>,
    occupied: CellRange,
}

impl CellIndex {
    pub(crate) fn from_sorted(keyed: &[([i64; 3], u32)]) -> Self {
        let mut index = Self {
            columns: Vec::new(),
            column_cell_starts: Vec::new(),
            cell_layers: Vec::new(),
            cell_point_starts: Vec::new(),
            occupied: CellRange::EMPTY,
        };
        let mut position = 0;
        for cell in keyed.chunk_by(|a, b| a.0 == b.0) {
            let [x, y, z] = cell[0].0;
            index.occupied = index.occupied.union(CellRange::at([x, y, z]));
            if index.columns.last() != Some(&[x, y]) {
                index.columns.push([x, y]);
                index
                    .column_cell_starts
                    .push(index.cell_layers.len() as u32);
            }
            index.cell_layers.push(z);
            index.cell_point_starts.push(position as u32);
            position += cell.len();
        }
        index
            .column_cell_starts
            .push(index.cell_layers.len() as u32);
        index.cell_point_starts.push(position as u32);
        index
    }

    pub(crate) fn cells(&self) -> impl Iterator<Item = Range<usize>> {
        self.cell_point_starts
            .windows(2)
            .map(|bounds| bounds[0] as usize..bounds[1] as usize)
    }

    pub(crate) fn occupied(&self) -> &CellRange {
        &self.occupied
    }

    fn scans_all_columns(&self, cells: &CellRange) -> bool {
        cells.row_count() > self.columns.len() as u128
    }

    pub(crate) fn for_each_run(&self, cells: &CellRange, mut visit: impl FnMut(Range<usize>)) {
        if self.scans_all_columns(cells) {
            for (column, &key) in self.columns.iter().enumerate() {
                if cells.contains_column(key) {
                    self.visit_column(column, cells, &mut visit);
                }
            }
            return;
        }
        let mut first = 0;
        for row in cells.rows() {
            let wanted = cells.row_columns(row);
            first += self.columns[first..].partition_point(|key| key < wanted.start());
            let count = self.columns[first..]
                .iter()
                .take_while(|key| *key <= wanted.end())
                .count();
            for column in first..first + count {
                self.visit_column(column, cells, &mut visit);
            }
            first += count;
        }
    }

    fn visit_column(&self, column: usize, cells: &CellRange, visit: &mut impl FnMut(Range<usize>)) {
        let cell_start = self.column_cell_starts[column] as usize;
        let cell_end = self.column_cell_starts[column + 1] as usize;
        let layers = &self.cell_layers[cell_start..cell_end];
        let wanted = cells.layers();
        let skipped = layers.partition_point(|z| z < wanted.start());
        let count = layers[skipped..]
            .iter()
            .take_while(|z| *z <= wanted.end())
            .count();
        let first = cell_start + skipped;
        let end = first + count;
        if count > 0 {
            let points_start = self.cell_point_starts[first] as usize;
            let points_end = self.cell_point_starts[end] as usize;
            visit(points_start..points_end);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn index_of(keys: &[[i64; 3]]) -> CellIndex {
        let mut keyed: Vec<_> = keys.iter().copied().zip(0u32..).collect();
        keyed.sort_unstable();
        CellIndex::from_sorted(&keyed)
    }

    fn runs(index: &CellIndex, cells: &CellRange) -> Vec<Range<usize>> {
        let mut runs = Vec::new();
        index.for_each_run(cells, |run| runs.push(run));
        runs
    }

    #[test]
    fn groups_points_into_columns_and_layers() {
        let index = index_of(&[[0, 0, 1], [0, 0, 0], [0, 0, 1], [5, -2, 9], [0, 1, 0]]);
        assert_eq!(index.columns, vec![[0, 0], [0, 1], [5, -2]]);
        assert_eq!(index.column_cell_starts, vec![0, 2, 3, 4]);
        assert_eq!(index.cell_layers, vec![0, 1, 0, 9]);
        assert_eq!(index.cell_point_starts, vec![0, 1, 3, 4, 5]);
        assert_eq!(
            index.occupied,
            CellRange::at([0, -2, 0]).union(CellRange::at([5, 1, 9]))
        );
    }

    #[test]
    fn sparse_keys_use_memory_proportional_to_points() {
        let keys: Vec<[i64; 3]> = (0..1000)
            .map(|i| [i * 1_000_003, -i * 7_919, i * i])
            .collect();
        let index = index_of(&keys);
        assert_eq!(index.columns.len(), 1000);
        assert_eq!(index.cell_layers.len(), 1000);
    }

    #[test]
    fn row_lookup_and_column_scan_cover_the_requested_cells() {
        let keys: Vec<[i64; 3]> = (-3..3)
            .flat_map(|x| (-3..3).flat_map(move |y| (-3..3).map(move |z| [x, y, z])))
            .collect();
        let index = index_of(&keys);
        let near = CellRange::around([0.5, 0.5, 0.5], 1.0, 1.0);
        assert!(!index.scans_all_columns(&near));
        let wide = CellRange::around([0.5, 0.5, 0.5], 100.0, 1.0);
        assert!(index.scans_all_columns(&wide));
        let all: usize = runs(&index, &wide).iter().map(ExactSizeIterator::len).sum();
        assert_eq!(all, keys.len());
        let nearby: usize = runs(&index, &near).iter().map(ExactSizeIterator::len).sum();
        assert_eq!(nearby, 3 * 3 * 3);
    }

    #[test]
    fn empty_index_visits_nothing() {
        let index = index_of(&[]);
        assert!(runs(&index, &CellRange::UNBOUNDED).is_empty());
        assert_eq!(index.occupied, CellRange::EMPTY);
    }
}

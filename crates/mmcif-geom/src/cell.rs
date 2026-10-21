use std::ops::RangeInclusive;

// The distance test can accept a point a few ulps, or one underflowed square, beyond the radius;
// padding the reach by these slacks keeps every accepted point inside the visited cells.
const RELATIVE_SLACK: f64 = 8.0 * f64::EPSILON;
const UNDERFLOW_SLACK: f64 = 1e-150;

pub(crate) fn distance_squared(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

pub(crate) fn cell_key(point: [f64; 3], cell_size: f64) -> [i64; 3] {
    point.map(|value| cell_coordinate(value, cell_size))
}

fn cell_coordinate(value: f64, cell_size: f64) -> i64 {
    (value / cell_size).floor() as i64
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CellRange {
    lower: [i64; 3],
    upper: [i64; 3],
}

impl CellRange {
    pub(crate) const UNBOUNDED: Self = Self {
        lower: [i64::MIN; 3],
        upper: [i64::MAX; 3],
    };

    pub(crate) const EMPTY: Self = Self {
        lower: [i64::MAX; 3],
        upper: [i64::MIN; 3],
    };

    pub(crate) fn at(key: [i64; 3]) -> Self {
        Self {
            lower: key,
            upper: key,
        }
    }

    pub(crate) fn around(point: [f64; 3], radius: f64, cell_size: f64) -> Self {
        let magnitude = point
            .iter()
            .fold(0.0_f64, |largest, value| largest.max(value.abs()));
        let reach = radius + (radius + magnitude) * RELATIVE_SLACK + UNDERFLOW_SLACK;
        Self {
            lower: point.map(|value| cell_coordinate(value - reach, cell_size)),
            upper: point.map(|value| cell_coordinate(value + reach, cell_size)),
        }
    }

    pub(crate) fn union(self, other: Self) -> Self {
        Self {
            lower: [0, 1, 2].map(|axis| self.lower[axis].min(other.lower[axis])),
            upper: [0, 1, 2].map(|axis| self.upper[axis].max(other.upper[axis])),
        }
    }

    pub(crate) fn contains(&self, other: &Self) -> bool {
        (0..3).all(|axis| {
            self.lower[axis] <= other.lower[axis] && other.upper[axis] <= self.upper[axis]
        })
    }

    pub(crate) fn row_count(&self) -> u128 {
        (i128::from(self.upper[0]) - i128::from(self.lower[0]) + 1) as u128
    }

    pub(crate) fn rows(&self) -> RangeInclusive<i64> {
        self.lower[0]..=self.upper[0]
    }

    pub(crate) fn row_columns(&self, row: i64) -> RangeInclusive<[i64; 2]> {
        [row, self.lower[1]]..=[row, self.upper[1]]
    }

    pub(crate) fn contains_column(&self, [x, y]: [i64; 2]) -> bool {
        (self.lower[0]..=self.upper[0]).contains(&x) && (self.lower[1]..=self.upper[1]).contains(&y)
    }

    pub(crate) fn layers(&self) -> RangeInclusive<i64> {
        self.lower[2]..=self.upper[2]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cell_key_floors_negative_coordinates() {
        assert_eq!(cell_key([-0.5, 0.0, 3.99], 2.0), [-1, 0, 1]);
        assert_eq!(cell_key([-2.0, 2.0, 4.0], 2.0), [-1, 1, 2]);
    }

    #[test]
    fn range_spans_every_cell_within_reach() {
        let range = CellRange::around([5.0, -5.0, 0.6], 2.5, 1.0);
        assert_eq!(range.rows(), 2..=7);
        assert_eq!(range.row_columns(3), [3, -8]..=[3, -3]);
        assert_eq!(range.layers(), -2..=3);
        assert_eq!(range.row_count(), 6);
    }

    #[test]
    fn range_reaching_a_cell_boundary_includes_the_cell_beyond() {
        let range = CellRange::around([0.5, 0.5, 0.5], 2.5, 1.0);
        assert_eq!(range.layers(), -3..=3);
    }

    #[test]
    fn range_keeps_points_lost_to_an_underflowed_square() {
        let range = CellRange::around([0.0; 3], 0.0, 1.0);
        assert_eq!(range.rows(), -1..=0);
        assert_eq!(distance_squared([-1e-170, 0.0, 0.0], [0.0; 3]), 0.0);
    }

    #[test]
    fn empty_range_is_the_identity_of_union_and_inside_every_range() {
        let cell = CellRange::at([3, -4, 5]);
        assert_eq!(CellRange::EMPTY.union(cell), cell);
        assert!(cell.contains(&CellRange::EMPTY));
        assert!(!CellRange::EMPTY.contains(&cell));
        assert!(CellRange::UNBOUNDED.contains(&cell));
        assert!(!cell.contains(&cell.union(CellRange::at([3, -4, 6]))));
    }

    #[test]
    fn unbounded_range_counts_every_row() {
        assert_eq!(CellRange::UNBOUNDED.row_count(), 1 << 64);
        assert!(CellRange::UNBOUNDED.contains_column([i64::MIN, i64::MAX]));
    }
}

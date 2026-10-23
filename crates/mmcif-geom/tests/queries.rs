use mmcif_geom::{NeighbourGrid, Radius};

const ORIGIN: [f64; 3] = [0.0; 3];

fn grid(points: &[[f64; 3]], cell_size: f64) -> NeighbourGrid {
    NeighbourGrid::new(points, cell_size).unwrap()
}

fn r(value: f64) -> Radius {
    Radius::new(value).unwrap()
}

fn lattice(extent: i32) -> Vec<[f64; 3]> {
    let steps = || (-extent..=extent).map(f64::from);
    steps()
        .flat_map(|x| steps().flat_map(move |y| steps().map(move |z| [x, y, z])))
        .collect()
}

fn brute_within(points: &[[f64; 3]], point: [f64; 3], radius: f64) -> Vec<u32> {
    (0u32..)
        .zip(points)
        .filter(|&(_, &other)| {
            let [dx, dy, dz] = [0, 1, 2].map(|axis| point[axis] - other[axis]);
            dx * dx + dy * dy + dz * dz <= radius * radius
        })
        .map(|(index, _)| index)
        .collect()
}

#[test]
fn within_includes_points_exactly_at_the_radius() {
    let points = [
        [3.0, 4.0, 0.0],
        [3.0, 4.0, 0.001],
        ORIGIN,
        [-3.0, -4.0, 0.0],
        [0.0, 0.0, 5.0],
        [0.0, 0.0, -5.000_000_1],
    ];
    assert_eq!(grid(&points, 1.0).within(ORIGIN, r(5.0)), vec![0, 2, 3, 4]);
}

#[test]
fn within_stays_exact_when_the_radius_spans_many_cells() {
    let points = lattice(8);
    let query = [0.3, -0.2, 0.1];
    for (cell_size, radius) in [(0.7, 6.5), (1.0, 3.0), (0.25, 2.0), (50.0, 4.0)] {
        assert_eq!(
            grid(&points, cell_size).within(query, r(radius)),
            brute_within(&points, query, radius)
        );
    }
}

#[test]
fn within_zero_radius_finds_only_coincident_points() {
    let points = [[1.0, 1.0, 1.0], [1.0, 1.0, 1.000_000_1], [1.0, 1.0, 1.0]];
    assert_eq!(grid(&points, 0.5).within([1.0; 3], r(0.0)), vec![0, 2]);
}

#[test]
fn within_keeps_matches_that_rounding_places_past_the_reach() {
    let (query, radius) = (-3.4, 3.6158235594456634);
    let above = 0.21582355944566348;
    assert!(above > query + radius);
    let grid_above = grid(&[[above, 0.0, 0.0]], above);
    assert_eq!(grid_above.within([query, 0.0, 0.0], r(radius)), vec![0]);

    let (query, radius) = (6.025727701281269, 6.4);
    let below = -0.3742722987187311;
    assert!(below < query - radius);
    let grid_below = grid(&[[below, 0.0, 0.0]], radius - query);
    assert_eq!(grid_below.within([query, 0.0, 0.0], r(radius)), vec![0]);
}

#[test]
fn within_zero_radius_keeps_points_whose_square_underflows() {
    let grid = grid(&[[-1e-170, 0.0, 0.0], [1e-100, 0.0, 0.0]], 1.0);
    assert_eq!(grid.within(ORIGIN, r(0.0)), vec![0]);
}

#[test]
fn within_huge_radius_returns_every_point() {
    let points = [
        [1e6, -1e6, 0.0],
        ORIGIN,
        [-3.0, 7.0, 2.0],
        [1e300, 0.0, -1e300],
    ];
    let fine = grid(&points, 1e-3);
    assert_eq!(fine.within(ORIGIN, r(2e6)), vec![0, 1, 2]);
    assert_eq!(fine.within(ORIGIN, Radius::INFINITE), vec![0, 1, 2, 3]);
}

#[test]
fn radius_whose_square_overflows_matches_every_point() {
    let grid = grid(&[[f64::MAX, 0.0, 0.0], ORIGIN], 1e300);
    assert_eq!(grid.within(ORIGIN, r(1e155)), vec![0, 1]);
    let mut visited = Vec::new();
    grid.for_each_within(ORIGIN, r(1e155), |index, distance| {
        visited.push((index, distance))
    });
    visited.sort_by_key(|&(index, _)| index);
    assert_eq!(visited, vec![(0, f64::INFINITY), (1, 0.0)]);
}

#[test]
fn for_each_within_reports_plain_distances() {
    let points = [
        [3.0, 4.0, 0.0],
        [0.0, 0.0, 1.0],
        [6.0, 0.0, 0.0],
        [0.0, -2.0, 0.0],
    ];
    let mut visited = Vec::new();
    grid(&points, 1.0).for_each_within(ORIGIN, r(5.0), |index, distance| {
        visited.push((index, distance))
    });
    visited.sort_by_key(|&(index, _)| index);
    assert_eq!(visited, vec![(0, 5.0), (1, 1.0), (3, 2.0)]);
}

#[test]
fn nearest_returns_the_closest_point_and_its_distance() {
    let points = [[10.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
    assert_eq!(grid(&points, 1.0).nearest(ORIGIN, r(100.0)), Some((1, 2.0)));
}

#[test]
fn nearest_breaks_ties_by_lower_index() {
    let points = [
        [5.0, 5.0, 5.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    assert_eq!(grid(&points, 0.5).nearest(ORIGIN, r(10.0)), Some((1, 1.0)));
}

#[test]
fn nearest_includes_a_point_exactly_at_the_max_radius() {
    let grid = grid(&[[3.0, 0.0, 4.0]], 1.0);
    assert_eq!(grid.nearest(ORIGIN, r(5.0)), Some((0, 5.0)));
    assert_eq!(grid.nearest(ORIGIN, r(4.999)), None);
}

#[test]
fn nearest_searches_beyond_the_neighbouring_cells() {
    let grid = grid(&[[1000.0, 0.0, 0.0], [-1500.0, 0.0, 0.0]], 1.0);
    assert_eq!(grid.nearest(ORIGIN, Radius::INFINITE), Some((0, 1000.0)));
    assert_eq!(grid.nearest(ORIGIN, r(999.0)), None);
}

#[test]
fn nearest_with_zero_max_radius_needs_a_coincident_point() {
    let grid = grid(&[[1.0, 0.0, 0.0], ORIGIN, ORIGIN], 1.0);
    assert_eq!(grid.nearest(ORIGIN, r(0.0)), Some((1, 0.0)));
    assert_eq!(grid.nearest([0.5, 0.0, 0.0], r(0.0)), None);
}

#[test]
fn pairs_within_lists_each_pair_once_in_order() {
    let points = [ORIGIN, [0.0, 0.0, 1.0], [0.0, 0.0, 2.5], [0.0, 0.0, 1.0]];
    assert_eq!(
        grid(&points, 1.0).pairs_within(r(1.5)),
        vec![(0, 1), (0, 3), (1, 2), (1, 3), (2, 3)]
    );
}

#[test]
fn pairs_within_stays_exact_when_the_radius_spans_many_cells() {
    let points = lattice(3);
    let radius = 2.0;
    let expected: Vec<(u32, u32)> = (0u32..)
        .zip(&points)
        .flat_map(|(first, &point)| {
            brute_within(&points, point, radius)
                .into_iter()
                .filter(move |&second| second > first)
                .map(move |second| (first, second))
        })
        .collect();
    assert_eq!(grid(&points, 0.3).pairs_within(r(radius)), expected);
}

#[test]
fn pairs_within_zero_radius_pairs_duplicates() {
    let points = [[2.0; 3], [1.0; 3], [2.0; 3], [2.0; 3]];
    assert_eq!(
        grid(&points, 1.0).pairs_within(r(0.0)),
        vec![(0, 2), (0, 3), (2, 3)]
    );
}

#[test]
fn contacts_pairs_other_index_with_grid_index() {
    let grid = grid(&[ORIGIN, [4.0, 0.0, 0.0]], 1.0);
    let others = [[2.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    assert_eq!(
        grid.contacts(others, r(2.0)).unwrap(),
        vec![(0, 0), (0, 1), (2, 0)]
    );
    assert_eq!(
        grid.contacts(std::iter::empty::<[f64; 3]>(), r(2.0))
            .unwrap(),
        Vec::new()
    );
}

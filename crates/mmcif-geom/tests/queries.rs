use mmcif_geom::{GeomError, NeighbourGrid};

const ORIGIN: [f64; 3] = [0.0; 3];
const INVALID_RADII: [f64; 3] = [-1.0, f64::NEG_INFINITY, f64::NAN];
const NON_FINITE_QUERIES: [[f64; 3]; 3] = [
    [f64::NAN, 0.0, 0.0],
    [0.0, f64::INFINITY, 0.0],
    [0.0, 0.0, f64::NEG_INFINITY],
];

fn grid(points: &[[f64; 3]], cell_size: f64) -> NeighbourGrid {
    NeighbourGrid::new(points, cell_size).unwrap()
}

fn is_invalid_radius(error: GeomError, radius: f64) -> bool {
    match error {
        GeomError::InvalidRadius(reported) => reported.total_cmp(&radius).is_eq(),
        _ => false,
    }
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
    assert_eq!(
        grid(&points, 1.0).within(ORIGIN, 5.0).unwrap(),
        vec![0, 2, 3, 4]
    );
}

#[test]
fn within_stays_exact_when_the_radius_spans_many_cells() {
    let points = lattice(8);
    let query = [0.3, -0.2, 0.1];
    for (cell_size, radius) in [(0.7, 6.5), (1.0, 3.0), (0.25, 2.0), (50.0, 4.0)] {
        assert_eq!(
            grid(&points, cell_size).within(query, radius).unwrap(),
            brute_within(&points, query, radius)
        );
    }
}

#[test]
fn within_zero_radius_finds_only_coincident_points() {
    let points = [[1.0, 1.0, 1.0], [1.0, 1.0, 1.000_000_1], [1.0, 1.0, 1.0]];
    assert_eq!(
        grid(&points, 0.5).within([1.0; 3], 0.0).unwrap(),
        vec![0, 2]
    );
}

#[test]
fn within_keeps_matches_that_rounding_places_past_the_reach() {
    let (query, radius) = (-3.4, 3.6158235594456634);
    let above = 0.21582355944566348;
    assert!(above > query + radius);
    let grid_above = grid(&[[above, 0.0, 0.0]], above);
    assert_eq!(
        grid_above.within([query, 0.0, 0.0], radius).unwrap(),
        vec![0]
    );

    let (query, radius) = (6.025727701281269, 6.4);
    let below = -0.3742722987187311;
    assert!(below < query - radius);
    let grid_below = grid(&[[below, 0.0, 0.0]], radius - query);
    assert_eq!(
        grid_below.within([query, 0.0, 0.0], radius).unwrap(),
        vec![0]
    );
}

#[test]
fn within_zero_radius_keeps_points_whose_square_underflows() {
    let grid = grid(&[[-1e-170, 0.0, 0.0], [1e-100, 0.0, 0.0]], 1.0);
    assert_eq!(grid.within(ORIGIN, 0.0).unwrap(), vec![0]);
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
    assert_eq!(fine.within(ORIGIN, 2e6).unwrap(), vec![0, 1, 2]);
    assert_eq!(
        fine.within(ORIGIN, f64::INFINITY).unwrap(),
        vec![0, 1, 2, 3]
    );
}

#[test]
fn radius_whose_square_overflows_matches_every_point() {
    let grid = grid(&[[f64::MAX, 0.0, 0.0], ORIGIN], 1e300);
    assert_eq!(grid.within(ORIGIN, 1e155).unwrap(), vec![0, 1]);
    let mut visited = Vec::new();
    grid.for_each_within(ORIGIN, 1e155, |index, distance| {
        visited.push((index, distance))
    })
    .unwrap();
    visited.sort_by_key(|&(index, _)| index);
    assert_eq!(visited, vec![(0, f64::INFINITY), (1, 0.0)]);
}

#[test]
fn within_rejects_invalid_radius_and_query() {
    let grid = grid(&[ORIGIN], 1.0);
    for radius in INVALID_RADII {
        assert!(is_invalid_radius(
            grid.within(ORIGIN, radius).unwrap_err(),
            radius
        ));
    }
    for query in NON_FINITE_QUERIES {
        assert_eq!(
            grid.within(query, 1.0).unwrap_err(),
            GeomError::NonFiniteQuery
        );
    }
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
    grid(&points, 1.0)
        .for_each_within(ORIGIN, 5.0, |index, distance| {
            visited.push((index, distance))
        })
        .unwrap();
    visited.sort_by_key(|&(index, _)| index);
    assert_eq!(visited, vec![(0, 5.0), (1, 1.0), (3, 2.0)]);
}

#[test]
fn for_each_within_rejects_invalid_radius_and_query() {
    let grid = grid(&[ORIGIN], 1.0);
    let mut calls = 0;
    for radius in INVALID_RADII {
        let error = grid
            .for_each_within(ORIGIN, radius, |_, _| calls += 1)
            .unwrap_err();
        assert!(is_invalid_radius(error, radius));
    }
    for query in NON_FINITE_QUERIES {
        let error = grid
            .for_each_within(query, 1.0, |_, _| calls += 1)
            .unwrap_err();
        assert_eq!(error, GeomError::NonFiniteQuery);
    }
    assert_eq!(calls, 0);
}

#[test]
fn nearest_returns_the_closest_point_and_its_distance() {
    let points = [[10.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
    assert_eq!(
        grid(&points, 1.0).nearest(ORIGIN, 100.0).unwrap(),
        Some((1, 2.0))
    );
}

#[test]
fn nearest_breaks_ties_by_lower_index() {
    let points = [
        [5.0, 5.0, 5.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    assert_eq!(
        grid(&points, 0.5).nearest(ORIGIN, 10.0).unwrap(),
        Some((1, 1.0))
    );
}

#[test]
fn nearest_includes_a_point_exactly_at_the_max_radius() {
    let grid = grid(&[[3.0, 0.0, 4.0]], 1.0);
    assert_eq!(grid.nearest(ORIGIN, 5.0).unwrap(), Some((0, 5.0)));
    assert_eq!(grid.nearest(ORIGIN, 4.999).unwrap(), None);
}

#[test]
fn nearest_searches_beyond_the_neighbouring_cells() {
    let grid = grid(&[[1000.0, 0.0, 0.0], [-1500.0, 0.0, 0.0]], 1.0);
    assert_eq!(
        grid.nearest(ORIGIN, f64::INFINITY).unwrap(),
        Some((0, 1000.0))
    );
    assert_eq!(grid.nearest(ORIGIN, 999.0).unwrap(), None);
}

#[test]
fn nearest_with_zero_max_radius_needs_a_coincident_point() {
    let grid = grid(&[[1.0, 0.0, 0.0], ORIGIN, ORIGIN], 1.0);
    assert_eq!(grid.nearest(ORIGIN, 0.0).unwrap(), Some((1, 0.0)));
    assert_eq!(grid.nearest([0.5, 0.0, 0.0], 0.0).unwrap(), None);
}

#[test]
fn nearest_rejects_invalid_radius_and_query() {
    let grid = grid(&[ORIGIN], 1.0);
    for radius in INVALID_RADII {
        assert!(is_invalid_radius(
            grid.nearest(ORIGIN, radius).unwrap_err(),
            radius
        ));
    }
    for query in NON_FINITE_QUERIES {
        assert_eq!(
            grid.nearest(query, 1.0).unwrap_err(),
            GeomError::NonFiniteQuery
        );
    }
}

#[test]
fn pairs_within_lists_each_pair_once_in_order() {
    let points = [ORIGIN, [0.0, 0.0, 1.0], [0.0, 0.0, 2.5], [0.0, 0.0, 1.0]];
    assert_eq!(
        grid(&points, 1.0).pairs_within(1.5).unwrap(),
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
    assert_eq!(grid(&points, 0.3).pairs_within(radius).unwrap(), expected);
}

#[test]
fn pairs_within_zero_radius_pairs_duplicates() {
    let points = [[2.0; 3], [1.0; 3], [2.0; 3], [2.0; 3]];
    assert_eq!(
        grid(&points, 1.0).pairs_within(0.0).unwrap(),
        vec![(0, 2), (0, 3), (2, 3)]
    );
}

#[test]
fn pairs_within_rejects_invalid_radius() {
    let grid = grid(&[ORIGIN], 1.0);
    for radius in INVALID_RADII {
        assert!(is_invalid_radius(
            grid.pairs_within(radius).unwrap_err(),
            radius
        ));
    }
}

#[test]
fn contacts_pairs_other_index_with_grid_index() {
    let grid = grid(&[ORIGIN, [4.0, 0.0, 0.0]], 1.0);
    let others = [[2.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    assert_eq!(
        grid.contacts(&others, 2.0).unwrap(),
        vec![(0, 0), (0, 1), (2, 0)]
    );
    assert_eq!(grid.contacts(&[], 2.0).unwrap(), Vec::new());
}

#[test]
fn contacts_rejects_invalid_radius_and_non_finite_others() {
    let grid = grid(&[ORIGIN], 1.0);
    for radius in INVALID_RADII {
        assert!(is_invalid_radius(
            grid.contacts(&[ORIGIN], radius).unwrap_err(),
            radius
        ));
    }
    for query in NON_FINITE_QUERIES {
        assert_eq!(
            grid.contacts(&[ORIGIN, ORIGIN, query], 1.0).unwrap_err(),
            GeomError::NonFiniteCoordinate { index: 2 }
        );
    }
}

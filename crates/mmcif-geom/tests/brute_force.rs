use mmcif_geom::NeighbourGrid;
use proptest::array::uniform3;
use proptest::collection::vec;
use proptest::prelude::*;
use proptest::sample::{Index, select};

const QUERIES_FROM_POINTS: usize = 4;
const SMALLEST_SUBNORMAL: f64 = f64::from_bits(1);
const EXTREME_VALUES: [f64; 10] = [
    f64::MAX,
    -f64::MAX,
    1e300,
    -1e300,
    1e154,
    -3.5,
    f64::MIN_POSITIVE,
    SMALLEST_SUBNORMAL,
    0.0,
    -0.0,
];

fn distance_squared(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

fn brute_within(points: &[[f64; 3]], point: [f64; 3], radius: f64) -> Vec<(u32, f64)> {
    (0u32..)
        .zip(points)
        .map(|(index, &other)| (index, distance_squared(point, other)))
        .filter(|&(_, squared)| squared <= radius * radius)
        .collect()
}

fn brute_nearest(points: &[[f64; 3]], point: [f64; 3], max_radius: f64) -> Option<(u32, f64)> {
    brute_within(points, point, max_radius)
        .into_iter()
        .min_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)))
        .map(|(index, squared)| (index, squared.sqrt()))
}

fn brute_contacts(grid_points: &[[f64; 3]], others: &[[f64; 3]], radius: f64) -> Vec<(u32, u32)> {
    (0u32..)
        .zip(others)
        .flat_map(|(other, &point)| {
            brute_within(grid_points, point, radius)
                .into_iter()
                .map(move |(index, _)| (other, index))
        })
        .collect()
}

fn brute_pairs(points: &[[f64; 3]], radius: f64) -> Vec<(u32, u32)> {
    brute_contacts(points, points, radius)
        .into_iter()
        .filter(|&(first, second)| first < second)
        .collect()
}

fn lattice_coordinate() -> impl Strategy<Value = f64> {
    (-6i32..=6).prop_map(|step| f64::from(step) * 0.5)
}

fn extreme_coordinate() -> impl Strategy<Value = f64> {
    select(EXTREME_VALUES.to_vec())
}

fn cloud() -> impl Strategy<Value = Vec<[f64; 3]>> {
    prop_oneof![
        vec(uniform3(extreme_coordinate()), 0..40),
        vec(uniform3(-30.0f64..30.0), 0..100),
        vec(uniform3(-5000.0f64..5000.0), 0..60),
        (
            uniform3(-1000.0f64..1000.0),
            vec(uniform3(-0.75f64..0.75), 0..100)
        )
            .prop_map(|(centre, offsets)| {
                offsets
                    .into_iter()
                    .map(|offset| [0, 1, 2].map(|axis| centre[axis] + offset[axis]))
                    .collect()
            }),
        vec(uniform3(lattice_coordinate()), 0..100),
        (
            vec(uniform3(-4.0f64..4.0), 1..5),
            vec(any::<Index>(), 0..80)
        )
            .prop_map(|(distinct, picks)| picks.iter().map(|pick| *pick.get(&distinct)).collect()),
    ]
}

fn cell_size() -> impl Strategy<Value = f64> {
    prop_oneof![
        3 => 0.05f64..20.0,
        2 => select(vec![0.5, 1.0, 1.5, 2.0, 4.0]),
        1 => 1e3f64..1e7,
        1 => 1e-4f64..0.05,
        1 => select(vec![SMALLEST_SUBNORMAL, 1e-300, 1e150, 1e300, f64::MAX]),
    ]
}

fn radius() -> impl Strategy<Value = f64> {
    prop_oneof![
        4 => 0.0f64..12.0,
        2 => (0u8..=8).prop_map(|step| f64::from(step) * 0.5),
        1 => 12.0f64..3000.0,
        1 => Just(0.0),
        1 => Just(f64::INFINITY),
        1 => select(vec![1e154, 1.35e154, 1e155, 1e300, f64::MAX]),
    ]
}

fn query() -> impl Strategy<Value = [f64; 3]> {
    prop_oneof![
        uniform3(-40.0f64..40.0),
        uniform3(lattice_coordinate()),
        uniform3(-6000.0f64..6000.0),
        uniform3(extreme_coordinate()),
    ]
}

fn shifted_by_ulps(value: f64, ulps: i64) -> f64 {
    let bits = value.to_bits() as i64;
    let away_from_zero = if value.is_sign_negative() {
        -ulps
    } else {
        ulps
    };
    f64::from_bits((bits + away_from_zero) as u64)
}

fn tenths(range: std::ops::Range<i32>) -> impl Strategy<Value = f64> {
    range.prop_map(|tenth| f64::from(tenth) / 10.0)
}

fn queries_for(points: &[[f64; 3]], generated: &[[f64; 3]]) -> Vec<[f64; 3]> {
    generated
        .iter()
        .chain(points.iter().take(QUERIES_FROM_POINTS))
        .copied()
        .collect()
}

proptest! {
    #[test]
    fn within_matches_brute_force(
        points in cloud(),
        cell_size in cell_size(),
        radius in radius(),
        generated in vec(query(), 1..6),
    ) {
        let grid = NeighbourGrid::new(&points, cell_size).unwrap();
        for point in queries_for(&points, &generated) {
            let expected = brute_within(&points, point, radius);
            let indices: Vec<u32> = expected.iter().map(|&(index, _)| index).collect();
            prop_assert_eq!(grid.within(point, radius).unwrap(), indices);

            let mut visited = Vec::new();
            grid.for_each_within(point, radius, |index, distance| visited.push((index, distance)))
                .unwrap();
            visited.sort_by_key(|&(index, _)| index);
            let distances: Vec<(u32, f64)> = expected
                .iter()
                .map(|&(index, squared)| (index, squared.sqrt()))
                .collect();
            prop_assert_eq!(visited, distances);
        }
    }

    #[test]
    fn points_within_ulps_of_the_reach_match_brute_force(
        centre in prop_oneof![-100.0f64..100.0, tenths(-1000..1000)],
        radius in prop_oneof![0.0f64..20.0, tenths(0..200)],
        above in any::<bool>(),
        ulps in -3i64..=3,
        axis in 0usize..3,
    ) {
        let edge = if above { centre + radius } else { centre - radius };
        let coordinate = shifted_by_ulps(edge, ulps);
        let boundary = if above { coordinate.max(edge) } else { edge.max(coordinate) };
        prop_assume!(coordinate.is_finite() && boundary != 0.0);
        let query = [centre; 3];
        let mut point = query;
        point[axis] = coordinate;
        let grid = NeighbourGrid::new(&[point], boundary.abs()).unwrap();
        let expected: Vec<u32> = brute_within(&[point], query, radius)
            .into_iter()
            .map(|(index, _)| index)
            .collect();
        prop_assert_eq!(grid.within(query, radius).unwrap(), expected);
    }

    #[test]
    fn nearest_matches_brute_force(
        points in cloud(),
        cell_size in cell_size(),
        max_radius in radius(),
        generated in vec(query(), 1..6),
    ) {
        let grid = NeighbourGrid::new(&points, cell_size).unwrap();
        for point in queries_for(&points, &generated) {
            prop_assert_eq!(
                grid.nearest(point, max_radius).unwrap(),
                brute_nearest(&points, point, max_radius)
            );
        }
    }

    #[test]
    fn pairs_within_matches_brute_force(
        points in cloud(),
        cell_size in cell_size(),
        radius in radius(),
    ) {
        let grid = NeighbourGrid::new(&points, cell_size).unwrap();
        prop_assert_eq!(grid.pairs_within(radius).unwrap(), brute_pairs(&points, radius));
    }

    #[test]
    fn contacts_matches_brute_force(
        points in cloud(),
        others in cloud(),
        cell_size in cell_size(),
        radius in radius(),
    ) {
        let grid = NeighbourGrid::new(&points, cell_size).unwrap();
        prop_assert_eq!(grid.contacts(&others, radius).unwrap(), brute_contacts(&points, &others, radius));
    }
}

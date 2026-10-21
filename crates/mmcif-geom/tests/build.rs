use mmcif_geom::{GeomError, NeighbourGrid};

const ORIGIN: [f64; 3] = [0.0; 3];

#[test]
fn keeps_a_copy_of_the_points_in_input_order() {
    let points = [[1.0, 2.0, 3.0], [-4.0, 5.5, 0.0], [1.0, 2.0, 3.0]];
    let grid = NeighbourGrid::new(&points, 2.0).unwrap();
    assert_eq!(grid.len(), 3);
    assert!(!grid.is_empty());
    assert_eq!(grid.cell_size(), 2.0);
    assert_eq!(grid.points(), &points);
}

#[test]
fn empty_grid_answers_every_query_with_nothing() {
    let grid = NeighbourGrid::new(&[], 1.0).unwrap();
    assert!(grid.is_empty());
    assert_eq!(grid.len(), 0);
    assert!(grid.points().is_empty());
    assert_eq!(
        grid.within(ORIGIN, f64::INFINITY).unwrap(),
        Vec::<u32>::new()
    );
    assert_eq!(grid.nearest(ORIGIN, f64::INFINITY).unwrap(), None);
    assert_eq!(grid.pairs_within(f64::INFINITY).unwrap(), Vec::new());
    assert_eq!(grid.contacts(&[ORIGIN], f64::INFINITY).unwrap(), Vec::new());
}

#[test]
fn single_point_grid_has_no_pairs() {
    let grid = NeighbourGrid::new(&[[1.0, 1.0, 1.0]], 1.0).unwrap();
    assert_eq!(grid.pairs_within(f64::INFINITY).unwrap(), Vec::new());
    assert_eq!(
        grid.nearest(ORIGIN, 2.0).unwrap(),
        Some((0, 3.0_f64.sqrt()))
    );
}

#[test]
fn rejects_a_non_finite_coordinate_with_its_index() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let points = [ORIGIN, [1.0, 1.0, 1.0], [0.0, bad, 0.0], [bad; 3]];
        assert_eq!(
            NeighbourGrid::new(&points, 1.0).unwrap_err(),
            GeomError::NonFiniteCoordinate { index: 2 }
        );
    }
}

#[test]
fn rejects_a_cell_size_that_is_not_positive_and_finite() {
    for cell_size in [0.0, -0.0, -1.5, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            NeighbourGrid::new(&[ORIGIN], cell_size).unwrap_err(),
            GeomError::InvalidCellSize(cell_size)
        );
    }
    assert!(matches!(
        NeighbourGrid::new(&[ORIGIN], f64::NAN),
        Err(GeomError::InvalidCellSize(size)) if size.is_nan()
    ));
}

#[test]
fn accepts_extreme_but_finite_coordinates() {
    let points = [
        [f64::MAX, f64::MIN, 0.0],
        [f64::MIN_POSITIVE, -0.0, 1e-300],
        ORIGIN,
    ];
    let grid = NeighbourGrid::new(&points, 1e-300).unwrap();
    assert_eq!(grid.within(ORIGIN, 1e-290).unwrap(), vec![1, 2]);
    assert_eq!(grid.within(ORIGIN, f64::INFINITY).unwrap(), vec![0, 1, 2]);
    assert_eq!(
        grid.nearest([f64::MAX, f64::MIN, 1.0], f64::INFINITY)
            .unwrap(),
        Some((0, 1.0))
    );
}

#[test]
fn sparse_capsid_shell_builds_and_queries_exactly() {
    let shell = capsid_shell(20_000, 900.0);
    let grid = NeighbourGrid::new(&shell, 4.0).unwrap();
    let neighbours = grid.within(shell[0], 30.0).unwrap();
    assert!(neighbours.contains(&0));
    assert!(neighbours.len() < shell.len() / 100);

    let squared = |point: [f64; 3]| point.iter().map(|value| value * value).sum::<f64>();
    let closest_to_centre = (0u32..)
        .zip(&shell)
        .min_by(|a, b| squared(*a.1).total_cmp(&squared(*b.1)).then(a.0.cmp(&b.0)))
        .map(|(index, _)| index)
        .unwrap();
    let (index, distance) = grid.nearest(ORIGIN, 1000.0).unwrap().unwrap();
    assert_eq!(index, closest_to_centre);
    assert!((distance - 900.0).abs() < 1e-9);
    assert_eq!(grid.nearest(ORIGIN, 800.0).unwrap(), None);
}

fn capsid_shell(count: usize, radius: f64) -> Vec<[f64; 3]> {
    let golden_angle = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    (0..count)
        .map(|i| {
            let height = 1.0 - 2.0 * (i as f64 + 0.5) / count as f64;
            let ring = (1.0 - height * height).sqrt();
            let angle = golden_angle * i as f64;
            [
                radius * ring * angle.cos(),
                radius * ring * angle.sin(),
                radius * height,
            ]
        })
        .collect()
}

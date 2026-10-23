use mmcif_geom::{GeomError, NeighbourGrid, Radius};

#[test]
fn a_grid_is_built_from_any_iterator_of_points() {
    let points = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let from_slice = NeighbourGrid::new(&points, 4.0).unwrap();
    let from_iterator =
        NeighbourGrid::new((0..2).map(|i| [2.0 * f64::from(i), 0.0, 0.0]), 4.0).unwrap();
    let radius = Radius::new(2.0).unwrap();
    assert_eq!(
        from_slice.within([0.0; 3], radius),
        from_iterator.within([0.0; 3], radius)
    );
    assert_eq!(from_iterator.len(), 2);
    assert_eq!(
        NeighbourGrid::new([[0.0, 0.0, 0.0], [0.0, f64::NAN, 0.0]], 4.0).unwrap_err(),
        GeomError::NonFiniteCoordinate { index: 1 }
    );
}

#[test]
fn a_radius_is_checked_once() {
    assert_eq!(Radius::new(2.5).unwrap().get(), 2.5);
    assert_eq!(Radius::new(0.0).unwrap().get(), 0.0);
    assert_eq!(Radius::new(-1.0), Err(GeomError::InvalidRadius(-1.0)));
    assert_eq!(
        Radius::new(f64::NEG_INFINITY),
        Err(GeomError::InvalidRadius(f64::NEG_INFINITY))
    );
    assert!(matches!(Radius::new(f64::NAN), Err(GeomError::InvalidRadius(r)) if r.is_nan()));
    assert_eq!(Radius::INFINITE.get(), f64::INFINITY);
    assert_eq!(Radius::try_from(1.0), Radius::new(1.0));
}

#[test]
fn queries_do_not_fail() {
    let grid = NeighbourGrid::new([[0.0; 3], [1.0, 0.0, 0.0], [5.0, 0.0, 0.0]], 4.0).unwrap();
    let radius = Radius::new(1.5).unwrap();
    assert_eq!(grid.within([0.0; 3], radius), [0, 1]);
    assert_eq!(grid.pairs_within(radius), [(0, 1)]);
    assert_eq!(
        grid.nearest([4.0, 0.0, 0.0], Radius::INFINITE),
        Some((2, 1.0))
    );
    let mut seen = Vec::new();
    grid.for_each_within([0.0; 3], radius, |index, _| seen.push(index));
    seen.sort_unstable();
    assert_eq!(seen, [0, 1]);
}

#[test]
fn a_query_point_that_is_not_finite_finds_nothing() {
    let grid = NeighbourGrid::new([[0.0; 3], [1.0, 0.0, 0.0]], 4.0).unwrap();
    for point in [
        [f64::NAN, 0.0, 0.0],
        [f64::INFINITY, 0.0, 0.0],
        [0.0, f64::NEG_INFINITY, 0.0],
    ] {
        assert!(grid.within(point, Radius::INFINITE).is_empty());
        assert_eq!(grid.nearest(point, Radius::INFINITE), None);
        let mut called = false;
        grid.for_each_within(point, Radius::INFINITE, |_, _| called = true);
        assert!(!called);
    }
    let contacts = grid
        .contacts(
            [[f64::NAN, 0.0, 0.0], [0.5, 0.0, 0.0]],
            Radius::new(1.0).unwrap(),
        )
        .unwrap();
    assert_eq!(contacts, [(1, 0), (1, 1)]);
}

use std::hint::black_box;
use std::time::Duration;

use criterion::{Criterion, SamplingMode, criterion_group, criterion_main};
use mmcif_geom::NeighbourGrid;

const POINT_COUNT: usize = 100_000;
const BOX_EDGE: f64 = 100.0;
const RADIUS: f64 = 4.0;
const CELL_SIZE: f64 = 4.0;
const SEED: u64 = 0x6d6d_6369_6667_656f;

struct SplitMix64(u64);

impl SplitMix64 {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut mixed = self.0;
        mixed = (mixed ^ (mixed >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        mixed ^ (mixed >> 31)
    }

    fn next_unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn random_points() -> Vec<[f64; 3]> {
    let mut generator = SplitMix64(SEED);
    (0..POINT_COUNT)
        .map(|_| [(); 3].map(|()| generator.next_unit() * BOX_EDGE))
        .collect()
}

fn grid(c: &mut Criterion) {
    let points = random_points();
    let grid = NeighbourGrid::new(&points, CELL_SIZE).expect("generated points are finite");
    let mut group = c.benchmark_group("grid");
    group
        .sampling_mode(SamplingMode::Flat)
        .sample_size(10)
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(2));
    group.bench_function("build_100k", |b| {
        b.iter(|| NeighbourGrid::new(black_box(&points), CELL_SIZE))
    });
    group.bench_function("pairs_within_4A_100k", |b| {
        b.iter(|| grid.pairs_within(black_box(RADIUS)))
    });
    group.bench_function("within_4A_100k", |b| {
        b.iter(|| {
            points
                .iter()
                .map(|&point| {
                    grid.within(black_box(point), RADIUS)
                        .map_or(0, |found| found.len())
                })
                .sum::<usize>()
        })
    });
    group.finish();
}

criterion_group!(benches, grid);
criterion_main!(benches);

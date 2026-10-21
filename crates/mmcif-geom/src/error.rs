#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum GeomError {
    #[error("point {index} has a non-finite coordinate")]
    NonFiniteCoordinate { index: usize },
    #[error("query point has a non-finite coordinate")]
    NonFiniteQuery,
    #[error("cell size must be positive and finite, got {0}")]
    InvalidCellSize(f64),
    #[error("radius must be a non-negative number, got {0}")]
    InvalidRadius(f64),
    #[error("{0} points exceed the supported maximum of {max}", max = u32::MAX)]
    TooManyPoints(usize),
}

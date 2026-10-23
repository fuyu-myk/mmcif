/// An error from building or querying a [`NeighbourGrid`](crate::NeighbourGrid).
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum GeomError {
    /// A point passed to the grid has a NaN or infinite coordinate.
    #[error("point {index} has a non-finite coordinate")]
    NonFiniteCoordinate {
        /// The position of the point in the input.
        index: usize,
    },
    /// The cell size is zero, negative, NaN or infinite.
    #[error("cell size must be positive and finite, got {0}")]
    InvalidCellSize(f64),
    /// The radius is negative or NaN.
    #[error("radius must be a non-negative number, got {0}")]
    InvalidRadius(f64),
    /// More than `u32::MAX` points were given.
    #[error("{0} points exceed the supported maximum of {max}", max = u32::MAX)]
    TooManyPoints(usize),
}

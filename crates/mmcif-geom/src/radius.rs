use crate::check;
use crate::error::GeomError;

/// A query radius: zero, positive or infinite, but never negative or NaN.
///
/// A radius is checked once when it is made, so the queries of a
/// [`NeighbourGrid`](crate::NeighbourGrid) that take one cannot fail.
///
/// # Examples
///
/// ```
/// use mmcif_geom::{GeomError, Radius};
///
/// let radius = Radius::new(4.5)?;
/// assert_eq!(radius.get(), 4.5);
///
/// assert_eq!(Radius::new(-1.0), Err(GeomError::InvalidRadius(-1.0)));
/// assert_eq!(Radius::INFINITE.get(), f64::INFINITY);
/// # Ok::<(), GeomError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Radius(f64);

impl Radius {
    /// A radius that every point lies within.
    pub const INFINITE: Radius = Radius(f64::INFINITY);

    /// Creates a radius.
    ///
    /// # Errors
    ///
    /// Returns [`GeomError::InvalidRadius`] if `radius` is negative or NaN.
    pub fn new(radius: f64) -> Result<Radius, GeomError> {
        check::radius(radius).map(Radius)
    }

    /// Returns the radius.
    pub fn get(self) -> f64 {
        self.0
    }
}

impl TryFrom<f64> for Radius {
    type Error = GeomError;

    fn try_from(radius: f64) -> Result<Radius, GeomError> {
        Radius::new(radius)
    }
}

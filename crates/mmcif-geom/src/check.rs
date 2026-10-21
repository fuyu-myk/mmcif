use crate::error::GeomError;

pub(crate) fn cell_size(cell_size: f64) -> Result<f64, GeomError> {
    if cell_size > 0.0 && cell_size.is_finite() {
        Ok(cell_size)
    } else {
        Err(GeomError::InvalidCellSize(cell_size))
    }
}

pub(crate) fn radius(radius: f64) -> Result<f64, GeomError> {
    if radius >= 0.0 {
        Ok(radius)
    } else {
        Err(GeomError::InvalidRadius(radius))
    }
}

pub(crate) fn point_count(count: usize) -> Result<(), GeomError> {
    match u32::try_from(count) {
        Ok(_) => Ok(()),
        Err(_) => Err(GeomError::TooManyPoints(count)),
    }
}

pub(crate) fn finite_points(points: &[[f64; 3]]) -> Result<(), GeomError> {
    match points.iter().position(|&point| !is_finite(point)) {
        Some(index) => Err(GeomError::NonFiniteCoordinate { index }),
        None => Ok(()),
    }
}

pub(crate) fn query_point(point: [f64; 3]) -> Result<[f64; 3], GeomError> {
    if is_finite(point) {
        Ok(point)
    } else {
        Err(GeomError::NonFiniteQuery)
    }
}

fn is_finite(point: [f64; 3]) -> bool {
    point.iter().all(|value| value.is_finite())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_count_accepts_up_to_u32_max() {
        assert_eq!(point_count(0), Ok(()));
        assert_eq!(point_count(u32::MAX as usize), Ok(()));
    }

    #[test]
    fn point_count_rejects_more_than_u32_max() {
        let count = u32::MAX as usize + 1;
        assert_eq!(point_count(count), Err(GeomError::TooManyPoints(count)));
    }

    #[test]
    fn negative_zero_radius_is_accepted() {
        assert_eq!(radius(-0.0), Ok(-0.0));
    }
}

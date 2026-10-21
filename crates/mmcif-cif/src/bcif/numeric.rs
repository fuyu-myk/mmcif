use super::array::{Array, FloatWidth};
use super::{BcifError, Budget};

const MAX_DECIMALS: usize = 18;

fn invalid(encoding: &'static str, parameter: &'static str) -> BcifError {
    BcifError::InvalidParameter {
        encoding,
        parameter,
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Sentinels {
    upper: i64,
    lower: i64,
}

impl Sentinels {
    pub fn new(byte_count: i64, unsigned: bool) -> Result<Self, BcifError> {
        let (upper, lower) = match (byte_count, unsigned) {
            (1, false) => (i64::from(i8::MAX), i64::from(i8::MIN)),
            (2, false) => (i64::from(i16::MAX), i64::from(i16::MIN)),
            (1, true) => (i64::from(u8::MAX), i64::from(u8::MAX)),
            (2, true) => (i64::from(u16::MAX), i64::from(u16::MAX)),
            _ => return Err(invalid("IntegerPacking", "byteCount")),
        };
        Ok(Sentinels { upper, lower })
    }

    fn contains(self, value: i64) -> bool {
        value == self.upper || value == self.lower
    }
}

fn power_of_ten(factor: f64) -> Option<usize> {
    std::iter::successors(Some(1.0), |scale| Some(scale * 10.0))
        .take(MAX_DECIMALS + 1)
        .position(|scale| scale == factor)
}

pub(crate) fn fixed_point<'v>(
    values: Vec<i64>,
    factor: f64,
    width: FloatWidth,
) -> Result<Array<'v>, BcifError> {
    if !(factor.is_finite() && factor > 0.0) {
        return Err(invalid("FixedPoint", "factor"));
    }
    Ok(match power_of_ten(factor) {
        Some(decimals) => Array::Decimals(values, decimals),
        None => Array::Floats(
            values
                .into_iter()
                .map(|value| width.round(value as f64 / factor))
                .collect(),
            width,
        ),
    })
}

pub(crate) fn interval_quantization<'v>(
    values: Vec<i64>,
    min: f64,
    max: f64,
    steps: i64,
    width: FloatWidth,
) -> Result<Array<'v>, BcifError> {
    if !(min.is_finite() && max.is_finite()) {
        return Err(invalid("IntervalQuantization", "min"));
    }
    if steps < 2 {
        return Err(invalid("IntervalQuantization", "numSteps"));
    }
    let step = (max - min) / (steps - 1) as f64;
    Ok(Array::Floats(
        values
            .into_iter()
            .map(|value| width.round(min + step * value as f64))
            .collect(),
        width,
    ))
}

fn run_count(count: i64) -> Result<usize, BcifError> {
    usize::try_from(count).map_err(|_| invalid("RunLength", "count"))
}

pub(crate) fn run_length<'v>(
    pairs: Vec<i64>,
    size: usize,
    budget: &mut Budget,
) -> Result<Array<'v>, BcifError> {
    if pairs.len() % 2 != 0 {
        return Err(BcifError::OddRunLength(pairs.len()));
    }
    let mut total: usize = 0;
    for pair in pairs.chunks_exact(2) {
        total = total
            .checked_add(run_count(pair[1])?)
            .ok_or(BcifError::Overflow("RunLength"))?;
    }
    if total != size {
        return Err(BcifError::LengthMismatch {
            what: "RunLength output",
            expected: size,
            actual: total,
        });
    }
    budget.take(size)?;
    let mut output = Vec::with_capacity(size);
    for pair in pairs.chunks_exact(2) {
        output.extend(std::iter::repeat_n(pair[0], run_count(pair[1])?));
    }
    Ok(Array::Ints(output))
}

pub(crate) fn delta<'v>(mut values: Vec<i64>, origin: i64) -> Result<Array<'v>, BcifError> {
    let mut running = origin;
    for value in &mut values {
        running = running
            .checked_add(*value)
            .ok_or(BcifError::Overflow("Delta"))?;
        *value = running;
    }
    Ok(Array::Ints(values))
}

pub(crate) fn integer_packing<'v>(
    packed: Vec<i64>,
    sentinels: Sentinels,
    size: usize,
    budget: &mut Budget,
) -> Result<Array<'v>, BcifError> {
    if packed.len() == size {
        return Ok(Array::Ints(packed));
    }
    let count = packed
        .iter()
        .filter(|&&value| !sentinels.contains(value))
        .count();
    if count != size {
        return Err(BcifError::LengthMismatch {
            what: "IntegerPacking output",
            expected: size,
            actual: count,
        });
    }
    if packed
        .last()
        .is_some_and(|&value| sentinels.contains(value))
    {
        return Err(BcifError::UnterminatedPacking);
    }
    budget.take(size)?;
    let mut output = Vec::with_capacity(size);
    let mut sum: i64 = 0;
    for value in packed {
        sum = sum
            .checked_add(value)
            .ok_or(BcifError::Overflow("IntegerPacking"))?;
        if !sentinels.contains(value) {
            output.push(sum);
            sum = 0;
        }
    }
    Ok(Array::Ints(output))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcif::BcifLimits;

    fn budget() -> Budget {
        Budget::new(BcifLimits::default())
    }

    fn ints(array: Result<Array<'_>, BcifError>) -> Vec<i64> {
        match array.unwrap() {
            Array::Ints(values) => values,
            other => panic!("expected integers, got {other:?}"),
        }
    }

    #[test]
    fn scales_fixed_point_values() {
        assert_eq!(
            fixed_point(vec![12345, -500], 1000.0, FloatWidth::Single),
            Ok(Array::Decimals(vec![12345, -500], 3))
        );
        assert_eq!(
            fixed_point(vec![7], 1.0, FloatWidth::Double),
            Ok(Array::Decimals(vec![7], 0))
        );
        assert_eq!(
            fixed_point(vec![1, -2], 4.0, FloatWidth::Double),
            Ok(Array::Floats(vec![0.25, -0.5], FloatWidth::Double))
        );
        assert_eq!(
            fixed_point(vec![1], 3.0, FloatWidth::Single),
            Ok(Array::Floats(
                vec![f64::from(1.0f32 / 3.0)],
                FloatWidth::Single
            ))
        );
        for factor in [0.0, -10.0, f64::NAN, f64::INFINITY] {
            assert!(fixed_point(vec![1], factor, FloatWidth::Double).is_err());
        }
    }

    #[test]
    fn quantizes_intervals() {
        assert_eq!(
            interval_quantization(vec![0, 1, 4], 0.0, 1.0, 5, FloatWidth::Double),
            Ok(Array::Floats(vec![0.0, 0.25, 1.0], FloatWidth::Double))
        );
        assert_eq!(
            interval_quantization(vec![2], -1.0, 1.0, 3, FloatWidth::Single),
            Ok(Array::Floats(vec![1.0], FloatWidth::Single))
        );
        assert!(interval_quantization(vec![0], 0.0, 1.0, 1, FloatWidth::Double).is_err());
        assert!(interval_quantization(vec![0], f64::NAN, 1.0, 2, FloatWidth::Double).is_err());
    }

    #[test]
    fn expands_runs() {
        assert_eq!(
            ints(run_length(vec![7, 2, -1, 3, 4, 0], 5, &mut budget())),
            [7, 7, -1, -1, -1]
        );
        assert!(ints(run_length(Vec::new(), 0, &mut budget())).is_empty());
        assert_eq!(
            run_length(vec![7, 2], 3, &mut budget()),
            Err(BcifError::LengthMismatch {
                what: "RunLength output",
                expected: 3,
                actual: 2
            })
        );
        assert_eq!(
            run_length(vec![7, 2, 1], 2, &mut budget()),
            Err(BcifError::OddRunLength(3))
        );
        assert!(matches!(
            run_length(vec![7, -1], 0, &mut budget()),
            Err(BcifError::InvalidParameter {
                parameter: "count",
                ..
            })
        ));
        assert_eq!(
            run_length(vec![1, i64::MAX, 1, i64::MAX, 1, 2], 0, &mut budget()),
            Err(BcifError::Overflow("RunLength"))
        );
    }

    #[test]
    fn checks_run_length_size_against_the_limit_before_allocating() {
        let mut small = Budget::new(BcifLimits {
            max_elements: 1000,
            ..BcifLimits::default()
        });
        assert_eq!(
            run_length(vec![1, 1 << 40], 1 << 40, &mut small),
            Err(BcifError::LimitExceeded { limit: 1000 })
        );
    }

    #[test]
    fn accumulates_deltas_from_the_origin() {
        assert_eq!(ints(delta(vec![0, 1, 1, -3], 100)), [100, 101, 102, 99]);
        assert_eq!(ints(delta(vec![5], -5)), [0]);
        assert!(ints(delta(Vec::new(), 9)).is_empty());
        assert_eq!(delta(vec![1], i64::MAX), Err(BcifError::Overflow("Delta")));
    }

    fn unpack(
        packed: Vec<i64>,
        byte_count: i64,
        unsigned: bool,
        size: usize,
    ) -> Result<Array<'static>, BcifError> {
        let sentinels = Sentinels::new(byte_count, unsigned)?;
        integer_packing(packed, sentinels, size, &mut budget())
    }

    #[test]
    fn unpacks_signed_integers_with_sentinels() {
        let packed = vec![127, 127, 46, -128, -128, -44, 127, 0, -128, 0, 5];
        assert_eq!(ints(unpack(packed, 1, false, 5)), [300, -300, 127, -128, 5]);
        let packed = vec![32767, 1, -32768, -1, 7];
        assert_eq!(ints(unpack(packed, 2, false, 3)), [32768, -32769, 7]);
    }

    #[test]
    fn unpacks_unsigned_integers_with_sentinels() {
        assert_eq!(
            ints(unpack(vec![255, 255, 10, 255, 0, 3], 1, true, 3)),
            [520, 255, 3]
        );
        assert_eq!(ints(unpack(vec![65535, 5, 9], 2, true, 2)), [65540, 9]);
        assert_eq!(ints(unpack(vec![254, 255, 0], 1, true, 2)), [254, 255]);
    }

    #[test]
    fn keeps_data_that_already_has_the_declared_size() {
        assert_eq!(ints(unpack(vec![127, -128], 1, false, 2)), [127, -128]);
    }

    #[test]
    fn rejects_malformed_packing() {
        assert_eq!(
            unpack(vec![1, 2, 3], 1, false, 2),
            Err(BcifError::LengthMismatch {
                what: "IntegerPacking output",
                expected: 2,
                actual: 3
            })
        );
        assert_eq!(
            unpack(vec![5, 127, 127], 1, false, 1),
            Err(BcifError::UnterminatedPacking)
        );
        assert!(matches!(
            unpack(vec![1], 4, false, 1),
            Err(BcifError::InvalidParameter {
                parameter: "byteCount",
                ..
            })
        ));
    }
}

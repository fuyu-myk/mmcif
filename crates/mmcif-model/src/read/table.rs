use std::str::FromStr;

use mmcif_cif::{Block, Category, Column, TypedValue, Value};

use crate::error::ReadError;

pub(crate) struct Table<'a> {
    name: &'static str,
    category: Category<'a>,
}

impl<'a> Table<'a> {
    pub fn optional(block: &Block<'a>, name: &'static str) -> Option<Table<'a>> {
        block
            .category(name)
            .map(|category| Table { name, category })
    }

    pub fn required(block: &Block<'a>, name: &'static str) -> Result<Table<'a>, ReadError> {
        Table::optional(block, name).ok_or(ReadError::MissingCategory(name))
    }

    pub fn rows(&self) -> usize {
        self.category.row_count()
    }

    pub fn column(&self, item: &'static str) -> Field<'a> {
        Field {
            category: self.name,
            item,
            column: self.category.column(item),
        }
    }

    pub fn required_column(&self, item: &'static str) -> Result<Field<'a>, ReadError> {
        let field = self.column(item);
        if field.column.is_none() {
            return Err(ReadError::MissingItem {
                category: self.name,
                item,
            });
        }
        Ok(field)
    }
}

#[derive(Clone, Copy)]
pub(crate) struct Field<'a> {
    category: &'static str,
    item: &'static str,
    column: Option<Column<'a>>,
}

impl<'a> Field<'a> {
    pub fn is_present(&self) -> bool {
        self.column.is_some()
    }

    fn value(&self, row: usize) -> Option<Value<'a>> {
        self.column.and_then(|column| column.get(row))
    }

    pub fn invalid(&self, row: usize, value: &str, problem: &'static str) -> ReadError {
        ReadError::InvalidValue {
            category: self.category,
            item: self.item,
            row,
            value: value.to_owned(),
            problem,
        }
    }

    pub fn text(&self, row: usize) -> Option<&'a str> {
        self.value(row).and_then(|value| value.as_str())
    }

    pub fn required_text(&self, row: usize) -> Result<&'a str, ReadError> {
        self.text(row).ok_or_else(|| {
            let raw = self.value(row).map_or("", |value| value.raw());
            self.invalid(row, raw, "missing value")
        })
    }

    fn typed<T: Number>(&self, row: usize) -> Option<Result<Option<T>, ReadError>> {
        let typed = self.column?.typed(row)?;
        let invalid = |shown: String| self.invalid(row, &shown, "not a valid number");
        Some(match typed {
            TypedValue::Inapplicable | TypedValue::Unknown => Ok(None),
            TypedValue::Int(value) => T::from_int(value)
                .map(Some)
                .ok_or_else(|| invalid(value.to_string())),
            TypedValue::Float(value) => T::from_float(value)
                .map(Some)
                .ok_or_else(|| invalid(value.to_string())),
        })
    }

    pub fn parse<T: Number>(&self, row: usize) -> Result<Option<T>, ReadError> {
        if let Some(typed) = self.typed(row) {
            return typed;
        }
        self.text(row)
            .map(|text| {
                text.parse()
                    .map_err(|_| self.invalid(row, text, "not a valid number"))
            })
            .transpose()
    }

    pub fn required_parse<T: Number>(&self, row: usize) -> Result<T, ReadError> {
        match self.typed(row) {
            Some(Ok(Some(value))) => Ok(value),
            Some(Ok(None)) => Err(self.invalid(row, "", "missing value")),
            Some(Err(error)) => Err(error),
            None => {
                let text = self.required_text(row)?;
                text.parse()
                    .map_err(|_| self.invalid(row, text, "not a valid number"))
            }
        }
    }

    pub fn single_char(&self, row: usize) -> Result<Option<char>, ReadError> {
        let Some(text) = self.text(row) else {
            return Ok(None);
        };
        let mut chars = text.chars();
        match (chars.next(), chars.next()) {
            (Some(code), None) => Ok(Some(code)),
            _ => Err(self.invalid(row, text, "expected a single character")),
        }
    }
}

pub(crate) trait Number: FromStr + Sized {
    fn from_int(value: i64) -> Option<Self>;
    fn from_float(value: f64) -> Option<Self>;
}

fn integral(value: f64) -> Option<i128> {
    (value.fract() == 0.0 && value.abs() < 1e30).then_some(value as i128)
}

macro_rules! integer_number {
    ($($ty:ty),*) => {
        $(impl Number for $ty {
            fn from_int(value: i64) -> Option<Self> {
                Self::try_from(value).ok()
            }

            fn from_float(value: f64) -> Option<Self> {
                if Self::MIN == 0 && value.is_sign_negative() {
                    return None;
                }
                integral(value).and_then(|value| Self::try_from(value).ok())
            }
        })*
    };
}

integer_number!(i8, i32, u32, i64);

impl Number for f64 {
    fn from_int(value: i64) -> Option<Self> {
        Some(value as f64)
    }

    fn from_float(value: f64) -> Option<Self> {
        Some(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn text_parse<T: FromStr>(shown: &str) -> Option<T> {
        shown.parse().ok()
    }

    #[test]
    fn typed_conversions_agree_with_parsing_the_text() {
        for value in [
            5.0,
            -0.0,
            0.5,
            1e20,
            f64::NAN,
            f64::INFINITY,
            -3.0,
            4294967295.0,
            4294967296.0,
        ] {
            let shown = value.to_string();
            assert_eq!(i32::from_float(value), text_parse::<i32>(&shown), "{shown}");
            assert_eq!(u32::from_float(value), text_parse::<u32>(&shown), "{shown}");
            assert_eq!(i64::from_float(value), text_parse::<i64>(&shown), "{shown}");
        }
        for value in [0i64, -1, 300, i64::MAX, 1 << 53, (1 << 53) + 1] {
            let shown = value.to_string();
            assert_eq!(i8::from_int(value), text_parse::<i8>(&shown), "{shown}");
            assert_eq!(u32::from_int(value), text_parse::<u32>(&shown), "{shown}");
            assert_eq!(
                f64::from_int(value).map(f64::to_bits),
                text_parse::<f64>(&shown).map(f64::to_bits)
            );
        }
    }
}

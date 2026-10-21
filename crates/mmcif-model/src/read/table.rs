use std::str::FromStr;

use mmcif_cif::{Block, Category, Column, Value};

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

    pub fn parse<T: FromStr>(&self, row: usize) -> Result<Option<T>, ReadError> {
        self.text(row)
            .map(|text| {
                text.parse()
                    .map_err(|_| self.invalid(row, text, "not a valid number"))
            })
            .transpose()
    }

    pub fn required_parse<T: FromStr>(&self, row: usize) -> Result<T, ReadError> {
        let text = self.required_text(row)?;
        text.parse()
            .map_err(|_| self.invalid(row, text, "not a valid number"))
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

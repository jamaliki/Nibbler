//! Borrowed, case-normalized category rows over one generic CIF block.

use std::collections::BTreeMap;
use std::rc::Rc;

use crate::cif::numeric::{parse_float, parse_integer};
use crate::cif::{CifBlock, CifEntry, CifRow, CifValueRef};

pub(crate) struct Row<'a> {
    values: RowValues<'a>,
}

enum RowValues<'a> {
    Scalar(BTreeMap<String, CifValueRef<'a>>),
    Loop {
        row: CifRow<'a>,
        columns: Rc<BTreeMap<String, usize>>,
    },
}

impl<'a> Row<'a> {
    fn value(&self, item: &str) -> Option<CifValueRef<'a>> {
        match &self.values {
            RowValues::Scalar(values) => values.get(item).copied(),
            RowValues::Loop { row, columns } => {
                columns.get(item).and_then(|column| row.get(*column))
            }
        }
    }

    pub fn text(&self, item: &str) -> Option<String> {
        match self.value(item)? {
            CifValueRef::Text(text) => Some(text.as_str().to_owned()),
            CifValueRef::Integer(number, _) => Some(number.to_string()),
            CifValueRef::Float(number, _, _) => Some(number.to_string()),
            CifValueRef::Unknown | CifValueRef::NotApplicable => None,
        }
    }

    pub fn integer(&self, item: &str) -> Result<Option<i64>, String> {
        let Some(text) = self.text(item) else {
            return Ok(None);
        };
        parse_integer(&text).map(Some).ok_or(text)
    }

    pub fn float(&self, item: &str) -> Result<Option<f64>, String> {
        let Some(text) = self.text(item) else {
            return Ok(None);
        };
        parse_float(&text).map(Some).ok_or(text)
    }
}

pub(crate) struct CategoryIndex<'a> {
    rows: BTreeMap<String, Vec<Row<'a>>>,
}

impl<'a> CategoryIndex<'a> {
    pub fn new(block: &'a CifBlock) -> Self {
        let mut rows: BTreeMap<String, Vec<Row<'a>>> = BTreeMap::new();
        let mut scalars: BTreeMap<String, BTreeMap<String, CifValueRef<'a>>> = BTreeMap::new();
        for entry in block.entries() {
            match entry {
                CifEntry::Item(item) => {
                    if let Some((category, name)) = split_tag(item.tag()) {
                        scalars
                            .entry(case_key(category))
                            .or_default()
                            .insert(case_key(name), item.value().as_ref());
                    }
                }
                CifEntry::Loop(cif_loop) => {
                    let mut columns: BTreeMap<String, Vec<(usize, String)>> = BTreeMap::new();
                    for (column, tag) in cif_loop.tags().iter().enumerate() {
                        if let Some((category, name)) = split_tag(tag) {
                            columns
                                .entry(case_key(category))
                                .or_default()
                                .push((column, case_key(name)));
                        }
                    }
                    for (category, category_columns) in columns {
                        let category_rows = rows.entry(category).or_default();
                        let columns = Rc::new(
                            category_columns
                                .into_iter()
                                .map(|(column, item)| (item, column))
                                .collect(),
                        );
                        for row_index in 0..cif_loop.row_count() {
                            if let Some(row) = cif_loop.row(row_index) {
                                category_rows.push(Row {
                                    values: RowValues::Loop {
                                        row,
                                        columns: Rc::clone(&columns),
                                    },
                                });
                            }
                        }
                    }
                }
                CifEntry::Frame(_) => {}
            }
        }
        for (category, values) in scalars {
            rows.entry(category).or_default().push(Row {
                values: RowValues::Scalar(values),
            });
        }
        Self { rows }
    }

    pub fn rows(&self, category: &str) -> &[Row<'a>] {
        self.rows.get(category).map_or(&[], Vec::as_slice)
    }
}

pub(crate) fn split_tag(tag: &str) -> Option<(&str, &str)> {
    tag.strip_prefix('_')?.split_once('.')
}

pub(crate) fn case_key(value: &str) -> String {
    value.to_ascii_lowercase()
}

pub(crate) fn category_name(tag: &str) -> Option<&str> {
    split_tag(tag).map(|(category, _)| category)
}

pub(crate) fn item_name(tag: &str) -> Option<&str> {
    split_tag(tag).map(|(_, item)| item)
}

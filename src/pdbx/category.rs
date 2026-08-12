//! A single-pass category view over one generic CIF block.

use std::collections::BTreeMap;

use crate::cif::{CifBlock, CifEntry, CifValueRef};

pub(crate) struct Row<'a> {
    pub values: BTreeMap<String, CifValueRef<'a>>,
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
                            .entry(fold(category))
                            .or_default()
                            .insert(fold(name), item.value().as_ref());
                    }
                }
                CifEntry::Loop(cif_loop) => {
                    let mut columns: BTreeMap<String, Vec<(usize, String)>> = BTreeMap::new();
                    for (column, tag) in cif_loop.tags().iter().enumerate() {
                        if let Some((category, name)) = split_tag(tag) {
                            columns
                                .entry(fold(category))
                                .or_default()
                                .push((column, fold(name)));
                        }
                    }
                    for (category, category_columns) in columns {
                        let category_rows = rows.entry(category).or_default();
                        for row_index in 0..cif_loop.row_count() {
                            let Some(values) = cif_loop.row(row_index) else {
                                continue;
                            };
                            let values = category_columns
                                .iter()
                                .filter_map(|(column, name)| {
                                    values.get(*column).map(|value| (name.clone(), value))
                                })
                                .collect();
                            category_rows.push(Row { values });
                        }
                    }
                }
                CifEntry::Frame(_) => {}
            }
        }
        for (category, values) in scalars {
            rows.entry(category).or_default().push(Row { values });
        }
        Self { rows }
    }

    pub fn rows(&self, category: &str) -> &[Row<'a>] {
        self.rows.get(&fold(category)).map_or(&[], Vec::as_slice)
    }
}

fn split_tag(tag: &str) -> Option<(&str, &str)> {
    tag.strip_prefix('_')?.split_once('.')
}

fn fold(value: &str) -> String {
    value.to_ascii_lowercase()
}

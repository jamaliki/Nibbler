//! Missing-aware extraction of DDL2 records from a logical CIF document.

use super::dictionary::{DictionaryError, ItemRange};
use super::document::{CifEntry, CifValueRef};

pub(super) fn range_records(entries: &[CifEntry]) -> Vec<ItemRange> {
    let mut output = Vec::new();
    for entry in entries {
        let CifEntry::Loop(cif_loop) = entry else {
            continue;
        };
        let minimum = cif_loop
            .tags()
            .iter()
            .position(|tag| tag.eq_ignore_ascii_case("_item_range.minimum"));
        let maximum = cif_loop
            .tags()
            .iter()
            .position(|tag| tag.eq_ignore_ascii_case("_item_range.maximum"));
        let (Some(minimum), Some(maximum)) = (minimum, maximum) else {
            continue;
        };
        for row_index in 0..cif_loop.row_count() {
            let Some(row) = cif_loop.row(row_index) else {
                continue;
            };
            output.push(ItemRange {
                minimum: row.get(minimum).and_then(range_bound),
                maximum: row.get(maximum).and_then(range_bound),
            });
        }
    }
    if !output.is_empty() {
        return output;
    }

    let minimum = entries.iter().find_map(|entry| match entry {
        CifEntry::Item(item) if item.tag().eq_ignore_ascii_case("_item_range.minimum") => {
            Some(range_bound(item.value().as_ref()))
        }
        _ => None,
    });
    let maximum = entries.iter().find_map(|entry| match entry {
        CifEntry::Item(item) if item.tag().eq_ignore_ascii_case("_item_range.maximum") => {
            Some(range_bound(item.value().as_ref()))
        }
        _ => None,
    });
    match (minimum, maximum) {
        (Some(minimum), Some(maximum)) => vec![ItemRange { minimum, maximum }],
        _ => Vec::new(),
    }
}

fn range_bound(value: CifValueRef<'_>) -> Option<String> {
    match value {
        CifValueRef::Text(value) => Some(value.as_str().to_owned()),
        CifValueRef::Integer(value, _) => Some(value.to_string()),
        CifValueRef::Float(value, _, _) => Some(value.to_string()),
        CifValueRef::Unknown | CifValueRef::NotApplicable => None,
    }
}

pub(super) fn required_scalar<'a>(
    entries: &'a [CifEntry],
    tag: &str,
) -> Result<&'a str, DictionaryError> {
    scalar(entries, tag)
        .ok_or_else(|| DictionaryError::new(format!("missing required DDL2 item {tag}")))
}

pub(super) fn scalar<'a>(entries: &'a [CifEntry], tag: &str) -> Option<&'a str> {
    entries.iter().find_map(|entry| match entry {
        CifEntry::Item(item) if item.tag().eq_ignore_ascii_case(tag) => text(item.value().as_ref()),
        _ => None,
    })
}

pub(super) fn values<'a>(entries: &'a [CifEntry], tag: &str) -> Vec<&'a str> {
    let mut output = Vec::new();
    for entry in entries {
        match entry {
            CifEntry::Item(item) if item.tag().eq_ignore_ascii_case(tag) => {
                if let Some(value) = text(item.value().as_ref()) {
                    output.push(value);
                }
            }
            CifEntry::Loop(cif_loop) => {
                let Some(column) = cif_loop
                    .tags()
                    .iter()
                    .position(|candidate| candidate.eq_ignore_ascii_case(tag))
                else {
                    continue;
                };
                for row_index in 0..cif_loop.row_count() {
                    if let Some(value) = cif_loop
                        .row(row_index)
                        .and_then(|row| row.get(column))
                        .and_then(text)
                    {
                        output.push(value);
                    }
                }
            }
            _ => {}
        }
    }
    output
}

pub(super) fn item_records(entries: &[CifEntry]) -> Vec<Vec<&str>> {
    let mut output = Vec::new();
    for entry in entries {
        let CifEntry::Loop(cif_loop) = entry else {
            continue;
        };
        let Some(name_column) = cif_loop
            .tags()
            .iter()
            .position(|tag| tag.eq_ignore_ascii_case("_item.name"))
        else {
            continue;
        };
        let category_column = cif_loop
            .tags()
            .iter()
            .position(|tag| tag.eq_ignore_ascii_case("_item.category_id"));
        let mandatory_column = cif_loop
            .tags()
            .iter()
            .position(|tag| tag.eq_ignore_ascii_case("_item.mandatory_code"));
        for row_index in 0..cif_loop.row_count() {
            let Some(row) = cif_loop.row(row_index) else {
                continue;
            };
            let Some(name) = row.get(name_column).and_then(text) else {
                continue;
            };
            let category = category_column
                .and_then(|column| row.get(column))
                .and_then(text)
                .or_else(|| item_category(name));
            let mandatory = mandatory_column
                .and_then(|column| row.get(column))
                .and_then(text)
                .unwrap_or("no");
            if let Some(category) = category {
                output.push(vec![name, category, mandatory]);
            }
        }
    }
    if !output.is_empty() {
        return output;
    }

    let Some(name) = scalar(entries, "_item.name") else {
        return output;
    };
    if let Some(category) = scalar(entries, "_item.category_id").or_else(|| item_category(name)) {
        output.push(vec![
            name,
            category,
            scalar(entries, "_item.mandatory_code").unwrap_or("no"),
        ]);
    }
    output
}

fn item_category(name: &str) -> Option<&str> {
    name.strip_prefix('_')?
        .split_once('.')
        .map(|(category, _)| category)
}

pub(super) fn loop_records<'a>(entries: &'a [CifEntry], tags: &[&str]) -> Vec<Vec<&'a str>> {
    let mut output = Vec::new();
    for entry in entries {
        let CifEntry::Loop(cif_loop) = entry else {
            continue;
        };
        let Some(columns) = tags
            .iter()
            .map(|tag| {
                cif_loop
                    .tags()
                    .iter()
                    .position(|candidate| candidate.eq_ignore_ascii_case(tag))
            })
            .collect::<Option<Vec<_>>>()
        else {
            continue;
        };
        for row_index in 0..cif_loop.row_count() {
            let Some(row) = cif_loop.row(row_index) else {
                continue;
            };
            if let Some(record) = columns
                .iter()
                .map(|column| row.get(*column).and_then(text))
                .collect::<Option<Vec<_>>>()
            {
                output.push(record);
            }
        }
    }
    output
}

fn text(value: CifValueRef<'_>) -> Option<&str> {
    value.as_text()
}

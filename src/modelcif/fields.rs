//! Typed field decoding at the generic-CIF to ModelCIF trust boundary.

use std::borrow::Cow;

use crate::cif::CifValueRef;
use crate::cif::numeric::{parse_float, parse_integer};
use crate::pdbx::category::Row;

use super::source::ModelCifError;

pub(super) fn required_text(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<String, ModelCifError> {
    optional_text(row, item).ok_or_else(|| {
        ModelCifError::new(
            "MODELCIF_ITEM_REQUIRED",
            format!("required item _{category}.{item} is absent or missing"),
            row_context(category, item, row_index),
        )
    })
}

pub(super) fn optional_text(row: &Row<'_>, item: &str) -> Option<String> {
    row.values
        .get(&item.to_ascii_lowercase())
        .and_then(|value| present_text(value))
        .map(Cow::into_owned)
}

pub(super) fn required_integer(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<i64, ModelCifError> {
    let text = required_text(row, category, item, row_index)?;
    parse_integer(&text).ok_or_else(|| typed(category, item, row_index, &text, "integer"))
}

pub(super) fn optional_integer(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<Option<i64>, ModelCifError> {
    let Some(text) = optional_text(row, item) else {
        return Ok(None);
    };
    parse_integer(&text)
        .map(Some)
        .ok_or_else(|| typed(category, item, row_index, &text, "integer"))
}

pub(super) fn required_float(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<f64, ModelCifError> {
    let text = required_text(row, category, item, row_index)?;
    parse_float(&text).ok_or_else(|| typed(category, item, row_index, &text, "finite float"))
}

pub(super) fn optional_float(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<Option<f64>, ModelCifError> {
    let Some(text) = optional_text(row, item) else {
        return Ok(None);
    };
    parse_float(&text)
        .map(Some)
        .ok_or_else(|| typed(category, item, row_index, &text, "finite float"))
}

fn row_context(category: &str, item: &str, row_index: usize) -> Vec<String> {
    vec![
        format!("category={category}"),
        format!("item=_{category}.{item}"),
        format!("row={}", row_index + 1),
    ]
}

fn present_text<'a>(value: &CifValueRef<'a>) -> Option<Cow<'a, str>> {
    match *value {
        CifValueRef::Text(text) => Some(Cow::Borrowed(text.as_str())),
        CifValueRef::Integer(number, _) => Some(Cow::Owned(number.to_string())),
        CifValueRef::Float(number, _, _) => Some(Cow::Owned(number.to_string())),
        CifValueRef::Unknown | CifValueRef::NotApplicable => None,
    }
}

fn typed(
    category: &str,
    item: &str,
    row_index: usize,
    text: &str,
    expected: &str,
) -> ModelCifError {
    ModelCifError::new(
        "MODELCIF_ITEM_TYPE",
        format!("value {text:?} for _{category}.{item} is not a {expected}"),
        row_context(category, item, row_index),
    )
}

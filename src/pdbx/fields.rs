//! Typed field decoding at the generic-CIF to PDBx trust boundary.

use std::borrow::Cow;

use crate::cif::CifValueRef;
use crate::cif::numeric::{parse_float, parse_integer};

use super::category::Row;
use super::source::ModelError;

pub(super) fn required_single(
    rows: &[Row<'_>],
    category: &str,
    item: &str,
) -> Result<String, ModelError> {
    match rows {
        [row] => required_text(row, category, item, 0),
        [] => Err(ModelError::new(
            "PDBX_CATEGORY_REQUIRED",
            format!("required category {category:?} is absent"),
            vec![format!("category={category}")],
        )),
        _ => Err(ModelError::new(
            "PDBX_CATEGORY_SINGLETON",
            format!("category {category:?} must contain exactly one row"),
            vec![format!("category={category}")],
        )),
    }
}

pub(super) fn required_text(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<String, ModelError> {
    row.values
        .get(&item.to_ascii_lowercase())
        .and_then(|value| present_text(value))
        .map(Cow::into_owned)
        .ok_or_else(|| {
            ModelError::new(
                "PDBX_ITEM_REQUIRED",
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
) -> Result<i64, ModelError> {
    let text = required_text(row, category, item, row_index)?;
    parse_integer(&text).ok_or_else(|| typed(category, item, row_index, &text, "integer"))
}

pub(super) fn optional_integer(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<Option<i64>, ModelError> {
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
) -> Result<f64, ModelError> {
    let text = required_text(row, category, item, row_index)?;
    parse_float(&text).ok_or_else(|| typed(category, item, row_index, &text, "finite float"))
}

pub(super) fn optional_float(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<Option<f64>, ModelError> {
    let Some(text) = optional_text(row, item) else {
        return Ok(None);
    };
    parse_float(&text)
        .map(Some)
        .ok_or_else(|| typed(category, item, row_index, &text, "finite float"))
}

pub(super) fn optional_charge(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<Option<i32>, ModelError> {
    let Some(text) = optional_text(row, item) else {
        return Ok(None);
    };
    parse_charge(&text)
        .map(Some)
        .ok_or_else(|| typed(category, item, row_index, &text, "formal charge"))
}

pub(super) fn optional_bool(row: &Row<'_>, item: &str) -> Option<bool> {
    optional_text(row, item).and_then(|value| match value.to_ascii_lowercase().as_str() {
        "y" | "yes" => Some(true),
        "n" | "no" => Some(false),
        _ => None,
    })
}

pub(super) fn duplicate(category: &str, item: &str, value: &str, row_index: usize) -> ModelError {
    ModelError::new(
        "PDBX_ID_DUPLICATE",
        format!("value {value:?} for _{category}.{item} is not unique"),
        row_context(category, item, row_index),
    )
}

pub(super) fn row_context(category: &str, item: &str, row_index: usize) -> Vec<String> {
    vec![
        format!("category={category}"),
        format!("item=_{category}.{item}"),
        format!("row={}", row_index + 1),
    ]
}

fn parse_charge(text: &str) -> Option<i32> {
    let normalized = match text.as_bytes() {
        [digits @ .., b'+'] => format!("+{}", std::str::from_utf8(digits).ok()?),
        [digits @ .., b'-'] => format!("-{}", std::str::from_utf8(digits).ok()?),
        _ => text.to_owned(),
    };
    normalized.parse().ok()
}

fn present_text<'a>(value: &CifValueRef<'a>) -> Option<Cow<'a, str>> {
    match *value {
        CifValueRef::Text(text) => Some(Cow::Borrowed(text.as_str())),
        CifValueRef::Integer(number, _) => Some(Cow::Owned(number.to_string())),
        CifValueRef::Float(number, _, _) => Some(Cow::Owned(number.to_string())),
        CifValueRef::Unknown | CifValueRef::NotApplicable => None,
    }
}

fn typed(category: &str, item: &str, row_index: usize, text: &str, expected: &str) -> ModelError {
    ModelError::new(
        "PDBX_ITEM_TYPE",
        format!("value {text:?} for _{category}.{item} is not a {expected}"),
        row_context(category, item, row_index),
    )
}

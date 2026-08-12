//! Typed field decoding at the generic-CIF to ModelCIF trust boundary.

use crate::pdbx::SemanticError;
use crate::pdbx::category::Row;

pub(super) fn required_text(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<String, SemanticError> {
    row.text(item)
        .ok_or_else(|| missing(category, item, row_index))
}

pub(super) fn optional_text(row: &Row<'_>, item: &str) -> Option<String> {
    row.text(item)
}

pub(super) fn required_integer(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<i64, SemanticError> {
    match row.integer(item) {
        Ok(Some(value)) => Ok(value),
        Ok(None) => Err(missing(category, item, row_index)),
        Err(text) => Err(typed(category, item, row_index, &text, "integer")),
    }
}

pub(super) fn optional_integer(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<Option<i64>, SemanticError> {
    row.integer(item)
        .map_err(|text| typed(category, item, row_index, &text, "integer"))
}

pub(super) fn required_float(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<f64, SemanticError> {
    match row.float(item) {
        Ok(Some(value)) => Ok(value),
        Ok(None) => Err(missing(category, item, row_index)),
        Err(text) => Err(typed(category, item, row_index, &text, "finite float")),
    }
}

pub(super) fn optional_float(
    row: &Row<'_>,
    category: &str,
    item: &str,
    row_index: usize,
) -> Result<Option<f64>, SemanticError> {
    row.float(item)
        .map_err(|text| typed(category, item, row_index, &text, "finite float"))
}

fn row_context(category: &str, item: &str, row_index: usize) -> Vec<String> {
    vec![
        format!("category={category}"),
        format!("item=_{category}.{item}"),
        format!("row={}", row_index + 1),
    ]
}

fn typed(
    category: &str,
    item: &str,
    row_index: usize,
    text: &str,
    expected: &str,
) -> SemanticError {
    SemanticError::new(
        "MODELCIF_ITEM_TYPE",
        format!("value {text:?} for _{category}.{item} is not a {expected}"),
        row_context(category, item, row_index),
    )
}

fn missing(category: &str, item: &str, row_index: usize) -> SemanticError {
    SemanticError::new(
        "MODELCIF_ITEM_REQUIRED",
        format!("required item _{category}.{item} is absent or missing"),
        row_context(category, item, row_index),
    )
}

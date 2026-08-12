//! Typed field decoding at the generic-CIF to ModelCIF trust boundary.

use crate::pdbx::SemanticError;
use crate::pdbx::fields::Row;

pub(super) fn required_text<const N: usize>(
    row: &Row<'_, N>,
    category: &str,
    item: &str,
    field: usize,
) -> Result<String, SemanticError> {
    row.text(field)
        .ok_or_else(|| missing(category, item, row.row_index()))
}

pub(super) fn optional_text<const N: usize>(row: &Row<'_, N>, field: usize) -> Option<String> {
    row.text(field)
}

pub(super) fn required_integer<const N: usize>(
    row: &Row<'_, N>,
    category: &str,
    item: &str,
    field: usize,
) -> Result<i64, SemanticError> {
    match row.integer(field) {
        Ok(Some(value)) => Ok(value),
        Ok(None) => Err(missing(category, item, row.row_index())),
        Err(text) => Err(typed(category, item, row.row_index(), &text, "integer")),
    }
}

pub(super) fn optional_integer<const N: usize>(
    row: &Row<'_, N>,
    category: &str,
    item: &str,
    field: usize,
) -> Result<Option<i64>, SemanticError> {
    row.integer(field)
        .map_err(|text| typed(category, item, row.row_index(), &text, "integer"))
}

pub(super) fn required_float<const N: usize>(
    row: &Row<'_, N>,
    category: &str,
    item: &str,
    field: usize,
) -> Result<f64, SemanticError> {
    match row.float(field) {
        Ok(Some(value)) => Ok(value),
        Ok(None) => Err(missing(category, item, row.row_index())),
        Err(text) => Err(typed(
            category,
            item,
            row.row_index(),
            &text,
            "finite float",
        )),
    }
}

pub(super) fn optional_float<const N: usize>(
    row: &Row<'_, N>,
    category: &str,
    item: &str,
    field: usize,
) -> Result<Option<f64>, SemanticError> {
    row.float(field)
        .map_err(|text| typed(category, item, row.row_index(), &text, "finite float"))
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

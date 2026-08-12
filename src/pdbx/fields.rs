//! Lazy category rows and typed decoding at the generic-CIF trust boundary.

use crate::cif::numeric::{parse_float, parse_integer};
use crate::cif::{CategoryColumn, CategoryOccurrence, CategoryView, CifValueRef, split_tag};

use super::error::SemanticError;

#[derive(Clone, Copy)]
enum Occurrence<'a> {
    Scalars(CategoryView<'a>),
    Loop(CategoryOccurrence<'a>),
}

impl<'a> Occurrence<'a> {
    fn row_count(self) -> usize {
        match self {
            Self::Scalars(_) => 1,
            Self::Loop(occurrence) => occurrence.row_count(),
        }
    }

    fn column(self, category: &str, item: &str) -> Option<CategoryColumn<'a>> {
        let matches = |tag: &str| {
            split_tag(tag).is_some_and(|(tag_category, tag_item)| {
                tag_category.eq_ignore_ascii_case(category) && tag_item.eq_ignore_ascii_case(item)
            })
        };
        match self {
            Self::Scalars(view) => view
                .occurrences()
                .filter(|occurrence| occurrence.is_scalar())
                .flat_map(CategoryOccurrence::columns)
                .find(|column| matches(column.tag())),
            Self::Loop(occurrence) => occurrence.columns().find(|column| matches(column.tag())),
        }
    }
}

/// One logical category row with its requested columns resolved once per occurrence.
#[derive(Clone, Copy)]
pub(crate) struct Row<'a, const N: usize> {
    columns: [Option<CategoryColumn<'a>>; N],
    source_row: usize,
    logical_row: usize,
}

impl<const N: usize> Row<'_, N> {
    pub(crate) const fn row_index(&self) -> usize {
        self.logical_row
    }

    pub(crate) fn text(&self, field: usize) -> Option<String> {
        match self
            .columns
            .get(field)
            .copied()
            .flatten()?
            .value(self.source_row)?
        {
            CifValueRef::Text(text) => Some(text.as_str().to_owned()),
            CifValueRef::Integer(number, _) => Some(number.to_string()),
            CifValueRef::Float(number, _, _) => Some(number.to_string()),
            CifValueRef::Unknown | CifValueRef::NotApplicable => None,
        }
    }

    pub(crate) fn integer(&self, field: usize) -> Result<Option<i64>, String> {
        let Some(text) = self.text(field) else {
            return Ok(None);
        };
        parse_integer(&text).map(Some).ok_or(text)
    }

    pub(crate) fn float(&self, field: usize) -> Result<Option<f64>, String> {
        let Some(text) = self.text(field) else {
            return Ok(None);
        };
        parse_float(&text).map(Some).ok_or(text)
    }
}

pub(crate) fn category_rows<'a, const N: usize>(
    view: CategoryView<'a>,
    category: &'a str,
    items: [&'a str; N],
) -> impl Iterator<Item = Row<'a, N>> {
    let has_category = move |occurrence: CategoryOccurrence<'a>| {
        occurrence.columns().any(|column| {
            split_tag(column.tag()).is_some_and(|(name, _)| name.eq_ignore_ascii_case(category))
        })
    };
    view.occurrences()
        .scan(false, move |scalar_emitted, occurrence| {
            let selected = if !has_category(occurrence) {
                None
            } else if occurrence.is_scalar() {
                (!std::mem::replace(scalar_emitted, true)).then_some(Occurrence::Scalars(view))
            } else {
                Some(Occurrence::Loop(occurrence))
            };
            Some(selected)
        })
        .flatten()
        .flat_map(move |occurrence| {
            let columns = items.map(|item| occurrence.column(category, item));
            (0..occurrence.row_count()).map(move |source_row| Row {
                columns,
                source_row,
                logical_row: 0,
            })
        })
        .enumerate()
        .map(|(logical_row, mut row)| {
            row.logical_row = logical_row;
            row
        })
}

pub(super) fn required_single(
    view: CategoryView<'_>,
    category: &str,
    item: &str,
) -> Result<String, SemanticError> {
    let mut rows = category_rows(view, category, [item]);
    match (rows.next(), rows.next()) {
        (Some(row), None) => required_text(&row, category, item, 0),
        (None, _) => Err(SemanticError::new(
            "PDBX_CATEGORY_REQUIRED",
            format!("required category {category:?} is absent"),
            vec![format!("category={category}")],
        )),
        (Some(_), Some(_)) => Err(SemanticError::new(
            "PDBX_CATEGORY_SINGLETON",
            format!("category {category:?} must contain exactly one row"),
            vec![format!("category={category}")],
        )),
    }
}

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

pub(super) fn optional_charge<const N: usize>(
    row: &Row<'_, N>,
    category: &str,
    item: &str,
    field: usize,
) -> Result<Option<i32>, SemanticError> {
    let Some(text) = optional_text(row, field) else {
        return Ok(None);
    };
    parse_charge(&text)
        .map(Some)
        .ok_or_else(|| typed(category, item, row.row_index(), &text, "formal charge"))
}

pub(super) fn optional_bool<const N: usize>(row: &Row<'_, N>, field: usize) -> Option<bool> {
    optional_text(row, field).and_then(|value| match value.to_ascii_lowercase().as_str() {
        "y" | "yes" => Some(true),
        "n" | "no" => Some(false),
        _ => None,
    })
}

pub(super) fn duplicate(
    category: &str,
    item: &str,
    value: &str,
    row_index: usize,
) -> SemanticError {
    SemanticError::new(
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

pub(crate) fn case_key(value: &str) -> String {
    value.to_ascii_lowercase()
}

fn parse_charge(text: &str) -> Option<i32> {
    let normalized = match text.as_bytes() {
        [digits @ .., b'+'] => format!("+{}", std::str::from_utf8(digits).ok()?),
        [digits @ .., b'-'] => format!("-{}", std::str::from_utf8(digits).ok()?),
        _ => text.to_owned(),
    };
    normalized.parse().ok()
}

fn typed(
    category: &str,
    item: &str,
    row_index: usize,
    text: &str,
    expected: &str,
) -> SemanticError {
    SemanticError::new(
        "PDBX_ITEM_TYPE",
        format!("value {text:?} for _{category}.{item} is not a {expected}"),
        row_context(category, item, row_index),
    )
}

fn missing(category: &str, item: &str, row_index: usize) -> SemanticError {
    SemanticError::new(
        "PDBX_ITEM_REQUIRED",
        format!("required item _{category}.{item} is absent or missing"),
        row_context(category, item, row_index),
    )
}

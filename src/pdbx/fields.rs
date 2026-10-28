//! Lazy category rows and typed decoding at the generic-CIF trust boundary.

use std::borrow::Cow;
use std::rc::Rc;

use crate::cif::numeric::{parse_float, parse_integer};
use crate::cif::{
    CategoryColumn, CategoryOccurrence, CategoryView, CifValueRef, TextColumn, split_tag,
};

use super::error::SemanticError;
use super::short_text::ShortText;

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

/// One requested column. Text columns of columnar loops are read directly.
#[derive(Clone, Copy)]
enum Field<'a> {
    Absent,
    Text(TextColumn<'a>),
    Column(CategoryColumn<'a>),
}

impl<'a> Field<'a> {
    fn new(column: Option<CategoryColumn<'a>>) -> Self {
        match column {
            None => Self::Absent,
            Some(column) => column
                .text_column()
                .map_or(Self::Column(column), Self::Text),
        }
    }

    #[inline]
    fn value(self, row: usize) -> Option<CifValueRef<'a>> {
        match self {
            Self::Absent => None,
            Self::Text(column) => column.value(row),
            Self::Column(column) => column.value(row),
        }
    }

    /// Return a key that two rows share only if their values are equal, or `None` when
    /// the storage gives none.
    #[inline]
    fn key(self, row: usize) -> Option<u64> {
        match self {
            Self::Text(column) => column.key(row),
            Self::Absent | Self::Column(_) => None,
        }
    }
}

/// Return present text, borrowed when the document stores the value as text.
fn value_str(value: Option<CifValueRef<'_>>) -> Option<Cow<'_, str>> {
    match value? {
        CifValueRef::Text(text) => Some(Cow::Borrowed(text.as_str())),
        CifValueRef::Integer(number, _) => Some(Cow::Owned(number.to_string())),
        CifValueRef::Float(number, _, _) => Some(Cow::Owned(number.to_string())),
        CifValueRef::Unknown | CifValueRef::NotApplicable => None,
    }
}

fn value_integer(value: Option<CifValueRef<'_>>) -> Result<Option<i64>, String> {
    // A typed integer prints and parses back to itself, so it skips the text round trip.
    if let Some(CifValueRef::Integer(number, _)) = value {
        return Ok(Some(number));
    }
    let Some(text) = value_str(value) else {
        return Ok(None);
    };
    parse_integer(&text)
        .map(Some)
        .ok_or_else(|| text.into_owned())
}

fn value_float(value: Option<CifValueRef<'_>>) -> Result<Option<f64>, String> {
    // `f64` display is the shortest text that parses back to the same value, so a finite
    // typed float skips the text round trip.
    if let Some(CifValueRef::Float(number, _, _)) = value
        && number.is_finite()
    {
        return Ok(Some(number));
    }
    let Some(text) = value_str(value) else {
        return Ok(None);
    };
    parse_float(&text)
        .map(Some)
        .ok_or_else(|| text.into_owned())
}

/// One logical category row with its requested columns resolved once per occurrence.
///
/// Rows of one occurrence share their column array, so advancing to the next row copies
/// a pointer instead of `N` column descriptions.
#[derive(Clone)]
pub(crate) struct Row<'a, const N: usize> {
    columns: Rc<[Field<'a>; N]>,
    source_row: usize,
    logical_row: usize,
}

impl<'a, const N: usize> Row<'a, N> {
    pub(crate) const fn row_index(&self) -> usize {
        self.logical_row
    }

    fn value(&self, field: usize) -> Option<CifValueRef<'a>> {
        self.columns.get(field)?.value(self.source_row)
    }

    /// Return present text, borrowed when the document stores the value as text.
    pub(crate) fn str(&self, field: usize) -> Option<Cow<'a, str>> {
        value_str(self.value(field))
    }

    pub(crate) fn text(&self, field: usize) -> Option<String> {
        self.str(field).map(Cow::into_owned)
    }

    pub(crate) fn integer(&self, field: usize) -> Result<Option<i64>, String> {
        value_integer(self.value(field))
    }

    pub(crate) fn float(&self, field: usize) -> Result<Option<f64>, String> {
        value_float(self.value(field))
    }
}

/// Return the occurrences that hold `category`: each loop, and the block's scalar items
/// once, at the position of the first scalar occurrence.
fn occurrences<'a>(
    view: CategoryView<'a>,
    category: &'a str,
) -> impl Iterator<Item = Occurrence<'a>> {
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
}

pub(crate) fn category_rows<'a, const N: usize>(
    view: CategoryView<'a>,
    category: &'a str,
    items: [&'a str; N],
) -> impl Iterator<Item = Row<'a, N>> {
    occurrences(view, category)
        .flat_map(move |occurrence| {
            let columns = Rc::new(items.map(|item| Field::new(occurrence.column(category, item))));
            (0..occurrence.row_count()).map(move |source_row| Row {
                columns: Rc::clone(&columns),
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

/// One occurrence of a category, read column by column.
pub(crate) struct Table<'a, const N: usize> {
    fields: [Field<'a>; N],
    rows: usize,
}

impl<'a, const N: usize> Table<'a, N> {
    pub(crate) const fn rows(&self) -> usize {
        self.rows
    }

    /// Return a reader for requested column `field`.
    pub(crate) fn column<T>(&self, field: usize) -> Column<'a, T> {
        Column {
            field: self.fields[field],
            last: None,
        }
    }
}

/// Return the occurrences of `category` in the order [`category_rows`] reads them.
pub(crate) fn category_tables<'a, const N: usize>(
    view: CategoryView<'a>,
    category: &'a str,
    items: [&'a str; N],
) -> Vec<Table<'a, N>> {
    occurrences(view, category)
        .map(|occurrence| Table {
            fields: items.map(|item| Field::new(occurrence.column(category, item))),
            rows: occurrence.row_count(),
        })
        .collect()
}

/// Decodes one column of a [`Table`] row by row. When the storage shows that a row holds
/// the same value as the row decoded last, the reader returns a copy of that result.
pub(crate) struct Column<'a, T> {
    field: Field<'a>,
    last: Option<(u64, T)>,
}

impl<'a, T: Clone> Column<'a, T> {
    // Most rows repeat the previous value, so the check stays inline and decoding does not.
    #[inline(always)]
    fn get<E>(
        &mut self,
        row: usize,
        decode: impl FnOnce(Option<CifValueRef<'a>>) -> Result<T, E>,
    ) -> Result<T, E> {
        let key = self.field.key(row);
        if let (Some(key), Some((last, value))) = (key, &self.last)
            && *last == key
        {
            return Ok(value.clone());
        }
        self.decode(row, key, decode)
    }

    #[inline(never)]
    fn decode<E>(
        &mut self,
        row: usize,
        key: Option<u64>,
        decode: impl FnOnce(Option<CifValueRef<'a>>) -> Result<T, E>,
    ) -> Result<T, E> {
        let value = decode(self.field.value(row))?;
        if let Some(key) = key {
            self.last = Some((key, value.clone()));
        }
        Ok(value)
    }
}

impl Column<'_, Option<ShortText>> {
    #[inline(always)]
    pub(crate) fn optional_short(&mut self, row: usize) -> Option<ShortText> {
        self.get::<()>(row, |value| {
            Ok(value_str(value).map(|text| ShortText::new(&text)))
        })
        .unwrap_or_default()
    }

    /// Return the text in `row`; `logical_row` numbers the row in error contexts.
    #[inline(always)]
    pub(crate) fn required_short(
        &mut self,
        row: usize,
        logical_row: usize,
        category: &str,
        item: &str,
    ) -> Result<ShortText, SemanticError> {
        self.optional_short(row)
            .ok_or_else(|| missing(category, item, logical_row))
    }
}

impl Column<'_, Option<i64>> {
    #[inline(always)]
    pub(crate) fn optional_integer(
        &mut self,
        row: usize,
        logical_row: usize,
        category: &str,
        item: &str,
    ) -> Result<Option<i64>, SemanticError> {
        self.get(row, value_integer)
            .map_err(|text| typed(category, item, logical_row, &text, "integer"))
    }

    #[inline(always)]
    pub(crate) fn required_integer(
        &mut self,
        row: usize,
        logical_row: usize,
        category: &str,
        item: &str,
    ) -> Result<i64, SemanticError> {
        self.optional_integer(row, logical_row, category, item)?
            .ok_or_else(|| missing(category, item, logical_row))
    }
}

impl Column<'_, Option<f64>> {
    #[inline(always)]
    pub(crate) fn optional_float(
        &mut self,
        row: usize,
        logical_row: usize,
        category: &str,
        item: &str,
    ) -> Result<Option<f64>, SemanticError> {
        self.get(row, value_float)
            .map_err(|text| typed(category, item, logical_row, &text, "finite float"))
    }

    #[inline(always)]
    pub(crate) fn required_float(
        &mut self,
        row: usize,
        logical_row: usize,
        category: &str,
        item: &str,
    ) -> Result<f64, SemanticError> {
        self.optional_float(row, logical_row, category, item)?
            .ok_or_else(|| missing(category, item, logical_row))
    }
}

impl Column<'_, Option<i32>> {
    #[inline(always)]
    pub(crate) fn optional_charge(
        &mut self,
        row: usize,
        logical_row: usize,
        category: &str,
        item: &str,
    ) -> Result<Option<i32>, SemanticError> {
        self.get(row, |value| {
            let Some(text) = value_str(value) else {
                return Ok(None);
            };
            parse_charge(&text)
                .map(Some)
                .ok_or_else(|| typed(category, item, logical_row, &text, "formal charge"))
        })
    }
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
    let Some(text) = row.str(field) else {
        return Ok(None);
    };
    parse_charge(&text)
        .map(Some)
        .ok_or_else(|| typed(category, item, row.row_index(), &text, "formal charge"))
}

pub(super) fn optional_bool<const N: usize>(row: &Row<'_, N>, field: usize) -> Option<bool> {
    let value = row.str(field)?;
    if value.eq_ignore_ascii_case("y") || value.eq_ignore_ascii_case("yes") {
        Some(true)
    } else if value.eq_ignore_ascii_case("n") || value.eq_ignore_ascii_case("no") {
        Some(false)
    } else {
        None
    }
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
        _ => return text.parse().ok(),
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

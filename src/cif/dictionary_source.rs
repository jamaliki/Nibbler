//! Borrowed column views over generic CIF entries.

use super::dictionary::{DictionaryError, ItemRange};
use super::document::{CifEntry, CifItem, CifLoop, CifValueRef};

/// A zero-copy view of the scalar and loop occurrences in one block or save frame.
#[derive(Clone, Copy)]
pub(crate) struct CategoryView<'a>(&'a [CifEntry]);

impl<'a> CategoryView<'a> {
    pub(crate) const fn new(entries: &'a [CifEntry]) -> Self {
        Self(entries)
    }

    pub(crate) fn occurrences(self) -> impl Iterator<Item = CategoryOccurrence<'a>> {
        self.0.iter().filter_map(|entry| match entry {
            CifEntry::Item(item) => Some(CategoryOccurrence::Scalar(item)),
            CifEntry::Loop(cif_loop) => Some(CategoryOccurrence::Loop(cif_loop)),
            CifEntry::Frame(_) => None,
        })
    }

    pub(crate) fn scalar_text(self, tag: &str) -> Option<&'a str> {
        self.occurrences()
            .filter(|occurrence| occurrence.is_scalar())
            .find_map(|occurrence| occurrence.column(tag)?.value(0)?.as_text())
    }

    pub(crate) fn text_values(self, tag: &str) -> Vec<&'a str> {
        self.occurrences()
            .filter_map(|occurrence| occurrence.column(tag))
            .flat_map(CategoryColumn::values)
            .filter_map(CifValueRef::as_text)
            .collect()
    }
}

/// One scalar item or loop, before dictionary aliases assign canonical categories.
#[derive(Clone, Copy)]
pub(crate) enum CategoryOccurrence<'a> {
    Scalar(&'a CifItem),
    Loop(&'a CifLoop),
}

impl<'a> CategoryOccurrence<'a> {
    pub(crate) const fn is_scalar(self) -> bool {
        matches!(self, Self::Scalar { .. })
    }

    pub(crate) fn row_count(self) -> usize {
        match self {
            Self::Scalar(_) => 1,
            Self::Loop(cif_loop) => cif_loop.row_count(),
        }
    }

    pub(crate) const fn columns(self) -> CategoryColumns<'a> {
        CategoryColumns {
            occurrence: self,
            next: 0,
        }
    }

    pub(crate) fn column(self, tag: &str) -> Option<CategoryColumn<'a>> {
        self.columns()
            .find(|column| column.tag().eq_ignore_ascii_case(tag))
    }
}

/// One source column whose values are read lazily from the logical document.
#[derive(Clone, Copy)]
pub(crate) enum CategoryColumn<'a> {
    Scalar(&'a CifItem),
    Loop {
        tag: &'a str,
        cif_loop: &'a CifLoop,
        column: usize,
    },
}

impl<'a> CategoryColumn<'a> {
    pub(crate) fn tag(self) -> &'a str {
        match self {
            Self::Scalar(item) => item.tag(),
            Self::Loop { tag, .. } => tag,
        }
    }

    pub(crate) fn row_count(self) -> usize {
        match self {
            Self::Scalar(_) => 1,
            Self::Loop { cif_loop, .. } => cif_loop.row_count(),
        }
    }

    pub(crate) fn value(self, row: usize) -> Option<CifValueRef<'a>> {
        match self {
            Self::Scalar(item) => (row == 0).then(|| item.value().as_ref()),
            Self::Loop {
                cif_loop, column, ..
            } => cif_loop.value(row, column),
        }
    }

    pub(crate) fn values(self) -> impl Iterator<Item = CifValueRef<'a>> {
        (0..self.row_count()).filter_map(move |row| self.value(row))
    }
}

pub(crate) struct CategoryColumns<'a> {
    occurrence: CategoryOccurrence<'a>,
    next: usize,
}

impl<'a> Iterator for CategoryColumns<'a> {
    type Item = CategoryColumn<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let column = match self.occurrence {
            CategoryOccurrence::Scalar(item) => {
                (self.next == 0).then_some(CategoryColumn::Scalar(item))
            }
            CategoryOccurrence::Loop(cif_loop) => Some(CategoryColumn::Loop {
                tag: cif_loop.tags().get(self.next)?,
                cif_loop,
                column: self.next,
            }),
        };
        self.next += usize::from(column.is_some());
        column
    }
}

pub(crate) fn split_tag(tag: &str) -> Option<(&str, &str)> {
    tag.strip_prefix('_')?.split_once('.')
}

pub(super) fn required_scalar<'a>(
    entries: &'a [CifEntry],
    tag: &str,
) -> Result<&'a str, DictionaryError> {
    scalar(entries, tag)
        .ok_or_else(|| DictionaryError::new(format!("missing required DDL2 item {tag}")))
}

pub(super) fn scalar<'a>(entries: &'a [CifEntry], tag: &str) -> Option<&'a str> {
    CategoryView::new(entries).scalar_text(tag)
}

pub(super) fn values<'a>(entries: &'a [CifEntry], tag: &str) -> Vec<&'a str> {
    CategoryView::new(entries).text_values(tag)
}

pub(super) fn loop_records<'a>(entries: &'a [CifEntry], tags: &[&str]) -> Vec<Vec<&'a str>> {
    CategoryView::new(entries)
        .occurrences()
        .filter(|occurrence| !occurrence.is_scalar())
        .flat_map(|occurrence| {
            let columns = tags
                .iter()
                .map(|tag| occurrence.column(tag))
                .collect::<Option<Vec<_>>>();
            (0..occurrence.row_count()).filter_map(move |row| {
                columns
                    .as_ref()?
                    .iter()
                    .map(|column| column.value(row)?.as_text())
                    .collect()
            })
        })
        .collect()
}

pub(super) fn item_records(entries: &[CifEntry]) -> Vec<Vec<&str>> {
    let view = CategoryView::new(entries);
    let mut output = Vec::new();
    for occurrence in view
        .occurrences()
        .filter(|occurrence| !occurrence.is_scalar())
    {
        let Some(names) = occurrence.column("_item.name") else {
            continue;
        };
        let categories = occurrence.column("_item.category_id");
        let mandatory = occurrence.column("_item.mandatory_code");
        for row in 0..occurrence.row_count() {
            let Some(name) = names.value(row).and_then(CifValueRef::as_text) else {
                continue;
            };
            let category = categories
                .and_then(|column| column.value(row))
                .and_then(CifValueRef::as_text)
                .or_else(|| split_tag(name).map(|(category, _)| category));
            let mandatory = mandatory
                .and_then(|column| column.value(row))
                .and_then(CifValueRef::as_text)
                .unwrap_or("no");
            if let Some(category) = category {
                output.push(vec![name, category, mandatory]);
            }
        }
    }
    if output.is_empty() {
        let Some(name) = view.scalar_text("_item.name") else {
            return output;
        };
        let Some(category) = view
            .scalar_text("_item.category_id")
            .or_else(|| split_tag(name).map(|(category, _)| category))
        else {
            return output;
        };
        output.push(vec![
            name,
            category,
            view.scalar_text("_item.mandatory_code").unwrap_or("no"),
        ]);
    }
    output
}

pub(super) fn range_records(entries: &[CifEntry]) -> Vec<ItemRange> {
    let view = CategoryView::new(entries);
    let mut output = Vec::new();
    for occurrence in view
        .occurrences()
        .filter(|occurrence| !occurrence.is_scalar())
    {
        let (Some(minimum), Some(maximum)) = (
            occurrence.column("_item_range.minimum"),
            occurrence.column("_item_range.maximum"),
        ) else {
            continue;
        };
        output.extend((0..occurrence.row_count()).map(|row| ItemRange {
            minimum: minimum.value(row).and_then(range_bound),
            maximum: maximum.value(row).and_then(range_bound),
        }));
    }
    if output.is_empty()
        && let (Some(minimum), Some(maximum)) = (
            view.occurrences()
                .filter(|occurrence| occurrence.is_scalar())
                .find_map(|occurrence| occurrence.column("_item_range.minimum"))
                .map(|column| column.value(0).and_then(range_bound)),
            view.occurrences()
                .filter(|occurrence| occurrence.is_scalar())
                .find_map(|occurrence| occurrence.column("_item_range.maximum"))
                .map(|column| column.value(0).and_then(range_bound)),
        )
    {
        output.push(ItemRange { minimum, maximum });
    }
    output
}

fn range_bound(value: CifValueRef<'_>) -> Option<String> {
    match value {
        CifValueRef::Text(value) => Some(value.as_str().to_owned()),
        CifValueRef::Integer(value, _) => Some(value.to_string()),
        CifValueRef::Float(value, _, _) => Some(value.to_string()),
        CifValueRef::Unknown | CifValueRef::NotApplicable => None,
    }
}

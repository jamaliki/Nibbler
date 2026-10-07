//! One category of a data block read as a table.

use std::borrow::Cow;

use super::dictionary_source::split_tag;
use super::document::{CifEntry, CifItem, CifLoop, CifValueRef};

/// Whether `tag` belongs to `category`, ignoring ASCII case.
pub(crate) fn in_category(tag: &str, category: &str) -> bool {
    split_tag(tag).is_some_and(|(name, _)| name.eq_ignore_ascii_case(category))
}

/// Whether `tag` names `item` (in any category), ignoring ASCII case.
fn item_is(tag: &str, item: &str) -> bool {
    split_tag(tag).is_some_and(|(_, name)| name.eq_ignore_ascii_case(item))
}

/// One category of a block: its first loop, or its scalar items as one row in source
/// order. Values are borrowed from the block.
#[derive(Clone, Debug)]
pub(crate) enum BlockCategory<'a> {
    Loop(&'a CifLoop),
    Items(Vec<&'a CifItem>),
}

impl<'a> BlockCategory<'a> {
    /// Return `category` of the block whose entries are `entries`, or `None` when the
    /// block has none of its items.
    pub(crate) fn find(entries: &'a [CifEntry], category: &str) -> Option<Self> {
        let mut items = Vec::new();
        for entry in entries {
            match entry {
                CifEntry::Loop(cif_loop)
                    if cif_loop
                        .tags()
                        .first()
                        .is_some_and(|tag| in_category(tag, category)) =>
                {
                    return Some(Self::Loop(cif_loop));
                }
                CifEntry::Item(item) if in_category(item.tag(), category) => items.push(item),
                CifEntry::Loop(_) | CifEntry::Item(_) | CifEntry::Frame(_) => {}
            }
        }
        (!items.is_empty()).then_some(Self::Items(items))
    }

    pub(crate) fn row_count(&self) -> usize {
        match self {
            Self::Loop(cif_loop) => cif_loop.row_count(),
            Self::Items(_) => 1,
        }
    }

    #[cfg(feature = "reduce3")]
    pub(crate) const fn is_loop(&self) -> bool {
        matches!(self, Self::Loop(_))
    }

    /// Return the full tags in column order.
    pub(crate) fn tags(&self) -> Vec<&'a str> {
        match self {
            Self::Loop(cif_loop) => cif_loop.tags().iter().map(String::as_str).collect(),
            Self::Items(items) => items.iter().map(|item| item.tag()).collect(),
        }
    }

    /// Return the column index of `item`, ignoring ASCII case.
    pub(crate) fn column(&self, item: &str) -> Option<usize> {
        match self {
            Self::Loop(cif_loop) => cif_loop.tags().iter().position(|tag| item_is(tag, item)),
            Self::Items(items) => items.iter().position(|value| item_is(value.tag(), item)),
        }
    }

    pub(crate) fn value(&self, row: usize, column: usize) -> Option<CifValueRef<'a>> {
        match self {
            Self::Loop(cif_loop) => cif_loop.value(row, column),
            Self::Items(items) => (row == 0)
                .then(|| items.get(column).map(|item| item.value().as_ref()))
                .flatten(),
        }
    }

    /// Return a present value as text (a number as its source lexeme when it has one), or
    /// `None` for `?`, `.`, and absent values.
    pub(crate) fn text(&self, row: usize, column: usize) -> Option<Cow<'a, str>> {
        match self.value(row, column)? {
            CifValueRef::Text(text) => Some(Cow::Borrowed(text.as_str())),
            CifValueRef::Integer(number, original) => Some(original.map_or_else(
                || Cow::Owned(number.to_string()),
                |text| Cow::Borrowed(text.as_str()),
            )),
            CifValueRef::Float(number, _, original) => Some(original.map_or_else(
                || Cow::Owned(number.to_string()),
                |text| Cow::Borrowed(text.as_str()),
            )),
            CifValueRef::Unknown | CifValueRef::NotApplicable => None,
        }
    }
}

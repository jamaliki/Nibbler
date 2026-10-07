//! Loop columns built row by row from text, in the layout BinaryCIF decoding uses.

use super::document::{ColumnValues, LoopColumn, StringColumn};

/// One value appended to a [`TextColumnBuilder`].
#[derive(Clone, Copy, Debug)]
pub(crate) enum TextCell<'a> {
    Text(&'a str),
    NotApplicable,
    Unknown,
}

/// A string column under construction: one shared text buffer, offsets into it, one index
/// per row, and a missing-state mask created on the first `?` or `.`. Consecutive equal
/// values share one buffer entry.
#[derive(Debug)]
pub(crate) struct TextColumnBuilder {
    data: String,
    offsets: Vec<u32>,
    indices: Vec<u32>,
    mask: Option<Vec<u8>>,
}

impl TextColumnBuilder {
    pub(crate) fn with_capacity(rows: usize) -> Self {
        let mut offsets = Vec::with_capacity(rows.min(1 << 20) + 1);
        offsets.push(0);
        Self {
            data: String::new(),
            offsets,
            indices: Vec::with_capacity(rows),
            mask: None,
        }
    }

    fn last_entry(&self) -> Option<&str> {
        let [.., start, end] = self.offsets[..] else {
            return None;
        };
        self.data.get(start as usize..end as usize)
    }

    /// Append one value; `None` when the column outgrows 32-bit offsets.
    pub(crate) fn push(&mut self, value: TextCell<'_>) -> Option<()> {
        let row = self.indices.len();
        let (index, kind) = match value {
            TextCell::Text("") => (0, 0),
            TextCell::Text(text) => {
                if self.last_entry() != Some(text) {
                    self.data.push_str(text);
                    self.offsets.push(u32::try_from(self.data.len()).ok()?);
                }
                (u32::try_from(self.offsets.len() - 1).ok()?, 0)
            }
            TextCell::NotApplicable => (0, 1),
            TextCell::Unknown => (0, 2),
        };
        if kind != 0 || self.mask.is_some() {
            self.mask.get_or_insert_with(|| vec![0; row]).push(kind);
        }
        self.indices.push(index);
        Some(())
    }

    pub(crate) fn finish(self) -> LoopColumn {
        LoopColumn {
            values: ColumnValues::Strings(StringColumn::new(self.data, self.offsets, self.indices)),
            mask: self.mask,
        }
    }
}

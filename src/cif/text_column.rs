//! Loop columns built row by row from text, in the layout BinaryCIF decoding uses.

use super::document::{ColumnValues, LoopColumn, StringColumn, TextColumn};

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
    /// Whole UTF-8 values, one after another.
    data: Vec<u8>,
    offsets: Vec<u32>,
    indices: Vec<u32>,
    mask: Option<Vec<u8>>,
    /// Where the last entry starts in `data`. Entries are never empty, so it runs to the
    /// end of `data`, and `data` is empty until the first entry.
    last: usize,
}

impl TextColumnBuilder {
    pub(crate) fn with_capacity(rows: usize) -> Self {
        let mut offsets = Vec::with_capacity(rows.min(1 << 20) + 1);
        offsets.push(0);
        Self {
            data: Vec::new(),
            offsets,
            indices: Vec::with_capacity(rows),
            mask: None,
            last: 0,
        }
    }

    /// Append one value; `None` when the column outgrows 32-bit offsets.
    #[inline]
    pub(crate) fn push(&mut self, value: TextCell<'_>) -> Option<()> {
        let row = self.indices.len();
        let (index, kind) = match value {
            TextCell::Text("") => (0, 0),
            TextCell::Text(text) => {
                let last = self.data.get(self.last..).unwrap_or_default();
                if !same_bytes(last, text.as_bytes()) {
                    self.last = self.data.len();
                    append(&mut self.data, text.as_bytes());
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
        // `data` holds whole values copied from `str`s, so it is valid UTF-8.
        let data = String::from_utf8(self.data)
            .unwrap_or_else(|error| String::from_utf8_lossy(error.as_bytes()).into_owned());
        LoopColumn {
            values: ColumnValues::Strings(StringColumn::new(data, self.offsets, self.indices)),
            mask: self.mask,
        }
    }
}

/// Append `text` to `data`. A value of up to 16 bytes goes through a fixed-size buffer and
/// a truncation, which avoids a call to `memcpy` for each of them.
#[inline]
fn append(data: &mut Vec<u8>, text: &[u8]) {
    let length = text.len();
    let end = data.len() + length;
    match length {
        0..=3 => data.extend(text.iter().copied()),
        4..=8 => {
            let mut piece = [0; 8];
            piece[..4].copy_from_slice(&word::<4>(text, 0));
            piece[length - 4..length].copy_from_slice(&word::<4>(text, length - 4));
            data.extend_from_slice(&piece);
            data.truncate(end);
        }
        9..=16 => {
            let mut piece = [0; 16];
            piece[..8].copy_from_slice(&word::<8>(text, 0));
            piece[length - 8..length].copy_from_slice(&word::<8>(text, length - 8));
            data.extend_from_slice(&piece);
            data.truncate(end);
        }
        _ => data.extend_from_slice(text),
    }
}

/// Whether `a` and `b` are equal. Values of up to 16 bytes are compared with fixed-size
/// loads that may overlap, which avoids a call to `memcmp` for each of them.
#[inline]
fn same_bytes(a: &[u8], b: &[u8]) -> bool {
    let length = a.len();
    if length != b.len() {
        return false;
    }
    match length {
        0 => true,
        1..=3 => a[0] == b[0] && a[length / 2] == b[length / 2] && a[length - 1] == b[length - 1],
        4..=7 => {
            word::<4>(a, 0) == word::<4>(b, 0)
                && word::<4>(a, length - 4) == word::<4>(b, length - 4)
        }
        8..=16 => {
            word::<8>(a, 0) == word::<8>(b, 0)
                && word::<8>(a, length - 8) == word::<8>(b, length - 8)
        }
        _ => a == b,
    }
}

/// The `N` bytes of `bytes` from `start`, or zeros when `bytes` is too short.
#[inline]
fn word<const N: usize>(bytes: &[u8], start: usize) -> [u8; N] {
    bytes
        .get(start..start + N)
        .and_then(|piece| piece.try_into().ok())
        .unwrap_or([0; N])
}

/// A string column built from a source text column: a row takes a source row's value
/// by sharing its buffer entry, or new text, which goes after the source entries.
/// Consecutive equal new texts share one entry.
#[derive(Debug)]
pub(crate) struct SourceColumnBuilder<'a> {
    source_indices: &'a [u32],
    source_mask: Option<&'a [u8]>,
    data: String,
    offsets: Vec<u32>,
    indices: Vec<u32>,
    mask: Option<Vec<u8>>,
    /// The index of the last new text.
    last: Option<u32>,
}

impl<'a> SourceColumnBuilder<'a> {
    /// Start from a copy of `source`'s text, with room for `rows` rows and `entries` new
    /// entries of `text` bytes in all; `None` when the source entries do not start at 0
    /// and end at the end of its text.
    pub(crate) fn new(source: TextColumn<'a>, rows: usize, entries: usize, text: usize) -> Option<Self> {
        let (data, offsets, indices, mask, source_rows) = source.parts();
        let end = *offsets.last()? as usize;
        if offsets.first() != Some(&0) || data.len() != end || indices.len() < source_rows {
            return None;
        }
        Some(Self {
            source_indices: &indices[..source_rows],
            source_mask: mask.map(|mask| &mask[..source_rows.min(mask.len())]),
            data: {
                let mut copy = String::with_capacity(data.len().saturating_add(text));
                copy.push_str(data);
                copy
            },
            offsets: {
                let mut copy = Vec::with_capacity(offsets.len().saturating_add(entries));
                copy.extend_from_slice(offsets);
                copy
            },
            indices: Vec::with_capacity(rows),
            mask: None,
            last: None,
        })
    }

    /// Append the value of source row `row`, or `?` when the source has no such row.
    #[inline]
    pub(crate) fn push_source(&mut self, row: usize) {
        let state = self.source_mask.map_or(Some(0), |mask| mask.get(row).copied());
        match (state, self.source_indices.get(row)) {
            (Some(1), Some(_)) => self.push_index(0, 1),
            (Some(2), Some(_)) | (None, _) | (_, None) => self.push_index(0, 2),
            (Some(_), Some(&index)) => self.push_index(index, 0),
        }
    }

    /// Append new text or a missing state; `None` when the column outgrows 32-bit
    /// offsets.
    pub(crate) fn push(&mut self, value: TextCell<'_>) -> Option<()> {
        match value {
            TextCell::Text("") => self.push_index(0, 0),
            TextCell::Text(text) => {
                let last = self.last.and_then(|index| {
                    let start = *self.offsets.get((index as usize).checked_sub(1)?)? as usize;
                    let end = *self.offsets.get(index as usize)? as usize;
                    (self.data.get(start..end)? == text).then_some(index)
                });
                let index = match last {
                    Some(index) => index,
                    None => {
                        self.data.push_str(text);
                        self.offsets.push(u32::try_from(self.data.len()).ok()?);
                        u32::try_from(self.offsets.len() - 1).ok()?
                    }
                };
                self.last = Some(index);
                self.push_index(index, 0);
            }
            TextCell::NotApplicable => self.push_index(0, 1),
            TextCell::Unknown => self.push_index(0, 2),
        }
        Some(())
    }

    #[inline]
    fn push_index(&mut self, index: u32, state: u8) {
        let row = self.indices.len();
        if state != 0 || self.mask.is_some() {
            self.mask.get_or_insert_with(|| vec![0; row]).push(state);
        }
        self.indices.push(index);
    }

    pub(crate) fn finish(self) -> LoopColumn {
        LoopColumn {
            values: ColumnValues::Strings(StringColumn::new(self.data, self.offsets, self.indices)),
            mask: self.mask,
        }
    }
}

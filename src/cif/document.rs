use std::fmt::{self, Debug, Formatter};
use std::ops::Range;
use std::sync::Arc;

use super::source::SourceBuffer;

/// The lexical form used by a parsed text value.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum QuoteStyle {
    /// A whitespace-delimited unquoted value.
    Unquoted,
    /// A value delimited by single quotes.
    Single,
    /// A value delimited by double quotes.
    Double,
    /// A semicolon-delimited text field.
    TextField,
}

#[derive(Clone)]
enum TextStorage {
    Source {
        source: SourceBuffer,
        range: Range<usize>,
    },
    Owned(Arc<str>),
}

/// An immutable CIF string value.
///
/// Parsed values borrow logically from a reference-counted [`SourceBuffer`]; values
/// constructed with [`TextValue::new`] own one string allocation. Equality compares
/// decoded text, not the original quoting style.
#[derive(Clone)]
pub struct TextValue {
    storage: TextStorage,
    quote_style: QuoteStyle,
}

impl TextValue {
    /// Construct an owned text value.
    #[must_use]
    pub fn new(text: impl Into<String>) -> Self {
        Self {
            storage: TextStorage::Owned(Arc::from(text.into())),
            quote_style: QuoteStyle::Unquoted,
        }
    }

    pub(crate) fn from_source(
        source: SourceBuffer,
        range: Range<usize>,
        quote_style: QuoteStyle,
    ) -> Self {
        Self {
            storage: TextStorage::Source { source, range },
            quote_style,
        }
    }

    /// Return the decoded value without CIF delimiters.
    #[must_use]
    pub fn as_str(&self) -> &str {
        match &self.storage {
            TextStorage::Source { source, range } => source.slice(range.clone()),
            TextStorage::Owned(text) => text,
        }
    }

    /// Return the original lexical style, or [`QuoteStyle::Unquoted`] for a newly
    /// constructed value.
    #[must_use]
    pub const fn quote_style(&self) -> QuoteStyle {
        self.quote_style
    }
}

impl Debug for TextValue {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("TextValue")
            .field("text", &self.as_str())
            .field("quote_style", &self.quote_style)
            .finish()
    }
}

impl PartialEq for TextValue {
    fn eq(&self, other: &Self) -> bool {
        self.as_str() == other.as_str()
    }
}

impl Eq for TextValue {}

/// A retained numeric spelling for future preserving output.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OriginalLexeme(TextValue);

impl OriginalLexeme {
    /// Construct an owned original spelling.
    #[must_use]
    pub fn new(text: impl Into<String>) -> Self {
        Self(TextValue::new(text))
    }

    /// Return the original spelling.
    #[must_use]
    pub fn as_str(&self) -> &str {
        self.0.as_str()
    }
}

/// Parenthesized uncertainty digits associated with a CIF floating-point value.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct StandardUncertainty(u64);

impl StandardUncertainty {
    /// Construct uncertainty digits to be written in parentheses.
    #[must_use]
    pub const fn new(digits: u64) -> Self {
        Self(digits)
    }

    /// Return the parenthesized uncertainty digits.
    #[must_use]
    pub const fn digits(self) -> u64 {
        self.0
    }
}

/// A logical CIF value with distinct missing states.
///
/// Schema-less parsing keeps ordinary values as [`CifValue::Text`]. Integer and float
/// variants are available for later schema-guided parsing and constructed documents.
#[derive(Clone, Debug, PartialEq)]
pub enum CifValue {
    /// Present text, including an explicitly quoted empty string.
    Text(TextValue),
    /// A typed integer and optional original spelling.
    Integer(i64, Option<OriginalLexeme>),
    /// A typed finite float, optional standard uncertainty, and original spelling.
    Float(f64, Option<StandardUncertainty>, Option<OriginalLexeme>),
    /// CIF unknown (`?`).
    Unknown,
    /// CIF not-applicable (`.`).
    NotApplicable,
}

impl CifValue {
    /// Construct a present text value.
    #[must_use]
    pub fn text(text: impl Into<String>) -> Self {
        Self::Text(TextValue::new(text))
    }

    /// Return present text, or `None` for typed and missing values.
    #[must_use]
    pub fn as_text(&self) -> Option<&str> {
        match self {
            Self::Text(text) => Some(text.as_str()),
            Self::Integer(..) | Self::Float(..) | Self::Unknown | Self::NotApplicable => None,
        }
    }

    /// Borrow this owned value without allocating.
    #[must_use]
    pub fn as_ref(&self) -> CifValueRef<'_> {
        match self {
            Self::Text(text) => CifValueRef::Text(TextValueRef {
                text: text.as_str(),
                quote_style: text.quote_style(),
            }),
            Self::Integer(number, original) => CifValueRef::Integer(*number, original.as_ref()),
            Self::Float(number, uncertainty, original) => {
                CifValueRef::Float(*number, *uncertainty, original.as_ref())
            }
            Self::Unknown => CifValueRef::Unknown,
            Self::NotApplicable => CifValueRef::NotApplicable,
        }
    }
}

/// A borrowed CIF text value.
#[derive(Clone, Copy, Debug)]
pub struct TextValueRef<'a> {
    text: &'a str,
    quote_style: QuoteStyle,
}

impl<'a> TextValueRef<'a> {
    /// Return the decoded value without CIF delimiters.
    #[must_use]
    pub const fn as_str(self) -> &'a str {
        self.text
    }

    /// Return the source quoting style.
    #[must_use]
    pub const fn quote_style(self) -> QuoteStyle {
        self.quote_style
    }

    fn to_owned(self) -> TextValue {
        TextValue {
            storage: TextStorage::Owned(Arc::from(self.text)),
            quote_style: self.quote_style,
        }
    }
}

impl PartialEq for TextValueRef<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.text == other.text
    }
}

impl Eq for TextValueRef<'_> {}

/// A non-owning view of one logical CIF value.
///
/// Parsed loop storage returns this compact view instead of allocating an owned
/// [`CifValue`] for every cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CifValueRef<'a> {
    /// Present text, including an explicitly quoted empty string.
    Text(TextValueRef<'a>),
    /// A typed integer and optional original spelling.
    Integer(i64, Option<&'a OriginalLexeme>),
    /// A typed finite float, optional standard uncertainty, and original spelling.
    Float(f64, Option<StandardUncertainty>, Option<&'a OriginalLexeme>),
    /// CIF unknown (`?`).
    Unknown,
    /// CIF not-applicable (`.`).
    NotApplicable,
}

impl<'a> CifValueRef<'a> {
    /// Return present text, or `None` for typed and missing values.
    #[must_use]
    pub const fn as_text(self) -> Option<&'a str> {
        match self {
            Self::Text(text) => Some(text.as_str()),
            Self::Integer(..) | Self::Float(..) | Self::Unknown | Self::NotApplicable => None,
        }
    }

    /// Allocate an independent owned value.
    #[must_use]
    pub fn to_owned(self) -> CifValue {
        match self {
            Self::Text(text) => CifValue::Text(text.to_owned()),
            Self::Integer(number, original) => CifValue::Integer(number, original.cloned()),
            Self::Float(number, uncertainty, original) => {
                CifValue::Float(number, uncertainty, original.cloned())
            }
            Self::Unknown => CifValue::Unknown,
            Self::NotApplicable => CifValue::NotApplicable,
        }
    }
}

impl PartialEq<CifValue> for CifValueRef<'_> {
    fn eq(&self, other: &CifValue) -> bool {
        *self == other.as_ref()
    }
}

impl PartialEq<CifValueRef<'_>> for CifValue {
    fn eq(&self, other: &CifValueRef<'_>) -> bool {
        self.as_ref() == *other
    }
}

/// A scalar tag-value pair.
#[derive(Clone, Debug, PartialEq)]
pub struct CifItem {
    tag: String,
    value: CifValue,
}

impl CifItem {
    pub(crate) fn new(tag: String, value: CifValue) -> Self {
        Self { tag, value }
    }

    /// Return the tag with its leading underscore.
    #[must_use]
    pub fn tag(&self) -> &str {
        &self.tag
    }

    /// Return the scalar value.
    #[must_use]
    pub const fn value(&self) -> &CifValue {
        &self.value
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct SourceCell {
    span_start: u32,
    span_end: u32,
}

#[derive(Clone, Debug)]
struct SourceCells {
    chunks: Vec<Vec<SourceCell>>,
    offsets: Vec<usize>,
    len: usize,
}

impl SourceCells {
    fn chunked(chunks: Vec<Vec<SourceCell>>) -> Self {
        let mut offsets = Vec::with_capacity(chunks.len());
        let mut len = 0;
        for chunk in &chunks {
            offsets.push(len);
            len += chunk.len();
        }
        Self {
            chunks,
            offsets,
            len,
        }
    }

    fn get(&self, index: usize) -> Option<SourceCell> {
        if index >= self.len {
            return None;
        }
        let chunk_index = self.offsets.partition_point(|offset| *offset <= index) - 1;
        self.chunks[chunk_index]
            .get(index - self.offsets[chunk_index])
            .copied()
    }
}

impl SourceCell {
    pub(crate) fn from_token(token: &super::token::Token) -> Self {
        debug_assert!(token.span.start <= token.span.end);
        debug_assert!(u32::try_from(token.span.end).is_ok());
        // Full-document parsing rejects sources larger than `u32::MAX` before
        // tokenization. Lexer ranges are bounded by that source, so these casts cannot
        // truncate for any accepted document.
        Self {
            span_start: token.span.start as u32,
            span_end: token.span.end as u32,
        }
    }

    fn value<'a>(self, source: &'a SourceBuffer) -> CifValueRef<'a> {
        let bytes = source.as_str().as_bytes();
        let start = self.span_start as usize;
        let end = self.span_end as usize;
        let (content, quote_style) = match bytes[start] {
            b'\'' => (start + 1..end - 1, QuoteStyle::Single),
            b'"' => (start + 1..end - 1, QuoteStyle::Double),
            b';' if start == 0 || matches!(bytes[start - 1], b'\n' | b'\r') => (
                start + 1..preceding_line_ending_start(bytes, end - 1),
                QuoteStyle::TextField,
            ),
            _ => (start..end, QuoteStyle::Unquoted),
        };
        let text = source.slice(content);
        if quote_style == QuoteStyle::Unquoted {
            if text == "?" {
                return CifValueRef::Unknown;
            }
            if text == "." {
                return CifValueRef::NotApplicable;
            }
        }
        CifValueRef::Text(TextValueRef { text, quote_style })
    }
}

fn preceding_line_ending_start(bytes: &[u8], byte_offset: usize) -> usize {
    let mut content_end = byte_offset;
    if content_end > 0 && bytes[content_end - 1] == b'\n' {
        content_end -= 1;
        if content_end > 0 && bytes[content_end - 1] == b'\r' {
            content_end -= 1;
        }
    } else if content_end > 0 && bytes[content_end - 1] == b'\r' {
        content_end -= 1;
    }
    content_end
}

#[derive(Clone, Debug)]
pub(crate) struct StringColumn {
    data: String,
    offsets: Vec<u32>,
    indices: Vec<u32>,
}

impl StringColumn {
    pub(crate) fn new(data: String, offsets: Vec<u32>, indices: Vec<u32>) -> Self {
        Self {
            data,
            offsets,
            indices,
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.indices.len()
    }

    pub(crate) fn value(&self, row_index: usize) -> &str {
        let encoded = self.indices[row_index];
        if encoded == 0 {
            return "";
        }
        let index = encoded as usize - 1;
        let start = self.offsets[index] as usize;
        let end = self.offsets[index + 1] as usize;
        &self.data[start..end]
    }
}

#[derive(Clone, Debug)]
pub(crate) enum ColumnValues {
    Integers(Vec<i64>),
    Floats(Vec<f64>),
    Strings(StringColumn),
}

impl ColumnValues {
    pub(crate) fn len(&self) -> usize {
        match self {
            Self::Integers(values) => values.len(),
            Self::Floats(values) => values.len(),
            Self::Strings(values) => values.len(),
        }
    }

    fn value(&self, row_index: usize) -> CifValueRef<'_> {
        match self {
            Self::Integers(values) => CifValueRef::Integer(values[row_index], None),
            Self::Floats(values) => CifValueRef::Float(values[row_index], None, None),
            Self::Strings(values) => CifValueRef::Text(TextValueRef {
                text: values.value(row_index),
                quote_style: QuoteStyle::Unquoted,
            }),
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct LoopColumn {
    pub(crate) values: ColumnValues,
    pub(crate) mask: Option<Vec<u8>>,
}

impl LoopColumn {
    fn value(&self, row_index: usize) -> CifValueRef<'_> {
        match self.mask.as_ref().map(|mask| mask[row_index]).unwrap_or(0) {
            1 => CifValueRef::NotApplicable,
            2 => CifValueRef::Unknown,
            _ => self.values.value(row_index),
        }
    }
}

#[derive(Clone, Debug)]
enum LoopStorage {
    Owned(Vec<CifValue>),
    Source {
        source: SourceBuffer,
        cells: SourceCells,
    },
    Columns {
        columns: Vec<LoopColumn>,
        row_count: usize,
    },
}

/// One CIF loop with its original tag order and immediately addressable values.
#[derive(Clone, Debug)]
pub struct CifLoop {
    tags: Vec<String>,
    storage: LoopStorage,
}

impl CifLoop {
    pub(crate) fn new(tags: Vec<String>, values: Vec<CifValue>) -> Self {
        Self {
            tags,
            storage: LoopStorage::Owned(values),
        }
    }

    pub(crate) fn from_source_chunks(
        tags: Vec<String>,
        source: SourceBuffer,
        cells: Vec<Vec<SourceCell>>,
    ) -> Self {
        Self {
            tags,
            storage: LoopStorage::Source {
                source,
                cells: SourceCells::chunked(cells),
            },
        }
    }

    pub(crate) fn from_columns(
        tags: Vec<String>,
        columns: Vec<LoopColumn>,
        row_count: usize,
    ) -> Self {
        Self {
            tags,
            storage: LoopStorage::Columns { columns, row_count },
        }
    }

    /// Return loop tags in source order.
    #[must_use]
    pub fn tags(&self) -> &[String] {
        &self.tags
    }

    /// Iterate over values in row-major order.
    #[must_use]
    pub fn values(&self) -> CifValues<'_> {
        CifValues {
            cif_loop: self,
            next_index: 0,
            source_chunk: 0,
            source_cell: 0,
        }
    }

    /// Return the number of logical values.
    #[must_use]
    pub fn value_count(&self) -> usize {
        match &self.storage {
            LoopStorage::Owned(values) => values.len(),
            LoopStorage::Source { cells, .. } => cells.len,
            LoopStorage::Columns { columns, row_count } => columns.len().saturating_mul(*row_count),
        }
    }

    /// Return the number of columns.
    #[must_use]
    pub fn column_count(&self) -> usize {
        self.tags.len()
    }

    /// Return the number of complete rows.
    #[must_use]
    pub fn row_count(&self) -> usize {
        self.value_count().checked_div(self.tags.len()).unwrap_or(0)
    }

    /// Return one complete row, or `None` if `row_index` is out of range.
    #[must_use]
    pub fn row(&self, row_index: usize) -> Option<CifRow<'_>> {
        (row_index < self.row_count()).then_some(CifRow {
            cif_loop: self,
            row_index,
        })
    }

    /// Return one value by row and column index.
    #[must_use]
    pub fn value(&self, row_index: usize, column_index: usize) -> Option<CifValueRef<'_>> {
        if column_index >= self.column_count() || row_index >= self.row_count() {
            return None;
        }
        match &self.storage {
            LoopStorage::Owned(values) => values
                .get(row_index * self.column_count() + column_index)
                .map(CifValue::as_ref),
            LoopStorage::Source { source, cells } => cells
                .get(row_index * self.column_count() + column_index)
                .map(|cell| cell.value(source)),
            LoopStorage::Columns { columns, .. } => columns
                .get(column_index)
                .map(|column| column.value(row_index)),
        }
    }
}

impl PartialEq for CifLoop {
    fn eq(&self, other: &Self) -> bool {
        self.tags == other.tags
            && self.value_count() == other.value_count()
            && self.values().eq(other.values())
    }
}

/// A borrowed row of a [`CifLoop`].
#[derive(Clone, Copy)]
pub struct CifRow<'a> {
    cif_loop: &'a CifLoop,
    row_index: usize,
}

impl<'a> CifRow<'a> {
    /// Return the number of columns in the row.
    #[must_use]
    pub fn len(self) -> usize {
        self.cif_loop.column_count()
    }

    /// Return whether the row has no columns.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.len() == 0
    }

    /// Return one value by column index.
    #[must_use]
    pub fn get(self, column_index: usize) -> Option<CifValueRef<'a>> {
        self.cif_loop.value(self.row_index, column_index)
    }

    /// Iterate over the row in column order.
    #[must_use]
    pub fn iter(self) -> CifRowValues<'a> {
        CifRowValues {
            row: self,
            next_column: 0,
        }
    }
}

impl Debug for CifRow<'_> {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter.debug_list().entries(self.iter()).finish()
    }
}

/// Iterator over one loop row.
pub struct CifRowValues<'a> {
    row: CifRow<'a>,
    next_column: usize,
}

impl<'a> Iterator for CifRowValues<'a> {
    type Item = CifValueRef<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let value = self.row.get(self.next_column)?;
        self.next_column += 1;
        Some(value)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.row.len().saturating_sub(self.next_column);
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for CifRowValues<'_> {}

/// Iterator over loop values in row-major order.
pub struct CifValues<'a> {
    cif_loop: &'a CifLoop,
    next_index: usize,
    source_chunk: usize,
    source_cell: usize,
}

impl<'a> Iterator for CifValues<'a> {
    type Item = CifValueRef<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next_index >= self.cif_loop.value_count() {
            return None;
        }
        let value = match &self.cif_loop.storage {
            LoopStorage::Owned(values) => values.get(self.next_index).map(CifValue::as_ref),
            LoopStorage::Source { source, cells } => {
                while cells
                    .chunks
                    .get(self.source_chunk)
                    .is_some_and(|chunk| self.source_cell == chunk.len())
                {
                    self.source_chunk += 1;
                    self.source_cell = 0;
                }
                let value = cells
                    .chunks
                    .get(self.source_chunk)
                    .and_then(|chunk| chunk.get(self.source_cell))
                    .map(|cell| cell.value(source));
                self.source_cell += usize::from(value.is_some());
                value
            }
            LoopStorage::Columns { columns, row_count } => {
                let column_count = columns.len();
                let row_index = self.next_index.checked_div(column_count)?;
                let column_index = self.next_index % column_count;
                (row_index < *row_count).then(|| columns[column_index].value(row_index))
            }
        };
        self.next_index += 1;
        value
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.cif_loop.value_count().saturating_sub(self.next_index);
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for CifValues<'_> {}

/// A scalar item, distinct loop, or save frame in source order.
#[derive(Clone, Debug, PartialEq)]
pub enum CifEntry {
    /// One scalar item.
    Item(CifItem),
    /// One loop; loops with equal category prefixes remain distinct.
    Loop(CifLoop),
    /// One non-nested save frame.
    Frame(CifFrame),
}

/// A named CIF save frame.
#[derive(Clone, Debug, PartialEq)]
pub struct CifFrame {
    code: String,
    entries: Vec<CifEntry>,
}

impl CifFrame {
    pub(crate) fn new(code: String, entries: Vec<CifEntry>) -> Self {
        Self { code, entries }
    }

    /// Return the frame code without the `save_` prefix.
    #[must_use]
    pub fn code(&self) -> &str {
        &self.code
    }

    /// Return frame entries in source order.
    #[must_use]
    pub fn entries(&self) -> &[CifEntry] {
        &self.entries
    }
}

/// The kind of top-level CIF block.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BlockKind {
    /// A named `data_` block.
    Data,
    /// The unnamed `global_` block.
    Global,
}

/// One top-level CIF data or global block.
#[derive(Clone, Debug, PartialEq)]
pub struct CifBlock {
    kind: BlockKind,
    code: Option<String>,
    entries: Vec<CifEntry>,
}

impl CifBlock {
    pub(crate) fn data(code: String, entries: Vec<CifEntry>) -> Self {
        Self {
            kind: BlockKind::Data,
            code: Some(code),
            entries,
        }
    }

    pub(crate) fn global(entries: Vec<CifEntry>) -> Self {
        Self {
            kind: BlockKind::Global,
            code: None,
            entries,
        }
    }

    /// Return whether this is a data or global block.
    #[must_use]
    pub const fn kind(&self) -> BlockKind {
        self.kind
    }

    /// Return the code of a data block, or `None` for a global block.
    #[must_use]
    pub fn code(&self) -> Option<&str> {
        self.code.as_deref()
    }

    /// Return block entries in source order.
    #[must_use]
    pub fn entries(&self) -> &[CifEntry] {
        &self.entries
    }
}

/// An order-preserving generic CIF document.
#[derive(Clone, Debug, PartialEq)]
pub struct CifDocument {
    blocks: Arc<[CifBlock]>,
}

impl CifDocument {
    pub(crate) fn new(blocks: Vec<CifBlock>) -> Self {
        Self {
            blocks: blocks.into(),
        }
    }

    /// Return blocks in source order.
    #[must_use]
    pub fn blocks(&self) -> &[CifBlock] {
        &self.blocks
    }
}

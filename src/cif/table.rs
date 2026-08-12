use std::sync::Arc;

use arrow_array::builder::{Float64Builder, Int64Builder, StringBuilder, UInt8Builder};
use arrow_array::{Array, ArrayRef, Float64Array, Int64Array, StringArray, UInt8Array};
use arrow_schema::DataType;

use super::numeric::{parse_float, parse_integer};
use super::sink::ParsedValue;

/// Physical Arrow type selected by a compiled CIF dictionary.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum ColumnType {
    /// Schema-less or character data.
    #[default]
    Text,
    /// DDL2 integer data.
    Integer,
    /// DDL2 floating-point data.
    Float,
}

/// The reason a CIF value is absent.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MissingKind {
    /// CIF unknown (`?`).
    Unknown,
    /// CIF not-applicable (`.`).
    NotApplicable,
}

/// One schema-less projected value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CifCell {
    /// Present text, including quoted empty text.
    Text(String),
    /// CIF unknown (`?`).
    Unknown,
    /// CIF not-applicable (`.`).
    NotApplicable,
}

/// One borrowed schema-less projected value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CifCellRef<'a> {
    /// Present text, including quoted empty text.
    Text(&'a str),
    /// A present schema-typed integer.
    Integer(i64),
    /// A present schema-typed finite float.
    Float(f64),
    /// CIF unknown (`?`).
    Unknown,
    /// CIF not-applicable (`.`).
    NotApplicable,
}

impl CifCell {
    /// Return present text, or `None` for either missing state.
    #[must_use]
    pub fn as_text(&self) -> Option<&str> {
        match self {
            Self::Text(text) => Some(text),
            Self::Unknown | Self::NotApplicable => None,
        }
    }

    /// Return the precise missing state, or `None` for present text.
    #[must_use]
    pub const fn missing_kind(&self) -> Option<MissingKind> {
        match self {
            Self::Text(_) => None,
            Self::Unknown => Some(MissingKind::Unknown),
            Self::NotApplicable => Some(MissingKind::NotApplicable),
        }
    }
}

/// One column in a projected CIF table.
#[derive(Clone, Debug)]
pub struct CifColumn {
    name: String,
    values: Vec<ColumnData>,
    missing: Vec<UInt8Array>,
}

#[derive(Clone, Debug)]
enum ColumnData {
    Text(StringArray),
    Integer(Int64Array),
    Float(Float64Array),
}

impl ColumnData {
    fn len(&self) -> usize {
        match self {
            Self::Text(values) => values.len(),
            Self::Integer(values) => values.len(),
            Self::Float(values) => values.len(),
        }
    }
}

impl CifColumn {
    /// Return the item name without its category prefix.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Return column values in source row order.
    #[must_use]
    pub fn values(&self) -> impl ExactSizeIterator<Item = CifCellRef<'_>> {
        ColumnValues {
            column: self,
            chunk_index: 0,
            row_index: 0,
            remaining: self.len(),
        }
    }

    fn value_at(&self, chunk_index: usize, row_index: usize) -> CifCellRef<'_> {
        match self.missing[chunk_index].value(row_index) {
            0 => match &self.values[chunk_index] {
                ColumnData::Text(values) => CifCellRef::Text(values.value(row_index)),
                ColumnData::Integer(values) => CifCellRef::Integer(values.value(row_index)),
                ColumnData::Float(values) => CifCellRef::Float(values.value(row_index)),
            },
            1 => CifCellRef::Unknown,
            2 => CifCellRef::NotApplicable,
            kind => {
                debug_assert!(false, "invalid internal CIF missing kind {kind}");
                CifCellRef::Unknown
            }
        }
    }

    pub(crate) fn data_type(&self) -> DataType {
        match self.values.first() {
            Some(ColumnData::Integer(_)) => DataType::Int64,
            Some(ColumnData::Float(_)) => DataType::Float64,
            Some(ColumnData::Text(_)) | None => DataType::Utf8,
        }
    }

    pub(crate) fn array(&self, chunk_index: usize) -> ArrayRef {
        match &self.values[chunk_index] {
            ColumnData::Text(values) => Arc::new(values.clone()),
            ColumnData::Integer(values) => Arc::new(values.clone()),
            ColumnData::Float(values) => Arc::new(values.clone()),
        }
    }

    pub(crate) fn missing_array(&self, chunk_index: usize) -> &UInt8Array {
        &self.missing[chunk_index]
    }

    pub(crate) fn chunk_count(&self) -> usize {
        self.values.len()
    }

    fn len(&self) -> usize {
        self.values.iter().map(ColumnData::len).sum()
    }
}

struct ColumnValues<'a> {
    column: &'a CifColumn,
    chunk_index: usize,
    row_index: usize,
    remaining: usize,
}

impl<'a> Iterator for ColumnValues<'a> {
    type Item = CifCellRef<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        while self
            .column
            .values
            .get(self.chunk_index)
            .is_some_and(|chunk| self.row_index == chunk.len())
        {
            self.chunk_index += 1;
            self.row_index = 0;
        }
        if self.remaining == 0 {
            return None;
        }
        let value = self.column.value_at(self.chunk_index, self.row_index);
        self.row_index += 1;
        self.remaining -= 1;
        Some(value)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl ExactSizeIterator for ColumnValues<'_> {}

/// Source identity retained for every projected row.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RowProvenance {
    source_name: String,
    block_code: Option<String>,
    frame_code: Option<String>,
    row_count: usize,
}

impl RowProvenance {
    /// Return the caller-supplied source name.
    #[must_use]
    pub fn source_name(&self) -> &str {
        &self.source_name
    }

    /// Return the data-block code, or `None` for a global block.
    #[must_use]
    pub fn block_code(&self) -> Option<&str> {
        self.block_code.as_deref()
    }

    /// Return the save-frame code when the row came from a frame.
    #[must_use]
    pub fn frame_code(&self) -> Option<&str> {
        self.frame_code.as_deref()
    }

    /// Return the number of consecutive rows sharing this identity.
    #[must_use]
    pub const fn row_count(&self) -> usize {
        self.row_count
    }
}

/// A projected category with column-major values and lossless CIF missing states.
#[derive(Clone, Debug)]
pub struct CifTable {
    category: String,
    columns: Vec<CifColumn>,
    provenance: Vec<RowProvenance>,
    chunk_rows: Vec<usize>,
}

pub(super) struct TableBuilder {
    pub(super) category: String,
    pub(super) columns: Vec<ColumnBuilder>,
    pub(super) provenance: Vec<RowProvenance>,
}

pub(super) struct ColumnBuilder {
    pub(super) name: String,
    values: ColumnDataBuilder,
    missing: UInt8Builder,
}

enum ColumnDataBuilder {
    Text(StringBuilder),
    Integer(Int64Builder),
    Float(Float64Builder),
}

impl ColumnBuilder {
    pub(super) fn new(name: String, column_type: ColumnType) -> Self {
        Self {
            name,
            values: match column_type {
                ColumnType::Text => ColumnDataBuilder::Text(StringBuilder::new()),
                ColumnType::Integer => ColumnDataBuilder::Integer(Int64Builder::new()),
                ColumnType::Float => ColumnDataBuilder::Float(Float64Builder::new()),
            },
            missing: UInt8Builder::new(),
        }
    }

    pub(super) fn column_type(&self) -> ColumnType {
        match self.values {
            ColumnDataBuilder::Text(_) => ColumnType::Text,
            ColumnDataBuilder::Integer(_) => ColumnType::Integer,
            ColumnDataBuilder::Float(_) => ColumnType::Float,
        }
    }

    pub(super) fn append(&mut self, value: CifCell) -> Result<(), String> {
        match value {
            CifCell::Text(text) => self.append_text_value(&text)?,
            CifCell::Unknown => self.append_missing(MissingKind::Unknown),
            CifCell::NotApplicable => self.append_missing(MissingKind::NotApplicable),
        }
        Ok(())
    }

    pub(super) fn append_text_value(&mut self, text: &str) -> Result<(), String> {
        self.append_text(text)?;
        self.missing.append_value(0);
        Ok(())
    }

    pub(super) fn append_integer_value(&mut self, number: i64) {
        match &mut self.values {
            ColumnDataBuilder::Text(values) => values.append_value(number.to_string()),
            ColumnDataBuilder::Integer(values) => values.append_value(number),
            ColumnDataBuilder::Float(values) => values.append_value(number as f64),
        }
        self.missing.append_value(0);
    }

    pub(super) fn append_float_value(&mut self, number: f64) -> Result<(), String> {
        match &mut self.values {
            ColumnDataBuilder::Text(values) => values.append_value(number.to_string()),
            ColumnDataBuilder::Integer(_) => {
                return Err(format!(
                    "floating-point value {number} is not a DDL2 integer"
                ));
            }
            ColumnDataBuilder::Float(values) => values.append_value(number),
        }
        self.missing.append_value(0);
        Ok(())
    }

    pub(super) fn append_missing(&mut self, kind: MissingKind) {
        self.append_null();
        self.missing.append_value(match kind {
            MissingKind::Unknown => 1,
            MissingKind::NotApplicable => 2,
        });
    }

    pub(super) fn append_parsed(&mut self, value: &ParsedValue<'_>) -> Result<(), String> {
        match (value.is_unquoted(), value.as_str()) {
            (true, "?") => {
                self.append_null();
                self.missing.append_value(1);
            }
            (true, ".") => {
                self.append_null();
                self.missing.append_value(2);
            }
            _ => {
                self.append_text(value.as_str())?;
                self.missing.append_value(0);
            }
        }
        Ok(())
    }

    fn append_text(&mut self, text: &str) -> Result<(), String> {
        match &mut self.values {
            ColumnDataBuilder::Text(values) => values.append_value(text),
            ColumnDataBuilder::Integer(values) => {
                let Some(number) = parse_integer(text) else {
                    return Err(format!("value {text:?} is not a DDL2 integer"));
                };
                values.append_value(number);
            }
            ColumnDataBuilder::Float(values) => {
                let Some(number) = parse_float(text) else {
                    return Err(format!("value {text:?} is not a finite DDL2 float"));
                };
                values.append_value(number);
            }
        }
        Ok(())
    }

    fn append_null(&mut self) {
        match &mut self.values {
            ColumnDataBuilder::Text(values) => values.append_null(),
            ColumnDataBuilder::Integer(values) => values.append_null(),
            ColumnDataBuilder::Float(values) => values.append_null(),
        }
    }

    fn finish(self) -> CifColumn {
        let values = match self.values {
            ColumnDataBuilder::Text(mut values) => ColumnData::Text(values.finish()),
            ColumnDataBuilder::Integer(mut values) => ColumnData::Integer(values.finish()),
            ColumnDataBuilder::Float(mut values) => ColumnData::Float(values.finish()),
        };
        let mut missing = self.missing;
        CifColumn {
            name: self.name,
            values: vec![values],
            missing: vec![missing.finish()],
        }
    }
}

impl TableBuilder {
    pub(super) fn empty_like(&self) -> Self {
        Self {
            category: self.category.clone(),
            columns: self
                .columns
                .iter()
                .map(|column| ColumnBuilder::new(column.name.clone(), column.column_type()))
                .collect(),
            provenance: Vec::new(),
        }
    }
    pub(super) fn append_provenance(
        &mut self,
        source_name: &str,
        block_code: Option<&str>,
        frame_code: Option<&str>,
    ) {
        self.append_provenance_run(source_name, block_code, frame_code, 1);
    }

    pub(super) fn append_provenance_run(
        &mut self,
        source_name: &str,
        block_code: Option<&str>,
        frame_code: Option<&str>,
        row_count: usize,
    ) {
        if row_count == 0 {
            return;
        }
        if let Some(run) = self.provenance.last_mut()
            && run.source_name == source_name
            && run.block_code.as_deref() == block_code
            && run.frame_code.as_deref() == frame_code
        {
            run.row_count += row_count;
            return;
        }
        self.provenance.push(RowProvenance {
            source_name: source_name.to_owned(),
            block_code: block_code.map(str::to_owned),
            frame_code: frame_code.map(str::to_owned),
            row_count,
        });
    }

    pub(super) fn finish(self) -> CifTable {
        let row_count = self.provenance.iter().map(RowProvenance::row_count).sum();
        CifTable {
            category: self.category,
            columns: self
                .columns
                .into_iter()
                .map(ColumnBuilder::finish)
                .collect(),
            provenance: self.provenance,
            chunk_rows: vec![row_count],
        }
    }
}

impl CifTable {
    /// Return the category name without a leading underscore.
    #[must_use]
    pub fn category(&self) -> &str {
        &self.category
    }

    /// Return projected columns in requested or source order.
    #[must_use]
    pub fn columns(&self) -> &[CifColumn] {
        &self.columns
    }

    /// Find a column by ASCII case-insensitive item name.
    #[must_use]
    pub fn column(&self, name: &str) -> Option<&CifColumn> {
        self.columns
            .iter()
            .find(|column| column.name.eq_ignore_ascii_case(name))
    }

    /// Return the number of projected rows.
    #[must_use]
    pub fn row_count(&self) -> usize {
        self.provenance.iter().map(RowProvenance::row_count).sum()
    }

    pub(super) fn append(&mut self, mut other: Self) -> Result<(), String> {
        if self.category != other.category
            || self.columns.len() != other.columns.len()
            || self
                .columns
                .iter()
                .zip(&other.columns)
                .any(|(left, right)| {
                    !left.name.eq_ignore_ascii_case(&right.name)
                        || left.data_type() != right.data_type()
                })
        {
            return Err("parallel projection chunks have incompatible columns".to_owned());
        }
        if self.chunk_rows == [0] {
            self.chunk_rows.clear();
            for column in &mut self.columns {
                column.values.clear();
                column.missing.clear();
            }
        }
        for (column, mut incoming) in self.columns.iter_mut().zip(other.columns) {
            column.values.append(&mut incoming.values);
            column.missing.append(&mut incoming.missing);
        }
        self.chunk_rows.append(&mut other.chunk_rows);
        for run in other.provenance {
            if let Some(previous) = self.provenance.last_mut()
                && previous.source_name == run.source_name
                && previous.block_code == run.block_code
                && previous.frame_code == run.frame_code
            {
                previous.row_count += run.row_count;
            } else {
                self.provenance.push(run);
            }
        }
        Ok(())
    }

    pub(crate) fn chunk_rows(&self) -> &[usize] {
        &self.chunk_rows
    }

    /// Return run-length encoded source, block, and frame identity.
    ///
    /// Each run's [`RowProvenance::row_count`] gives its consecutive row count.
    #[must_use]
    pub fn provenance(&self) -> &[RowProvenance] {
        &self.provenance
    }

    /// Return source identity for one row.
    #[must_use]
    pub fn row_provenance(&self, row_index: usize) -> Option<&RowProvenance> {
        let mut remaining = row_index;
        self.provenance.iter().find(|run| {
            if remaining < run.row_count {
                true
            } else {
                remaining -= run.row_count;
                false
            }
        })
    }
}

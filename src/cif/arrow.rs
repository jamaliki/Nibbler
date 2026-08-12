use std::collections::HashMap;
use std::error::Error;
use std::fmt::{self, Display, Formatter};
use std::sync::Arc;

use arrow_array::ffi_stream::FFI_ArrowArrayStream;
use arrow_array::{ArrayRef, RecordBatch, RecordBatchIterator, StringArray, StructArray};
use arrow_schema::extension::{EXTENSION_TYPE_METADATA_KEY, EXTENSION_TYPE_NAME_KEY};
use arrow_schema::{ArrowError as ApacheArrowError, DataType, Field, Fields, Schema};

use super::table::{CifColumn, CifTable};

/// How lossless CIF missing states cross the Arrow interface.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum ArrowMissingPolicy {
    /// Map both CIF missing states to Arrow null and record counts in field metadata.
    #[default]
    Collapse,
    /// Add a `__missing_kind` UInt8 companion for every projected value column.
    Columns,
    /// Emit struct storage annotated as the `nibbler.cif_missing` extension type.
    Extension,
}

/// A failure while constructing an Arrow record batch or stream.
#[derive(Debug)]
pub struct ArrowExportError(ApacheArrowError);

impl Display for ArrowExportError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(formatter, "Arrow export failed: {}", self.0)
    }
}

impl Error for ArrowExportError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.0)
    }
}

impl From<ApacheArrowError> for ArrowExportError {
    fn from(error: ApacheArrowError) -> Self {
        Self(error)
    }
}

/// Export a table as one owned Arrow C stream.
///
/// Consumers assume ownership through the Arrow C Stream release callback. Source
/// provenance columns are included when `include_provenance` is true.
///
/// # Errors
///
/// Returns [`ArrowExportError`] if Arrow rejects the constructed schema or arrays.
pub fn export_arrow_stream(
    table: &CifTable,
    missing: ArrowMissingPolicy,
    include_provenance: bool,
) -> Result<FFI_ArrowArrayStream, ArrowExportError> {
    let batches = record_batches(table, missing, include_provenance)?;
    let schema = batches[0].schema();
    let reader = RecordBatchIterator::new(batches.into_iter().map(Ok), schema);
    Ok(FFI_ArrowArrayStream::new(Box::new(reader)))
}

fn record_batches(
    table: &CifTable,
    missing: ArrowMissingPolicy,
    include_provenance: bool,
) -> Result<Vec<RecordBatch>, ArrowExportError> {
    let mut batches = Vec::with_capacity(table.chunk_rows().len());
    let missing_counts = table
        .columns()
        .iter()
        .map(count_missing)
        .collect::<Vec<_>>();
    let mut row_start = 0;
    for (chunk_index, row_count) in table.chunk_rows().iter().copied().enumerate() {
        let mut fields = Vec::new();
        let mut arrays = Vec::new();
        for (column, counts) in table.columns().iter().zip(&missing_counts) {
            append_column(
                column,
                chunk_index,
                *counts,
                missing,
                &mut fields,
                &mut arrays,
            );
        }
        if include_provenance {
            append_provenance(table, row_start, row_count, &mut fields, &mut arrays);
        }
        let mut metadata = HashMap::new();
        metadata.insert(
            "nibbler:cif_category".to_owned(),
            table.category().to_owned(),
        );
        let schema = Arc::new(Schema::new_with_metadata(fields, metadata));
        batches.push(RecordBatch::try_new(schema, arrays)?);
        row_start += row_count;
    }
    Ok(batches)
}

fn append_column(
    column: &CifColumn,
    chunk_index: usize,
    counts: (usize, usize),
    missing: ArrowMissingPolicy,
    fields: &mut Vec<Field>,
    arrays: &mut Vec<ArrayRef>,
) {
    let values = column.array(chunk_index);
    let data_type = column.data_type();
    let kinds = Arc::new(column.missing_array(chunk_index).clone());
    let mut metadata = HashMap::new();
    metadata.insert("nibbler:cif_unknown_count".to_owned(), counts.0.to_string());
    metadata.insert(
        "nibbler:cif_not_applicable_count".to_owned(),
        counts.1.to_string(),
    );

    match missing {
        ArrowMissingPolicy::Collapse => {
            fields.push(Field::new(column.name(), data_type.clone(), true).with_metadata(metadata));
            arrays.push(values);
        }
        ArrowMissingPolicy::Columns => {
            fields.push(Field::new(column.name(), data_type.clone(), true).with_metadata(metadata));
            arrays.push(values);
            fields.push(Field::new(
                format!("{}__missing_kind", column.name()),
                DataType::UInt8,
                false,
            ));
            arrays.push(kinds);
        }
        ArrowMissingPolicy::Extension => {
            metadata.insert(
                EXTENSION_TYPE_NAME_KEY.to_owned(),
                "nibbler.cif_missing".to_owned(),
            );
            metadata.insert(
                EXTENSION_TYPE_METADATA_KEY.to_owned(),
                r#"{"unknown":1,"not_applicable":2}"#.to_owned(),
            );
            let storage_fields = Fields::from(vec![
                Field::new("value", data_type, true),
                Field::new("kind", DataType::UInt8, false),
            ]);
            let storage = StructArray::new(storage_fields.clone(), vec![values, kinds], None);
            fields.push(
                Field::new(column.name(), DataType::Struct(storage_fields), false)
                    .with_metadata(metadata),
            );
            arrays.push(Arc::new(storage));
        }
    }
}

fn append_provenance(
    table: &CifTable,
    row_start: usize,
    row_count: usize,
    fields: &mut Vec<Field>,
    arrays: &mut Vec<ArrayRef>,
) {
    fields.extend([
        Field::new("_nibbler_source", DataType::Utf8, false),
        Field::new("_nibbler_block", DataType::Utf8, true),
        Field::new("_nibbler_frame", DataType::Utf8, true),
    ]);
    let selected = provenance_runs(table, row_start, row_count);
    arrays.extend([
        Arc::new(
            selected
                .iter()
                .flat_map(|(row, count)| std::iter::repeat_n(Some(row.source_name()), *count))
                .collect::<StringArray>(),
        ) as ArrayRef,
        Arc::new(
            selected
                .iter()
                .flat_map(|(row, count)| std::iter::repeat_n(row.block_code(), *count))
                .collect::<StringArray>(),
        ),
        Arc::new(
            selected
                .iter()
                .flat_map(|(row, count)| std::iter::repeat_n(row.frame_code(), *count))
                .collect::<StringArray>(),
        ),
    ]);
}

fn provenance_runs(
    table: &CifTable,
    row_start: usize,
    row_count: usize,
) -> Vec<(&super::table::RowProvenance, usize)> {
    let mut skip = row_start;
    let mut remaining = row_count;
    let mut selected = Vec::new();
    for run in table.provenance() {
        if skip >= run.row_count() {
            skip -= run.row_count();
            continue;
        }
        let count = (run.row_count() - skip).min(remaining);
        selected.push((run, count));
        remaining -= count;
        skip = 0;
        if remaining == 0 {
            break;
        }
    }
    selected
}

fn count_missing(column: &CifColumn) -> (usize, usize) {
    (0..column.chunk_count()).fold((0, 0), |counts, chunk_index| {
        (0..column.missing_array(chunk_index).len()).fold(
            counts,
            |(unknown, not_applicable), row_index| match column
                .missing_array(chunk_index)
                .value(row_index)
            {
                1 => (unknown + 1, not_applicable),
                2 => (unknown, not_applicable + 1),
                _ => (unknown, not_applicable),
            },
        )
    })
}

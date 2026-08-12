//! BinaryCIF container and logical-document decoding.

use crate::cif::document::{ColumnValues, LoopColumn};
use crate::cif::{CifBlock, CifDocument, CifEntry, CifLoop, Limits};

use super::codec::decode_data;
use super::error::{BinaryCifError, BinaryCifErrorCode};
use super::model::{BinaryCategory, BinaryFile, Encoding};

const MAX_BINARY_EXPANSION_RATIO: usize = 1_000;

pub(super) fn decode_file(bytes: &[u8]) -> Result<BinaryFile<'_>, BinaryCifError> {
    let file: BinaryFile<'_> = rmp_serde::from_slice(bytes).map_err(|error| {
        BinaryCifError::new(
            BinaryCifErrorCode::Container,
            format!("invalid BinaryCIF MessagePack container: {error}"),
        )
    })?;
    if !file.version.starts_with("0.3.") {
        return Err(BinaryCifError::new(
            BinaryCifErrorCode::Version,
            format!(
                "unsupported BinaryCIF version {:?}; expected 0.3.x",
                file.version
            ),
        ));
    }
    if file.data_blocks.is_empty() {
        return Err(BinaryCifError::new(
            BinaryCifErrorCode::Shape,
            "BinaryCIF contains no data blocks",
        ));
    }
    validate_file_shape(&file, bytes.len())?;
    Ok(file)
}

fn validate_file_shape(file: &BinaryFile<'_>, input_bytes: usize) -> Result<(), BinaryCifError> {
    let limits = Limits::default();
    let expanded_values = input_bytes.saturating_mul(MAX_BINARY_EXPANSION_RATIO);
    let value_limit = limits.values.min(expanded_values);
    let row_limit = limits.rows.min(expanded_values);
    if file.data_blocks.len() > limits.blocks {
        return Err(shape_error(
            "BinaryCIF block count exceeds the resource limit",
        ));
    }
    let mut categories = 0_usize;
    let mut tags = 0_usize;
    let mut rows = 0_usize;
    let mut values = 0_usize;
    for block in &file.data_blocks {
        categories = categories
            .checked_add(block.categories.len())
            .ok_or_else(|| shape_error("BinaryCIF category count overflows usize"))?;
        if categories > limits.loops {
            return Err(shape_error(
                "BinaryCIF category count exceeds the resource limit",
            ));
        }
        for category in &block.categories {
            if category.columns.len() > limits.loop_columns {
                return Err(shape_error(
                    "BinaryCIF category width exceeds the resource limit",
                ));
            }
            tags = tags
                .checked_add(category.columns.len())
                .ok_or_else(|| shape_error("BinaryCIF column count overflows usize"))?;
            rows = rows
                .checked_add(category.row_count)
                .ok_or_else(|| shape_error("BinaryCIF row count overflows usize"))?;
            let category_values = category
                .row_count
                .checked_mul(category.columns.len())
                .ok_or_else(|| shape_error("BinaryCIF category size overflows usize"))?;
            values = values
                .checked_add(category_values)
                .ok_or_else(|| shape_error("BinaryCIF value count overflows usize"))?;
            if tags > limits.tags || rows > row_limit || values > value_limit {
                return Err(shape_error(
                    "BinaryCIF logical size exceeds the resource limit",
                ));
            }
            for column in &category.columns {
                validate_encoding_sizes(&column.data.encoding, value_limit)?;
                if let Some(mask) = &column.mask {
                    validate_encoding_sizes(&mask.encoding, value_limit)?;
                }
            }
        }
    }
    Ok(())
}

fn validate_encoding_sizes(
    encodings: &[Encoding<'_>],
    value_limit: usize,
) -> Result<(), BinaryCifError> {
    for encoding in encodings {
        match encoding {
            Encoding::RunLength { src_size, .. } | Encoding::IntegerPacking { src_size, .. }
                if *src_size > value_limit =>
            {
                return Err(shape_error(
                    "BinaryCIF encoding expansion exceeds the resource limit",
                ));
            }
            Encoding::StringArray {
                data_encoding,
                offset_encoding,
                ..
            } => {
                validate_encoding_sizes(data_encoding, value_limit)?;
                validate_encoding_sizes(offset_encoding, value_limit)?;
            }
            _ => {}
        }
    }
    Ok(())
}

/// Decode a BinaryCIF 0.3 MessagePack document into the shared logical model.
///
/// # Errors
///
/// Returns a structured error for malformed containers, unsupported versions or
/// encoding chains, invalid string dictionaries, and inconsistent row counts.
pub fn decode_binary(bytes: &[u8]) -> Result<CifDocument, BinaryCifError> {
    let file = decode_file(bytes)?;
    let mut blocks = Vec::with_capacity(file.data_blocks.len());
    for block in file.data_blocks {
        if block.header.is_empty() {
            return Err(shape_error("BinaryCIF data block has an empty header"));
        }
        let mut entries = Vec::with_capacity(block.categories.len());
        for category in block.categories {
            if let Some(cif_loop) = decode_category(category)? {
                entries.push(CifEntry::Loop(cif_loop));
            }
        }
        blocks.push(CifBlock::data(block.header, entries));
    }
    Ok(CifDocument::new(blocks))
}

fn decode_category(category: BinaryCategory<'_>) -> Result<Option<CifLoop>, BinaryCifError> {
    let category_name = normalize_category_name(&category.name)?;
    if category.row_count == 0 {
        return Ok(None);
    }
    if category.columns.is_empty() {
        return Err(shape_error(format!(
            "category {category_name:?} has rows but no columns"
        )));
    }
    let mut tags = Vec::with_capacity(category.columns.len());
    let mut columns = Vec::with_capacity(category.columns.len());
    for column in category.columns {
        if column.name.is_empty() || column.name.contains('.') {
            return Err(shape_error(format!(
                "category {category_name:?} has invalid column name {:?}",
                column.name
            )));
        }
        tags.push(format!("_{category_name}.{}", column.name));
        let values = decode_data(column.data)?;
        if values.len() != category.row_count {
            return Err(shape_error(format!(
                "_{category_name}.{} decodes {} rows, expected {}",
                column.name,
                values.len(),
                category.row_count
            )));
        }
        let mask = column.mask.map(decode_mask).transpose()?;
        if mask
            .as_ref()
            .is_some_and(|mask| mask.len() != category.row_count)
        {
            return Err(shape_error(format!(
                "_{category_name}.{} mask decodes {} rows, expected {}",
                column.name,
                mask.as_ref().map_or(0, Vec::len),
                category.row_count
            )));
        }
        if let ColumnValues::Floats(numbers) = &values
            && numbers.iter().enumerate().any(|(row_index, number)| {
                mask.as_ref().map_or(0, |mask| mask[row_index]) == 0 && !number.is_finite()
            })
        {
            return Err(BinaryCifError::new(
                BinaryCifErrorCode::Encoding,
                "BinaryCIF present floating-point values must be finite",
            ));
        }
        columns.push(LoopColumn { values, mask });
    }
    Ok(Some(CifLoop::from_columns(
        tags,
        columns,
        category.row_count,
    )))
}

pub(super) fn decode_mask(data: super::model::BinaryData) -> Result<Vec<u8>, BinaryCifError> {
    let ColumnValues::Integers(values) = decode_data(data)? else {
        return Err(shape_error("BinaryCIF mask must decode to integers"));
    };
    values
        .into_iter()
        .map(|value| match value {
            0..=2 => Ok(value as u8),
            _ => Err(shape_error(format!(
                "BinaryCIF mask contains invalid value {value}"
            ))),
        })
        .collect()
}

pub(super) fn normalize_category_name(name: &str) -> Result<&str, BinaryCifError> {
    let name = name.strip_prefix('_').unwrap_or(name);
    if name.is_empty() || name.contains('.') {
        return Err(shape_error(format!(
            "invalid BinaryCIF category name {name:?}"
        )));
    }
    Ok(name)
}

pub(super) fn shape_error(message: impl Into<String>) -> BinaryCifError {
    BinaryCifError::new(BinaryCifErrorCode::Shape, message)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::expect_used)]

    use std::borrow::Cow;

    use super::{decode_category, validate_encoding_sizes};
    use crate::cif::CifValueRef;
    use crate::cif::binary::model::{BinaryCategory, BinaryColumn, BinaryData, Encoding};

    #[test]
    fn rejects_declared_encoding_expansion_before_allocation() {
        let encodings = [Encoding::RunLength {
            src_type: 3,
            src_size: usize::MAX,
        }];
        assert!(validate_encoding_sizes(&encodings, 1_000).is_err());
    }

    #[test]
    fn ignores_non_finite_binary_placeholders_for_missing_values() {
        let mut data = Vec::new();
        data.extend_from_slice(&f64::NAN.to_le_bytes());
        data.extend_from_slice(&1.25_f64.to_le_bytes());
        let category = BinaryCategory {
            name: "_measurement".to_owned(),
            row_count: 2,
            columns: vec![BinaryColumn {
                name: "value".to_owned(),
                data: BinaryData {
                    encoding: vec![Encoding::ByteArray { data_type: 33 }],
                    data: Cow::Owned(data),
                },
                mask: Some(BinaryData {
                    encoding: vec![Encoding::ByteArray { data_type: 4 }],
                    data: Cow::Owned(vec![1, 0]),
                }),
            }],
        };

        let cif_loop = decode_category(category)
            .expect("masked non-finite placeholder is valid")
            .expect("non-empty category produces a loop");
        assert_eq!(cif_loop.value(0, 0), Some(CifValueRef::NotApplicable));
        assert_eq!(
            cif_loop.value(1, 0),
            Some(CifValueRef::Float(1.25, None, None))
        );
    }
}

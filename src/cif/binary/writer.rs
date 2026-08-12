//! Deterministic BinaryCIF serialization.

use std::borrow::Cow;
use std::collections::HashMap;

use crate::cif::{BlockKind, CifDocument, CifEntry, CifValue};

use super::document::normalize_category_name;
use super::error::{BinaryCifError, BinaryCifErrorCode};
use super::model::{BinaryBlock, BinaryCategory, BinaryColumn, BinaryData, BinaryFile, Encoding};

struct CategoryValues {
    name: String,
    row_count: usize,
    columns: Vec<(String, Vec<CifValue>)>,
    scalar: bool,
}

/// Encode the shared logical model as deterministic BinaryCIF 0.3 MessagePack.
///
/// # Errors
///
/// Returns an error for global blocks, save frames, mixed scalar/loop forms of one
/// category, cross-category loops, non-finite numbers, or container serialization
/// failures.
pub fn encode_binary(document: &CifDocument) -> Result<Vec<u8>, BinaryCifError> {
    if document.blocks().is_empty() {
        return Err(unrepresentable("cannot encode an empty CIF document"));
    }
    let mut data_blocks = Vec::with_capacity(document.blocks().len());
    for block in document.blocks() {
        if block.kind() != BlockKind::Data {
            return Err(unrepresentable(
                "BinaryCIF 0.3 cannot represent a global_ block",
            ));
        }
        let header = block
            .code()
            .filter(|value| !value.is_empty())
            .ok_or_else(|| unrepresentable("data block has no header"))?
            .to_owned();
        let categories = collect_categories(block.entries())?
            .into_iter()
            .map(|category| {
                let columns = category
                    .columns
                    .into_iter()
                    .map(|(name, values)| encode_column(name, values))
                    .collect::<Result<Vec<_>, _>>()?;
                Ok(BinaryCategory {
                    name: format!("_{}", category.name),
                    row_count: category.row_count,
                    columns,
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        data_blocks.push(BinaryBlock { header, categories });
    }
    let file = BinaryFile {
        version: "0.3.0".to_owned(),
        encoder: format!("Nibbler {}", env!("CARGO_PKG_VERSION")),
        data_blocks,
    };
    rmp_serde::to_vec_named(&file).map_err(|error| {
        BinaryCifError::new(
            BinaryCifErrorCode::Container,
            format!("failed to encode BinaryCIF MessagePack: {error}"),
        )
    })
}

fn collect_categories(entries: &[CifEntry]) -> Result<Vec<CategoryValues>, BinaryCifError> {
    let mut categories: Vec<CategoryValues> = Vec::new();
    for entry in entries {
        match entry {
            CifEntry::Frame(_) => {
                return Err(unrepresentable(
                    "BinaryCIF 0.3 cannot represent save frames",
                ));
            }
            CifEntry::Item(item) => {
                let (category, column) = split_tag(item.tag())?;
                if let Some(existing) = categories.iter_mut().find(|value| value.name == category) {
                    if !existing.scalar {
                        return Err(unrepresentable(format!(
                            "category {category:?} mixes scalar and loop forms"
                        )));
                    }
                    existing
                        .columns
                        .push((column.to_owned(), vec![item.value().clone()]));
                } else {
                    categories.push(CategoryValues {
                        name: category.to_owned(),
                        row_count: 1,
                        columns: vec![(column.to_owned(), vec![item.value().clone()])],
                        scalar: true,
                    });
                }
            }
            CifEntry::Loop(cif_loop) => {
                let Some(first_tag) = cif_loop.tags().first() else {
                    return Err(unrepresentable("cannot encode a loop without tags"));
                };
                let (category, _) = split_tag(first_tag)?;
                if categories.iter().any(|value| value.name == category) {
                    return Err(unrepresentable(format!(
                        "category {category:?} occurs in multiple CIF entries"
                    )));
                }
                if cif_loop.row_count() == 0 {
                    return Err(unrepresentable("cannot encode an empty CIF loop"));
                }
                let mut columns = Vec::with_capacity(cif_loop.column_count());
                for (column_index, tag) in cif_loop.tags().iter().enumerate() {
                    let (candidate, column) = split_tag(tag)?;
                    if candidate != category {
                        return Err(unrepresentable(format!(
                            "one loop mixes categories {category:?} and {candidate:?}"
                        )));
                    }
                    let values = (0..cif_loop.row_count())
                        .filter_map(|row_index| cif_loop.value(row_index, column_index))
                        .map(|value| value.to_owned())
                        .collect::<Vec<_>>();
                    if values.len() != cif_loop.row_count() {
                        return Err(unrepresentable("cannot encode an incomplete CIF loop"));
                    }
                    columns.push((column.to_owned(), values));
                }
                categories.push(CategoryValues {
                    name: category.to_owned(),
                    row_count: cif_loop.row_count(),
                    columns,
                    scalar: false,
                });
            }
        }
    }
    Ok(categories)
}

fn encode_column(
    name: String,
    values: Vec<CifValue>,
) -> Result<BinaryColumn<'static>, BinaryCifError> {
    let mask_values = values
        .iter()
        .map(|value| match value {
            CifValue::NotApplicable => 1,
            CifValue::Unknown => 2,
            _ => 0,
        })
        .collect::<Vec<_>>();
    let mask = mask_values
        .iter()
        .any(|value| *value != 0)
        .then(|| BinaryData {
            encoding: vec![Encoding::ByteArray { data_type: 4 }],
            data: Cow::Owned(mask_values),
        });
    let present = values
        .iter()
        .filter(|value| !matches!(value, CifValue::Unknown | CifValue::NotApplicable))
        .collect::<Vec<_>>();
    let data = if present.iter().all(
        |value| matches!(value, CifValue::Integer(number, _) if i32::try_from(*number).is_ok()),
    ) {
        let mut data = Vec::with_capacity(values.len() * 4);
        for value in &values {
            let number = match value {
                CifValue::Integer(number, _) => i32::try_from(*number)
                    .map_err(|_| unrepresentable("integer exceeds BinaryCIF Int32"))?,
                CifValue::Unknown | CifValue::NotApplicable => 0,
                _ => return Err(unrepresentable("column has mixed numeric value types")),
            };
            data.extend_from_slice(&number.to_le_bytes());
        }
        BinaryData {
            encoding: vec![Encoding::ByteArray { data_type: 3 }],
            data: Cow::Owned(data),
        }
    } else if present
        .iter()
        .all(|value| matches!(value, CifValue::Float(number, None, _) if number.is_finite()))
    {
        let mut data = Vec::with_capacity(values.len() * 8);
        for value in &values {
            let number = match value {
                CifValue::Float(number, None, _) if number.is_finite() => *number,
                CifValue::Unknown | CifValue::NotApplicable => 0.0,
                CifValue::Float(..) => {
                    return Err(unrepresentable(
                        "BinaryCIF present floating-point values must be finite",
                    ));
                }
                _ => return Err(unrepresentable("column has mixed numeric value types")),
            };
            data.extend_from_slice(&number.to_le_bytes());
        }
        BinaryData {
            encoding: vec![Encoding::ByteArray { data_type: 33 }],
            data: Cow::Owned(data),
        }
    } else {
        encode_strings(&values)?
    };
    Ok(BinaryColumn { name, data, mask })
}

fn encode_strings(values: &[CifValue]) -> Result<BinaryData<'static>, BinaryCifError> {
    let mut dictionary = HashMap::new();
    let mut string_data = String::new();
    let mut offsets = vec![0_i32];
    let mut indices = Vec::with_capacity(values.len());
    for value in values {
        let text = match value {
            CifValue::Unknown | CifValue::NotApplicable => String::new(),
            CifValue::Text(text) => text.as_str().to_owned(),
            _ => super::super::writer::format_value(value).map_err(|error| {
                unrepresentable(format!("cannot format BinaryCIF string value: {error}"))
            })?,
        };
        let index = if let Some(index) = dictionary.get(&text) {
            *index
        } else {
            let index = i32::try_from(dictionary.len())
                .map_err(|_| unrepresentable("StringArray dictionary exceeds Int32"))?;
            string_data.push_str(&text);
            offsets.push(
                i32::try_from(string_data.len())
                    .map_err(|_| unrepresentable("StringArray data exceeds Int32"))?,
            );
            dictionary.insert(text, index);
            index
        };
        indices.push(index);
    }
    let mut data = Vec::with_capacity(indices.len() * 4);
    for index in indices {
        data.extend_from_slice(&index.to_le_bytes());
    }
    let mut offset_bytes = Vec::with_capacity(offsets.len() * 4);
    for offset in offsets {
        offset_bytes.extend_from_slice(&offset.to_le_bytes());
    }
    Ok(BinaryData {
        encoding: vec![Encoding::StringArray {
            data_encoding: vec![Encoding::ByteArray { data_type: 3 }],
            string_data,
            offset_encoding: vec![Encoding::ByteArray { data_type: 3 }],
            offsets: Cow::Owned(offset_bytes),
        }],
        data: Cow::Owned(data),
    })
}

fn split_tag(tag: &str) -> Result<(&str, &str), BinaryCifError> {
    let Some((category, column)) = tag
        .strip_prefix('_')
        .and_then(|value| value.split_once('.'))
    else {
        return Err(unrepresentable(format!("invalid CIF tag {tag:?}")));
    };
    normalize_category_name(category)?;
    if column.is_empty() || column.contains('.') {
        return Err(unrepresentable(format!("invalid CIF tag {tag:?}")));
    }
    Ok((category, column))
}

fn unrepresentable(message: impl Into<String>) -> BinaryCifError {
    BinaryCifError::new(BinaryCifErrorCode::Unrepresentable, message)
}

#[cfg(test)]
mod tests {
    use super::{Encoding, encode_column};
    use crate::cif::{CifValue, StandardUncertainty};

    #[test]
    fn standard_uncertainty_uses_lossless_string_encoding() {
        let result = encode_column(
            "measurement".to_owned(),
            vec![CifValue::Float(
                1.25,
                Some(StandardUncertainty::new(3)),
                None,
            )],
        );
        assert!(matches!(
            result,
            Ok(column) if matches!(column.data.encoding.as_slice(), [Encoding::StringArray { .. }])
        ));
    }
}

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
    let mut block_headers = Vec::with_capacity(file.data_blocks.len());
    let mut category_names = Vec::new();
    let mut column_names = Vec::new();
    for block in &file.data_blocks {
        require_text_name(&block.header, "data block header")?;
        block_headers.push(block.header.as_str());
        categories = categories
            .checked_add(block.categories.len())
            .ok_or_else(|| shape_error("BinaryCIF category count overflows usize"))?;
        if categories > limits.loops {
            return Err(shape_error(
                "BinaryCIF category count exceeds the resource limit",
            ));
        }
        category_names.clear();
        for category in &block.categories {
            require_text_name(&category.name, "category name")?;
            let category_name = normalize_category_name(&category.name)?;
            category_names.push(category_name);
            if category.row_count > 0 && category.columns.is_empty() {
                return Err(shape_error(format!(
                    "category {category_name:?} has rows but no columns"
                )));
            }
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
            column_names.clear();
            // Zero-row categories are checked too: the document decoder drops them, but
            // projection takes its column layout from the first matching category.
            for column in &category.columns {
                require_column_name(category_name, &column.name)?;
                column_names.push(column.name.as_str());
                validate_encoding_sizes(&column.data.encoding, value_limit)?;
                if let Some(mask) = &column.mask {
                    validate_encoding_sizes(&mask.encoding, value_limit)?;
                }
            }
            require_unique_names(&mut column_names, "column name")?;
        }
        require_unique_names(&mut category_names, "category name")?;
    }
    require_unique_names(&mut block_headers, "data block header")
}

/// Require a BinaryCIF name to have a text CIF spelling as part of one bare token.
///
/// The text lexer ends a token at space, tab, CR, or LF and rejects every other control
/// character, so a name containing either would decode to a document whose canonical
/// text does not parse. Strict text CIF also rejects an empty block code, and an empty
/// column name would leave its tag without an item name. Checking here covers both
/// document decoding and projection.
fn require_text_name(name: &str, kind: &str) -> Result<(), BinaryCifError> {
    if name.is_empty() {
        return Err(shape_error(format!("BinaryCIF {kind} is empty")));
    }
    if name.chars().any(breaks_bare_token) {
        return Err(shape_error(format!(
            "BinaryCIF {kind} {name:?} contains whitespace or a control character"
        )));
    }
    Ok(())
}

/// Require a BinaryCIF column name to have a text CIF spelling that splits back out of
/// its `_category.column` tag.
///
/// Text CIF splits a tag at its first '.', so the name must not contain one. That test
/// shares one scan with [`require_text_name`]'s because shape validation runs on every
/// read: a second scan of each column name took validation of PDB 1CRN from 14.4 µs to
/// 16.8 µs, while one scan stays within noise of 14.7 µs.
fn require_column_name(category_name: &str, name: &str) -> Result<(), BinaryCifError> {
    if name.is_empty()
        || name
            .chars()
            .any(|character| character == '.' || breaks_bare_token(character))
    {
        // Report a missing text spelling as such; otherwise the name contains '.'.
        require_text_name(name, "column name")?;
        return Err(shape_error(format!(
            "category {category_name:?} has invalid column name {name:?}"
        )));
    }
    Ok(())
}

fn breaks_bare_token(character: char) -> bool {
    character == ' ' || character.is_control()
}

/// Require the names of one scope to differ under ASCII case folding.
///
/// Text CIF folds block codes and tags this way, so a repeated block header, category
/// within one block, or column within one category would decode to a document whose
/// canonical text fails with `DuplicateBlock` or `DuplicateTag`.
///
/// Every BinaryCIF read runs this. Sorting the borrowed names in place adds about 4 µs to
/// shape validation of PDB 1CRN, against about 26 µs for a set of folded copies. Folding
/// preserves byte length, so ordering by length first settles most comparisons without
/// reading bytes and still places folded equals next to each other.
fn require_unique_names(names: &mut [&str], kind: &str) -> Result<(), BinaryCifError> {
    names.sort_unstable_by(|left, right| {
        left.len().cmp(&right.len()).then_with(|| {
            let left = left.bytes().map(|byte| byte.to_ascii_lowercase());
            left.cmp(right.bytes().map(|byte| byte.to_ascii_lowercase()))
        })
    });
    for pair in names.windows(2) {
        if let [left, right] = pair
            && left.eq_ignore_ascii_case(right)
        {
            return Err(shape_error(format!(
                "duplicate BinaryCIF {kind} {left:?} and {right:?}"
            )));
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
/// Returns a structured error for malformed containers, names without a text CIF
/// spelling (including empty names and column names containing `.`), block headers,
/// categories within a block, or columns within a category repeated under ASCII case
/// folding, categories with rows but no columns, unsupported versions or encoding
/// chains, invalid string dictionaries, and inconsistent row counts. Categories with no
/// rows are checked like any other and then omitted from the document.
pub fn decode_binary(bytes: &[u8]) -> Result<CifDocument, BinaryCifError> {
    let file = decode_file(bytes)?;
    let mut blocks = Vec::with_capacity(file.data_blocks.len());
    for block in file.data_blocks {
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
    let mut tags = Vec::with_capacity(category.columns.len());
    let mut columns = Vec::with_capacity(category.columns.len());
    for column in category.columns {
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

    use super::{BinaryCifErrorCode, decode_binary, decode_category, validate_encoding_sizes};
    use crate::cif::binary::model::{
        BinaryBlock, BinaryCategory, BinaryColumn, BinaryData, BinaryFile, Encoding,
    };
    use crate::cif::{CifValueRef, ProjectionPlan, parse, project_binary, write_canonical};

    /// Serialize one single-row integer column per `(header, category, column)` triple.
    ///
    /// A triple joins the previous block or category only when it spells the same name
    /// exactly, so case variants and repeats of non-adjacent names stay separate.
    fn single_row_container(columns: &[(&str, &str, &str)]) -> Vec<u8> {
        let mut data_blocks: Vec<BinaryBlock<'_>> = Vec::new();
        for &(header, category, column) in columns {
            if data_blocks
                .last()
                .is_none_or(|block| block.header != header)
            {
                data_blocks.push(BinaryBlock {
                    header: header.to_owned(),
                    categories: Vec::new(),
                });
            }
            let block = data_blocks.last_mut().expect("a block was just ensured");
            if block
                .categories
                .last()
                .is_none_or(|existing| existing.name != category)
            {
                block.categories.push(BinaryCategory {
                    name: category.to_owned(),
                    row_count: 1,
                    columns: Vec::new(),
                });
            }
            let category = block.categories.last_mut().expect("a category was ensured");
            category.columns.push(uint8_column(column, vec![6]));
        }
        serialize(data_blocks)
    }

    /// Serialize `categories` as the only block of a container.
    fn single_block_container(categories: Vec<BinaryCategory<'static>>) -> Vec<u8> {
        serialize(vec![BinaryBlock {
            header: "glycan".to_owned(),
            categories,
        }])
    }

    fn uint8_column(name: &str, values: Vec<u8>) -> BinaryColumn<'static> {
        BinaryColumn {
            name: name.to_owned(),
            data: BinaryData {
                encoding: vec![Encoding::ByteArray { data_type: 4 }],
                data: Cow::Owned(values),
            },
            mask: None,
        }
    }

    fn serialize(data_blocks: Vec<BinaryBlock<'_>>) -> Vec<u8> {
        let file = BinaryFile {
            version: "0.3.0".to_owned(),
            encoder: "test".to_owned(),
            data_blocks,
        };
        rmp_serde::to_vec_named(&file).expect("test container serializes")
    }

    fn assert_shape_error_in_documents_and_projections(bytes: &[u8]) {
        assert_eq!(
            decode_binary(bytes).map_err(|error| error.code()).err(),
            Some(BinaryCifErrorCode::Shape)
        );
        let plan = ProjectionPlan::new("entry").expect("static category is valid");
        assert_eq!(
            project_binary(bytes, plan)
                .map_err(|error| error.code())
                .err(),
            Some(BinaryCifErrorCode::Shape)
        );
    }

    #[test]
    fn rejects_names_without_a_text_spelling_in_documents_and_projections() {
        // Minimized from a robustness seed: one mutated byte turned `_atom_type` into
        // `\0atom_type`, which decoded and projected but wrote text the lexer rejects.
        for column in [
            ("glycan", "\0atom_type", "symbol"),
            ("glycan", "_atom_type", "type symbol"),
            ("gly\u{85}can", "_atom_type", "symbol"),
        ] {
            assert_shape_error_in_documents_and_projections(&single_row_container(&[column]));
        }

        let document = decode_binary(&single_row_container(&[("glycan", "_atom_typé", "symbol")]))
            .expect("printable non-ASCII names have a text spelling");
        let text = write_canonical(&document).expect("decoded document is writable");
        assert!(parse(text.as_bytes()).is_ok());
    }

    #[test]
    fn rejects_names_that_repeat_under_text_case_folding() {
        // Minimized from a robustness seed: one mutated byte turned `Cartn_y` into
        // `Cartn_x`, which decoded and projected but wrote text with a duplicate tag.
        for columns in [
            [
                ("glycan", "_atom_type", "symbol"),
                ("GLYCAN", "_atom_type", "symbol"),
            ],
            [
                ("glycan", "_atom_type", "symbol"),
                ("glycan", "ATOM_TYPE", "symbol"),
            ],
            [
                ("glycan", "_atom_type", "symbol"),
                ("glycan", "_atom_type", "Symbol"),
            ],
        ] {
            assert_shape_error_in_documents_and_projections(&single_row_container(&columns));
        }

        let repeated_across_blocks = single_row_container(&[
            ("glycan", "_atom_type", "symbol"),
            ("ligand", "_ATOM_TYPE", "symbol"),
        ]);
        let plan = ProjectionPlan::new("atom_type").expect("static category is valid");
        let table = project_binary(&repeated_across_blocks, plan)
            .expect("one category may repeat across distinct blocks");
        assert_eq!(table.row_count(), 2);
    }

    #[test]
    fn rejects_empty_block_headers_in_documents_and_projections() {
        // Text CIF rejects a bare `data_` with an empty block code.
        assert_shape_error_in_documents_and_projections(&single_row_container(&[(
            "",
            "_atom_type",
            "symbol",
        )]));
    }

    #[test]
    fn rejects_column_names_that_do_not_split_out_of_a_tag_in_documents_and_projections() {
        for column in [("glycan", "_entry", ""), ("glycan", "_entry", "a.b")] {
            assert_shape_error_in_documents_and_projections(&single_row_container(&[column]));
        }
    }

    #[test]
    fn rejects_rows_without_columns_in_documents_and_projections() {
        assert_shape_error_in_documents_and_projections(&single_block_container(vec![
            BinaryCategory {
                name: "_entry".to_owned(),
                row_count: 1,
                columns: Vec::new(),
            },
        ]));
    }

    #[test]
    fn checks_column_names_of_zero_row_categories() {
        // A zero-row category decodes to no loop, but projection still takes its column
        // layout from it, so its names must form text tags like any other category's.
        let zero_rows = |column| {
            single_block_container(vec![BinaryCategory {
                name: "_entry".to_owned(),
                row_count: 0,
                columns: vec![uint8_column(column, Vec::new())],
            }])
        };
        assert_shape_error_in_documents_and_projections(&zero_rows("a.b"));

        let accepted = zero_rows("id");
        let document = decode_binary(&accepted).expect("a zero-row category is accepted");
        assert!(document.blocks()[0].entries().is_empty());
        let plan = ProjectionPlan::new("entry").expect("static category is valid");
        let table = project_binary(&accepted, plan).expect("a zero-row category projects");
        assert_eq!(table.row_count(), 0);
    }

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

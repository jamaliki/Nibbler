use std::borrow::Cow;

use super::codec::{DecodedColumn, decode_data};
use super::model::{BinaryData, Encoding};

fn bytes(data: Vec<u8>, encoding: Vec<Encoding<'static>>) -> BinaryData<'static> {
    BinaryData {
        encoding,
        data: Cow::Owned(data),
    }
}

fn int32(values: &[i32]) -> Vec<u8> {
    values
        .iter()
        .flat_map(|value| value.to_le_bytes())
        .collect()
}

#[test]
fn decodes_numeric_encoding_chains_in_reverse_order() {
    let fixed = decode_data(bytes(
        int32(&[12, -25]),
        vec![
            Encoding::FixedPoint {
                factor: 10.0,
                src_type: 33,
            },
            Encoding::ByteArray { data_type: 3 },
        ],
    ));
    assert!(matches!(fixed, Ok(DecodedColumn::Floats(values)) if values == [1.2, -2.5]));

    let interval = decode_data(bytes(
        int32(&[0, 2, 4]),
        vec![
            Encoding::IntervalQuantization {
                min: -1.0,
                max: 1.0,
                num_steps: 5,
                src_type: 33,
            },
            Encoding::ByteArray { data_type: 3 },
        ],
    ));
    assert!(matches!(interval, Ok(DecodedColumn::Floats(values)) if values == [-1.0, 0.0, 1.0]));

    let packed = decode_data(bytes(
        vec![127, 3, 128, 254],
        vec![
            Encoding::IntegerPacking {
                byte_count: 1,
                src_size: 2,
                is_unsigned: false,
            },
            Encoding::ByteArray { data_type: 1 },
        ],
    ));
    assert!(matches!(packed, Ok(DecodedColumn::Integers(values)) if values == [130, -130]));
}

#[test]
fn decodes_run_length_delta_and_compact_string_arrays() {
    let run_length = decode_data(bytes(
        int32(&[7, 2, 9, 3]),
        vec![
            Encoding::RunLength {
                src_type: 3,
                src_size: 5,
            },
            Encoding::ByteArray { data_type: 3 },
        ],
    ));
    assert!(matches!(run_length, Ok(DecodedColumn::Integers(values)) if values == [7, 7, 9, 9, 9]));

    let delta = decode_data(bytes(
        int32(&[2, -1, 4]),
        vec![
            Encoding::Delta {
                origin: 10,
                src_type: 3,
            },
            Encoding::ByteArray { data_type: 3 },
        ],
    ));
    assert!(matches!(delta, Ok(DecodedColumn::Integers(values)) if values == [12, 11, 15]));

    let strings = decode_data(bytes(
        int32(&[1, 0, -1, 1]),
        vec![Encoding::StringArray {
            data_encoding: vec![Encoding::ByteArray { data_type: 3 }],
            string_data: "ALAATP".to_owned(),
            offset_encoding: vec![Encoding::ByteArray { data_type: 3 }],
            offsets: Cow::Owned(int32(&[0, 3, 6])),
        }],
    ));
    let Ok(DecodedColumn::Strings(values)) = strings else {
        unreachable!("StringArray fixture must decode to strings");
    };
    assert_eq!(values.value(0), "ATP");
    assert_eq!(values.value(1), "ALA");
    assert_eq!(values.value(2), "");
    assert_eq!(values.value(3), "ATP");
}

#[test]
fn rejects_truncated_arrays_and_out_of_range_string_indices() {
    let truncated = decode_data(bytes(
        vec![1, 2, 3],
        vec![Encoding::ByteArray { data_type: 3 }],
    ));
    assert!(truncated.is_err());

    let invalid_index = decode_data(bytes(
        int32(&[1]),
        vec![Encoding::StringArray {
            data_encoding: vec![Encoding::ByteArray { data_type: 3 }],
            string_data: "ALA".to_owned(),
            offset_encoding: vec![Encoding::ByteArray { data_type: 3 }],
            offsets: Cow::Owned(int32(&[0, 3])),
        }],
    ));
    assert!(invalid_index.is_err());

    let empty_dictionary = decode_data(bytes(
        int32(&[-1, -1]),
        vec![Encoding::StringArray {
            data_encoding: vec![Encoding::ByteArray { data_type: 3 }],
            string_data: String::new(),
            offset_encoding: vec![Encoding::ByteArray { data_type: 3 }],
            offsets: Cow::Owned(int32(&[0])),
        }],
    ));
    assert!(matches!(
        empty_dictionary,
        Ok(DecodedColumn::Strings(values)) if values.value(0).is_empty()
    ));
}

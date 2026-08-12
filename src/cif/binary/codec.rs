//! BinaryCIF encoding-chain inversion.

use std::borrow::Cow;

use crate::cif::document::{ColumnValues, StringColumn};

use super::error::{BinaryCifError, BinaryCifErrorCode};
use super::model::{BinaryData, Encoding};

enum Decoded<'a> {
    Bytes(Cow<'a, [u8]>),
    Integers(Vec<i64>),
    Floats(Vec<f64>),
    Strings(StringColumn),
}

#[derive(Clone, Copy)]
enum IntegerEncoding {
    I8,
    I16,
    I32,
    U8,
    U16,
    U32,
}

impl IntegerEncoding {
    const fn width(self) -> usize {
        match self {
            Self::I8 | Self::U8 => 1,
            Self::I16 | Self::U16 => 2,
            Self::I32 | Self::U32 => 4,
        }
    }

    fn decode(self, chunk: &[u8]) -> i64 {
        match self {
            Self::I8 => i64::from(i8::from_le_bytes([chunk[0]])),
            Self::I16 => i64::from(i16::from_le_bytes([chunk[0], chunk[1]])),
            Self::I32 => i64::from(i32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]])),
            Self::U8 => i64::from(chunk[0]),
            Self::U16 => i64::from(u16::from_le_bytes([chunk[0], chunk[1]])),
            Self::U32 => i64::from(u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]])),
        }
    }
}

pub(super) fn decode_data(data: BinaryData<'_>) -> Result<ColumnValues, BinaryCifError> {
    match decode_encoded(data.data, data.encoding)? {
        Decoded::Integers(values) => Ok(ColumnValues::Integers(values)),
        Decoded::Floats(values) => Ok(ColumnValues::Floats(values)),
        Decoded::Strings(values) => Ok(ColumnValues::Strings(values)),
        Decoded::Bytes(_) => Err(encoding_error(
            "encoding chain leaves raw bytes without a ByteArray descriptor",
        )),
    }
}

fn decode_encoded<'a>(
    bytes: Cow<'a, [u8]>,
    encodings: Vec<Encoding<'a>>,
) -> Result<Decoded<'a>, BinaryCifError> {
    let mut decoded = Decoded::Bytes(bytes);
    for encoding in encodings.into_iter().rev() {
        decoded = decode_step(decoded, encoding)?;
    }
    Ok(decoded)
}

fn decode_step<'a>(
    decoded: Decoded<'a>,
    encoding: Encoding<'a>,
) -> Result<Decoded<'a>, BinaryCifError> {
    match encoding {
        Encoding::ByteArray { data_type } => decode_byte_array(decoded, data_type),
        Encoding::IntegerPacking {
            byte_count,
            src_size,
            is_unsigned,
        } => decode_integer_packing(decoded, byte_count, src_size, is_unsigned),
        Encoding::RunLength { src_type, src_size } => {
            decode_run_length(decoded, src_type, src_size)
        }
        Encoding::Delta { origin, src_type } => decode_delta(decoded, origin, src_type),
        Encoding::FixedPoint { factor, src_type } => decode_fixed_point(decoded, factor, src_type),
        Encoding::IntervalQuantization {
            min,
            max,
            num_steps,
            src_type,
        } => decode_interval(decoded, min, max, num_steps, src_type),
        Encoding::StringArray {
            data_encoding,
            string_data,
            offset_encoding,
            offsets,
        } => decode_string_array(
            decoded,
            data_encoding,
            &string_data,
            offset_encoding,
            offsets,
        ),
    }
}

fn decode_byte_array(decoded: Decoded<'_>, data_type: i32) -> Result<Decoded<'_>, BinaryCifError> {
    let Decoded::Bytes(bytes) = decoded else {
        return Err(encoding_error(
            "ByteArray must be the last applied encoding",
        ));
    };
    match data_type {
        1 => decode_integer_bytes(&bytes, IntegerEncoding::I8),
        2 => decode_integer_bytes(&bytes, IntegerEncoding::I16),
        3 => decode_integer_bytes(&bytes, IntegerEncoding::I32),
        4 => decode_integer_bytes(&bytes, IntegerEncoding::U8),
        5 => decode_integer_bytes(&bytes, IntegerEncoding::U16),
        6 => decode_integer_bytes(&bytes, IntegerEncoding::U32),
        32 => {
            if bytes.len() % 4 != 0 {
                return Err(encoding_error("Float32 ByteArray has a partial value"));
            }
            Ok(Decoded::Floats(
                bytes
                    .chunks_exact(4)
                    .map(|chunk| {
                        f64::from(f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
                    })
                    .collect(),
            ))
        }
        33 => {
            if bytes.len() % 8 != 0 {
                return Err(encoding_error("Float64 ByteArray has a partial value"));
            }
            Ok(Decoded::Floats(
                bytes
                    .chunks_exact(8)
                    .map(|chunk| {
                        f64::from_le_bytes([
                            chunk[0], chunk[1], chunk[2], chunk[3], chunk[4], chunk[5], chunk[6],
                            chunk[7],
                        ])
                    })
                    .collect(),
            ))
        }
        value => Err(encoding_error(format!(
            "unsupported ByteArray type code {value}"
        ))),
    }
}

fn decode_integer_bytes(
    bytes: &[u8],
    encoding: IntegerEncoding,
) -> Result<Decoded<'static>, BinaryCifError> {
    let width = encoding.width();
    if bytes.len() % width != 0 {
        return Err(encoding_error(format!(
            "ByteArray length {} is not divisible by element width {width}",
            bytes.len()
        )));
    }
    Ok(Decoded::Integers(
        bytes
            .chunks_exact(width)
            .map(|chunk| encoding.decode(chunk))
            .collect(),
    ))
}

fn decode_integer_packing(
    decoded: Decoded<'_>,
    byte_count: usize,
    src_size: usize,
    is_unsigned: bool,
) -> Result<Decoded<'_>, BinaryCifError> {
    let Decoded::Integers(values) = decoded else {
        return Err(encoding_error("IntegerPacking requires integer input"));
    };
    let (upper, lower) = match (byte_count, is_unsigned) {
        (1, true) => (i64::from(u8::MAX), 0),
        (2, true) => (i64::from(u16::MAX), 0),
        (1, false) => (i64::from(i8::MAX), i64::from(i8::MIN)),
        (2, false) => (i64::from(i16::MAX), i64::from(i16::MIN)),
        _ => return Err(encoding_error("IntegerPacking byteCount must be 1 or 2")),
    };
    let mut output = Vec::with_capacity(src_size);
    let mut accumulated = 0_i64;
    let mut continuation = false;
    for value in values {
        accumulated = accumulated
            .checked_add(value)
            .ok_or_else(|| encoding_error("IntegerPacking value overflows i64"))?;
        continuation = value == upper || (!is_unsigned && value == lower);
        if !continuation {
            output.push(accumulated);
            accumulated = 0;
        }
    }
    if continuation || output.len() != src_size {
        return Err(encoding_error(format!(
            "IntegerPacking decodes {} values, expected {src_size}",
            output.len()
        )));
    }
    Ok(Decoded::Integers(output))
}

fn decode_run_length(
    decoded: Decoded<'_>,
    src_type: i32,
    src_size: usize,
) -> Result<Decoded<'_>, BinaryCifError> {
    require_integer_type(src_type, "RunLength")?;
    let Decoded::Integers(values) = decoded else {
        return Err(encoding_error("RunLength requires integer input"));
    };
    if values.len() % 2 != 0 {
        return Err(encoding_error(
            "RunLength input must contain value/count pairs",
        ));
    }
    let mut output = Vec::with_capacity(src_size);
    for pair in values.chunks_exact(2) {
        let count = usize::try_from(pair[1])
            .map_err(|_| encoding_error("RunLength count must be non-negative"))?;
        if output.len().saturating_add(count) > src_size {
            return Err(encoding_error("RunLength output exceeds declared srcSize"));
        }
        output.extend(std::iter::repeat_n(pair[0], count));
    }
    if output.len() != src_size {
        return Err(encoding_error(format!(
            "RunLength decodes {} values, expected {src_size}",
            output.len()
        )));
    }
    Ok(Decoded::Integers(output))
}

fn decode_delta(
    decoded: Decoded<'_>,
    origin: i64,
    src_type: i32,
) -> Result<Decoded<'_>, BinaryCifError> {
    require_integer_type(src_type, "Delta")?;
    let Decoded::Integers(values) = decoded else {
        return Err(encoding_error("Delta requires integer input"));
    };
    let mut current = origin;
    let mut output = Vec::with_capacity(values.len());
    for difference in values {
        current = current
            .checked_add(difference)
            .ok_or_else(|| encoding_error("Delta value overflows i64"))?;
        output.push(current);
    }
    Ok(Decoded::Integers(output))
}

fn decode_fixed_point(
    decoded: Decoded<'_>,
    factor: f64,
    src_type: i32,
) -> Result<Decoded<'_>, BinaryCifError> {
    require_float_type(src_type, "FixedPoint")?;
    if !factor.is_finite() || factor == 0.0 {
        return Err(encoding_error(
            "FixedPoint factor must be finite and non-zero",
        ));
    }
    let Decoded::Integers(values) = decoded else {
        return Err(encoding_error("FixedPoint requires integer input"));
    };
    Ok(Decoded::Floats(
        values
            .into_iter()
            .map(|value| value as f64 / factor)
            .collect(),
    ))
}

fn decode_interval(
    decoded: Decoded<'_>,
    min: f64,
    max: f64,
    num_steps: i64,
    src_type: i32,
) -> Result<Decoded<'_>, BinaryCifError> {
    require_float_type(src_type, "IntervalQuantization")?;
    if !min.is_finite() || !max.is_finite() || max < min || num_steps < 2 {
        return Err(encoding_error(
            "IntervalQuantization requires a finite ordered interval and at least two steps",
        ));
    }
    let Decoded::Integers(values) = decoded else {
        return Err(encoding_error(
            "IntervalQuantization requires integer input",
        ));
    };
    let delta = (max - min) / (num_steps - 1) as f64;
    Ok(Decoded::Floats(
        values
            .into_iter()
            .map(|value| min + value as f64 * delta)
            .collect(),
    ))
}

fn decode_string_array<'a>(
    decoded: Decoded<'a>,
    data_encoding: Vec<Encoding<'a>>,
    string_data: &str,
    offset_encoding: Vec<Encoding<'a>>,
    offsets: Cow<'a, [u8]>,
) -> Result<Decoded<'a>, BinaryCifError> {
    let Decoded::Bytes(index_bytes) = decoded else {
        return Err(encoding_error("StringArray must receive raw bytes"));
    };
    let Decoded::Integers(indices) = decode_encoded(index_bytes, data_encoding)? else {
        return Err(encoding_error(
            "StringArray indices must decode to integers",
        ));
    };
    let Decoded::Integers(offsets) = decode_encoded(offsets, offset_encoding)? else {
        return Err(encoding_error(
            "StringArray offsets must decode to integers",
        ));
    };
    if offsets.is_empty() {
        return Err(string_error("StringArray requires at least one offset"));
    }
    let offsets = offsets
        .into_iter()
        .map(|offset| {
            u32::try_from(offset)
                .map_err(|_| string_error("StringArray offsets must be non-negative"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    if offsets.first() != Some(&0)
        || offsets.last().copied().map(|offset| offset as usize) != Some(string_data.len())
        || offsets.windows(2).any(|pair| pair[0] > pair[1])
        || offsets
            .iter()
            .any(|offset| !string_data.is_char_boundary(*offset as usize))
    {
        return Err(string_error(
            "StringArray offsets must be ordered UTF-8 boundaries spanning stringData",
        ));
    }
    let mut decoded_indices = Vec::with_capacity(indices.len());
    for index in indices {
        if index == -1 {
            decoded_indices.push(0);
            continue;
        }
        let index = u32::try_from(index)
            .map_err(|_| string_error("StringArray index must be -1 or non-negative"))?;
        let encoded = index
            .checked_add(1)
            .ok_or_else(|| string_error("StringArray dictionary exceeds UInt32"))?;
        let index = index as usize;
        if offsets.get(index..=index + 1).is_none() {
            return Err(string_error("StringArray index exceeds its dictionary"));
        }
        decoded_indices.push(encoded);
    }
    Ok(Decoded::Strings(StringColumn::new(
        string_data.to_owned(),
        offsets,
        decoded_indices,
    )))
}

fn require_integer_type(value: i32, encoding: &str) -> Result<(), BinaryCifError> {
    if (1..=6).contains(&value) {
        Ok(())
    } else {
        Err(encoding_error(format!(
            "{encoding} has non-integer srcType {value}"
        )))
    }
}

fn require_float_type(value: i32, encoding: &str) -> Result<(), BinaryCifError> {
    if matches!(value, 32 | 33) {
        Ok(())
    } else {
        Err(encoding_error(format!(
            "{encoding} has non-floating srcType {value}"
        )))
    }
}

fn encoding_error(message: impl Into<String>) -> BinaryCifError {
    BinaryCifError::new(BinaryCifErrorCode::Encoding, message)
}

fn string_error(message: impl Into<String>) -> BinaryCifError {
    BinaryCifError::new(BinaryCifErrorCode::StringData, message)
}

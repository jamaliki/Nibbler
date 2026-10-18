//! Owned native-endian buffers that Python wraps with `numpy.frombuffer`.

/// One NumPy array as its dtype string, shape and bytes.
pub(crate) struct Array {
    pub dtype: &'static str,
    pub shape: Vec<usize>,
    pub bytes: Vec<u8>,
}

pub(super) fn int64(values: &[i64], shape: Vec<usize>) -> Array {
    Array { dtype: "=i8", shape, bytes: values.iter().flat_map(|value| value.to_ne_bytes()).collect() }
}

pub(super) fn float32(values: &[f32], shape: Vec<usize>) -> Array {
    Array { dtype: "=f4", shape, bytes: values.iter().flat_map(|value| value.to_ne_bytes()).collect() }
}

pub(super) fn int16(values: &[i16], shape: Vec<usize>) -> Array {
    Array { dtype: "=i2", shape, bytes: values.iter().flat_map(|value| value.to_ne_bytes()).collect() }
}

pub(super) fn int8(values: &[i8], shape: Vec<usize>) -> Array {
    Array { dtype: "i1", shape, bytes: values.iter().flat_map(|value| value.to_ne_bytes()).collect() }
}

pub(super) fn boolean(values: &[bool], shape: Vec<usize>) -> Array {
    Array { dtype: "?", shape, bytes: values.iter().map(|&value| u8::from(value)).collect() }
}

pub(super) fn uint8(values: Vec<u8>, shape: Vec<usize>) -> Array {
    Array { dtype: "u1", shape, bytes: values }
}

pub(super) fn text(value: String) -> Array {
    let count = value.len();
    uint8(value.into_bytes(), vec![count])
}

pub(super) fn points(values: &[[f32; 3]], shape: Vec<usize>) -> Array {
    Array { dtype: "=f4", shape, bytes: values.iter().flatten().flat_map(|value| value.to_ne_bytes()).collect() }
}

/// Columns of index pairs, as an `(2, n)` array.
pub(super) fn columns(values: &[[i64; 2]]) -> Array {
    let data: Vec<i64> = values.iter().map(|pair| pair[0]).chain(values.iter().map(|pair| pair[1])).collect();
    int64(&data, vec![2, values.len()])
}

/// Rows of index pairs, as an `(n, 2)` array.
pub(super) fn pairs(values: &[[i64; 2]]) -> Array {
    let data: Vec<i64> = values.iter().flatten().copied().collect();
    int64(&data, vec![values.len(), 2])
}

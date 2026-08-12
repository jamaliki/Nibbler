//! Serde representation of the BinaryCIF MessagePack container.

use std::borrow::Cow;

use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub(super) struct BinaryFile<'a> {
    pub(super) version: String,
    pub(super) encoder: String,
    #[serde(borrow)]
    pub(super) data_blocks: Vec<BinaryBlock<'a>>,
}

#[derive(Debug, Deserialize, Serialize)]
pub(super) struct BinaryBlock<'a> {
    pub(super) header: String,
    #[serde(borrow)]
    pub(super) categories: Vec<BinaryCategory<'a>>,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub(super) struct BinaryCategory<'a> {
    pub(super) name: String,
    pub(super) row_count: usize,
    #[serde(borrow)]
    pub(super) columns: Vec<BinaryColumn<'a>>,
}

#[derive(Debug, Deserialize, Serialize)]
pub(super) struct BinaryColumn<'a> {
    pub(super) name: String,
    #[serde(borrow)]
    pub(super) data: BinaryData<'a>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(borrow)]
    pub(super) mask: Option<BinaryData<'a>>,
}

#[derive(Debug, Deserialize, Serialize)]
pub(super) struct BinaryData<'a> {
    #[serde(borrow)]
    pub(super) encoding: Vec<Encoding<'a>>,
    #[serde(borrow, with = "serde_bytes")]
    pub(super) data: Cow<'a, [u8]>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "kind")]
pub(super) enum Encoding<'a> {
    #[serde(rename = "ByteArray")]
    ByteArray {
        #[serde(rename = "type")]
        data_type: i32,
    },
    #[serde(rename = "FixedPoint")]
    FixedPoint {
        factor: f64,
        #[serde(rename = "srcType")]
        src_type: i32,
    },
    #[serde(rename = "IntervalQuantization")]
    IntervalQuantization {
        min: f64,
        max: f64,
        #[serde(rename = "numSteps")]
        num_steps: i64,
        #[serde(rename = "srcType")]
        src_type: i32,
    },
    #[serde(rename = "RunLength")]
    RunLength {
        #[serde(rename = "srcType")]
        src_type: i32,
        #[serde(rename = "srcSize")]
        src_size: usize,
    },
    #[serde(rename = "Delta")]
    Delta {
        origin: i64,
        #[serde(rename = "srcType")]
        src_type: i32,
    },
    #[serde(rename = "IntegerPacking")]
    IntegerPacking {
        #[serde(rename = "byteCount")]
        byte_count: usize,
        #[serde(rename = "srcSize")]
        src_size: usize,
        #[serde(rename = "isUnsigned")]
        is_unsigned: bool,
    },
    #[serde(rename = "StringArray")]
    StringArray {
        #[serde(rename = "dataEncoding")]
        #[serde(borrow)]
        data_encoding: Vec<Encoding<'a>>,
        #[serde(rename = "stringData")]
        string_data: String,
        #[serde(rename = "offsetEncoding")]
        #[serde(borrow)]
        offset_encoding: Vec<Encoding<'a>>,
        #[serde(borrow, with = "serde_bytes")]
        offsets: Cow<'a, [u8]>,
    },
}

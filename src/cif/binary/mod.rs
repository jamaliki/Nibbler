//! BinaryCIF 0.3 decoding, encoding, and direct projection.

mod codec;
#[cfg(test)]
mod codec_tests;
mod document;
mod error;
mod model;
mod projection;
mod writer;

pub use document::decode_binary;
pub use error::{BinaryCifError, BinaryCifErrorCode};
pub use projection::project_binary;
#[cfg(feature = "python")]
pub(crate) use projection::project_binary_named;
pub use writer::encode_binary;

#[cfg(feature = "python")]
pub(crate) fn is_binary_container(bytes: &[u8]) -> bool {
    matches!(bytes.first(), Some(0x80..=0x8f | 0xde | 0xdf))
}

//! Stable BinaryCIF diagnostics.

use std::error::Error;
use std::fmt::{self, Display, Formatter};

/// Stable categories for BinaryCIF failures.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BinaryCifErrorCode {
    /// The MessagePack container or BinaryCIF object shape is malformed.
    Container,
    /// The BinaryCIF version is unsupported.
    Version,
    /// An encoding descriptor or encoding chain is unsupported or malformed.
    Encoding,
    /// Decoded column or mask lengths do not match the declared row count.
    Shape,
    /// A decoded string table or UTF-8 value is invalid.
    StringData,
    /// The CIF document cannot be represented by BinaryCIF without data loss.
    Unrepresentable,
}

impl BinaryCifErrorCode {
    /// Return the stable machine-readable diagnostic code.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Container => "CIF_BINARY_CONTAINER",
            Self::Version => "CIF_BINARY_VERSION",
            Self::Encoding => "CIF_BINARY_ENCODING",
            Self::Shape => "CIF_BINARY_SHAPE",
            Self::StringData => "CIF_BINARY_STRING_DATA",
            Self::Unrepresentable => "CIF_BINARY_UNREPRESENTABLE",
        }
    }
}

/// A structured BinaryCIF decode, projection, or encode failure.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BinaryCifError {
    code: BinaryCifErrorCode,
    message: String,
}

impl BinaryCifError {
    pub(super) fn new(code: BinaryCifErrorCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }

    /// Return the stable error category.
    #[must_use]
    pub const fn code(&self) -> BinaryCifErrorCode {
        self.code
    }

    /// Return the human-readable diagnostic detail.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }
}

impl Display for BinaryCifError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {}", self.code.as_str(), self.message)
    }
}

impl Error for BinaryCifError {}

//! Errors raised while constructing typed semantic models.

use std::error::Error;
use std::fmt::{self, Display, Formatter};

/// A structural or typed failure while constructing a semantic model.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SemanticError {
    code: &'static str,
    message: String,
    context: Vec<String>,
}

impl SemanticError {
    pub(crate) fn new(
        code: &'static str,
        message: impl Into<String>,
        context: Vec<String>,
    ) -> Self {
        Self {
            code,
            message: message.into(),
            context,
        }
    }

    /// Return the stable machine-readable error code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.code
    }

    /// Return the human-readable failure detail.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }

    /// Return category, item, and row context when available.
    #[must_use]
    pub fn context(&self) -> &[String] {
        &self.context
    }
}

impl Display for SemanticError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {}", self.code, self.message)
    }
}

impl Error for SemanticError {}

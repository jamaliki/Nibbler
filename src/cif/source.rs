use std::ops::Range;
use std::sync::{Arc, OnceLock};

use super::error::{ParseError, ParseErrorCode, SourceSpan};

#[derive(Debug)]
struct SourceInner {
    name: Box<str>,
    text: String,
    line_starts: OnceLock<Vec<usize>>,
}

/// Immutable, reference-counted UTF-8 source bytes.
///
/// Clones share the underlying allocation. Line starts are indexed only if a
/// diagnostic asks for a display location.
#[derive(Clone, Debug)]
pub struct SourceBuffer(Arc<SourceInner>);

impl SourceBuffer {
    /// Copy and validate source bytes.
    ///
    /// # Errors
    ///
    /// Returns [`ParseErrorCode::InvalidUtf8`] at the first invalid byte.
    pub fn from_bytes(name: impl Into<String>, bytes: &[u8]) -> Result<Self, ParseError> {
        let source_name = name.into();
        let text = match std::str::from_utf8(bytes) {
            Ok(text) => text,
            Err(error) => {
                let byte_offset = error.valid_up_to();
                return Err(ParseError::new(
                    ParseErrorCode::InvalidUtf8,
                    "CIF input must be valid UTF-8",
                    source_name,
                    raw_span(bytes, byte_offset),
                ));
            }
        };
        Ok(Self::from_valid_text(source_name, text.to_owned()))
    }

    pub(crate) fn from_owned_bytes(name: String, bytes: Vec<u8>) -> Result<Self, ParseError> {
        match String::from_utf8(bytes) {
            Ok(text) => Ok(Self::from_valid_text(name, text)),
            Err(error) => {
                let byte_offset = error.utf8_error().valid_up_to();
                let bytes = error.into_bytes();
                Err(ParseError::new(
                    ParseErrorCode::InvalidUtf8,
                    "CIF input must be valid UTF-8",
                    name,
                    raw_span(&bytes, byte_offset),
                ))
            }
        }
    }

    /// Copy an already-valid UTF-8 string into a shared source buffer.
    #[must_use]
    pub fn from_text(name: impl Into<String>, text: &str) -> Self {
        Self::from_valid_text(name.into(), text.to_owned())
    }

    fn from_valid_text(name: String, text: String) -> Self {
        Self(Arc::new(SourceInner {
            name: name.into_boxed_str(),
            text,
            line_starts: OnceLock::new(),
        }))
    }

    /// Return the caller-supplied source name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.0.name
    }

    /// Return the entire source as UTF-8 text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0.text
    }

    /// Return the source length in bytes.
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.text.len()
    }

    /// Return whether the source is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.text.is_empty()
    }

    pub(crate) fn slice(&self, range: Range<usize>) -> &str {
        // Lexer-produced ranges begin and end at ASCII syntax bytes, hence at UTF-8
        // character boundaries. No file content can violate this invariant.
        &self.0.text[range]
    }

    pub(crate) fn span(&self, byte_start: usize, byte_end: usize) -> SourceSpan {
        let bounded_start = byte_start.min(self.len());
        let line_starts = self
            .0
            .line_starts
            .get_or_init(|| collect_line_starts(self.0.text.as_bytes()));
        let line_index = line_starts.partition_point(|start| *start <= bounded_start) - 1;
        let line_start = line_starts[line_index];
        let column = self.0.text[line_start..bounded_start].chars().count() + 1;
        SourceSpan::new(
            bounded_start,
            byte_end.min(self.len()),
            line_index + 1,
            column,
        )
    }

    pub(crate) fn error(
        &self,
        code: ParseErrorCode,
        message: impl Into<String>,
        byte_start: usize,
        byte_end: usize,
    ) -> ParseError {
        ParseError::new(code, message, self.name(), self.span(byte_start, byte_end))
    }
}

fn raw_span(bytes: &[u8], byte_offset: usize) -> SourceSpan {
    let prefix = &bytes[..byte_offset.min(bytes.len())];
    let line_starts = collect_line_starts(prefix);
    let line = line_starts.len();
    let line_start = line_starts[line - 1];
    let column = std::str::from_utf8(&prefix[line_start..])
        .map_or(prefix.len() - line_start + 1, |text| {
            text.chars().count() + 1
        });
    SourceSpan::new(
        byte_offset,
        byte_offset.saturating_add(1).min(bytes.len()),
        line,
        column,
    )
}

fn collect_line_starts(bytes: &[u8]) -> Vec<usize> {
    let mut starts = vec![0];
    let mut byte_offset = 0;
    while byte_offset < bytes.len() {
        match bytes[byte_offset] {
            b'\r' if bytes.get(byte_offset + 1) == Some(&b'\n') => {
                byte_offset += 2;
                starts.push(byte_offset);
            }
            b'\r' | b'\n' => {
                byte_offset += 1;
                starts.push(byte_offset);
            }
            _ => byte_offset += 1,
        }
    }
    starts
}

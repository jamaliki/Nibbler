use std::error::Error;
use std::fmt::{self, Display, Formatter};

/// A half-open byte span with a one-based display location.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SourceSpan {
    byte_start: usize,
    byte_end: usize,
    line: usize,
    column: usize,
}

impl SourceSpan {
    pub(crate) const fn new(
        byte_start: usize,
        byte_end: usize,
        line: usize,
        column: usize,
    ) -> Self {
        Self {
            byte_start,
            byte_end,
            line,
            column,
        }
    }

    /// Return the first byte offset of the span.
    #[must_use]
    pub const fn byte_start(self) -> usize {
        self.byte_start
    }

    /// Return the exclusive final byte offset of the span.
    #[must_use]
    pub const fn byte_end(self) -> usize {
        self.byte_end
    }

    /// Return the one-based line of the first byte.
    #[must_use]
    pub const fn line(self) -> usize {
        self.line
    }

    /// Return the one-based character column of the first byte.
    #[must_use]
    pub const fn column(self) -> usize {
        self.column
    }
}

/// Stable categories for strict CIF parse failures.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ParseErrorCode {
    /// The source exceeds a configured hard limit.
    ResourceLimit,
    /// The source is not valid UTF-8.
    InvalidUtf8,
    /// A forbidden control character occurs in the source.
    InvalidCharacter,
    /// A quoted value has no closing quote.
    UnterminatedQuotedValue,
    /// A semicolon-delimited value has no closing delimiter.
    UnterminatedTextField,
    /// A text-field closing semicolon is followed by non-whitespace on its line.
    InvalidTextFieldTerminator,
    /// A token appears before the first block header.
    ExpectedBlock,
    /// A `data_` header has no block code.
    EmptyBlockCode,
    /// A tag consists only of its leading underscore.
    InvalidTag,
    /// Two block codes differ only in case or are identical.
    DuplicateBlock,
    /// A tag is repeated in one block or save frame.
    DuplicateTag,
    /// A scalar tag is not followed by a value.
    MissingItemValue,
    /// A loop contains no tags.
    MissingLoopTag,
    /// A loop contains no values.
    MissingLoopValue,
    /// A loop value count is not divisible by its tag count.
    LoopValueCount,
    /// A save frame was opened inside another save frame.
    NestedSaveFrame,
    /// A closing `save_` occurs outside a save frame.
    UnexpectedSaveEnd,
    /// A save frame reaches a block boundary or end of input before `save_`.
    UnterminatedSaveFrame,
    /// A `stop_` occurs outside a loop.
    UnexpectedStop,
    /// A bare value occurs where a tag or control word is required.
    UnexpectedValue,
    /// A control word is not legal in the current parser state.
    UnexpectedControl,
}

impl ParseErrorCode {
    /// Return the stable machine-readable diagnostic code.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ResourceLimit => "CIF_RESOURCE_LIMIT",
            Self::InvalidUtf8 => "CIF_INVALID_UTF8",
            Self::InvalidCharacter => "CIF_INVALID_CHARACTER",
            Self::UnterminatedQuotedValue => "CIF_UNTERMINATED_QUOTED_VALUE",
            Self::UnterminatedTextField => "CIF_UNTERMINATED_TEXT_FIELD",
            Self::InvalidTextFieldTerminator => "CIF_TEXT_FIELD_TERMINATOR_INVALID",
            Self::ExpectedBlock => "CIF_EXPECTED_BLOCK",
            Self::EmptyBlockCode => "CIF_BLOCK_CODE_EMPTY",
            Self::InvalidTag => "CIF_TAG_INVALID",
            Self::DuplicateBlock => "CIF_DUPLICATE_BLOCK",
            Self::DuplicateTag => "CIF_DUPLICATE_TAG",
            Self::MissingItemValue => "CIF_ITEM_VALUE_MISSING",
            Self::MissingLoopTag => "CIF_LOOP_TAG_MISSING",
            Self::MissingLoopValue => "CIF_LOOP_VALUE_MISSING",
            Self::LoopValueCount => "CIF_LOOP_VALUE_COUNT",
            Self::NestedSaveFrame => "CIF_SAVE_FRAME_NESTED",
            Self::UnexpectedSaveEnd => "CIF_SAVE_FRAME_END_UNEXPECTED",
            Self::UnterminatedSaveFrame => "CIF_SAVE_FRAME_UNTERMINATED",
            Self::UnexpectedStop => "CIF_STOP_UNEXPECTED",
            Self::UnexpectedValue => "CIF_VALUE_UNEXPECTED",
            Self::UnexpectedControl => "CIF_CONTROL_UNEXPECTED",
        }
    }
}

/// A located, structured failure produced by the CIF lexer or parser.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseError {
    code: ParseErrorCode,
    message: String,
    source_name: String,
    span: SourceSpan,
}

impl ParseError {
    pub(crate) fn new(
        code: ParseErrorCode,
        message: impl Into<String>,
        source_name: impl Into<String>,
        span: SourceSpan,
    ) -> Self {
        Self {
            code,
            message: message.into(),
            source_name: source_name.into(),
            span,
        }
    }

    /// Return the stable error category.
    #[must_use]
    pub const fn code(&self) -> ParseErrorCode {
        self.code
    }

    /// Return the stable machine-readable diagnostic code.
    #[must_use]
    pub const fn code_str(&self) -> &'static str {
        self.code.as_str()
    }

    /// Return the human-readable error detail.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }

    /// Return the source name supplied by the caller.
    #[must_use]
    pub fn source_name(&self) -> &str {
        &self.source_name
    }

    /// Return the byte and display location of the failure.
    #[must_use]
    pub const fn span(&self) -> SourceSpan {
        self.span
    }
}

impl Display for ParseError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{}:{}:{}: {}: {}",
            self.source_name,
            self.span.line(),
            self.span.column(),
            self.code.as_str(),
            self.message
        )
    }
}

impl Error for ParseError {}

/// Stable categories for canonical serialization failures.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum WriteErrorCode {
    /// The logical document violates an internal structural invariant.
    InvalidDocument,
    /// A text value cannot be represented losslessly as CIF 1.1.
    UnrepresentableText,
    /// A floating-point value is NaN or infinite.
    NonFiniteFloat,
}

impl WriteErrorCode {
    /// Return the stable machine-readable diagnostic code.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::InvalidDocument => "CIF_WRITE_INVALID_DOCUMENT",
            Self::UnrepresentableText => "CIF_WRITE_UNREPRESENTABLE_TEXT",
            Self::NonFiniteFloat => "CIF_WRITE_NON_FINITE_FLOAT",
        }
    }
}

/// A structured canonical serialization failure.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct WriteError {
    code: WriteErrorCode,
    message: String,
}

impl WriteError {
    pub(crate) fn new(code: WriteErrorCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }

    /// Return the stable error category.
    #[must_use]
    pub const fn code(&self) -> WriteErrorCode {
        self.code
    }

    /// Return the stable machine-readable diagnostic code.
    #[must_use]
    pub const fn code_str(&self) -> &'static str {
        self.code.as_str()
    }

    /// Return the human-readable error detail.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }
}

impl Display for WriteError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {}", self.code.as_str(), self.message)
    }
}

impl Error for WriteError {}

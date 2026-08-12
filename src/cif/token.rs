use std::ops::Range;

use super::document::QuoteStyle;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TokenKind {
    Data,
    Global,
    Loop,
    Stop,
    SaveStart,
    SaveEnd,
    Tag,
    Value,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Token {
    pub(crate) kind: TokenKind,
    pub(crate) span: Range<usize>,
    pub(crate) content: Range<usize>,
    pub(crate) quote_style: QuoteStyle,
    pub(crate) line_start: bool,
}

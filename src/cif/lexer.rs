use memchr::memchr3;

use super::document::QuoteStyle;
use super::error::{ParseError, ParseErrorCode};
use super::source::SourceBuffer;
use super::token::{Token, TokenKind};

pub(crate) struct Lexer {
    source: SourceBuffer,
    byte_offset: usize,
    byte_end: usize,
    max_token_bytes: usize,
}

impl Lexer {
    pub(crate) fn new(source: SourceBuffer, max_token_bytes: usize) -> Self {
        let byte_end = source.len();
        Self {
            source,
            byte_offset: 0,
            byte_end,
            max_token_bytes,
        }
    }

    pub(crate) fn bounded(
        source: SourceBuffer,
        byte_offset: usize,
        byte_end: usize,
        max_token_bytes: usize,
    ) -> Self {
        debug_assert!(byte_offset <= byte_end);
        debug_assert!(byte_end <= source.len());
        Self {
            source,
            byte_offset,
            byte_end,
            max_token_bytes,
        }
    }

    pub(crate) fn set_offset(&mut self, byte_offset: usize) {
        debug_assert!(byte_offset <= self.byte_end);
        self.byte_offset = byte_offset;
    }

    pub(crate) fn next_token(&mut self) -> Result<Option<Token>, ParseError> {
        self.skip_trivia()?;
        let bytes = self.source.as_str().as_bytes();
        if self.byte_offset == self.byte_end {
            return Ok(None);
        }

        let line_start = is_line_start(bytes, self.byte_offset);
        match bytes[self.byte_offset] {
            b'\'' | b'"' => self.quoted_value(line_start).map(Some),
            b';' if line_start => self.text_field().map(Some),
            _ => self.unquoted_token(line_start).map(Some),
        }
    }

    fn skip_trivia(&mut self) -> Result<(), ParseError> {
        let bytes = self.source.as_str().as_bytes();
        while self.byte_offset < self.byte_end {
            let byte = bytes[self.byte_offset];
            if is_whitespace(byte) {
                self.byte_offset += 1;
                continue;
            }
            if byte == b'#' {
                self.byte_offset += 1;
                while self.byte_offset < self.byte_end
                    && !matches!(bytes[self.byte_offset], b'\n' | b'\r')
                {
                    self.reject_forbidden(bytes[self.byte_offset], self.byte_offset)?;
                    self.byte_offset += 1;
                }
                continue;
            }
            break;
        }
        Ok(())
    }

    fn quoted_value(&mut self, line_start: bool) -> Result<Token, ParseError> {
        let bytes = self.source.as_str().as_bytes();
        let start = self.byte_offset;
        let quote = bytes[start];
        let quote_style = if quote == b'\'' {
            QuoteStyle::Single
        } else {
            QuoteStyle::Double
        };
        let content_start = start + 1;
        self.byte_offset = content_start;

        while self.byte_offset < self.byte_end {
            let Some(relative) =
                memchr3(quote, b'\n', b'\r', &bytes[self.byte_offset..self.byte_end])
            else {
                break;
            };
            self.byte_offset += relative;
            let byte = bytes[self.byte_offset];
            if byte == quote
                && (self.byte_offset + 1 == bytes.len()
                    || is_whitespace(bytes[self.byte_offset + 1]))
            {
                let content_end = self.byte_offset;
                self.byte_offset += 1;
                self.validate_token(start, self.byte_offset)?;
                return Ok(Token {
                    kind: TokenKind::Value,
                    span: start..self.byte_offset,
                    content: content_start..content_end,
                    quote_style,
                    line_start,
                });
            }
            if matches!(byte, b'\n' | b'\r') {
                return Err(self.source.error(
                    ParseErrorCode::UnterminatedQuotedValue,
                    "quoted CIF values cannot cross a line boundary",
                    start,
                    self.byte_offset,
                ));
            }
            self.byte_offset += 1;
        }

        self.byte_offset = self.byte_end;
        self.validate_token(start, self.byte_offset)?;

        Err(self.source.error(
            ParseErrorCode::UnterminatedQuotedValue,
            "quoted CIF value reaches end of input before its closing quote",
            start,
            self.byte_end,
        ))
    }

    fn text_field(&mut self) -> Result<Token, ParseError> {
        let bytes = self.source.as_str().as_bytes();
        let start = self.byte_offset;
        let content_start = start + 1;
        self.byte_offset = content_start;

        while self.byte_offset < self.byte_end {
            let byte = bytes[self.byte_offset];
            if byte == b';' && is_line_start(bytes, self.byte_offset) {
                if self.byte_offset + 1 < bytes.len() && !is_whitespace(bytes[self.byte_offset + 1])
                {
                    return Err(self.source.error(
                        ParseErrorCode::InvalidTextFieldTerminator,
                        "a text-field closing semicolon must be followed by whitespace or end of input",
                        self.byte_offset,
                        self.byte_offset + 1,
                    ));
                }
                let content_end = preceding_line_ending_start(bytes, self.byte_offset);
                self.byte_offset += 1;
                self.check_token_size(start, self.byte_offset)?;
                return Ok(Token {
                    kind: TokenKind::Value,
                    span: start..self.byte_offset,
                    content: content_start..content_end,
                    quote_style: QuoteStyle::TextField,
                    line_start: true,
                });
            }
            self.reject_forbidden(byte, self.byte_offset)?;
            self.byte_offset += 1;
            self.check_token_size(start, self.byte_offset)?;
        }

        Err(self.source.error(
            ParseErrorCode::UnterminatedTextField,
            "semicolon-delimited CIF value has no closing delimiter",
            start,
            self.byte_end,
        ))
    }

    fn unquoted_token(&mut self, line_start: bool) -> Result<Token, ParseError> {
        let bytes = self.source.as_str().as_bytes();
        let start = self.byte_offset;
        while self.byte_offset < self.byte_end {
            let byte = bytes[self.byte_offset];
            if byte > b' ' && byte != 0x7f {
                if byte >= 0xc0 {
                    self.reject_forbidden(byte, self.byte_offset)?;
                }
                self.byte_offset += 1;
                continue;
            }
            if is_whitespace(byte) {
                break;
            }
            self.reject_forbidden(byte, self.byte_offset)?;
            self.byte_offset += 1;
        }
        self.check_token_size(start, self.byte_offset)?;

        let content = start..self.byte_offset;
        let text = self.source.slice(content.clone());
        let kind = classify_unquoted(text);
        Ok(Token {
            kind,
            span: content.clone(),
            content,
            quote_style: QuoteStyle::Unquoted,
            line_start,
        })
    }

    fn check_token_size(&self, start: usize, end: usize) -> Result<(), ParseError> {
        if end - start > self.max_token_bytes {
            return Err(self.source.error(
                ParseErrorCode::ResourceLimit,
                format!(
                    "token exceeds the configured limit of {} bytes",
                    self.max_token_bytes
                ),
                start,
                end,
            ));
        }
        Ok(())
    }

    fn validate_token(&self, start: usize, end: usize) -> Result<(), ParseError> {
        self.check_token_size(start, end)?;
        let text = self.source.slice(start..end);
        if !needs_scalar_validation(text.as_bytes()) {
            return Ok(());
        }
        if let Some(relative) = text
            .as_bytes()
            .iter()
            .position(|byte| (*byte < b' ' && !is_whitespace(*byte)) || *byte == 0x7f)
        {
            return self.reject_forbidden(text.as_bytes()[relative], start + relative);
        }
        if !text.is_ascii()
            && let Some((relative, character)) = text
                .char_indices()
                .find(|(_, character)| character.is_control())
        {
            return Err(self.source.error(
                ParseErrorCode::InvalidCharacter,
                format!("forbidden control character U+{:04X}", u32::from(character)),
                start + relative,
                start + relative + character.len_utf8(),
            ));
        }
        Ok(())
    }

    fn reject_forbidden(&self, byte: u8, byte_offset: usize) -> Result<(), ParseError> {
        if (byte < b' ' && !is_whitespace(byte)) || byte == 0x7f {
            return Err(self.source.error(
                ParseErrorCode::InvalidCharacter,
                format!("forbidden control byte 0x{byte:02x}"),
                byte_offset,
                byte_offset + 1,
            ));
        }
        if byte >= 0xc0 {
            // UTF-8 validity is established by SourceBuffer, and a non-continuation
            // byte is therefore a character boundary.
            if let Some(character) = self.source.as_str()[byte_offset..].chars().next()
                && character.is_control()
            {
                return Err(self.source.error(
                    ParseErrorCode::InvalidCharacter,
                    format!("forbidden control character U+{:04X}", u32::from(character)),
                    byte_offset,
                    byte_offset + character.len_utf8(),
                ));
            }
        }
        Ok(())
    }
}

const BYTE_HIGH_BITS: u64 = 0x8080_8080_8080_8080;
const BYTE_LOW_BITS: u64 = 0x0101_0101_0101_0101;
const BYTE_UPPER_BITS: u64 = 0xe0e0_e0e0_e0e0_e0e0;
const BYTE_DEL: u64 = 0x7f7f_7f7f_7f7f_7f7f;

fn needs_scalar_validation(bytes: &[u8]) -> bool {
    let mut chunks = bytes.chunks_exact(8);
    if chunks.any(|chunk| suspicious_word(word(chunk))) {
        return true;
    }
    let remainder = chunks.remainder();
    remainder.iter().any(|byte| *byte < b' ' || *byte >= 0x7f)
}

fn word(chunk: &[u8]) -> u64 {
    debug_assert_eq!(chunk.len(), 8);
    u64::from_ne_bytes([
        chunk[0], chunk[1], chunk[2], chunk[3], chunk[4], chunk[5], chunk[6], chunk[7],
    ])
}

fn suspicious_word(word: u64) -> bool {
    word & BYTE_HIGH_BITS != 0
        || has_zero_byte(word & BYTE_UPPER_BITS)
        || has_zero_byte(word ^ BYTE_DEL)
}

fn has_zero_byte(word: u64) -> bool {
    word.wrapping_sub(BYTE_LOW_BITS) & !word & BYTE_HIGH_BITS != 0
}

fn classify_unquoted(text: &str) -> TokenKind {
    if text.starts_with('_') {
        return TokenKind::Tag;
    }
    match text.as_bytes().first().map(u8::to_ascii_lowercase) {
        Some(b'l') if eq_ascii_case(text, "loop_") => TokenKind::Loop,
        Some(b'g') if eq_ascii_case(text, "global_") => TokenKind::Global,
        Some(b'd') if starts_ascii_case(text, "data_") => TokenKind::Data,
        Some(b's') if eq_ascii_case(text, "stop_") => TokenKind::Stop,
        Some(b's') if eq_ascii_case(text, "save_") => TokenKind::SaveEnd,
        Some(b's') if starts_ascii_case(text, "save_") => TokenKind::SaveStart,
        _ => TokenKind::Value,
    }
}

pub(crate) const fn is_whitespace(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r')
}

fn is_line_start(bytes: &[u8], byte_offset: usize) -> bool {
    byte_offset == 0 || matches!(bytes[byte_offset - 1], b'\n' | b'\r')
}

fn preceding_line_ending_start(bytes: &[u8], byte_offset: usize) -> usize {
    let mut content_end = byte_offset;
    if content_end > 0 && bytes[content_end - 1] == b'\n' {
        content_end -= 1;
        if content_end > 0 && bytes[content_end - 1] == b'\r' {
            content_end -= 1;
        }
    } else if content_end > 0 && bytes[content_end - 1] == b'\r' {
        content_end -= 1;
    }
    content_end
}

fn eq_ascii_case(left: &str, right: &str) -> bool {
    left.eq_ignore_ascii_case(right)
}

fn starts_ascii_case(text: &str, prefix: &str) -> bool {
    text.get(..prefix.len())
        .is_some_and(|candidate| candidate.eq_ignore_ascii_case(prefix))
}

#[cfg(feature = "fuzzing")]
/// Exercise only tokenization for `cargo fuzz`.
pub fn fuzz_lexer(bytes: &[u8]) {
    let Ok(source) = SourceBuffer::from_bytes("<fuzz>", bytes) else {
        return;
    };
    let mut lexer = Lexer::new(source, 64 * 1024 * 1024);
    loop {
        match lexer.next_token() {
            Ok(Some(_)) => {}
            Ok(None) | Err(_) => return,
        }
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use super::{Lexer, needs_scalar_validation};
    use crate::cif::{ParseErrorCode, SourceBuffer};

    #[test]
    fn bulk_control_precheck_has_no_forbidden_ascii_false_negatives() {
        for byte in 0_u8..=31 {
            assert!(needs_scalar_validation(&[
                b'A', byte, b'Z', b'A', b'A', b'A', b'A', b'A',
            ]));
        }
        assert!(needs_scalar_validation(&[
            b'A', 0x7f, b'Z', b'A', b'A', b'A', b'A', b'A',
        ]));
        assert!(needs_scalar_validation("A\u{0085}Z".as_bytes()));
        for byte in 32_u8..=126 {
            assert!(!needs_scalar_validation(&[
                b'A', byte, b'Z', b'A', b'A', b'A', b'A', b'A',
            ]));
        }
    }

    #[test]
    fn fused_unquoted_scan_rejects_ascii_and_unicode_controls() {
        for text in ["bad\0value", "bad\u{007f}value", "bad\u{0085}value"] {
            let source = SourceBuffer::from_text("<test>", text);
            let error = Lexer::new(source, usize::MAX)
                .next_token()
                .expect_err("unquoted control characters must be rejected");
            assert_eq!(error.code(), ParseErrorCode::InvalidCharacter);
        }
    }
}

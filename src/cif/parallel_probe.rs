use memchr::{memchr, memchr3};

use super::lexer::is_whitespace;
use super::source::SourceBuffer;

pub(super) fn loop_reaches(
    source: &SourceBuffer,
    start: usize,
    target: usize,
    max_token_bytes: usize,
) -> bool {
    let bytes = source.as_str().as_bytes();
    let mut cursor = start;
    while cursor < bytes.len() {
        while cursor < bytes.len() && is_whitespace(bytes[cursor]) {
            cursor += 1;
        }
        if cursor == bytes.len() {
            return false;
        }
        if bytes[cursor] == b'#' {
            let Some(relative) = memchr3(b'\n', b'\r', 0x7f, &bytes[cursor + 1..]) else {
                return false;
            };
            let terminator = bytes[cursor + 1 + relative];
            if terminator == 0x7f {
                return false;
            }
            cursor += relative + 1;
            continue;
        }

        let token_start = cursor;
        cursor = match bytes[cursor] {
            b'\'' | b'"' => match quoted_end(bytes, cursor) {
                Some(end) => end,
                None => return false,
            },
            b';' if is_line_start(bytes, cursor) => match text_end(bytes, cursor) {
                Some(end) => end,
                None => return false,
            },
            _ => match bare_end(bytes, cursor) {
                Some(end) if is_value(&bytes[cursor..end]) => end,
                _ => return false,
            },
        };
        if cursor - token_start > max_token_bytes {
            return false;
        }
        if cursor >= target {
            return true;
        }
    }
    false
}

fn starts_ascii_case(token: &[u8], prefix: &[u8]) -> bool {
    token
        .get(..prefix.len())
        .is_some_and(|candidate| candidate.eq_ignore_ascii_case(prefix))
}

fn is_line_start(bytes: &[u8], offset: usize) -> bool {
    offset == 0 || matches!(bytes[offset - 1], b'\n' | b'\r')
}

fn quoted_end(bytes: &[u8], start: usize) -> Option<usize> {
    let quote = bytes[start];
    let mut cursor = start + 1;
    while cursor < bytes.len() {
        let relative = memchr3(quote, b'\n', b'\r', &bytes[cursor..])?;
        cursor += relative;
        if bytes[cursor] == quote && (cursor + 1 == bytes.len() || is_whitespace(bytes[cursor + 1]))
        {
            return Some(cursor + 1);
        }
        if matches!(bytes[cursor], b'\n' | b'\r') {
            return None;
        }
        cursor += 1;
    }
    None
}

fn text_end(bytes: &[u8], start: usize) -> Option<usize> {
    let mut cursor = start + 1;
    while let Some(relative) = memchr(b';', &bytes[cursor..]) {
        cursor += relative;
        if is_line_start(bytes, cursor) {
            return (cursor + 1 == bytes.len() || is_whitespace(bytes[cursor + 1]))
                .then_some(cursor + 1);
        }
        cursor += 1;
    }
    None
}

fn bare_end(bytes: &[u8], start: usize) -> Option<usize> {
    let mut cursor = start;
    while cursor < bytes.len() {
        let byte = bytes[cursor];
        if byte > b' ' && byte != 0x7f {
            cursor += 1;
        } else if is_whitespace(byte) {
            break;
        } else {
            return None;
        }
    }
    Some(cursor)
}

fn is_value(token: &[u8]) -> bool {
    if token.starts_with(b"_") {
        return false;
    }
    match token.first().map(u8::to_ascii_lowercase) {
        Some(b'l') if token.eq_ignore_ascii_case(b"loop_") => false,
        Some(b'g') if token.eq_ignore_ascii_case(b"global_") => false,
        Some(b'd') if starts_ascii_case(token, b"data_") => false,
        Some(b's') if token.eq_ignore_ascii_case(b"stop_") => false,
        Some(b's') if starts_ascii_case(token, b"save_") => false,
        _ => true,
    }
}

#[cfg(test)]
mod tests {
    use super::loop_reaches;
    use crate::cif::{SourceBuffer, lexer::Lexer, token::TokenKind};

    fn generic(source: &SourceBuffer, start: usize, target: usize) -> bool {
        let mut lexer = Lexer::new(source.clone(), usize::MAX);
        lexer.set_offset(start);
        while let Ok(Some(token)) = lexer.next_token() {
            if token.kind != TokenKind::Value {
                return false;
            }
            if token.span.end >= target {
                return true;
            }
        }
        false
    }

    #[test]
    fn conservative_probe_matches_valid_loop_prefixes() {
        let cases = [
            "one two three stop_ tail",
            "one # comment\ntwo 'three four' \"five\" loop_",
            "one\n;text\nbody\n;\ntwo data_next",
            "one 'stop_' \"loop_\" _tag",
            "one\r\ntwo\r\nthree\r\n",
        ];
        for text in cases {
            let source = SourceBuffer::from_text("<probe-test>", text);
            for target in 1..=source.len() {
                assert_eq!(
                    loop_reaches(&source, 0, target, usize::MAX),
                    generic(&source, 0, target),
                    "text={text:?}, target={target}"
                );
            }
        }
    }
}

//! Text that matches Python's `json.dumps(sort_keys=True)` and `repr` byte for byte.

use std::fmt::Write;

/// Append `value` as a JSON string with `ensure_ascii` escaping.
pub(super) fn string(out: &mut String, value: &str) {
    out.push('"');
    if value.bytes().all(|byte| matches!(byte, b' '..=b'~') && byte != b'"' && byte != b'\\') {
        out.push_str(value);
        out.push('"');
        return;
    }
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            ' '..='~' => out.push(character),
            _ => {
                let mut units = [0_u16; 2];
                for unit in character.encode_utf16(&mut units) {
                    let _ = write!(out, "\\u{unit:04x}");
                }
            }
        }
    }
    out.push('"');
}

/// Append a list of strings in compact JSON.
pub(super) fn strings<'a>(out: &mut String, values: impl IntoIterator<Item = &'a str>) {
    out.push('[');
    for (index, value) in values.into_iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        string(out, value);
    }
    out.push(']');
}

/// Append Python's `repr` of a string.
pub(super) fn repr(out: &mut String, value: &str) {
    let quote = if value.contains('\'') && !value.contains('"') { '"' } else { '\'' };
    out.push(quote);
    for character in value.chars() {
        let code = character as u32;
        match character {
            _ if character == quote || character == '\\' => {
                out.push('\\');
                out.push(character);
            }
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            _ if code < 0x20 || (0x7f..0xa0).contains(&code) => {
                let _ = write!(out, "\\x{code:02x}");
            }
            _ => out.push(character),
        }
    }
    out.push(quote);
}

/// Python's `repr` of a list or tuple of strings.
pub(super) fn repr_sequence<'a>(values: impl IntoIterator<Item = &'a str>, tuple: bool) -> String {
    let mut out = String::from(if tuple { "(" } else { "[" });
    let mut count = 0;
    for value in values {
        if count > 0 {
            out.push_str(", ");
        }
        repr(&mut out, value);
        count += 1;
    }
    if tuple && count == 1 {
        out.push(',');
    }
    out.push(if tuple { ')' } else { ']' });
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escapes_like_python() {
        let mut out = String::new();
        string(&mut out, "C1'\"\\\u{7f}é😀\n");
        assert_eq!(out, r#""C1'\"\\\u007f\u00e9\ud83d\ude00\n""#);
        let mut out = String::new();
        repr(&mut out, "O3'");
        assert_eq!(out, "\"O3'\"");
        assert_eq!(repr_sequence(["A", "1", ""], false), "['A', '1', '']");
    }
}

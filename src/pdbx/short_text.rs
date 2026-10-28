//! Immutable text that stores short values inline.
//!
//! Atom-site identifiers are almost always a few bytes long. Storing them inline lets
//! model construction decode 190,000 atoms without one heap allocation per field.
//! Comparison, ordering and hashing match `str`.

use std::borrow::Borrow;
use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::ops::Deref;

const INLINE: usize = 22;

#[derive(Clone)]
enum Storage {
    Inline(u8, [u8; INLINE]),
    Heap(Box<str>),
}

/// Text of any length; values of up to 22 bytes need no heap allocation.
#[derive(Clone)]
pub struct ShortText(Storage);

impl ShortText {
    /// Copy `text`.
    #[must_use]
    pub fn new(text: &str) -> Self {
        let source = text.as_bytes();
        let length = source.len();
        if length > INLINE {
            return Self(Storage::Heap(text.into()));
        }
        // A copy of `length` bytes would call `memcpy`. Two fixed-size copies that may
        // overlap compile to plain loads and stores instead.
        let mut bytes = [0; INLINE];
        match length {
            0 => {}
            1..=3 => {
                bytes[0] = source[0];
                bytes[length / 2] = source[length / 2];
                bytes[length - 1] = source[length - 1];
            }
            4..=7 => {
                copy_fixed::<4>(&mut bytes, source, 0);
                copy_fixed::<4>(&mut bytes, source, length - 4);
            }
            8..=16 => {
                copy_fixed::<8>(&mut bytes, source, 0);
                copy_fixed::<8>(&mut bytes, source, length - 8);
            }
            _ => {
                copy_fixed::<16>(&mut bytes, source, 0);
                copy_fixed::<16>(&mut bytes, source, length - 16);
            }
        }
        // `length` is at most `INLINE`, so it fits in a byte.
        Self(Storage::Inline(length as u8, bytes))
    }

    /// Copy `text` with ASCII letters in upper case, as `str::to_ascii_uppercase` does.
    #[must_use]
    pub fn ascii_uppercase(text: &str) -> Self {
        let mut value = Self::new(text);
        match &mut value.0 {
            Storage::Inline(length, bytes) => bytes[..usize::from(*length)].make_ascii_uppercase(),
            Storage::Heap(text) => text.make_ascii_uppercase(),
        }
        value
    }

    /// Return the text.
    #[inline]
    #[must_use]
    pub fn as_str(&self) -> &str {
        match &self.0 {
            // Inline bytes always hold a complete UTF-8 string copied from a `str`.
            Storage::Inline(length, bytes) => {
                std::str::from_utf8(&bytes[..usize::from(*length)]).unwrap_or_default()
            }
            Storage::Heap(text) => text,
        }
    }

    /// Return the UTF-8 bytes of the text.
    #[inline]
    pub(crate) fn as_bytes(&self) -> &[u8] {
        match &self.0 {
            Storage::Inline(length, bytes) => &bytes[..usize::from(*length)],
            Storage::Heap(text) => text.as_bytes(),
        }
    }
}

/// Copy `N` bytes of `source` from `start` to the same place in `target`.
#[inline]
fn copy_fixed<const N: usize>(target: &mut [u8; INLINE], source: &[u8], start: usize) {
    let mut piece = [0; N];
    piece.copy_from_slice(&source[start..start + N]);
    target[start..start + N].copy_from_slice(&piece);
}

impl Deref for ShortText {
    type Target = str;

    #[inline]
    fn deref(&self) -> &str {
        self.as_str()
    }
}

impl AsRef<str> for ShortText {
    fn as_ref(&self) -> &str {
        self.as_str()
    }
}

impl Borrow<str> for ShortText {
    fn borrow(&self) -> &str {
        self.as_str()
    }
}

impl From<&str> for ShortText {
    fn from(text: &str) -> Self {
        Self::new(text)
    }
}

impl From<&ShortText> for String {
    fn from(text: &ShortText) -> Self {
        text.as_str().to_owned()
    }
}

impl PartialEq for ShortText {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        match (&self.0, &other.0) {
            // `new` zeroes the bytes past `length`, so equal texts have equal arrays,
            // and comparing the fixed-size arrays avoids a `memcmp` call.
            (Storage::Inline(left_length, left), Storage::Inline(right_length, right)) => {
                left_length == right_length && left == right
            }
            _ => self.as_bytes() == other.as_bytes(),
        }
    }
}

impl Eq for ShortText {}

impl PartialEq<str> for ShortText {
    #[inline]
    fn eq(&self, other: &str) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialEq<&str> for ShortText {
    #[inline]
    fn eq(&self, other: &&str) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialEq<String> for ShortText {
    #[inline]
    fn eq(&self, other: &String) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialEq<ShortText> for str {
    #[inline]
    fn eq(&self, other: &ShortText) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialEq<ShortText> for &str {
    #[inline]
    fn eq(&self, other: &ShortText) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialEq<ShortText> for String {
    #[inline]
    fn eq(&self, other: &ShortText) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialOrd for ShortText {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ShortText {
    // `str` orders by bytes, so byte order keeps every sorted map unchanged.
    fn cmp(&self, other: &Self) -> Ordering {
        self.as_bytes().cmp(other.as_bytes())
    }
}

impl Hash for ShortText {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.as_str().hash(state);
    }
}

impl fmt::Debug for ShortText {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self.as_str(), formatter)
    }
}

impl fmt::Display for ShortText {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self.as_str(), formatter)
    }
}

impl serde::Serialize for ShortText {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(self.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::ShortText;

    #[test]
    fn matches_str_for_inline_and_heap_values() {
        for text in [
            "",
            "C",
            "HETATM",
            "1234567890123456789012",
            "12345678901234567890123",
        ] {
            let value = ShortText::new(text);
            assert_eq!(value.as_str(), text);
            assert!(value == text);
            assert_eq!(format!("{value:?}"), format!("{text:?}"));
        }
        let alphabet = "abcdefghijklmnopqrstuvwxyz0123456789";
        for length in 0..=alphabet.len() {
            assert_eq!(
                ShortText::new(&alphabet[..length]).as_str(),
                &alphabet[..length]
            );
        }
        assert!(ShortText::ascii_uppercase("ca") == ShortText::new("CA"));
        assert!(ShortText::new("CA") != ShortText::new("C"));
        assert!(ShortText::new("AB") < ShortText::new("B"));
        assert!(ShortText::new("A") < ShortText::new("AA"));
        assert_eq!(std::mem::size_of::<ShortText>(), 24);
        assert_eq!(std::mem::size_of::<Option<ShortText>>(), 24);
    }
}

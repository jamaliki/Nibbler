#![no_main]

use _core::cif::{CifEntry, CifValue, format_text, parse};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|bytes: &[u8]| {
    let Ok(value) = std::str::from_utf8(bytes) else {
        return;
    };
    let Ok(formatted) = format_text(value) else {
        return;
    };
    let source = if formatted.starts_with(';') {
        format!("data_fuzz\n_fuzz.value\n{formatted}\n")
    } else {
        format!("data_fuzz\n_fuzz.value {formatted}\n")
    };
    let document = parse(source.as_bytes()).expect("formatted text must parse");
    let CifEntry::Item(item) = &document.blocks()[0].entries()[0] else {
        panic!("fuzz document must contain one scalar item");
    };
    let CifValue::Text(parsed) = item.value() else {
        panic!("formatted present text must remain present text");
    };
    assert_eq!(parsed.as_str(), value);
});

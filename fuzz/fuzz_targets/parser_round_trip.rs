#![no_main]

use _core::cif::{parse, write_canonical};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|bytes: &[u8]| {
    let Ok(document) = parse(bytes) else {
        return;
    };
    let canonical = write_canonical(&document)
        .expect("every successfully parsed document must be canonically representable");
    let reparsed = parse(canonical.as_bytes())
        .expect("canonical serialization of a parsed document must parse");
    assert_eq!(reparsed, document);
});

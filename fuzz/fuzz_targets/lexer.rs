#![no_main]

use _core::cif::fuzz_lexer;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|bytes: &[u8]| {
    fuzz_lexer(bytes);
});

//! Phase 1 integration tests for strict CIF syntax and logical round trips.
//!
//! Assertions intentionally use panic-based helpers, as permitted for tests by
//! `ENGINEERING.md`; production code retains the workspace-wide deny policy.

#![allow(clippy::expect_used, clippy::panic)]

use _core::cif::{
    BlockKind, CifDocument, CifEntry, CifValueRef, Limits, ParseErrorCode, ParseOptions,
    QuoteStyle, SourceBuffer, WriteErrorCode, format_text, parse, parse_source, parse_with_options,
    write_canonical,
};
use proptest::prelude::*;

const MINIMAL: &[u8] = include_bytes!("fixtures/syntax/minimal.cif");
const MISSING_VALUES: &[u8] = include_bytes!("fixtures/syntax/missing_values.cif");
const MALFORMED_LOOP: &[u8] = include_bytes!("fixtures/syntax/malformed_loop.cif");
const CHEMISTRY: &[u8] = include_bytes!("fixtures/chemistry/ligand_ion_water.cif");

fn only_block(document: &CifDocument) -> &_core::cif::CifBlock {
    assert_eq!(document.blocks().len(), 1);
    &document.blocks()[0]
}

fn text(value: CifValueRef<'_>) -> &str {
    value.as_text().expect("test expects a present text value")
}

#[test]
fn parses_scalar_loop_quotes_and_text_field() {
    let document = parse(MINIMAL).expect("committed minimal fixture must parse");
    let block = only_block(&document);
    assert_eq!(block.kind(), BlockKind::Data);
    assert_eq!(block.code(), Some("minimal"));
    assert_eq!(block.entries().len(), 2);

    let CifEntry::Item(item) = &block.entries()[0] else {
        panic!("first entry must be a scalar item");
    };
    assert_eq!(item.tag(), "_entry.id");
    assert_eq!(text(item.value().as_ref()), "minimal");

    let CifEntry::Loop(cif_loop) = &block.entries()[1] else {
        panic!("second entry must be a loop");
    };
    assert_eq!(
        cif_loop.tags(),
        ["_nibbler_example.id", "_nibbler_example.value"]
    );
    assert_eq!(cif_loop.row_count(), 2);
    assert_eq!(
        text(
            cif_loop
                .row(0)
                .expect("row exists")
                .get(1)
                .expect("value exists")
        ),
        "quoted value"
    );
    assert_eq!(
        text(
            cif_loop
                .row(1)
                .expect("row exists")
                .get(1)
                .expect("value exists")
        ),
        "a multiline value\nwhose delimiter starts in column one"
    );
    let CifValueRef::Text(multiline) = cif_loop
        .row(1)
        .expect("row exists")
        .get(1)
        .expect("value exists")
    else {
        panic!("multiline value must be present text");
    };
    assert_eq!(multiline.quote_style(), QuoteStyle::TextField);
}

#[test]
fn preserves_both_missing_kinds_and_quoted_empty_text() {
    let document = parse(MISSING_VALUES).expect("missing-value fixture must parse");
    let CifEntry::Loop(cif_loop) = &only_block(&document).entries()[0] else {
        panic!("fixture entry must be a loop");
    };
    assert_eq!(cif_loop.row_count(), 2);
    assert_eq!(
        cif_loop.row(0).expect("row exists").get(1),
        Some(CifValueRef::Unknown)
    );
    assert_eq!(
        cif_loop.row(0).expect("row exists").get(2),
        Some(CifValueRef::NotApplicable)
    );
    assert_eq!(
        text(
            cif_loop
                .row(0)
                .expect("row exists")
                .get(3)
                .expect("value exists")
        ),
        ""
    );
    assert_eq!(
        cif_loop.row(1).expect("row exists").get(1),
        Some(CifValueRef::NotApplicable)
    );
    assert_eq!(
        cif_loop.row(1).expect("row exists").get(2),
        Some(CifValueRef::Unknown)
    );
}

#[test]
fn malformed_loop_has_stable_code_and_location() {
    let error = parse(MALFORMED_LOOP).expect_err("incomplete final row must fail");
    assert_eq!(error.code(), ParseErrorCode::LoopValueCount);
    assert_eq!(error.code_str(), "CIF_LOOP_VALUE_COUNT");
    assert_eq!(error.span().line(), 6);
    assert_eq!(error.span().column(), 9);
}

#[test]
fn parses_global_blocks_multiple_data_blocks_and_save_frames() {
    let source = br#"
global_
_global.value 1
data_FIRST
_outside.value 2
save_example
_outside.value 3
loop_
_inside.id
_inside.value
1 'loop_'
stop_
save_
data_second
_second.value "data_not_a_block"
"#;
    let document = parse(source).expect("valid controls and frame must parse");
    assert_eq!(document.blocks().len(), 3);
    assert_eq!(document.blocks()[0].kind(), BlockKind::Global);
    assert_eq!(document.blocks()[1].code(), Some("FIRST"));
    assert_eq!(document.blocks()[2].code(), Some("second"));

    let CifEntry::Frame(frame) = &document.blocks()[1].entries()[1] else {
        panic!("second data-block entry must be the save frame");
    };
    assert_eq!(frame.code(), "example");
    assert_eq!(frame.entries().len(), 2);
}

#[test]
fn loops_with_equal_prefixes_remain_distinct() {
    let document = parse(b"data_x\nloop_\n_same.a\n1\nloop_\n_same.b\n2\n")
        .expect("distinct loops with unique tags are valid");
    assert!(matches!(
        only_block(&document).entries(),
        [CifEntry::Loop(_), CifEntry::Loop(_)]
    ));
}

#[test]
fn control_words_are_case_insensitive_only_when_unquoted() {
    let document = parse(b"DaTa_Mixed\n_item.one 'LoOp_'\n_item.two \"SaVe_frame\" # GLOBAL_\n")
        .expect("quoted controls are values and comments are ignored");
    let block = only_block(&document);
    assert_eq!(block.code(), Some("Mixed"));
    let CifEntry::Item(first) = &block.entries()[0] else {
        panic!("entry must be an item");
    };
    assert_eq!(text(first.value().as_ref()), "LoOp_");
}

#[test]
fn strict_grammar_failures_have_stable_codes() {
    let cases: &[(&[u8], ParseErrorCode)] = &[
        (b"# only a comment\n", ParseErrorCode::ExpectedBlock),
        (b"data_\n", ParseErrorCode::EmptyBlockCode),
        (b"data_x\n_ 1\n", ParseErrorCode::InvalidTag),
        (b"data_x\n_a 1\n_A 2\n", ParseErrorCode::DuplicateTag),
        (
            b"data_x\nloop_\n_a\n1\nloop_\n_A\n2\n",
            ParseErrorCode::DuplicateTag,
        ),
        (b"data_x\n_a\n", ParseErrorCode::MissingItemValue),
        (b"data_x\nloop_\n1\n", ParseErrorCode::MissingLoopTag),
        (
            b"data_x\nloop_\n_a\nstop_\n",
            ParseErrorCode::MissingLoopValue,
        ),
        (b"data_x\nsave_a\nsave_b\n", ParseErrorCode::NestedSaveFrame),
        (
            b"data_x\nsave_a\n_a 1\n",
            ParseErrorCode::UnterminatedSaveFrame,
        ),
        (b"data_x\nsave_\n", ParseErrorCode::UnexpectedSaveEnd),
        (b"data_x\nstop_\n", ParseErrorCode::UnexpectedStop),
        (b"data_x\nloose\n", ParseErrorCode::UnexpectedValue),
        (
            b"data_x\n_a 'unterminated\n",
            ParseErrorCode::UnterminatedQuotedValue,
        ),
        (
            b"data_x\n_a\n;unterminated\n",
            ParseErrorCode::UnterminatedTextField,
        ),
        (
            b"data_x\n_a\n;value\n;trailing\n",
            ParseErrorCode::InvalidTextFieldTerminator,
        ),
        (b"data_x\n_a bad\0value\n", ParseErrorCode::InvalidCharacter),
        (
            "data_x\n_a 'bad\u{0085}value'\n".as_bytes(),
            ParseErrorCode::InvalidCharacter,
        ),
    ];

    for (source, expected_code) in cases {
        let error = parse(source).expect_err("table case must be rejected");
        assert_eq!(error.code(), *expected_code, "source: {source:?}");
    }
}

#[test]
fn rejects_invalid_utf8_at_the_first_invalid_byte() {
    let error = parse(b"data_x\n_a \xff\n").expect_err("invalid UTF-8 must fail");
    assert_eq!(error.code(), ParseErrorCode::InvalidUtf8);
    assert_eq!(error.span().line(), 2);
    assert_eq!(error.span().column(), 4);
}

#[test]
fn source_name_and_unicode_display_column_survive_diagnostics() {
    let source = SourceBuffer::from_text("named.cif", "data_x\n_title 'é'\n_a 1\n_A 2\n");
    let error = parse_source(source).expect_err("duplicate tag must fail");
    assert_eq!(error.source_name(), "named.cif");
    assert_eq!(error.span().line(), 4);
    assert_eq!(error.span().column(), 1);
    assert!(error.to_string().starts_with("named.cif:4:1:"));
}

#[test]
fn diagnostics_count_cr_and_crlf_as_single_line_endings() {
    for source_text in ["data_x\r_a 1\r_A 2\r", "data_x\r\n_a 1\r\n_A 2\r\n"] {
        let error = parse(source_text.as_bytes()).expect_err("duplicate tag must fail");
        assert_eq!(error.code(), ParseErrorCode::DuplicateTag);
        assert_eq!(error.span().line(), 3);
        assert_eq!(error.span().column(), 1);
    }
}

#[test]
fn every_resource_limit_is_enforced() {
    type ConfigureLimit = fn(&mut Limits);
    let cases: &[(&[u8], ConfigureLimit)] = &[
        (b"data_x\n", |limits| limits.source_bytes = 2),
        (b"data_x\n", |limits| limits.token_bytes = 4),
        (b"data_x\n", |limits| limits.blocks = 0),
        (b"data_x\nsave_a\nsave_\n", |limits| limits.frames = 0),
        (b"data_x\nloop_\n_a\n1\n", |limits| limits.loops = 0),
        (b"data_x\n_a 1\n", |limits| limits.tags = 0),
        (b"data_x\nloop_\n_a\n1\n", |limits| limits.loop_columns = 0),
        (b"data_x\nloop_\n_a\n1\n", |limits| limits.rows = 0),
        (b"data_x\n_a 1\n", |limits| limits.values = 0),
    ];
    for (source, configure) in cases {
        let mut limits = Limits::default();
        configure(&mut limits);
        let error = parse_with_options(source, ParseOptions { limits })
            .expect_err("configured hard limit must fail");
        assert_eq!(error.code(), ParseErrorCode::ResourceLimit);
    }
}

#[test]
fn formatter_uses_one_deterministic_lossless_quoting_policy() {
    let cases = [
        ("bare", "bare"),
        ("", "''"),
        ("?", "'?'"),
        (".", "'.'"),
        ("hello world", "'hello world'"),
        ("it's", "\"it's\""),
        ("ADENOSINE-5'-TRIPHOSPHATE", "\"ADENOSINE-5'-TRIPHOSPHATE\""),
        ("both '\" quotes", ";both '\" quotes\n;"),
        ("line one\nline two", ";line one\nline two\n;"),
        ("carriage return\r", ";carriage return\r\r;"),
        ("DATA_value", "'DATA_value'"),
    ];
    for (value, expected) in cases {
        assert_eq!(format_text(value).expect("text is representable"), expected);
    }
    let delimiter_error = format_text("first\n;terminates")
        .expect_err("delimiter collision is not losslessly representable");
    assert_eq!(delimiter_error.code(), WriteErrorCode::UnrepresentableText);
    assert_eq!(
        format_text("forbidden\0control")
            .expect_err("forbidden control must fail")
            .code(),
        WriteErrorCode::UnrepresentableText
    );
}

#[test]
fn all_valid_fixtures_round_trip_logically_and_canonically() {
    for fixture in [MINIMAL, MISSING_VALUES, CHEMISTRY] {
        let parsed = parse(fixture).expect("valid fixture must parse");
        let canonical = write_canonical(&parsed).expect("valid document must serialize");
        assert!(canonical.ends_with("#\n"));
        assert!(!canonical.ends_with("\n\n"));
        let reparsed = parse(canonical.as_bytes()).expect("canonical output must reparse");
        assert_eq!(reparsed, parsed);
        assert_eq!(
            write_canonical(&reparsed).expect("reparsed document must serialize"),
            canonical,
            "canonical serialization must be idempotent"
        );
    }
}

fn one_text_value(document: &CifDocument) -> &str {
    let CifEntry::Item(item) = &only_block(document).entries()[0] else {
        panic!("generated document contains one scalar item");
    };
    text(item.value().as_ref())
}

proptest! {
    #[test]
    fn arbitrary_representable_text_survives_format_and_parse(value in any::<String>()) {
        let Ok(formatted) = format_text(&value) else {
            return Ok(());
        };
        let source = if formatted.starts_with(';') {
            format!("data_property\n_property.value\n{formatted}\n")
        } else {
            format!("data_property\n_property.value {formatted}\n")
        };
        let document = parse(source.as_bytes()).expect("formatter output must always parse");
        prop_assert_eq!(one_text_value(&document), value);
    }
}

use super::document::{
    BlockKind, CifDocument, CifEntry, CifFrame, CifItem, CifLoop, CifValue, CifValueRef,
    QuoteStyle, TextColumn,
};
use super::error::{WriteError, WriteErrorCode};
use super::numeric::{parse_float, parse_integer};
use std::borrow::Cow;

/// Formats one value; text written as it is stays borrowed from the document.
type ValueFormatter = for<'a> fn(CifValueRef<'a>) -> Result<Cow<'a, str>, WriteError>;

/// Serialize an order-preserving logical document using deterministic CIF 1.1 syntax.
///
/// The result uses `\n` line endings, `#` entry separators, and exactly one final
/// newline. Unknown and not-applicable values remain `?` and `.` respectively.
///
/// # Errors
///
/// Returns [`WriteErrorCode::InvalidDocument`] for an invalid internal structure,
/// [`WriteErrorCode::UnrepresentableText`] when text would collide with a
/// semicolon-field terminator, or [`WriteErrorCode::NonFiniteFloat`] for NaN or infinity.
pub fn write_canonical(document: &CifDocument) -> Result<String, WriteError> {
    write_document(document, format_value_ref)
}

/// Serialize a document while retaining source order and valid original value lexemes.
///
/// Layout is deterministic rather than byte-for-byte identical: comments and source
/// whitespace are not retained by the logical model. Parsed quote styles and numeric
/// spellings are reused when they still represent the value; constructed or changed
/// values use canonical formatting.
///
/// # Errors
///
/// Returns the same structured errors as [`write_canonical`].
pub fn write_preserving(document: &CifDocument) -> Result<String, WriteError> {
    write_document(document, format_preserving_value_ref)
}

fn write_document(document: &CifDocument, format: ValueFormatter) -> Result<String, WriteError> {
    if document.blocks().is_empty() {
        return Err(WriteError::new(
            WriteErrorCode::InvalidDocument,
            "a CIF document must contain at least one block",
        ));
    }
    let mut output = Output(Vec::new());
    for block in document.blocks() {
        match block.kind() {
            BlockKind::Data => {
                let Some(code) = block.code() else {
                    return Err(WriteError::new(
                        WriteErrorCode::InvalidDocument,
                        "a data block must have a block code",
                    ));
                };
                output.push_str("data_");
                output.push_str(code);
                output.push('\n');
            }
            BlockKind::Global => output.push_str("global_\n"),
        }
        for entry in block.entries() {
            output.push_str("#\n");
            write_entry(&mut output, entry, format)?;
        }
        output.push_str("#\n");
    }
    String::from_utf8(output.0).map_err(|_| {
        WriteError::new(
            WriteErrorCode::InvalidDocument,
            "the written document is not valid UTF-8",
        )
    })
}

/// The document text being written, held as bytes so that rows can copy stored text in
/// pieces of fixed size. Every byte comes from a `str`, so the text stays UTF-8.
struct Output(Vec<u8>);

impl Output {
    #[inline]
    fn push_str(&mut self, text: &str) {
        self.0.extend_from_slice(text.as_bytes());
    }

    #[inline]
    fn push(&mut self, character: char) {
        self.push_str(character.encode_utf8(&mut [0; 4]));
    }

    fn reserve(&mut self, additional: usize) {
        self.0.reserve(additional);
    }

    /// Append `data[start..end]`, or return `None` when the range is not in `data`.
    #[inline]
    fn push_stored(&mut self, data: &[u8], start: usize, end: usize) -> Option<()> {
        self.0.extend_from_slice(data.get(start..end)?);
        Some(())
    }
}

/// Format one present text value with the least intrusive deterministic quote style.
///
/// The returned semicolon form, if any, must begin in column one. Callers serializing
/// complete documents should normally use [`write_canonical`], which enforces this
/// layout rule.
///
/// # Errors
///
/// Returns [`WriteErrorCode::UnrepresentableText`] for forbidden control characters
/// or text containing a line that begins with `;`, which CIF 1.1 cannot represent in a
/// semicolon-delimited value without altering the content.
pub fn format_text(text: &str) -> Result<String, WriteError> {
    format_text_cow(text).map(Cow::into_owned)
}

fn format_text_cow(text: &str) -> Result<Cow<'_, str>, WriteError> {
    Ok(match text_form(text)? {
        TextForm::Bare => Cow::Borrowed(text),
        TextForm::Quoted(quote) => Cow::Owned(format!("{quote}{text}{quote}")),
        TextForm::Field(newline) => Cow::Owned(format!(";{text}{newline};")),
    })
}

/// How [`format_text`] writes a text value.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TextForm {
    /// As it is.
    Bare,
    /// Between two of the given quote character.
    Quoted(char),
    /// As a semicolon text field whose closing line starts after the given line ending.
    Field(char),
}

#[inline]
fn text_form(text: &str) -> Result<TextForm, WriteError> {
    // One pass decides text of printable ASCII and spaces, which needs none of the checks
    // in `other_text_form` for control characters, line endings and text-field delimiters.
    let (mut printable, mut space, mut single, mut double) = (true, false, false, false);
    for &byte in text.as_bytes() {
        printable &= matches!(byte, b' '..=b'~');
        space |= byte == b' ';
        single |= byte == b'\'';
        double |= byte == b'"';
    }
    if printable {
        return Ok(
            if !(space || single || double) && plain_start(text.as_bytes()) {
                TextForm::Bare
            } else if !single {
                TextForm::Quoted('\'')
            } else if !double {
                TextForm::Quoted('"')
            } else {
                TextForm::Field('\n')
            },
        );
    }
    other_text_form(text)
}

/// Return the form of text that holds a byte other than printable ASCII or a space.
fn other_text_form(text: &str) -> Result<TextForm, WriteError> {
    validate_text_characters(text)?;
    if is_bare_value(text) {
        return Ok(TextForm::Bare);
    }
    let is_single_line = !text.contains('\n') && !text.contains('\r');
    if is_single_line && !text.contains('\'') {
        return Ok(TextForm::Quoted('\''));
    }
    if is_single_line && !text.contains('"') {
        return Ok(TextForm::Quoted('"'));
    }
    if contains_text_field_delimiter(text) {
        return Err(WriteError::new(
            WriteErrorCode::UnrepresentableText,
            "text contains a line beginning with ';', which terminates a CIF 1.1 text field",
        ));
    }
    Ok(TextForm::Field(if text.ends_with('\r') {
        '\r'
    } else {
        '\n'
    }))
}

fn write_entry(
    output: &mut Output,
    entry: &CifEntry,
    format: ValueFormatter,
) -> Result<(), WriteError> {
    match entry {
        CifEntry::Item(item) => write_item(output, item, format),
        CifEntry::Loop(cif_loop) => write_loop(output, cif_loop, format),
        CifEntry::Frame(frame) => write_frame(output, frame, format),
    }
}

fn write_item(
    output: &mut Output,
    item: &CifItem,
    format: ValueFormatter,
) -> Result<(), WriteError> {
    let value = format(item.value().as_ref())?;
    output.push_str(item.tag());
    if is_text_field(&value) {
        output.push('\n');
    } else {
        output.push(' ');
    }
    output.push_str(&value);
    output.push('\n');
    Ok(())
}

fn write_loop(
    output: &mut Output,
    cif_loop: &CifLoop,
    format: ValueFormatter,
) -> Result<(), WriteError> {
    if cif_loop.column_count() == 0
        || cif_loop.value_count() == 0
        || !cif_loop
            .value_count()
            .is_multiple_of(cif_loop.column_count())
    {
        return Err(WriteError::new(
            WriteErrorCode::InvalidDocument,
            "a CIF loop must contain tags and complete non-empty rows",
        ));
    }
    output.push_str("loop_\n");
    for tag in cif_loop.tags() {
        output.push_str(tag);
        output.push('\n');
    }

    if write_text_rows(output, cif_loop)? {
        return Ok(());
    }
    for row_index in 0..cif_loop.row_count() {
        let Some(row) = cif_loop.row(row_index) else {
            return Err(WriteError::new(
                WriteErrorCode::InvalidDocument,
                "a CIF loop must contain complete rows",
            ));
        };
        let mut line_is_empty = true;
        for value in row.iter() {
            let value = format(value)?;
            if is_text_field(&value) {
                if !line_is_empty {
                    output.push('\n');
                }
                output.push_str(&value);
                output.push('\n');
                line_is_empty = true;
            } else {
                if !line_is_empty {
                    output.push(' ');
                }
                output.push_str(&value);
                line_is_empty = false;
            }
        }
        if !line_is_empty {
            output.push('\n');
        }
    }
    Ok(())
}

/// Write the rows of a loop whose columns are all stored as text, and return whether the
/// loop was one.
///
/// Column storage holds unquoted text, which both [`write_canonical`] and
/// [`write_preserving`] format with [`format_text`], so the rows need no formatter. A
/// column whose stored entries all pass [`plain_bytes`] and [`plain_start`] writes each
/// value as it is. Other columns choose a form once for each run of rows that share a
/// stored entry, and write quotes around the text instead of formatting a new string.
fn write_text_rows(output: &mut Output, cif_loop: &CifLoop) -> Result<bool, WriteError> {
    let rows = cif_loop.row_count();
    let mut columns = Vec::with_capacity(cif_loop.column_count());
    // A row holds each column's text and one separator after it.
    let mut row_bytes = cif_loop.column_count();
    for column_index in 0..cif_loop.column_count() {
        let Some(column) = cif_loop.text_column(column_index) else {
            return Ok(false);
        };
        let Some(column) = RowText::new(column, rows) else {
            return Ok(false);
        };
        row_bytes = row_bytes.saturating_add(column.mean_len());
        columns.push(column);
    }
    output.reserve(row_bytes.saturating_mul(rows));
    for row in 0..rows {
        let mut line_is_empty = true;
        for column in &mut columns {
            let (text, form) = match column.value(row)? {
                Cell::Stored(start, end) => {
                    if !line_is_empty {
                        output.push(' ');
                    }
                    if output
                        .push_stored(column.data.as_bytes(), start, end)
                        .is_none()
                    {
                        return Err(stored_entry_error());
                    }
                    line_is_empty = false;
                    continue;
                }
                Cell::Text(text, form) => (text, form),
            };
            if let TextForm::Field(newline) = form {
                if !line_is_empty {
                    output.push('\n');
                }
                output.push(';');
                output.push_str(text);
                output.push(newline);
                output.push_str(";\n");
                line_is_empty = true;
                continue;
            }
            if !line_is_empty {
                output.push(' ');
            }
            if let TextForm::Quoted(quote) = form {
                output.push(quote);
                output.push_str(text);
                output.push(quote);
            } else {
                output.push_str(text);
            }
            line_is_empty = false;
        }
        if !line_is_empty {
            output.push('\n');
        }
    }
    Ok(true)
}

/// One text column as [`write_text_rows`] reads it.
struct RowText<'a> {
    data: &'a str,
    offsets: &'a [u32],
    indices: &'a [u32],
    mask: Option<&'a [u8]>,
    /// Whether every stored entry is written as it is.
    plain: bool,
    /// The storage key and form of the last value given a form.
    last: Option<(u64, TextForm)>,
}

impl<'a> RowText<'a> {
    /// Read `column`, or return `None` when it has fewer than `rows` rows.
    fn new(column: TextColumn<'a>, rows: usize) -> Option<Self> {
        let (data, offsets, indices, mask, _) = column.parts();
        if indices.len() < rows || mask.is_some_and(|mask| mask.len() < rows) {
            return None;
        }
        let plain = plain_bytes(data.as_bytes())
            && offsets.windows(2).all(|pair| {
                data.get(pair[0] as usize..pair[1] as usize)
                    .is_some_and(|entry| plain_start(entry.as_bytes()))
            });
        Some(Self {
            data,
            offsets,
            indices,
            mask,
            plain,
            last: None,
        })
    }

    /// The mean length of a stored entry, or 1 when there are none.
    fn mean_len(&self) -> usize {
        self.data
            .len()
            .checked_div(self.offsets.len().saturating_sub(1))
            .unwrap_or(1)
    }

    /// Return the value in `row`. Missing states are written as `?` and `.`, as
    /// [`TextColumn::value`] reads them.
    #[inline]
    fn value(&mut self, row: usize) -> Result<Cell<'a>, WriteError> {
        let state = self
            .mask
            .and_then(|mask| mask.get(row).copied())
            .unwrap_or(0);
        match state {
            1 => return Ok(Cell::Text(".", TextForm::Bare)),
            2 => return Ok(Cell::Text("?", TextForm::Bare)),
            _ => {}
        }
        let index = self.indices.get(row).copied().unwrap_or(0);
        let range = match index as usize {
            0 => Some((0, 0)),
            entry => self
                .offsets
                .get(entry - 1)
                .zip(self.offsets.get(entry))
                .map(|(start, end)| (*start as usize, *end as usize)),
        };
        let Some((start, end)) = range else {
            return Err(stored_entry_error());
        };
        if self.plain && index != 0 {
            return Ok(Cell::Stored(start, end));
        }
        let Some(text) = self.data.get(start..end) else {
            return Err(stored_entry_error());
        };
        let key = (u64::from(state) << 32) | u64::from(index);
        let form = match self.last {
            Some((last, form)) if last == key => form,
            _ => text_form(text)?,
        };
        self.last = Some((key, form));
        Ok(Cell::Text(text, form))
    }
}

/// One value of a text column: a stored entry to write as it is, given by its byte range,
/// or text and the form to write it in.
enum Cell<'a> {
    Stored(usize, usize),
    Text(&'a str, TextForm),
}

fn stored_entry_error() -> WriteError {
    WriteError::new(
        WriteErrorCode::InvalidDocument,
        "a text column refers to an entry it does not hold",
    )
}

fn write_frame(
    output: &mut Output,
    frame: &CifFrame,
    format: ValueFormatter,
) -> Result<(), WriteError> {
    output.push_str("save_");
    output.push_str(frame.code());
    output.push('\n');
    for entry in frame.entries() {
        output.push_str("#\n");
        write_entry(output, entry, format)?;
    }
    output.push_str("#\nsave_\n");
    Ok(())
}

pub(super) fn format_value(value: &CifValue) -> Result<String, WriteError> {
    format_value_ref(value.as_ref()).map(Cow::into_owned)
}

fn format_value_ref(value: CifValueRef<'_>) -> Result<Cow<'_, str>, WriteError> {
    match value {
        CifValueRef::Text(text) => format_text_cow(text.as_str()),
        CifValueRef::Integer(number, _) => Ok(Cow::Owned(number.to_string())),
        CifValueRef::Float(number, uncertainty, _) => {
            if !number.is_finite() {
                return Err(WriteError::new(
                    WriteErrorCode::NonFiniteFloat,
                    "CIF missing states must be used instead of NaN or infinity",
                ));
            }
            let mut mantissa = number.to_string();
            let exponent = mantissa
                .find(['e', 'E'])
                .map(|byte_offset| mantissa.split_off(byte_offset));
            if !mantissa.contains('.') {
                mantissa.push_str(".0");
            }
            if let Some(uncertainty) = uncertainty {
                mantissa.push('(');
                mantissa.push_str(&uncertainty.digits().to_string());
                mantissa.push(')');
            }
            if let Some(exponent) = exponent {
                mantissa.push_str(&exponent);
            }
            Ok(Cow::Owned(mantissa))
        }
        CifValueRef::Unknown => Ok(Cow::Borrowed("?")),
        CifValueRef::NotApplicable => Ok(Cow::Borrowed(".")),
    }
}

fn format_preserving_value_ref(value: CifValueRef<'_>) -> Result<Cow<'_, str>, WriteError> {
    match value {
        CifValueRef::Text(text) => format_preserving_text(text.as_str(), text.quote_style()),
        CifValueRef::Integer(number, Some(original))
            if parse_integer(original.as_str()) == Some(number) =>
        {
            Ok(Cow::Borrowed(original.as_str()))
        }
        CifValueRef::Float(number, uncertainty, Some(original))
            if parse_float(original.as_str()).is_some_and(|parsed| {
                parsed.to_bits() == number.to_bits()
                    && original_uncertainty(original.as_str())
                        == Some(uncertainty.map(|value| value.digits()))
            }) =>
        {
            Ok(Cow::Borrowed(original.as_str()))
        }
        _ => format_value_ref(value),
    }
}

fn format_preserving_text(text: &str, style: QuoteStyle) -> Result<Cow<'_, str>, WriteError> {
    validate_text_characters(text)?;
    match style {
        QuoteStyle::Unquoted if is_bare_value(text) => Ok(Cow::Borrowed(text)),
        QuoteStyle::Single if valid_quoted_content(text, b'\'') => {
            Ok(Cow::Owned(format!("'{text}'")))
        }
        QuoteStyle::Double if valid_quoted_content(text, b'"') => {
            Ok(Cow::Owned(format!("\"{text}\"")))
        }
        QuoteStyle::TextField if !contains_text_field_delimiter(text) => {
            let delimiter_newline = if text.ends_with('\r') { '\r' } else { '\n' };
            Ok(Cow::Owned(format!(";{text}{delimiter_newline};")))
        }
        _ => format_text_cow(text),
    }
}

fn valid_quoted_content(text: &str, quote: u8) -> bool {
    !text.contains(['\n', '\r'])
        && !text.as_bytes().iter().enumerate().any(|(index, byte)| {
            *byte == quote
                && text
                    .as_bytes()
                    .get(index + 1)
                    .is_some_and(|next| super::lexer::is_whitespace(*next))
        })
}

fn original_uncertainty(text: &str) -> Option<Option<u64>> {
    let mantissa = text.split_once(['e', 'E']).map_or(text, |(head, _)| head);
    if !mantissa.ends_with(')') {
        return Some(None);
    }
    let open = mantissa.rfind('(')?;
    let digits = &mantissa[open + 1..mantissa.len() - 1];
    (!digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit()))
        .then(|| digits.parse::<u64>().ok())
        .flatten()
        .map(Some)
}

fn validate_text_characters(text: &str) -> Result<(), WriteError> {
    // ASCII text: the control characters are U+0000-U+001F and U+007F
    if text.is_ascii()
        && !text
            .bytes()
            .any(|byte| (byte < 0x20 || byte == 0x7f) && !matches!(byte, b'\t' | b'\n' | b'\r'))
    {
        return Ok(());
    }
    if let Some(character) = text
        .chars()
        .find(|character| character.is_control() && !matches!(*character, '\t' | '\n' | '\r'))
    {
        return Err(WriteError::new(
            WriteErrorCode::UnrepresentableText,
            format!(
                "text contains forbidden control character U+{:04X}",
                u32::from(character)
            ),
        ));
    }
    Ok(())
}

/// Return whether every byte is printable ASCII other than a space or quote.
#[inline]
fn plain_bytes(bytes: &[u8]) -> bool {
    bytes
        .iter()
        .all(|byte| matches!(byte, b'!' | b'#'..=b'&' | b'('..=b'~'))
}

/// Return whether text without spaces or quotes may be written bare: it is not empty, `?`
/// or `.`, and does not start with a reserved character or word.
#[inline]
fn plain_start(bytes: &[u8]) -> bool {
    let Some(&first) = bytes.first() else {
        return false;
    };
    !matches!(first, b'_' | b'$' | b'#' | b';' | b'[' | b']')
        && bytes != b"?"
        && bytes != b"."
        && !(matches!(first | 0x20, b'd' | b'g' | b'l' | b's')
            && [&b"data_"[..], b"global_", b"loop_", b"save_", b"stop_"]
                .iter()
                .any(|prefix| {
                    bytes
                        .get(..prefix.len())
                        .is_some_and(|start| start.eq_ignore_ascii_case(prefix))
                }))
}

fn is_bare_value(text: &str) -> bool {
    // in ASCII text the whitespace characters are U+0009-U+000D and U+0020
    let has_quote_or_space = if text.is_ascii() {
        text.bytes()
            .any(|byte| matches!(byte, b'\'' | b'"' | b' ' | b'\t'..=b'\r'))
    } else {
        text.contains(['\'', '"']) || text.chars().any(char::is_whitespace)
    };
    if text.is_empty() || text == "?" || text == "." || has_quote_or_space {
        return false;
    }
    if text
        .as_bytes()
        .first()
        .is_some_and(|byte| matches!(byte, b'_' | b'$' | b'#' | b';' | b'[' | b']' | b'\'' | b'"'))
    {
        return false;
    }
    !["data_", "global_", "loop_", "save_", "stop_"]
        .iter()
        .any(|prefix| starts_ascii_case(text, prefix))
}

fn contains_text_field_delimiter(text: &str) -> bool {
    text.as_bytes()
        .windows(2)
        .any(|pair| matches!(pair, [b'\n' | b'\r', b';']))
}

fn is_text_field(formatted: &str) -> bool {
    formatted.starts_with(';')
}

fn starts_ascii_case(text: &str, prefix: &str) -> bool {
    text.get(..prefix.len())
        .is_some_and(|candidate| candidate.eq_ignore_ascii_case(prefix))
}

#[cfg(test)]
mod tests {
    use super::{WriteErrorCode, write_canonical, write_preserving};
    use crate::cif::document::{
        CifBlock, CifDocument, CifEntry, CifItem, CifLoop, CifValue, StandardUncertainty,
    };

    fn document_with(values: Vec<CifValue>) -> CifDocument {
        let entries = values
            .into_iter()
            .enumerate()
            .map(|(index, value)| {
                CifEntry::Item(CifItem::new(format!("_number.value_{index}"), value))
            })
            .collect();
        CifDocument::new(vec![CifBlock::data("numbers".to_owned(), entries)])
    }

    #[test]
    fn canonicalizes_typed_numbers_without_original_lexemes() {
        let document = document_with(vec![
            CifValue::Integer(-7, None),
            CifValue::Float(2.0, None, None),
            CifValue::Float(1.25, Some(StandardUncertainty::new(3)), None),
        ]);
        assert_eq!(
            write_canonical(&document),
            Ok("data_numbers\n#\n_number.value_0 -7\n#\n_number.value_1 2.0\n#\n_number.value_2 1.25(3)\n#\n".to_owned())
        );
    }

    #[test]
    fn preserving_output_reuses_valid_source_lexemes() -> Result<(), Box<dyn std::error::Error>> {
        let source =
            b"data_styles\n_a.one 'alpha beta'\n_a.two \"gamma\"\n_a.three\n;delta\nepsilon\n;\n";
        let document = crate::cif::parse(source)?;
        let output = write_preserving(&document)?;
        assert!(output.contains("_a.one 'alpha beta'\n"));
        assert!(output.contains("_a.two \"gamma\"\n"));
        assert!(output.contains("_a.three\n;delta\nepsilon\n;\n"));
        assert_eq!(crate::cif::parse(output.as_bytes()), Ok(document));
        Ok(())
    }

    #[test]
    fn text_columns_write_as_owned_values() -> Result<(), Box<dyn std::error::Error>> {
        use crate::cif::{TextCell, TextColumnBuilder};
        use TextCell::{NotApplicable, Text, Unknown};
        let rows = [
            [Text("ATOM"), Text("O5'"), Text("1.5"), Unknown],
            [Text("ATOM"), Text("O5'"), Text("-2"), NotApplicable],
            [Text("HETATM"), Text(""), Text("3"), Text("a b")],
            [
                Text("HETATM"),
                Text("it's \"x\""),
                Text("4"),
                Text("line\nnext"),
            ],
            [Text("ATOM"), Text("data_x"), Text("5"), Text("end\r")],
            [Text("ATOM"), Text("_x"), Text("6"), Text("?")],
            [Text("ATOM"), Text("tab\there"), Text("7"), Text(".")],
            [Text("ATOM"), Text("\u{e9}"), Text("8"), Text("#x")],
        ];
        let mut builders: Vec<_> = (0..4)
            .map(|_| TextColumnBuilder::with_capacity(rows.len()))
            .collect();
        let mut values = Vec::new();
        for row in &rows {
            for (builder, cell) in builders.iter_mut().zip(row) {
                builder
                    .push(*cell)
                    .ok_or("the column outgrew its offsets")?;
                values.push(match cell {
                    Text(text) => CifValue::text(*text),
                    NotApplicable => CifValue::NotApplicable,
                    Unknown => CifValue::Unknown,
                });
            }
        }
        let tags: Vec<String> = ["a", "b", "c", "d"]
            .iter()
            .map(|name| format!("_t.{name}"))
            .collect();
        let document = |cif_loop| {
            CifDocument::new(vec![CifBlock::data(
                "t".to_owned(),
                vec![CifEntry::Loop(cif_loop)],
            )])
        };
        let columns = document(CifLoop::from_columns(
            tags.clone(),
            builders
                .into_iter()
                .map(TextColumnBuilder::finish)
                .collect(),
            rows.len(),
        ));
        let owned = document(CifLoop::new(tags, values));
        assert_eq!(write_canonical(&columns)?, write_canonical(&owned)?);
        assert_eq!(write_preserving(&columns)?, write_preserving(&owned)?);
        Ok(())
    }

    #[test]
    fn rejects_non_finite_floats() {
        for number in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let error = write_canonical(&document_with(vec![CifValue::Float(number, None, None)]));
            assert_eq!(
                error.map_err(|error| error.code()),
                Err(WriteErrorCode::NonFiniteFloat)
            );
        }
    }

    #[test]
    fn rejects_empty_documents_and_loops() {
        assert_eq!(
            write_canonical(&CifDocument::new(Vec::new())).map_err(|error| error.code()),
            Err(WriteErrorCode::InvalidDocument)
        );
        let document = CifDocument::new(vec![CifBlock::data(
            "invalid".to_owned(),
            vec![CifEntry::Loop(CifLoop::new(Vec::new(), Vec::new()))],
        )]);
        assert_eq!(
            write_canonical(&document).map_err(|error| error.code()),
            Err(WriteErrorCode::InvalidDocument)
        );
    }
}

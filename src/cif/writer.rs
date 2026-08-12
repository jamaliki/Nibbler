use super::document::{
    BlockKind, CifDocument, CifEntry, CifFrame, CifItem, CifLoop, CifValue, CifValueRef, QuoteStyle,
};
use super::error::{WriteError, WriteErrorCode};
use super::numeric::{parse_float, parse_integer};

type ValueFormatter = fn(CifValueRef<'_>) -> Result<String, WriteError>;

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
    let mut output = String::new();
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
    Ok(output)
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
    validate_text_characters(text)?;
    if is_bare_value(text) {
        return Ok(text.to_owned());
    }
    let is_single_line = !text.contains('\n') && !text.contains('\r');
    if is_single_line && !text.contains('\'') {
        return Ok(format!("'{text}'"));
    }
    if is_single_line && !text.contains('"') {
        return Ok(format!("\"{text}\""));
    }
    if contains_text_field_delimiter(text) {
        return Err(WriteError::new(
            WriteErrorCode::UnrepresentableText,
            "text contains a line beginning with ';', which terminates a CIF 1.1 text field",
        ));
    }
    let delimiter_newline = if text.ends_with('\r') { '\r' } else { '\n' };
    Ok(format!(";{text}{delimiter_newline};"))
}

fn write_entry(
    output: &mut String,
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
    output: &mut String,
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
    output: &mut String,
    cif_loop: &CifLoop,
    format: ValueFormatter,
) -> Result<(), WriteError> {
    if cif_loop.column_count() == 0
        || cif_loop.value_count() == 0
        || cif_loop.value_count() % cif_loop.column_count() != 0
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

fn write_frame(
    output: &mut String,
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
    format_value_ref(value.as_ref())
}

fn format_value_ref(value: CifValueRef<'_>) -> Result<String, WriteError> {
    match value {
        CifValueRef::Text(text) => format_text(text.as_str()),
        CifValueRef::Integer(number, _) => Ok(number.to_string()),
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
            Ok(mantissa)
        }
        CifValueRef::Unknown => Ok("?".to_owned()),
        CifValueRef::NotApplicable => Ok(".".to_owned()),
    }
}

fn format_preserving_value_ref(value: CifValueRef<'_>) -> Result<String, WriteError> {
    match value {
        CifValueRef::Text(text) => format_preserving_text(text.as_str(), text.quote_style()),
        CifValueRef::Integer(number, Some(original))
            if parse_integer(original.as_str()) == Some(number) =>
        {
            Ok(original.as_str().to_owned())
        }
        CifValueRef::Float(number, uncertainty, Some(original))
            if parse_float(original.as_str()).is_some_and(|parsed| {
                parsed.to_bits() == number.to_bits()
                    && original_uncertainty(original.as_str())
                        == Some(uncertainty.map(|value| value.digits()))
            }) =>
        {
            Ok(original.as_str().to_owned())
        }
        _ => format_value_ref(value),
    }
}

fn format_preserving_text(text: &str, style: QuoteStyle) -> Result<String, WriteError> {
    validate_text_characters(text)?;
    match style {
        QuoteStyle::Unquoted if is_bare_value(text) => Ok(text.to_owned()),
        QuoteStyle::Single if valid_quoted_content(text, b'\'') => Ok(format!("'{text}'")),
        QuoteStyle::Double if valid_quoted_content(text, b'"') => Ok(format!("\"{text}\"")),
        QuoteStyle::TextField if !contains_text_field_delimiter(text) => {
            let delimiter_newline = if text.ends_with('\r') { '\r' } else { '\n' };
            Ok(format!(";{text}{delimiter_newline};"))
        }
        _ => format_text(text),
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

fn is_bare_value(text: &str) -> bool {
    if text.is_empty()
        || text == "?"
        || text == "."
        || text.contains(['\'', '"'])
        || text.chars().any(char::is_whitespace)
    {
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

/// Parse a DDL2 integer without accepting decimal or exponent notation.
pub(crate) fn parse_integer(text: &str) -> Option<i64> {
    text.parse().ok()
}

/// Powers of ten that `f64` represents exactly.
const EXACT_POWERS: [f64; 23] = [
    1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16,
    1e17, 1e18, 1e19, 1e20, 1e21, 1e22,
];

/// Parse `[+-]digits[.digits]` when its digits form an integer below 2^53.
///
/// Both that integer and the power of ten are exact in `f64`, so one IEEE division
/// rounds the decimal value correctly and matches `str::parse`. Other spellings
/// return `None` and take the general path.
#[inline]
fn parse_simple_decimal(text: &str) -> Option<f64> {
    let bytes = text.as_bytes();
    let (negative, digits) = match bytes.first()? {
        b'-' => (true, &bytes[1..]),
        b'+' => (false, &bytes[1..]),
        _ => (false, bytes),
    };
    let mut mantissa: u64 = 0;
    let mut count = 0_usize;
    let mut fraction = None;
    for (index, &byte) in digits.iter().enumerate() {
        match byte {
            b'0'..=b'9' => {
                count += 1;
                if count > 16 {
                    return None;
                }
                mantissa = mantissa * 10 + u64::from(byte - b'0');
            }
            b'.' if fraction.is_none() => fraction = Some(digits.len() - index - 1),
            _ => return None,
        }
    }
    if count == 0 || mantissa >= 1 << 53 {
        return None;
    }
    let value = mantissa as f64 / EXACT_POWERS[fraction.unwrap_or(0)];
    Some(if negative { -value } else { value })
}

/// Parse a finite CIF float, including a parenthesized standard uncertainty.
pub(crate) fn parse_float(text: &str) -> Option<f64> {
    if let Some(value) = parse_simple_decimal(text) {
        return Some(value);
    }
    let exponent = text.find(['e', 'E']).unwrap_or(text.len());
    let (mantissa, suffix) = text.split_at(exponent);
    if !mantissa.ends_with(')') {
        let number = text.parse::<f64>().ok()?;
        return number.is_finite().then_some(number);
    }
    let uncertainty = mantissa.rfind('(')?;
    let uncertainty_digits = &mantissa[uncertainty + 1..mantissa.len() - 1];
    if uncertainty_digits.is_empty()
        || !uncertainty_digits.bytes().all(|byte| byte.is_ascii_digit())
    {
        return None;
    }
    let mantissa = &mantissa[..uncertainty];
    let number = format!("{mantissa}{suffix}").parse::<f64>().ok()?;
    number.is_finite().then_some(number)
}

#[cfg(test)]
mod tests {
    use super::{parse_float, parse_integer};

    #[test]
    fn parses_dictionary_numeric_forms_without_accepting_non_finite_values() {
        assert_eq!(parse_integer("-12"), Some(-12));
        assert_eq!(parse_integer("1.0"), None);
        assert_eq!(parse_float("1.25(3)e2"), Some(125.0));
        assert_eq!(parse_float("."), None);
        assert_eq!(parse_float("NaN"), None);
    }
}

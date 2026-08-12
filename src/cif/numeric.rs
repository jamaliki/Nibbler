/// Parse a DDL2 integer without accepting decimal or exponent notation.
pub(crate) fn parse_integer(text: &str) -> Option<i64> {
    text.parse().ok()
}

/// Parse a finite CIF float, including a parenthesized standard uncertainty.
pub(crate) fn parse_float(text: &str) -> Option<f64> {
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

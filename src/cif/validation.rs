use super::document::CifDocument;
use super::schema::{SchemaError, SchemaName, loaded_schema};
use super::validation_engine::validate;

const MAX_DIAGNOSTICS: usize = 10_000;
const COVERAGE: &[&str] = &[
    "known-categories-items",
    "ddl2-type-patterns",
    "mandatory-categories-items",
    "enumerations",
    "numeric-ranges",
    "category-keys",
    "parent-child-links",
];

/// Validation diagnostic severity.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Severity {
    /// Input violates a supported dictionary constraint.
    Error,
    /// Validation coverage or source metadata deserves attention.
    Warning,
}

impl Severity {
    /// Return the stable Python-facing value.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Error => "error",
            Self::Warning => "warning",
        }
    }
}

/// One stable dictionary-validation finding.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Diagnostic {
    pub(super) code: &'static str,
    pub(super) severity: Severity,
    pub(super) message: String,
    pub(super) context: Vec<String>,
}

impl Diagnostic {
    /// Return the stable machine-readable code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.code
    }

    /// Return the finding severity.
    #[must_use]
    pub const fn severity(&self) -> Severity {
        self.severity
    }

    /// Return the human-readable detail.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }

    /// Return block, frame, category, item, row, and dictionary context when known.
    #[must_use]
    pub fn context(&self) -> &[String] {
        &self.context
    }
}

/// An immutable result for one exact compiled dictionary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ValidationReport {
    schema_name: String,
    dictionary_version: String,
    coverage: Vec<String>,
    diagnostics: Vec<Diagnostic>,
}

impl ValidationReport {
    /// Return the selected schema name.
    #[must_use]
    pub fn schema_name(&self) -> &str {
        &self.schema_name
    }

    /// Return the exact pinned dictionary version.
    #[must_use]
    pub fn dictionary_version(&self) -> &str {
        &self.dictionary_version
    }

    /// Return the constraint families applied by this report.
    #[must_use]
    pub fn coverage(&self) -> &[String] {
        &self.coverage
    }

    /// Return findings in deterministic source/check order.
    #[must_use]
    pub fn diagnostics(&self) -> &[Diagnostic] {
        &self.diagnostics
    }

    /// Return whether no error diagnostic was produced.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        !self
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.severity == Severity::Error)
    }
}

/// Validate a logical CIF document against one built-in DDL2 schema.
///
/// Syntax is already guaranteed by construction of [`CifDocument`]. Dictionary checks
/// are bounded to 10,000 findings.
///
/// # Errors
///
/// Returns [`SchemaError`] if the embedded schema cannot be loaded or verified.
pub fn validate_document(
    document: &CifDocument,
    schema_name: SchemaName,
) -> Result<ValidationReport, SchemaError> {
    let schema = loaded_schema(schema_name)?;
    let metadata = schema.dictionary().metadata();
    let (mut diagnostics, truncated) = validate(document, schema);
    if truncated {
        diagnostics.push(Diagnostic {
            code: "CIF_SCHEMA_DIAGNOSTICS_TRUNCATED",
            severity: Severity::Warning,
            message: format!("validation stopped after {MAX_DIAGNOSTICS} findings"),
            context: vec![format!("dictionary={}", metadata.dictionary_name())],
        });
    }
    Ok(ValidationReport {
        schema_name: metadata.schema_name().to_owned(),
        dictionary_version: metadata.version().to_owned(),
        coverage: COVERAGE.iter().map(|value| (*value).to_owned()).collect(),
        diagnostics,
    })
}

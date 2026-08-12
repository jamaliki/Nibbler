//! Stable semantic diagnostics and validation reports.

use crate::cif::{Severity, ValidationReport};

pub(super) const MAX_DIAGNOSTICS: usize = 10_000;
const COVERAGE: &[&str] = &[
    "pdbx-ddl2-5.416",
    "entity-asym-coherence",
    "polymer-sequence-schemes",
    "nonpolymer-water-schemes",
    "branched-entity-graphs",
    "chemical-component-resolution",
    "component-atom-bonds",
    "atom-site-identity-coordinates",
    "struct-conn-endpoints",
];

/// Severity of one semantic-profile finding.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProfileSeverity {
    /// The document is not a valid strict PDBx coordinate model.
    Error,
    /// The document is valid but contains a profile-quality concern.
    Warning,
}

impl ProfileSeverity {
    /// Return the stable Python-facing value.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Error => "error",
            Self::Warning => "warning",
        }
    }
}

/// One stable dictionary or semantic-profile validation finding.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ProfileDiagnostic {
    code: String,
    severity: ProfileSeverity,
    message: String,
    context: Vec<String>,
}

impl ProfileDiagnostic {
    pub(crate) fn new(
        code: impl Into<String>,
        severity: ProfileSeverity,
        message: impl Into<String>,
        context: Vec<String>,
    ) -> Self {
        Self {
            code: code.into(),
            severity,
            message: message.into(),
            context,
        }
    }

    pub(crate) fn from_dictionary(finding: &crate::cif::Diagnostic) -> Self {
        Self::new(
            finding.code(),
            match finding.severity() {
                Severity::Error => ProfileSeverity::Error,
                Severity::Warning => ProfileSeverity::Warning,
            },
            finding.message(),
            finding.context().to_vec(),
        )
    }

    pub(crate) fn from_semantic(error: &super::SemanticError) -> Self {
        Self::new(
            error.code(),
            ProfileSeverity::Error,
            error.message(),
            error.context().to_vec(),
        )
    }

    /// Return the stable machine-readable code.
    #[must_use]
    pub fn code(&self) -> &str {
        &self.code
    }

    /// Return the finding severity.
    #[must_use]
    pub const fn severity(&self) -> ProfileSeverity {
        self.severity
    }

    /// Return the human-readable detail.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }

    /// Return category, item, row, and semantic identity context.
    #[must_use]
    pub fn context(&self) -> &[String] {
        &self.context
    }
}

/// Immutable result of dictionary and semantic-profile validation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ProfileValidationReport {
    dictionary_version: &'static str,
    coverage: &'static [&'static str],
    diagnostics: Vec<ProfileDiagnostic>,
}

impl ProfileValidationReport {
    pub(super) fn pdbx(diagnostics: Vec<ProfileDiagnostic>) -> Self {
        Self::new("5.416", COVERAGE, "PDBX_DIAGNOSTICS_TRUNCATED", diagnostics)
    }

    pub(crate) fn new(
        dictionary_version: &'static str,
        coverage: &'static [&'static str],
        truncation_code: &'static str,
        mut diagnostics: Vec<ProfileDiagnostic>,
    ) -> Self {
        if diagnostics.len() > MAX_DIAGNOSTICS {
            diagnostics.truncate(MAX_DIAGNOSTICS);
            diagnostics.push(ProfileDiagnostic::new(
                truncation_code,
                ProfileSeverity::Warning,
                format!("validation stopped after {MAX_DIAGNOSTICS} findings"),
                Vec::new(),
            ));
        }
        Self {
            dictionary_version,
            coverage,
            diagnostics,
        }
    }

    /// Return the exact dictionary version used by this profile.
    #[must_use]
    pub const fn dictionary_version(&self) -> &'static str {
        self.dictionary_version
    }

    /// Return the constraint families applied by this report.
    #[must_use]
    pub const fn coverage(&self) -> &'static [&'static str] {
        self.coverage
    }

    /// Return findings in deterministic check order.
    #[must_use]
    pub fn diagnostics(&self) -> &[ProfileDiagnostic] {
        &self.diagnostics
    }

    /// Return whether the model has no error diagnostic.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        !self
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.severity == ProfileSeverity::Error)
    }
}

pub(crate) fn dictionary_diagnostics(report: &ValidationReport) -> Vec<ProfileDiagnostic> {
    report
        .diagnostics()
        .iter()
        .map(ProfileDiagnostic::from_dictionary)
        .collect()
}

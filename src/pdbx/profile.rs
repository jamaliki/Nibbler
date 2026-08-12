//! Stable PDBx profile diagnostics and validation report.

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

/// Severity of one PDBx profile finding.
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

/// One stable dictionary or semantic PDBx validation finding.
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

/// Immutable result of dictionary and PDBx semantic validation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ProfileValidationReport {
    diagnostics: Vec<ProfileDiagnostic>,
}

impl ProfileValidationReport {
    pub(super) fn new(mut diagnostics: Vec<ProfileDiagnostic>) -> Self {
        if diagnostics.len() > MAX_DIAGNOSTICS {
            diagnostics.truncate(MAX_DIAGNOSTICS);
            diagnostics.push(ProfileDiagnostic::new(
                "PDBX_DIAGNOSTICS_TRUNCATED",
                ProfileSeverity::Warning,
                format!("validation stopped after {MAX_DIAGNOSTICS} findings"),
                Vec::new(),
            ));
        }
        Self { diagnostics }
    }

    /// Return the exact PDBx dictionary version used by this profile.
    #[must_use]
    pub const fn dictionary_version(&self) -> &'static str {
        "5.416"
    }

    /// Return the constraint families applied by this report.
    #[must_use]
    pub const fn coverage(&self) -> &'static [&'static str] {
        COVERAGE
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

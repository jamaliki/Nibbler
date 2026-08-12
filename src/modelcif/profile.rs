//! Stable ModelCIF validation report metadata.

use crate::pdbx::{ProfileDiagnostic, ProfileSeverity};

const MAX_DIAGNOSTICS: usize = 10_000;
const COVERAGE: &[&str] = &[
    "modelcif-ddl2-1.4.9",
    "pdbx-coordinate-semantics",
    "prediction-targets",
    "model-groups",
    "software-provenance",
    "protocol-data-flow",
    "template-mappings-alignments",
    "qa-definitions-values",
    "associated-files-archives",
    "confidence-b-factor-separation",
];

/// Immutable result of dictionary, PDBx, and ModelCIF semantic validation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ModelCifValidationReport {
    diagnostics: Vec<ProfileDiagnostic>,
}

impl ModelCifValidationReport {
    pub(super) fn new(mut diagnostics: Vec<ProfileDiagnostic>) -> Self {
        if diagnostics.len() > MAX_DIAGNOSTICS {
            diagnostics.truncate(MAX_DIAGNOSTICS);
            diagnostics.push(ProfileDiagnostic::new(
                "MODELCIF_DIAGNOSTICS_TRUNCATED",
                ProfileSeverity::Warning,
                format!("validation stopped after {MAX_DIAGNOSTICS} findings"),
                Vec::new(),
            ));
        }
        Self { diagnostics }
    }

    /// Return the exact pinned ModelCIF dictionary version.
    #[must_use]
    pub const fn dictionary_version(&self) -> &'static str {
        "1.4.9"
    }

    /// Return the semantic constraint families applied by this report.
    #[must_use]
    pub const fn coverage(&self) -> &'static [&'static str] {
        COVERAGE
    }

    /// Return findings in deterministic check order.
    #[must_use]
    pub fn diagnostics(&self) -> &[ProfileDiagnostic] {
        &self.diagnostics
    }

    /// Return whether no error diagnostic was emitted.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        !self
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.severity() == ProfileSeverity::Error)
    }
}

//! Strict ModelCIF validation orchestration and shared semantic indices.

use std::collections::{BTreeMap, BTreeSet};

use crate::cif::{CifDocument, SchemaError, SchemaName, validate_document as validate_cif};
use crate::pdbx::category::case_key;
use crate::pdbx::profile::dictionary_diagnostics;
use crate::pdbx::{
    ProfileDiagnostic, ProfileSeverity, ProfileValidationReport,
    validate_model as validate_pdbx_model,
};

use super::aggregate::ModelCifModel;
use super::model::TargetEntity;
use super::provenance::Software;
use super::qa::QaMetric;
use super::source::build_model;

/// Validate a document against the pinned dictionary and strict ModelCIF semantics.
///
/// # Errors
///
/// Returns an error only if the embedded dictionary cannot be loaded.
pub fn validate_document(document: &CifDocument) -> Result<ProfileValidationReport, SchemaError> {
    let dictionary = validate_cif(document, SchemaName::ModelCif)?;
    let mut diagnostics = dictionary_diagnostics(&dictionary);
    match build_model(document) {
        Ok(model) => diagnostics.extend(validate_semantics(&model)),
        Err(error) => diagnostics.push(ProfileDiagnostic::from_semantic(&error)),
    }
    Ok(report(diagnostics))
}

/// Validate an already constructed ModelCIF model.
#[must_use]
pub fn validate_model(model: &ModelCifModel) -> ProfileValidationReport {
    report(validate_semantics(model))
}

fn report(diagnostics: Vec<ProfileDiagnostic>) -> ProfileValidationReport {
    ProfileValidationReport::new(
        "1.4.9",
        &[
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
        ],
        "MODELCIF_DIAGNOSTICS_TRUNCATED",
        diagnostics,
    )
}

fn validate_semantics(model: &ModelCifModel) -> Vec<ProfileDiagnostic> {
    let mut validator = Validator::new(model);
    validator.audit_conform();
    validator.coordinates();
    validator.data_flow();
    validator.targets();
    validator.models();
    validator.provenance();
    validator.templates();
    validator.qa();
    validator.files();
    validator.diagnostics
}

pub(super) struct Validator<'a> {
    pub(super) model: &'a ModelCifModel,
    pub(super) diagnostics: Vec<ProfileDiagnostic>,
    pub(super) data: BTreeSet<i64>,
    pub(super) data_groups: BTreeSet<i64>,
    pub(super) models: BTreeSet<i64>,
    pub(super) groups: BTreeSet<i64>,
    pub(super) software: BTreeMap<i64, &'a Software>,
    pub(super) software_groups: BTreeSet<i64>,
    pub(super) target_entities: BTreeMap<String, &'a TargetEntity>,
    pub(super) target_instances: BTreeSet<String>,
    pub(super) metrics: BTreeMap<i64, &'a QaMetric>,
}

impl<'a> Validator<'a> {
    fn new(model: &'a ModelCifModel) -> Self {
        Self {
            model,
            diagnostics: Vec::new(),
            data: model.data.iter().map(|value| value.id).collect(),
            data_groups: model
                .data_groups
                .iter()
                .map(|value| value.group_id)
                .collect(),
            models: model.models.iter().map(|value| value.id).collect(),
            groups: model.model_groups.iter().map(|value| value.id).collect(),
            software: model
                .software
                .iter()
                .map(|value| (value.id, value))
                .collect(),
            software_groups: model
                .software_groups
                .iter()
                .map(|value| value.group_id)
                .collect(),
            target_entities: model
                .targets
                .iter()
                .map(|value| (case_key(&value.entity_id), value))
                .collect(),
            target_instances: model
                .target_instances
                .iter()
                .map(|value| case_key(&value.asym_id))
                .collect(),
            metrics: model
                .qa_metrics
                .iter()
                .map(|value| (value.id, value))
                .collect(),
        }
    }

    fn audit_conform(&mut self) {
        if !self
            .model
            .audit_conform()
            .iter()
            .any(|(name, version)| name.eq_ignore_ascii_case("mmcif_ma.dic") && version == "1.4.9")
        {
            self.error(
                "MODELCIF_AUDIT_CONFORM",
                "the ModelCIF profile requires an audit_conform row for mmcif_ma.dic 1.4.9",
                vec!["category=audit_conform".to_owned()],
            );
        }
    }

    fn coordinates(&mut self) {
        self.diagnostics.extend(
            validate_pdbx_model(&self.model.coordinates)
                .diagnostics()
                .iter()
                .cloned(),
        );
    }

    fn data_flow(&mut self) {
        self.require_nonempty("MODELCIF_DATA_REQUIRED", "ma_data", self.model.data.len());
        self.unique_i64(
            "MODELCIF_DATA_ID_DUPLICATE",
            "data",
            self.model.data.iter().map(|value| value.id),
        );
        self.unique_i64(
            "MODELCIF_DATA_GROUP_ORDINAL_DUPLICATE",
            "data-group ordinal",
            self.model.data_groups.iter().map(|value| value.ordinal_id),
        );
        for member in &self.model.data_groups {
            self.reference(
                "MODELCIF_DATA_GROUP_DATA",
                "data",
                member.data_id,
                self.data.contains(&member.data_id),
            );
        }
    }

    pub(super) fn require_nonempty(&mut self, code: &'static str, category: &str, count: usize) {
        if count == 0 {
            self.error(
                code,
                format!("strict ModelCIF requires _{category}"),
                vec![format!("category={category}")],
            );
        }
    }

    pub(super) fn unique_i64(
        &mut self,
        code: &'static str,
        label: &str,
        values: impl Iterator<Item = i64>,
    ) {
        let mut seen = BTreeSet::new();
        for value in values {
            if !seen.insert(value) {
                self.error(
                    code,
                    format!("{label} identifier {value} is not unique"),
                    vec![format!("{label}={value}")],
                );
            }
        }
    }

    pub(super) fn reference(&mut self, code: &'static str, label: &str, value: i64, present: bool) {
        if !present {
            self.error(
                code,
                format!("referenced {label} {value} does not exist"),
                vec![format!("{label}={value}")],
            );
        }
    }

    pub(super) fn error(
        &mut self,
        code: &'static str,
        message: impl Into<String>,
        context: Vec<String>,
    ) {
        self.diagnostics.push(ProfileDiagnostic::new(
            code,
            ProfileSeverity::Error,
            message,
            context,
        ));
    }
}

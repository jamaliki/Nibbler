//! PDBx profile validation above the pinned DDL2 dictionary layer.

use std::collections::{BTreeMap, BTreeSet};

use crate::cif::{
    CifDocument, SchemaError, SchemaName, Severity, validate_document as validate_cif,
};

use super::model::{
    AsymUnit, AtomSite, ComponentDefinition, ComponentResolution, Entity, PdbxModel,
};
use super::profile::{
    MAX_DIAGNOSTICS, ProfileDiagnostic, ProfileSeverity, ProfileValidationReport,
};
use super::source::build_model;

/// Validate a source document against both the pinned PDBx dictionary and semantic profile.
///
/// # Errors
///
/// Returns [`SchemaError`] only if the embedded PDBx dictionary cannot be loaded.
pub fn validate_document(document: &CifDocument) -> Result<ProfileValidationReport, SchemaError> {
    let dictionary = validate_cif(document, SchemaName::Pdbx)?;
    let mut diagnostics = dictionary
        .diagnostics()
        .iter()
        .map(|finding| {
            ProfileDiagnostic::new(
                finding.code(),
                match finding.severity() {
                    Severity::Error => ProfileSeverity::Error,
                    Severity::Warning => ProfileSeverity::Warning,
                },
                finding.message(),
                finding.context().to_vec(),
            )
        })
        .collect::<Vec<_>>();
    match build_model(document) {
        Ok(model) => diagnostics.extend(validate_semantics(&model)),
        Err(error) => diagnostics.push(ProfileDiagnostic::new(
            error.code(),
            ProfileSeverity::Error,
            error.message(),
            error.context().to_vec(),
        )),
    }
    Ok(ProfileValidationReport::new(diagnostics))
}

/// Validate one constructed PDBx model against semantic profile invariants.
#[must_use]
pub fn validate_model(model: &PdbxModel) -> ProfileValidationReport {
    ProfileValidationReport::new(validate_semantics(model))
}

fn validate_semantics(model: &PdbxModel) -> Vec<ProfileDiagnostic> {
    let mut validator = Validator {
        model,
        diagnostics: Vec::new(),
        entities: index(model.entities(), Entity::id),
        asym_units: index(model.asym_units(), AsymUnit::id),
        components: index(model.components(), ComponentDefinition::id),
    };
    validator.audit_conform();
    validator.identifiers();
    validator.entities();
    validator.asym_units();
    validator.components();
    validator.atom_sites();
    validator.connections();
    validator.diagnostics
}

pub(super) struct Validator<'a> {
    pub(super) model: &'a PdbxModel,
    pub(super) diagnostics: Vec<ProfileDiagnostic>,
    pub(super) entities: BTreeMap<String, &'a Entity>,
    pub(super) asym_units: BTreeMap<String, &'a AsymUnit>,
    pub(super) components: BTreeMap<String, &'a ComponentDefinition>,
}

impl Validator<'_> {
    fn audit_conform(&mut self) {
        match self.model.audit_dictionary.as_ref() {
            Some((name, version))
                if name.eq_ignore_ascii_case("mmcif_pdbx.dic") && version == "5.416" => {}
            Some((name, version)) => self.error(
                "PDBX_AUDIT_CONFORM",
                format!("expected mmcif_pdbx.dic 5.416, found {name} {version}"),
                vec!["category=audit_conform".to_owned()],
            ),
            None => self.error(
                "PDBX_AUDIT_CONFORM",
                "the PDBx profile requires _audit_conform for mmcif_pdbx.dic 5.416",
                vec!["category=audit_conform".to_owned()],
            ),
        }
    }

    fn identifiers(&mut self) {
        self.unique(
            "PDBX_ENTITY_ID_DUPLICATE",
            "entity",
            self.model.entities(),
            Entity::id,
        );
        self.unique(
            "PDBX_ASYM_ID_DUPLICATE",
            "asym",
            self.model.asym_units(),
            AsymUnit::id,
        );
        self.unique(
            "PDBX_ATOM_ID_DUPLICATE",
            "atom",
            self.model.atom_sites(),
            AtomSite::id,
        );
        let mut connections = BTreeSet::new();
        for connection in self.model.connections() {
            let (id, _) = connection.identity();
            if !connections.insert(fold(id)) {
                self.error(
                    "PDBX_CONNECTION_ID_DUPLICATE",
                    format!("connection identifier {id:?} is not unique"),
                    vec![format!("connection={id}")],
                );
            }
        }
    }

    fn components(&mut self) {
        for component in self.model.components() {
            if component.resolution() == ComponentResolution::Unresolved {
                self.error(
                    "PDBX_COMPONENT_UNRESOLVED",
                    format!("component {:?} has no embedded, caller, CCD, or minimal-registry definition", component.id()),
                    component_context(component),
                );
            }
            if component.component_type().is_none() {
                self.error(
                    "PDBX_COMPONENT_DEFINITION",
                    format!("component {:?} has no PDBx component type", component.id()),
                    component_context(component),
                );
            }
            let atoms = component
                .atoms()
                .iter()
                .map(|atom| fold(atom.atom_id()))
                .collect::<BTreeSet<_>>();
            if atoms.len() != component.atoms().len() {
                self.error(
                    "PDBX_COMPONENT_ATOM_DUPLICATE",
                    format!(
                        "component {:?} has duplicate atom identifiers",
                        component.id()
                    ),
                    component_context(component),
                );
            }
            for bond in component.bonds() {
                let (first, second) = bond.atom_ids();
                if !atoms.contains(&fold(first)) || !atoms.contains(&fold(second)) {
                    self.error(
                        "PDBX_COMPONENT_BOND_ENDPOINT",
                        format!(
                            "component {:?} bond {first:?}-{second:?} references an absent atom",
                            component.id()
                        ),
                        component_context(component),
                    );
                }
            }
        }
    }

    fn unique<T>(
        &mut self,
        code: &'static str,
        label: &str,
        values: &[T],
        identifier: impl Fn(&T) -> &str,
    ) {
        let mut seen = BTreeSet::new();
        for value in values {
            let id = identifier(value);
            if !seen.insert(fold(id)) {
                self.error(
                    code,
                    format!("{label} identifier {id:?} is not unique"),
                    vec![format!("{label}={id}")],
                );
            }
        }
    }

    pub(super) fn error(
        &mut self,
        code: &'static str,
        message: impl Into<String>,
        context: Vec<String>,
    ) {
        self.push(code, ProfileSeverity::Error, message, context);
    }

    pub(super) fn warning(
        &mut self,
        code: &'static str,
        message: impl Into<String>,
        context: Vec<String>,
    ) {
        self.push(code, ProfileSeverity::Warning, message, context);
    }

    fn push(
        &mut self,
        code: &'static str,
        severity: ProfileSeverity,
        message: impl Into<String>,
        context: Vec<String>,
    ) {
        if self.diagnostics.len() < MAX_DIAGNOSTICS {
            self.diagnostics
                .push(ProfileDiagnostic::new(code, severity, message, context));
        }
    }
}

fn index<T>(values: &[T], identifier: impl Fn(&T) -> &str) -> BTreeMap<String, &T> {
    values
        .iter()
        .map(|value| (fold(identifier(value)), value))
        .collect()
}

pub(super) fn fold(value: &str) -> String {
    value.to_ascii_lowercase()
}

pub(super) fn equal(left: &str, right: &str) -> bool {
    left.eq_ignore_ascii_case(right)
}

fn component_context(component: &ComponentDefinition) -> Vec<String> {
    vec![format!("component={}", component.id())]
}

pub(super) fn atom_context(atom: &AtomSite) -> Vec<String> {
    vec![format!("atom_site={}", atom.id())]
}

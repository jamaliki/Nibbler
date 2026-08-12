//! Core ModelCIF prediction records.

/// One typed `_ma_data` record used by the modeling data-flow graph.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Data {
    pub(super) id: i64,
    pub(super) name: String,
    pub(super) content_type: String,
}

impl Data {
    /// Return the stable data identifier.
    #[must_use]
    pub const fn id(&self) -> i64 {
        self.id
    }

    /// Return the caller-supplied data name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Return the ModelCIF content-type enumeration value.
    #[must_use]
    pub fn content_type(&self) -> &str {
        &self.content_type
    }
}

/// One modeled target entity and its source-data identity.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TargetEntity {
    pub(super) entity_id: String,
    pub(super) data_id: i64,
    pub(super) origin: String,
}

impl TargetEntity {
    /// Return the referenced `_entity.id`.
    #[must_use]
    pub fn entity_id(&self) -> &str {
        &self.entity_id
    }

    /// Return the target-sequence data identifier.
    #[must_use]
    pub const fn data_id(&self) -> i64 {
        self.data_id
    }

    /// Return the ModelCIF target-origin enumeration value.
    #[must_use]
    pub fn origin(&self) -> &str {
        &self.origin
    }
}

/// One label-space instance of a modeled target entity.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TargetInstance {
    pub(super) asym_id: String,
    pub(super) entity_id: String,
    pub(super) details: Option<String>,
}

impl TargetInstance {
    /// Return the label-space asymmetric-unit identifier.
    #[must_use]
    pub fn asym_id(&self) -> &str {
        &self.asym_id
    }

    /// Return the referenced target entity identifier.
    #[must_use]
    pub fn entity_id(&self) -> &str {
        &self.entity_id
    }

    /// Return optional free-text details about this target instance.
    #[must_use]
    pub fn details(&self) -> Option<&str> {
        self.details.as_deref()
    }
}

/// One deposited prediction model. Its `id` is `_ma_model_list.ordinal_id`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PredictionModel {
    pub(super) id: i64,
    pub(super) name: Option<String>,
    pub(super) assembly_id: Option<i64>,
    pub(super) model_type: String,
    pub(super) type_details: Option<String>,
    pub(super) data_id: i64,
}

impl PredictionModel {
    /// Return `_ma_model_list.ordinal_id`, the deposited model identity.
    #[must_use]
    pub const fn id(&self) -> i64 {
        self.id
    }

    /// Return the optional model name.
    #[must_use]
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// Return the referenced model-assembly identifier, when present.
    #[must_use]
    pub const fn assembly_id(&self) -> Option<i64> {
        self.assembly_id
    }

    /// Return the ModelCIF model-type enumeration value.
    #[must_use]
    pub fn model_type(&self) -> &str {
        &self.model_type
    }

    /// Return optional details refining the model type.
    #[must_use]
    pub fn type_details(&self) -> Option<&str> {
        self.type_details.as_deref()
    }

    /// Return the data identifier representing this model.
    #[must_use]
    pub const fn data_id(&self) -> i64 {
        self.data_id
    }
}

/// A named group of prediction models.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ModelGroup {
    pub(super) id: i64,
    pub(super) name: Option<String>,
    pub(super) details: Option<String>,
}

impl ModelGroup {
    /// Return the stable model-group identifier.
    #[must_use]
    pub const fn id(&self) -> i64 {
        self.id
    }

    /// Return the optional group name.
    #[must_use]
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// Return optional free-text group details.
    #[must_use]
    pub fn details(&self) -> Option<&str> {
        self.details.as_deref()
    }
}

/// One explicit membership edge between a model and a group.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ModelGroupLink {
    pub(super) model_id: i64,
    pub(super) group_id: i64,
}

impl ModelGroupLink {
    /// Return the member model identifier.
    #[must_use]
    pub const fn model_id(&self) -> i64 {
        self.model_id
    }

    /// Return the containing group identifier.
    #[must_use]
    pub const fn group_id(&self) -> i64 {
        self.group_id
    }
}

/// One declared representative of a model group.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ModelRepresentative {
    pub(super) id: i64,
    pub(super) group_id: i64,
    pub(super) model_id: i64,
    pub(super) selection_criteria: String,
}

impl ModelRepresentative {
    /// Return the stable representative-row identifier.
    #[must_use]
    pub const fn id(&self) -> i64 {
        self.id
    }

    /// Return the represented group identifier.
    #[must_use]
    pub const fn group_id(&self) -> i64 {
        self.group_id
    }

    /// Return the representative model identifier.
    #[must_use]
    pub const fn model_id(&self) -> i64 {
        self.model_id
    }

    /// Return the declared representative-selection criterion.
    #[must_use]
    pub fn selection_criteria(&self) -> &str {
        &self.selection_criteria
    }
}

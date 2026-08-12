//! Explicit ModelCIF quality-metric definitions and values.

/// One named quality metric. Its mode determines which value category may use it.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct QaMetric {
    pub(super) id: i64,
    pub(super) name: String,
    pub(super) description: Option<String>,
    pub(super) metric_type: String,
    pub(super) mode: String,
    pub(super) software_group_id: Option<i64>,
    pub(super) data_id: Option<i64>,
}

impl QaMetric {
    /// Return the stable metric identifier.
    #[must_use]
    pub const fn id(&self) -> i64 {
        self.id
    }
    /// Return the deposited metric name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Return the optional metric description.
    #[must_use]
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }
    /// Return the ModelCIF metric-type enumeration value.
    #[must_use]
    pub fn metric_type(&self) -> &str {
        &self.metric_type
    }
    /// Return the metric mode, such as `global`, `local`, or `local-pairwise`.
    #[must_use]
    pub fn mode(&self) -> &str {
        &self.mode
    }
    /// Return the producing software-group identifier, when present.
    #[must_use]
    pub const fn software_group_id(&self) -> Option<i64> {
        self.software_group_id
    }
    /// Return the data record containing this metric, when present.
    #[must_use]
    pub const fn data_id(&self) -> Option<i64> {
        self.data_id
    }
}

/// A label-space residue identity used by local and pairwise QA values.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct ResidueSite {
    pub(super) asym_id: String,
    pub(super) sequence_id: i64,
    pub(super) component_id: String,
}

impl ResidueSite {
    /// Return the label-space asymmetric-unit identifier.
    #[must_use]
    pub fn asym_id(&self) -> &str {
        &self.asym_id
    }
    /// Return the label-space polymer sequence identifier.
    #[must_use]
    pub const fn sequence_id(&self) -> i64 {
        self.sequence_id
    }
    /// Return the label-space component identifier.
    #[must_use]
    pub fn component_id(&self) -> &str {
        &self.component_id
    }
}

/// One global, local, or residue-pair ModelCIF QA value.
#[derive(Clone, Debug, PartialEq)]
pub enum QaValue {
    /// A model-level metric value.
    Global {
        /// Stable value-row identifier.
        ordinal_id: i64,
        /// Referenced prediction-model identifier.
        model_id: i64,
        /// Referenced metric-definition identifier.
        metric_id: i64,
        /// Finite numeric metric value.
        value: f64,
    },
    /// A residue-level metric value, such as pLDDT.
    Local {
        /// Stable value-row identifier.
        ordinal_id: i64,
        /// Referenced prediction-model identifier.
        model_id: i64,
        /// Label-space residue identity.
        site: ResidueSite,
        /// Referenced metric-definition identifier.
        metric_id: i64,
        /// Finite numeric metric value.
        value: f64,
    },
    /// A metric value over an ordered residue pair, such as PAE.
    Pairwise {
        /// Stable value-row identifier.
        ordinal_id: i64,
        /// Referenced prediction-model identifier.
        model_id: i64,
        /// First label-space residue identity.
        first: ResidueSite,
        /// Second label-space residue identity.
        second: ResidueSite,
        /// Referenced metric-definition identifier.
        metric_id: i64,
        /// Finite numeric metric value.
        value: f64,
    },
}

impl QaValue {
    /// Return the referenced model identifier.
    #[must_use]
    pub const fn model_id(&self) -> i64 {
        match self {
            Self::Global { model_id, .. }
            | Self::Local { model_id, .. }
            | Self::Pairwise { model_id, .. } => *model_id,
        }
    }

    /// Return the referenced metric definition identifier.
    #[must_use]
    pub const fn metric_id(&self) -> i64 {
        match self {
            Self::Global { metric_id, .. }
            | Self::Local { metric_id, .. }
            | Self::Pairwise { metric_id, .. } => *metric_id,
        }
    }

    /// Return the finite numeric metric value.
    #[must_use]
    pub const fn value(&self) -> f64 {
        match self {
            Self::Global { value, .. }
            | Self::Local { value, .. }
            | Self::Pairwise { value, .. } => *value,
        }
    }
}

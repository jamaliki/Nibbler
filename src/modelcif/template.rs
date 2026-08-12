//! Structural-template, mapping, and alignment records.

/// One structural template used by a comparative model.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Template {
    pub(super) ordinal_id: i64,
    pub(super) id: i64,
    pub(super) origin: String,
    pub(super) entity_type: String,
    pub(super) data_id: i64,
    pub(super) target_asym_id: String,
    pub(super) auth_asym_id: String,
    pub(super) label_asym_id: Option<String>,
    pub(super) label_entity_id: Option<String>,
    pub(super) model_number: i64,
    pub(super) transform_id: i64,
    pub(super) name: Option<String>,
}

impl Template {
    /// Return the stable template identifier.
    #[must_use]
    pub const fn id(&self) -> i64 {
        self.id
    }
    /// Return the ModelCIF template-origin enumeration value.
    #[must_use]
    pub fn origin(&self) -> &str {
        &self.origin
    }
    /// Return the template entity-type enumeration value.
    #[must_use]
    pub fn entity_type(&self) -> &str {
        &self.entity_type
    }
    /// Return the data record representing this template.
    #[must_use]
    pub const fn data_id(&self) -> i64 {
        self.data_id
    }
    /// Return the target asymmetric-unit identifier mapped to this template.
    #[must_use]
    pub fn target_asym_id(&self) -> &str {
        &self.target_asym_id
    }
    /// Return the template's author asymmetric-unit identifier.
    #[must_use]
    pub fn auth_asym_id(&self) -> &str {
        &self.auth_asym_id
    }
    /// Return the template coordinate model number.
    #[must_use]
    pub const fn model_number(&self) -> i64 {
        self.model_number
    }
}

/// One contiguous polymer segment in a template.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TemplateSegment {
    pub(super) id: i64,
    pub(super) template_id: i64,
    pub(super) sequence_begin: Option<i64>,
    pub(super) sequence_end: Option<i64>,
}

impl TemplateSegment {
    /// Return the stable template-segment identifier.
    #[must_use]
    pub const fn id(&self) -> i64 {
        self.id
    }
    /// Return the containing template identifier.
    #[must_use]
    pub const fn template_id(&self) -> i64 {
        self.template_id
    }
    /// Return the inclusive template sequence range.
    #[must_use]
    pub const fn sequence_range(&self) -> (Option<i64>, Option<i64>) {
        (self.sequence_begin, self.sequence_end)
    }
}

/// One target-to-template segment mapping.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TemplateMapping {
    pub(super) id: i64,
    pub(super) template_segment_id: i64,
    pub(super) target_asym_id: String,
    pub(super) target_begin: Option<i64>,
    pub(super) target_end: Option<i64>,
}

impl TemplateMapping {
    /// Return the stable mapping identifier.
    #[must_use]
    pub const fn id(&self) -> i64 {
        self.id
    }
    /// Return the referenced template-segment identifier.
    #[must_use]
    pub const fn template_segment_id(&self) -> i64 {
        self.template_segment_id
    }
    /// Return the mapped target asymmetric-unit identifier.
    #[must_use]
    pub fn target_asym_id(&self) -> &str {
        &self.target_asym_id
    }
    /// Return the inclusive target sequence range.
    #[must_use]
    pub const fn target_range(&self) -> (Option<i64>, Option<i64>) {
        (self.target_begin, self.target_end)
    }
}

/// One target-template alignment declaration.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Alignment {
    pub(super) id: i64,
    pub(super) data_id: i64,
    pub(super) software_group_id: Option<i64>,
    pub(super) length: Option<i64>,
    pub(super) alignment_type: Option<String>,
    pub(super) mode: Option<String>,
}

impl Alignment {
    /// Return the stable alignment identifier.
    #[must_use]
    pub const fn id(&self) -> i64 {
        self.id
    }
    /// Return the data record representing this alignment.
    #[must_use]
    pub const fn data_id(&self) -> i64 {
        self.data_id
    }
    /// Return the producing software-group identifier, when present.
    #[must_use]
    pub const fn software_group_id(&self) -> Option<i64> {
        self.software_group_id
    }
    /// Return the deposited alignment length, when present.
    #[must_use]
    pub const fn length(&self) -> Option<i64> {
        self.length
    }
    /// Return the alignment-type enumeration value, when present.
    #[must_use]
    pub fn alignment_type(&self) -> Option<&str> {
        self.alignment_type.as_deref()
    }
    /// Return the alignment-mode enumeration value, when present.
    #[must_use]
    pub fn mode(&self) -> Option<&str> {
        self.mode.as_deref()
    }
}

/// One target/template participant and score record for an alignment.
#[derive(Clone, Debug, PartialEq)]
pub struct AlignmentDetail {
    pub(super) ordinal_id: i64,
    pub(super) alignment_id: i64,
    pub(super) template_segment_id: i64,
    pub(super) target_asym_id: String,
    pub(super) score_type: Option<String>,
    pub(super) score_value: Option<f64>,
    pub(super) sequence_identity: Option<f64>,
}

impl AlignmentDetail {
    /// Return the containing alignment identifier.
    #[must_use]
    pub const fn alignment_id(&self) -> i64 {
        self.alignment_id
    }
    /// Return the referenced template-segment identifier.
    #[must_use]
    pub const fn template_segment_id(&self) -> i64 {
        self.template_segment_id
    }
    /// Return the participating target asymmetric-unit identifier.
    #[must_use]
    pub fn target_asym_id(&self) -> &str {
        &self.target_asym_id
    }
}

/// One target or template sequence participating in an alignment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentSequence {
    pub(super) ordinal_id: i64,
    pub(super) alignment_id: i64,
    pub(super) target_template_flag: String,
    pub(super) sequence: String,
}

impl AlignmentSequence {
    /// Return the containing alignment identifier.
    #[must_use]
    pub const fn alignment_id(&self) -> i64 {
        self.alignment_id
    }
    /// Return whether this is the target or template sequence.
    #[must_use]
    pub fn target_template_flag(&self) -> &str {
        &self.target_template_flag
    }
    /// Return the deposited aligned sequence, including gaps.
    #[must_use]
    pub fn sequence(&self) -> &str {
        &self.sequence
    }
}

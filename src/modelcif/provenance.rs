//! Modeling provenance, data-flow, and associated-file records.

/// One membership edge in a protocol data group.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DataGroupMember {
    pub(super) ordinal_id: i64,
    pub(super) group_id: i64,
    pub(super) data_id: i64,
}

impl DataGroupMember {
    /// Return the stable membership-row identifier.
    #[must_use]
    pub const fn ordinal_id(&self) -> i64 {
        self.ordinal_id
    }
    /// Return the containing data-group identifier.
    #[must_use]
    pub const fn group_id(&self) -> i64 {
        self.group_id
    }
    /// Return the member data identifier.
    #[must_use]
    pub const fn data_id(&self) -> i64 {
        self.data_id
    }
}

/// One software package named by the prediction provenance.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Software {
    pub(super) id: i64,
    pub(super) name: String,
    pub(super) classification: String,
    pub(super) version: Option<String>,
    pub(super) description: Option<String>,
    pub(super) location: Option<String>,
}

impl Software {
    /// Return the stable software identifier.
    #[must_use]
    pub const fn id(&self) -> i64 {
        self.id
    }
    /// Return the software name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Return the ModelCIF software-classification value.
    #[must_use]
    pub fn classification(&self) -> &str {
        &self.classification
    }
    /// Return the deposited software version, when present.
    #[must_use]
    pub fn version(&self) -> Option<&str> {
        self.version.as_deref()
    }
    /// Return the optional software description.
    #[must_use]
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }
    /// Return the optional software location.
    #[must_use]
    pub fn location(&self) -> Option<&str> {
        self.location.as_deref()
    }
}

/// One software membership edge in a ModelCIF software group.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SoftwareGroupMember {
    pub(super) ordinal_id: i64,
    pub(super) group_id: i64,
    pub(super) software_id: i64,
    pub(super) parameter_group_id: Option<i64>,
}

impl SoftwareGroupMember {
    /// Return the stable membership-row identifier.
    #[must_use]
    pub const fn ordinal_id(&self) -> i64 {
        self.ordinal_id
    }
    /// Return the containing software-group identifier.
    #[must_use]
    pub const fn group_id(&self) -> i64 {
        self.group_id
    }
    /// Return the member software identifier.
    #[must_use]
    pub const fn software_id(&self) -> i64 {
        self.software_id
    }
    /// Return the associated parameter-group identifier, when present.
    #[must_use]
    pub const fn parameter_group_id(&self) -> Option<i64> {
        self.parameter_group_id
    }
}

/// One explicit step in the prediction protocol and its data dependencies.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ProtocolStep {
    pub(super) ordinal_id: i64,
    pub(super) protocol_id: i64,
    pub(super) step_id: i64,
    pub(super) method_type: String,
    pub(super) name: Option<String>,
    pub(super) details: Option<String>,
    pub(super) software_group_id: Option<i64>,
    pub(super) input_data_group_id: Option<i64>,
    pub(super) output_data_group_id: Option<i64>,
}

impl ProtocolStep {
    /// Return the stable protocol-step row identifier.
    #[must_use]
    pub const fn ordinal_id(&self) -> i64 {
        self.ordinal_id
    }
    /// Return the containing protocol identifier.
    #[must_use]
    pub const fn protocol_id(&self) -> i64 {
        self.protocol_id
    }
    /// Return the step number within the protocol.
    #[must_use]
    pub const fn step_id(&self) -> i64 {
        self.step_id
    }
    /// Return the ModelCIF method-type enumeration value.
    #[must_use]
    pub fn method_type(&self) -> &str {
        &self.method_type
    }
    /// Return the optional step name.
    #[must_use]
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }
    /// Return optional free-text step details.
    #[must_use]
    pub fn details(&self) -> Option<&str> {
        self.details.as_deref()
    }
    /// Return the software group used by this step, when present.
    #[must_use]
    pub const fn software_group_id(&self) -> Option<i64> {
        self.software_group_id
    }
    /// Return the input data-group identifier, when present.
    #[must_use]
    pub const fn input_data_group_id(&self) -> Option<i64> {
        self.input_data_group_id
    }
    /// Return the output data-group identifier, when present.
    #[must_use]
    pub const fn output_data_group_id(&self) -> Option<i64> {
        self.output_data_group_id
    }
}

/// A file associated directly with the ModelCIF entry.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AssociatedFile {
    pub(super) id: i64,
    pub(super) entry_id: String,
    pub(super) file_url: String,
    pub(super) file_type: Option<String>,
    pub(super) file_format: Option<String>,
    pub(super) file_content: Option<String>,
    pub(super) details: Option<String>,
    pub(super) data_id: Option<i64>,
}

impl AssociatedFile {
    /// Return the stable associated-file identifier.
    #[must_use]
    pub const fn id(&self) -> i64 {
        self.id
    }
    /// Return the entry identifier owning this file.
    #[must_use]
    pub fn entry_id(&self) -> &str {
        &self.entry_id
    }
    /// Return the deposited file URL.
    #[must_use]
    pub fn file_url(&self) -> &str {
        &self.file_url
    }
    /// Return the ModelCIF file-type value, when present.
    #[must_use]
    pub fn file_type(&self) -> Option<&str> {
        self.file_type.as_deref()
    }
    /// Return the ModelCIF file-format value, when present.
    #[must_use]
    pub fn file_format(&self) -> Option<&str> {
        self.file_format.as_deref()
    }
    /// Return the ModelCIF file-content value, when present.
    #[must_use]
    pub fn file_content(&self) -> Option<&str> {
        self.file_content.as_deref()
    }
    /// Return optional free-text file details.
    #[must_use]
    pub fn details(&self) -> Option<&str> {
        self.details.as_deref()
    }
    /// Return the data record represented by this file, when present.
    #[must_use]
    pub const fn data_id(&self) -> Option<i64> {
        self.data_id
    }
}

/// One file contained in an associated archive.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArchiveMember {
    pub(super) id: i64,
    pub(super) archive_file_id: i64,
    pub(super) file_path: String,
    pub(super) file_format: Option<String>,
    pub(super) file_content: Option<String>,
    pub(super) description: Option<String>,
    pub(super) data_id: Option<i64>,
}

impl ArchiveMember {
    /// Return the stable archive-member identifier.
    #[must_use]
    pub const fn id(&self) -> i64 {
        self.id
    }
    /// Return the containing associated archive-file identifier.
    #[must_use]
    pub const fn archive_file_id(&self) -> i64 {
        self.archive_file_id
    }
    /// Return the member's path inside the archive.
    #[must_use]
    pub fn file_path(&self) -> &str {
        &self.file_path
    }
    /// Return the member file format, when present.
    #[must_use]
    pub fn file_format(&self) -> Option<&str> {
        self.file_format.as_deref()
    }
    /// Return the member content classification, when present.
    #[must_use]
    pub fn file_content(&self) -> Option<&str> {
        self.file_content.as_deref()
    }
    /// Return the optional member description.
    #[must_use]
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }
    /// Return the data record represented by this member, when present.
    #[must_use]
    pub const fn data_id(&self) -> Option<i64> {
        self.data_id
    }
}

//! Provenance and associated-file records and row decoding.

use crate::cif::CategoryView;
use crate::pdbx::SemanticError;
use crate::pdbx::fields::category_rows;

use super::fields::{optional_integer, optional_text, required_integer, required_text};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct DataGroupMember {
    pub ordinal_id: i64,
    pub group_id: i64,
    pub data_id: i64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Software {
    pub id: i64,
    pub name: String,
    pub classification: String,
    pub version: Option<String>,
    pub description: Option<String>,
    pub location: Option<String>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct SoftwareGroupMember {
    pub ordinal_id: i64,
    pub group_id: i64,
    pub software_id: i64,
    pub parameter_group_id: Option<i64>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ProtocolStep {
    pub ordinal_id: i64,
    pub protocol_id: i64,
    pub step_id: i64,
    pub method_type: String,
    pub name: Option<String>,
    pub details: Option<String>,
    pub software_group_id: Option<i64>,
    pub input_data_group_id: Option<i64>,
    pub output_data_group_id: Option<i64>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct AssociatedFile {
    pub id: i64,
    pub entry_id: String,
    pub file_url: String,
    pub file_type: Option<String>,
    pub file_format: Option<String>,
    pub file_content: Option<String>,
    pub details: Option<String>,
    pub data_id: Option<i64>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ArchiveMember {
    pub id: i64,
    pub archive_file_id: i64,
    pub file_path: String,
    pub file_format: Option<String>,
    pub file_content: Option<String>,
    pub description: Option<String>,
    pub data_id: Option<i64>,
}

pub(super) fn software(categories: CategoryView<'_>) -> Result<Vec<Software>, SemanticError> {
    category_rows(
        categories,
        "software",
        [
            "pdbx_ordinal",
            "name",
            "classification",
            "version",
            "description",
            "location",
        ],
    )
    .map(|row| {
        Ok(Software {
            id: required_integer(&row, "software", "pdbx_ordinal", 0)?,
            name: required_text(&row, "software", "name", 1)?,
            classification: required_text(&row, "software", "classification", 2)?,
            version: optional_text(&row, 3),
            description: optional_text(&row, 4),
            location: optional_text(&row, 5),
        })
    })
    .collect()
}

pub(super) fn software_groups(
    categories: CategoryView<'_>,
) -> Result<Vec<SoftwareGroupMember>, SemanticError> {
    category_rows(
        categories,
        "ma_software_group",
        [
            "ordinal_id",
            "group_id",
            "software_id",
            "parameter_group_id",
        ],
    )
    .map(|row| {
        Ok(SoftwareGroupMember {
            ordinal_id: required_integer(&row, "ma_software_group", "ordinal_id", 0)?,
            group_id: required_integer(&row, "ma_software_group", "group_id", 1)?,
            software_id: required_integer(&row, "ma_software_group", "software_id", 2)?,
            parameter_group_id: optional_integer(
                &row,
                "ma_software_group",
                "parameter_group_id",
                3,
            )?,
        })
    })
    .collect()
}

pub(super) fn protocol_steps(
    categories: CategoryView<'_>,
) -> Result<Vec<ProtocolStep>, SemanticError> {
    category_rows(
        categories,
        "ma_protocol_step",
        [
            "ordinal_id",
            "protocol_id",
            "step_id",
            "method_type",
            "step_name",
            "details",
            "software_group_id",
            "input_data_group_id",
            "output_data_group_id",
        ],
    )
    .map(|row| {
        Ok(ProtocolStep {
            ordinal_id: required_integer(&row, "ma_protocol_step", "ordinal_id", 0)?,
            protocol_id: required_integer(&row, "ma_protocol_step", "protocol_id", 1)?,
            step_id: required_integer(&row, "ma_protocol_step", "step_id", 2)?,
            method_type: required_text(&row, "ma_protocol_step", "method_type", 3)?,
            name: optional_text(&row, 4),
            details: optional_text(&row, 5),
            software_group_id: optional_integer(&row, "ma_protocol_step", "software_group_id", 6)?,
            input_data_group_id: optional_integer(
                &row,
                "ma_protocol_step",
                "input_data_group_id",
                7,
            )?,
            output_data_group_id: optional_integer(
                &row,
                "ma_protocol_step",
                "output_data_group_id",
                8,
            )?,
        })
    })
    .collect()
}

pub(super) fn associated_files(
    categories: CategoryView<'_>,
) -> Result<Vec<AssociatedFile>, SemanticError> {
    category_rows(
        categories,
        "ma_entry_associated_files",
        [
            "id",
            "entry_id",
            "file_url",
            "file_type",
            "file_format",
            "file_content",
            "details",
            "data_id",
        ],
    )
    .map(|row| {
        Ok(AssociatedFile {
            id: required_integer(&row, "ma_entry_associated_files", "id", 0)?,
            entry_id: required_text(&row, "ma_entry_associated_files", "entry_id", 1)?,
            file_url: required_text(&row, "ma_entry_associated_files", "file_url", 2)?,
            file_type: optional_text(&row, 3),
            file_format: optional_text(&row, 4),
            file_content: optional_text(&row, 5),
            details: optional_text(&row, 6),
            data_id: optional_integer(&row, "ma_entry_associated_files", "data_id", 7)?,
        })
    })
    .collect()
}

pub(super) fn archive_members(
    categories: CategoryView<'_>,
) -> Result<Vec<ArchiveMember>, SemanticError> {
    category_rows(
        categories,
        "ma_associated_archive_file_details",
        [
            "id",
            "archive_file_id",
            "file_path",
            "file_format",
            "file_content",
            "description",
            "data_id",
        ],
    )
    .map(|row| {
        Ok(ArchiveMember {
            id: required_integer(&row, "ma_associated_archive_file_details", "id", 0)?,
            archive_file_id: required_integer(
                &row,
                "ma_associated_archive_file_details",
                "archive_file_id",
                1,
            )?,
            file_path: required_text(&row, "ma_associated_archive_file_details", "file_path", 2)?,
            file_format: optional_text(&row, 3),
            file_content: optional_text(&row, 4),
            description: optional_text(&row, 5),
            data_id: optional_integer(&row, "ma_associated_archive_file_details", "data_id", 6)?,
        })
    })
    .collect()
}

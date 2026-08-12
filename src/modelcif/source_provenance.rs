//! Software, protocol, and associated-file row decoding.

use crate::pdbx::SemanticError;
use crate::pdbx::category::CategoryIndex;

use super::fields::{optional_integer, optional_text, required_integer, required_text};
use super::provenance::{
    ArchiveMember, AssociatedFile, ProtocolStep, Software, SoftwareGroupMember,
};

pub(super) fn software(categories: &CategoryIndex<'_>) -> Result<Vec<Software>, SemanticError> {
    categories
        .rows("software")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(Software {
                id: required_integer(row, "software", "pdbx_ordinal", row_index)?,
                name: required_text(row, "software", "name", row_index)?,
                classification: required_text(row, "software", "classification", row_index)?,
                version: optional_text(row, "version"),
                description: optional_text(row, "description"),
                location: optional_text(row, "location"),
            })
        })
        .collect()
}

pub(super) fn software_groups(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<SoftwareGroupMember>, SemanticError> {
    categories
        .rows("ma_software_group")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(SoftwareGroupMember {
                ordinal_id: required_integer(row, "ma_software_group", "ordinal_id", row_index)?,
                group_id: required_integer(row, "ma_software_group", "group_id", row_index)?,
                software_id: required_integer(row, "ma_software_group", "software_id", row_index)?,
                parameter_group_id: optional_integer(
                    row,
                    "ma_software_group",
                    "parameter_group_id",
                    row_index,
                )?,
            })
        })
        .collect()
}

pub(super) fn protocol_steps(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<ProtocolStep>, SemanticError> {
    categories
        .rows("ma_protocol_step")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(ProtocolStep {
                ordinal_id: required_integer(row, "ma_protocol_step", "ordinal_id", row_index)?,
                protocol_id: required_integer(row, "ma_protocol_step", "protocol_id", row_index)?,
                step_id: required_integer(row, "ma_protocol_step", "step_id", row_index)?,
                method_type: required_text(row, "ma_protocol_step", "method_type", row_index)?,
                name: optional_text(row, "step_name"),
                details: optional_text(row, "details"),
                software_group_id: optional_integer(
                    row,
                    "ma_protocol_step",
                    "software_group_id",
                    row_index,
                )?,
                input_data_group_id: optional_integer(
                    row,
                    "ma_protocol_step",
                    "input_data_group_id",
                    row_index,
                )?,
                output_data_group_id: optional_integer(
                    row,
                    "ma_protocol_step",
                    "output_data_group_id",
                    row_index,
                )?,
            })
        })
        .collect()
}

pub(super) fn associated_files(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<AssociatedFile>, SemanticError> {
    categories
        .rows("ma_entry_associated_files")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(AssociatedFile {
                id: required_integer(row, "ma_entry_associated_files", "id", row_index)?,
                entry_id: required_text(row, "ma_entry_associated_files", "entry_id", row_index)?,
                file_url: required_text(row, "ma_entry_associated_files", "file_url", row_index)?,
                file_type: optional_text(row, "file_type"),
                file_format: optional_text(row, "file_format"),
                file_content: optional_text(row, "file_content"),
                details: optional_text(row, "details"),
                data_id: optional_integer(row, "ma_entry_associated_files", "data_id", row_index)?,
            })
        })
        .collect()
}

pub(super) fn archive_members(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<ArchiveMember>, SemanticError> {
    categories
        .rows("ma_associated_archive_file_details")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(ArchiveMember {
                id: required_integer(row, "ma_associated_archive_file_details", "id", row_index)?,
                archive_file_id: required_integer(
                    row,
                    "ma_associated_archive_file_details",
                    "archive_file_id",
                    row_index,
                )?,
                file_path: required_text(
                    row,
                    "ma_associated_archive_file_details",
                    "file_path",
                    row_index,
                )?,
                file_format: optional_text(row, "file_format"),
                file_content: optional_text(row, "file_content"),
                description: optional_text(row, "description"),
                data_id: optional_integer(
                    row,
                    "ma_associated_archive_file_details",
                    "data_id",
                    row_index,
                )?,
            })
        })
        .collect()
}

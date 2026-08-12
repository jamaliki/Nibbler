//! Software, protocol, and associated-file semantic validation.

use std::collections::BTreeSet;

use super::validation::Validator;

impl Validator<'_> {
    pub(super) fn provenance(&mut self) {
        self.require_nonempty(
            "MODELCIF_SOFTWARE_REQUIRED",
            "software",
            self.model.software.len(),
        );
        self.unique_i64(
            "MODELCIF_SOFTWARE_ID_DUPLICATE",
            "software",
            self.model.software.iter().map(|value| value.id),
        );
        self.require_nonempty(
            "MODELCIF_SOFTWARE_GROUP_REQUIRED",
            "ma_software_group",
            self.model.software_groups.len(),
        );
        for member in &self.model.software_groups {
            if let Some(software) = self.software.get(&member.software_id) {
                if software.version.is_none() {
                    self.error(
                        "MODELCIF_SOFTWARE_VERSION",
                        format!("software {:?} has no supplied version", software.name),
                        vec![format!("software={}", software.id)],
                    );
                }
            } else {
                self.error(
                    "MODELCIF_SOFTWARE_GROUP_MEMBER",
                    format!(
                        "software group {} references absent software {}",
                        member.group_id, member.software_id
                    ),
                    vec![format!("software_group={}", member.group_id)],
                );
            }
        }
        self.require_nonempty(
            "MODELCIF_PROTOCOL_REQUIRED",
            "ma_protocol_step",
            self.model.protocol_steps.len(),
        );
        self.unique_i64(
            "MODELCIF_PROTOCOL_ORDINAL_DUPLICATE",
            "protocol ordinal",
            self.model
                .protocol_steps
                .iter()
                .map(|value| value.ordinal_id),
        );
        for step in &self.model.protocol_steps {
            match step.software_group_id {
                Some(id) => self.reference(
                    "MODELCIF_PROTOCOL_SOFTWARE",
                    "software group",
                    id,
                    self.software_groups.contains(&id),
                ),
                None => self.error(
                    "MODELCIF_PROTOCOL_SOFTWARE",
                    "a strict protocol step requires an explicit software group",
                    vec![format!("protocol={}:{}", step.protocol_id, step.step_id)],
                ),
            }
            match step.input_data_group_id {
                Some(id) => self.reference(
                    "MODELCIF_PROTOCOL_INPUT",
                    "data group",
                    id,
                    self.data_groups.contains(&id),
                ),
                None => self.error(
                    "MODELCIF_PROTOCOL_INPUT",
                    "a strict protocol step requires an input data group",
                    vec![format!("protocol={}:{}", step.protocol_id, step.step_id)],
                ),
            }
            match step.output_data_group_id {
                Some(id) => self.reference(
                    "MODELCIF_PROTOCOL_OUTPUT",
                    "data group",
                    id,
                    self.data_groups.contains(&id),
                ),
                None => self.error(
                    "MODELCIF_PROTOCOL_OUTPUT",
                    "a strict protocol step requires an output data group",
                    vec![format!("protocol={}:{}", step.protocol_id, step.step_id)],
                ),
            }
        }
    }

    pub(super) fn files(&mut self) {
        let files = self
            .model
            .associated_files
            .iter()
            .map(|value| value.id)
            .collect::<BTreeSet<_>>();
        for file in &self.model.associated_files {
            if !file
                .entry_id
                .eq_ignore_ascii_case(self.model.coordinates.entry_id())
            {
                self.error(
                    "MODELCIF_ASSOCIATED_ENTRY",
                    format!(
                        "associated file {} names entry {:?}, expected {:?}",
                        file.id,
                        file.entry_id,
                        self.model.coordinates.entry_id()
                    ),
                    vec![format!("file={}", file.id)],
                );
            }
            if let Some(data_id) = file.data_id {
                self.reference(
                    "MODELCIF_ASSOCIATED_DATA",
                    "data",
                    data_id,
                    self.data.contains(&data_id),
                );
            }
        }
        for member in &self.model.archive_members {
            self.reference(
                "MODELCIF_ARCHIVE_FILE",
                "associated file",
                member.archive_file_id,
                files.contains(&member.archive_file_id),
            );
            if let Some(data_id) = member.data_id {
                self.reference(
                    "MODELCIF_ARCHIVE_DATA",
                    "data",
                    data_id,
                    self.data.contains(&data_id),
                );
            }
        }
    }
}

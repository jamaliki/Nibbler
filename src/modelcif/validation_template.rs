//! Template, mapping, and alignment semantic validation.

use std::collections::BTreeSet;

use crate::pdbx::fields::case_key;

use super::validation::Validator;

impl Validator<'_> {
    pub(super) fn templates(&mut self) {
        let templates = self
            .model
            .templates
            .iter()
            .map(|value| value.id)
            .collect::<BTreeSet<_>>();
        let segments = self
            .model
            .template_segments
            .iter()
            .map(|value| value.id)
            .collect::<BTreeSet<_>>();
        let alignments = self
            .model
            .alignments
            .iter()
            .map(|value| value.id)
            .collect::<BTreeSet<_>>();
        for template in &self.model.templates {
            self.reference(
                "MODELCIF_TEMPLATE_DATA",
                "data",
                template.data_id,
                self.data.contains(&template.data_id),
            );
            if !self
                .target_instances
                .contains(&case_key(&template.target_asym_id))
            {
                self.error(
                    "MODELCIF_TEMPLATE_TARGET",
                    format!(
                        "template {} references absent target instance {:?}",
                        template.id, template.target_asym_id
                    ),
                    vec![format!("template={}", template.id)],
                );
            }
        }
        for segment in &self.model.template_segments {
            self.reference(
                "MODELCIF_TEMPLATE_SEGMENT",
                "template",
                segment.template_id,
                templates.contains(&segment.template_id),
            );
            if matches!((segment.sequence_begin, segment.sequence_end), (Some(begin), Some(end)) if begin > end)
            {
                self.error(
                    "MODELCIF_TEMPLATE_SEGMENT_RANGE",
                    format!("template segment {} begins after it ends", segment.id),
                    vec![format!("segment={}", segment.id)],
                );
            }
        }
        for mapping in &self.model.template_mappings {
            self.reference(
                "MODELCIF_TEMPLATE_MAPPING_SEGMENT",
                "template segment",
                mapping.template_segment_id,
                segments.contains(&mapping.template_segment_id),
            );
            if !self
                .target_instances
                .contains(&case_key(&mapping.target_asym_id))
            {
                self.error(
                    "MODELCIF_TEMPLATE_MAPPING_TARGET",
                    format!(
                        "mapping {} references absent target instance {:?}",
                        mapping.id, mapping.target_asym_id
                    ),
                    vec![format!("mapping={}", mapping.id)],
                );
            }
        }
        for alignment in &self.model.alignments {
            self.reference(
                "MODELCIF_ALIGNMENT_DATA",
                "data",
                alignment.data_id,
                self.data.contains(&alignment.data_id),
            );
            if let Some(group) = alignment.software_group_id {
                self.reference(
                    "MODELCIF_ALIGNMENT_SOFTWARE",
                    "software group",
                    group,
                    self.software_groups.contains(&group),
                );
            }
        }
        for detail in &self.model.alignment_details {
            self.reference(
                "MODELCIF_ALIGNMENT_DETAIL",
                "alignment",
                detail.alignment_id,
                alignments.contains(&detail.alignment_id),
            );
            self.reference(
                "MODELCIF_ALIGNMENT_SEGMENT",
                "template segment",
                detail.template_segment_id,
                segments.contains(&detail.template_segment_id),
            );
            if !self
                .target_instances
                .contains(&case_key(&detail.target_asym_id))
            {
                self.error(
                    "MODELCIF_ALIGNMENT_TARGET",
                    format!(
                        "alignment {} references absent target instance {:?}",
                        detail.alignment_id, detail.target_asym_id
                    ),
                    vec![format!("alignment={}", detail.alignment_id)],
                );
            }
        }
        for sequence in &self.model.alignment_sequences {
            self.reference(
                "MODELCIF_ALIGNMENT_SEQUENCE",
                "alignment",
                sequence.alignment_id,
                alignments.contains(&sequence.alignment_id),
            );
        }
    }
}

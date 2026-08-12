//! QA definition, mode, model, and residue-site validation.

use std::collections::BTreeSet;

use crate::pdbx::category::case_key;

use super::qa::{QaValue, ResidueSite};
use super::validation::Validator;

impl Validator<'_> {
    pub(super) fn qa(&mut self) {
        self.unique_i64(
            "MODELCIF_QA_METRIC_DUPLICATE",
            "QA metric",
            self.model.qa_metrics.iter().map(|value| value.id),
        );
        let sites = self
            .model
            .coordinates
            .atom_sites()
            .iter()
            .filter_map(|atom| {
                let (_, component_id, asym_id, _) = atom.label_ids();
                atom.label_seq_id().map(|sequence_id| {
                    (
                        atom.model_number(),
                        case_key(asym_id),
                        sequence_id,
                        case_key(component_id),
                    )
                })
            })
            .collect::<BTreeSet<_>>();
        for metric in &self.model.qa_metrics {
            if let Some(group) = metric.software_group_id {
                self.reference(
                    "MODELCIF_QA_SOFTWARE",
                    "software group",
                    group,
                    self.software_groups.contains(&group),
                );
            }
            if let Some(data_id) = metric.data_id {
                self.reference(
                    "MODELCIF_QA_DATA",
                    "data",
                    data_id,
                    self.data.contains(&data_id),
                );
            }
        }
        for value in &self.model.qa_values {
            self.reference(
                "MODELCIF_QA_MODEL",
                "model",
                value.model_id(),
                self.models.contains(&value.model_id()),
            );
            let Some(metric) = self.metrics.get(&value.metric_id()).copied() else {
                self.error(
                    "MODELCIF_QA_DEFINITION",
                    format!("QA value references absent metric {}", value.metric_id()),
                    vec![format!("metric={}", value.metric_id())],
                );
                continue;
            };
            let expected_mode = match value {
                QaValue::Global { .. } => "global",
                QaValue::Local { .. } => "local",
                QaValue::Pairwise { .. } => "local-pairwise",
            };
            if !metric.mode.eq_ignore_ascii_case(expected_mode) {
                self.error(
                    "MODELCIF_QA_MODE",
                    format!(
                        "metric {} has mode {:?} but is used as {expected_mode}",
                        metric.id, metric.mode
                    ),
                    vec![format!("metric={}", metric.id)],
                );
            }
            match value {
                QaValue::Local { model_id, site, .. } => self.qa_site(*model_id, site, &sites),
                QaValue::Pairwise {
                    model_id,
                    first,
                    second,
                    ..
                } => {
                    self.qa_site(*model_id, first, &sites);
                    self.qa_site(*model_id, second, &sites);
                }
                QaValue::Global { .. } => {}
            }
        }
    }

    fn qa_site(
        &mut self,
        model_id: i64,
        site: &ResidueSite,
        sites: &BTreeSet<(i64, String, i64, String)>,
    ) {
        let key = (
            model_id,
            case_key(&site.asym_id),
            site.sequence_id,
            case_key(&site.component_id),
        );
        if !sites.contains(&key) {
            self.error(
                "MODELCIF_QA_SITE",
                format!(
                    "QA site {}:{} {} is absent from model {model_id}",
                    site.asym_id, site.sequence_id, site.component_id
                ),
                vec![
                    format!("model={model_id}"),
                    format!("asym={}", site.asym_id),
                    format!("seq={}", site.sequence_id),
                ],
            );
        }
    }
}

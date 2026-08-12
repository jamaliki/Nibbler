//! Target and prediction-model relationship validation.

use std::collections::{BTreeMap, BTreeSet};

use crate::pdbx::category::case_key;

use super::validation::Validator;

impl Validator<'_> {
    pub(super) fn targets(&mut self) {
        self.require_nonempty(
            "MODELCIF_TARGET_REQUIRED",
            "ma_target_entity",
            self.model.targets.len(),
        );
        self.require_nonempty(
            "MODELCIF_TARGET_INSTANCE_REQUIRED",
            "ma_target_entity_instance",
            self.model.target_instances.len(),
        );
        let coordinate_entities = self
            .model
            .coordinates
            .entities()
            .iter()
            .map(|value| case_key(value.id()))
            .collect::<BTreeSet<_>>();
        let coordinate_asym = self
            .model
            .coordinates
            .asym_units()
            .iter()
            .map(|value| (case_key(value.id()), case_key(value.entity_id())))
            .collect::<BTreeMap<_, _>>();
        for target in &self.model.targets {
            if !coordinate_entities.contains(&case_key(&target.entity_id)) {
                self.error(
                    "MODELCIF_TARGET_ENTITY",
                    format!(
                        "target entity {:?} is absent from _entity",
                        target.entity_id
                    ),
                    vec![format!("entity={}", target.entity_id)],
                );
            }
            self.reference(
                "MODELCIF_TARGET_DATA",
                "data",
                target.data_id,
                self.data.contains(&target.data_id),
            );
        }
        for instance in &self.model.target_instances {
            match coordinate_asym.get(&case_key(&instance.asym_id)) {
                Some(entity_id) if *entity_id == case_key(&instance.entity_id) => {}
                Some(_) => self.error(
                    "MODELCIF_TARGET_INSTANCE_ENTITY",
                    format!(
                        "target instance {:?} disagrees with _struct_asym",
                        instance.asym_id
                    ),
                    vec![format!("asym={}", instance.asym_id)],
                ),
                None => self.error(
                    "MODELCIF_TARGET_INSTANCE",
                    format!(
                        "target instance {:?} is absent from _struct_asym",
                        instance.asym_id
                    ),
                    vec![format!("asym={}", instance.asym_id)],
                ),
            }
            if !self
                .target_entities
                .contains_key(&case_key(&instance.entity_id))
            {
                self.error(
                    "MODELCIF_TARGET_INSTANCE_TARGET",
                    format!(
                        "target instance {:?} references undeclared target entity {:?}",
                        instance.asym_id, instance.entity_id
                    ),
                    vec![format!("asym={}", instance.asym_id)],
                );
            }
        }
    }

    pub(super) fn models(&mut self) {
        self.require_nonempty(
            "MODELCIF_MODEL_REQUIRED",
            "ma_model_list",
            self.model.models.len(),
        );
        self.unique_i64(
            "MODELCIF_MODEL_ID_DUPLICATE",
            "model",
            self.model.models.iter().map(|value| value.id),
        );
        self.require_nonempty(
            "MODELCIF_MODEL_GROUP_REQUIRED",
            "ma_model_group",
            self.model.model_groups.len(),
        );
        self.unique_i64(
            "MODELCIF_MODEL_GROUP_DUPLICATE",
            "model group",
            self.model.model_groups.iter().map(|value| value.id),
        );
        let coordinate_models = self
            .model
            .coordinates
            .atom_sites()
            .iter()
            .map(|atom| atom.model_number())
            .collect::<BTreeSet<_>>();
        for prediction in &self.model.models {
            self.reference(
                "MODELCIF_MODEL_DATA",
                "data",
                prediction.data_id,
                self.data.contains(&prediction.data_id),
            );
            if !coordinate_models.contains(&prediction.id) {
                self.error(
                    "MODELCIF_MODEL_COORDINATES",
                    format!("model {} has no atom-site coordinates", prediction.id),
                    vec![format!("model={}", prediction.id)],
                );
            }
        }
        for model_id in coordinate_models {
            self.reference(
                "MODELCIF_ATOM_MODEL",
                "model",
                model_id,
                self.models.contains(&model_id),
            );
        }
        let mut linked_models = BTreeSet::new();
        for link in &self.model.model_group_links {
            self.reference(
                "MODELCIF_MODEL_GROUP_MODEL",
                "model",
                link.model_id,
                self.models.contains(&link.model_id),
            );
            self.reference(
                "MODELCIF_MODEL_GROUP_LINK",
                "model group",
                link.group_id,
                self.groups.contains(&link.group_id),
            );
            linked_models.insert(link.model_id);
        }
        for prediction in &self.model.models {
            if !linked_models.contains(&prediction.id) {
                self.error(
                    "MODELCIF_MODEL_UNGROUPED",
                    format!("model {} has no explicit group membership", prediction.id),
                    vec![format!("model={}", prediction.id)],
                );
            }
        }
        for representative in &self.model.representatives {
            self.reference(
                "MODELCIF_REPRESENTATIVE_MODEL",
                "model",
                representative.model_id,
                self.models.contains(&representative.model_id),
            );
            self.reference(
                "MODELCIF_REPRESENTATIVE_GROUP",
                "model group",
                representative.group_id,
                self.groups.contains(&representative.group_id),
            );
            if !self.model.model_group_links.iter().any(|link| {
                link.model_id == representative.model_id && link.group_id == representative.group_id
            }) {
                self.error(
                    "MODELCIF_REPRESENTATIVE_MEMBERSHIP",
                    format!(
                        "representative model {} is not in group {}",
                        representative.model_id, representative.group_id
                    ),
                    vec![format!("representative={}", representative.id)],
                );
            }
        }
    }
}

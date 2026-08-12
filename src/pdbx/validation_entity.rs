//! Entity, asym-unit, sequence, and instance-scheme checks.

use std::collections::{BTreeMap, BTreeSet};

use super::model::{Entity, EntityKind};
use super::validation::{Validator, equal, fold};

impl Validator<'_> {
    pub(super) fn entities(&mut self) {
        for entity in self.model.entities() {
            match entity.kind() {
                EntityKind::Polymer => self.polymer(entity),
                EntityKind::NonPolymer | EntityKind::Water => self.nonpolymer(entity),
                EntityKind::Branched => self.branched(entity),
            }
        }
    }

    fn polymer(&mut self, entity: &Entity) {
        if entity.polymer_type().is_none() || entity.sequence().is_empty() {
            self.error(
                "PDBX_POLYMER_DEFINITION",
                format!(
                    "polymer entity {:?} requires a type and complete sequence",
                    entity.id()
                ),
                entity_context(entity),
            );
        }
        let sequence = entity
            .sequence()
            .iter()
            .map(|monomer| (monomer.number(), fold(monomer.component_id())))
            .collect::<BTreeMap<_, _>>();
        for asym_id in self.asym_ids_for_entity(entity.id()) {
            for (number, component_id) in &sequence {
                let count = self
                    .model
                    .poly_scheme
                    .iter()
                    .filter(|row| {
                        equal(&row.asym_id, &asym_id)
                            && equal(&row.entity_id, entity.id())
                            && row.seq_id == *number
                            && fold(&row.component_id) == *component_id
                    })
                    .count();
                if count != 1 {
                    self.error(
                        "PDBX_POLYMER_SCHEME",
                        format!(
                            "asym {asym_id:?} requires exactly one scheme row for entity {:?} sequence {number}",
                            entity.id()
                        ),
                        vec![format!("asym={asym_id}"), format!("seq_id={number}")],
                    );
                }
            }
        }
    }

    fn nonpolymer(&mut self, entity: &Entity) {
        let Some(component_id) = entity.component_id() else {
            self.error(
                "PDBX_NONPOLY_COMPONENT",
                format!(
                    "{} entity {:?} requires exactly one component",
                    entity.kind().as_str(),
                    entity.id()
                ),
                entity_context(entity),
            );
            return;
        };
        for asym_id in self.asym_ids_for_entity(entity.id()) {
            let count = self
                .model
                .nonpoly_scheme
                .iter()
                .filter(|row| {
                    equal(&row.asym_id, &asym_id)
                        && equal(&row.entity_id, entity.id())
                        && equal(&row.component_id, component_id)
                        && !row.auth_seq_id.is_empty()
                })
                .count();
            if count != 1 {
                self.error(
                    "PDBX_NONPOLY_SCHEME",
                    format!("asym {asym_id:?} requires exactly one non-polymer scheme row"),
                    vec![
                        format!("asym={asym_id}"),
                        format!("component={component_id}"),
                    ],
                );
            }
        }
    }

    fn branched(&mut self, entity: &Entity) {
        if entity.branch_nodes().is_empty() {
            self.error(
                "PDBX_BRANCH_EMPTY",
                format!("branched entity {:?} requires component nodes", entity.id()),
                entity_context(entity),
            );
        }
        let nodes = entity
            .branch_nodes()
            .iter()
            .map(|node| node.number())
            .collect::<BTreeSet<_>>();
        if nodes.len() != entity.branch_nodes().len() {
            self.error(
                "PDBX_BRANCH_NODE_DUPLICATE",
                format!(
                    "branched entity {:?} has duplicate node numbers",
                    entity.id()
                ),
                entity_context(entity),
            );
        }
        for link in entity.branch_links() {
            let (first, second) = link.nodes();
            if !nodes.contains(&first) || !nodes.contains(&second) {
                self.error(
                    "PDBX_BRANCH_LINK_ENDPOINT",
                    format!(
                        "branched entity {:?} link {first}-{second} references an absent node",
                        entity.id()
                    ),
                    entity_context(entity),
                );
            }
        }
        for asym_id in self.asym_ids_for_entity(entity.id()) {
            for node in entity.branch_nodes() {
                let count = self
                    .model
                    .branch_scheme
                    .iter()
                    .filter(|row| {
                        equal(&row.asym_id, &asym_id)
                            && equal(&row.entity_id, entity.id())
                            && row.number == node.number()
                            && equal(&row.component_id, node.component_id())
                    })
                    .count();
                if count != 1 {
                    self.error(
                        "PDBX_BRANCH_SCHEME",
                        format!(
                            "asym {asym_id:?} requires exactly one branch scheme row for node {}",
                            node.number()
                        ),
                        vec![format!("asym={asym_id}"), format!("node={}", node.number())],
                    );
                }
            }
        }
    }

    pub(super) fn asym_units(&mut self) {
        for asym in self.model.asym_units() {
            if !self.entities.contains_key(&fold(asym.entity_id())) {
                self.error(
                    "PDBX_ASYM_ENTITY",
                    format!(
                        "asym {:?} references absent entity {:?}",
                        asym.id(),
                        asym.entity_id()
                    ),
                    vec![format!("asym={}", asym.id())],
                );
            }
        }
    }

    fn asym_ids_for_entity(&self, entity_id: &str) -> Vec<String> {
        self.model
            .asym_units()
            .iter()
            .filter(|asym| equal(asym.entity_id(), entity_id))
            .map(|asym| asym.id().to_owned())
            .collect()
    }
}

fn entity_context(entity: &Entity) -> Vec<String> {
    vec![format!("entity={}", entity.id())]
}

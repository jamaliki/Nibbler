//! Atom-site identity, numeric, chemistry, and connection checks.

use super::category::case_key;
use super::model::{AtomSite, ConnectionEndpoint, Entity, EntityKind};
use super::validation::{Validator, atom_context, equal};

impl Validator<'_> {
    pub(super) fn atom_sites(&mut self) {
        for atom in self.model.atom_sites() {
            self.atom_site(atom);
        }
    }

    fn atom_site(&mut self, atom: &AtomSite) {
        let (atom_id, component_id, asym_id, entity_id) = atom.label_ids();
        let Some(asym) = self.asym_units.get(&case_key(asym_id)).copied() else {
            self.error(
                "PDBX_ATOM_ASYM",
                format!("atom {:?} references absent asym {asym_id:?}", atom.id()),
                atom_context(atom),
            );
            return;
        };
        if !equal(asym.entity_id(), entity_id) {
            self.error(
                "PDBX_ATOM_ENTITY",
                format!(
                    "atom {:?} entity {entity_id:?} disagrees with asym {:?} entity {:?}",
                    atom.id(),
                    asym.id(),
                    asym.entity_id()
                ),
                atom_context(atom),
            );
            return;
        }
        let Some(entity) = self.entities.get(&case_key(entity_id)).copied() else {
            return;
        };
        self.atom_entity_identity(atom, entity, component_id);
        self.numeric_fields(atom);
        self.component_atom(atom, atom_id, component_id);
    }

    fn numeric_fields(&mut self, atom: &AtomSite) {
        let coordinates = atom.coordinates();
        let (occupancy, b_iso) = atom.occupancy_and_b_iso();
        if coordinates.iter().any(|value| !value.is_finite())
            || !occupancy.is_finite()
            || b_iso.is_some_and(|value| !value.is_finite())
        {
            self.error(
                "PDBX_ATOM_FINITE",
                format!("atom {:?} contains a non-finite numeric value", atom.id()),
                atom_context(atom),
            );
        }
        if !(0.0..=1.0).contains(&occupancy) {
            self.warning(
                "PDBX_OCCUPANCY_RANGE",
                format!(
                    "atom {:?} occupancy {occupancy} is outside [0, 1]",
                    atom.id()
                ),
                atom_context(atom),
            );
        }
        if atom.model_number() <= 0 {
            self.error(
                "PDBX_MODEL_NUMBER",
                format!(
                    "atom {:?} has non-positive model number {}",
                    atom.id(),
                    atom.model_number()
                ),
                atom_context(atom),
            );
        }
    }

    fn component_atom(&mut self, atom: &AtomSite, atom_id: &str, component_id: &str) {
        let Some(component) = self.components.get(&case_key(component_id)).copied() else {
            return;
        };
        if component.atoms().is_empty() {
            return;
        }
        match component
            .atoms()
            .iter()
            .find(|candidate| equal(candidate.atom_id(), atom_id))
        {
            Some(candidate) if !equal(candidate.element(), atom.element()) => self.error(
                "PDBX_ATOM_ELEMENT",
                format!(
                    "atom {:?} element {:?} disagrees with component atom {:?}",
                    atom.id(),
                    atom.element(),
                    candidate.element()
                ),
                atom_context(atom),
            ),
            Some(candidate)
                if candidate.formal_charge().is_some()
                    && atom.formal_charge().is_some()
                    && candidate.formal_charge() != atom.formal_charge() =>
            {
                self.error(
                    "PDBX_ATOM_CHARGE",
                    format!(
                        "atom {:?} formal charge disagrees with its component definition",
                        atom.id()
                    ),
                    atom_context(atom),
                );
            }
            Some(_) => {}
            None => self.error(
                "PDBX_COMPONENT_ATOM_MISSING",
                format!(
                    "atom {:?} does not resolve in component {component_id:?}",
                    atom.id()
                ),
                atom_context(atom),
            ),
        }
    }

    fn atom_entity_identity(&mut self, atom: &AtomSite, entity: &Entity, component_id: &str) {
        match entity.kind() {
            EntityKind::Polymer => match atom.label_seq_id() {
                Some(number)
                    if entity.sequence().iter().any(|monomer| {
                        monomer.number() == number && equal(monomer.component_id(), component_id)
                    }) => {}
                _ => self.error(
                    "PDBX_ATOM_POLYMER_IDENTITY",
                    format!(
                        "atom {:?} does not match entity {:?} sequence identity",
                        atom.id(),
                        entity.id()
                    ),
                    atom_context(atom),
                ),
            },
            EntityKind::NonPolymer | EntityKind::Water => {
                if atom.label_seq_id().is_some()
                    || entity
                        .component_id()
                        .is_none_or(|expected| !equal(expected, component_id))
                {
                    self.error(
                        "PDBX_ATOM_NONPOLY_IDENTITY",
                        format!(
                            "atom {:?} must use its entity component and label_seq_id=.",
                            atom.id()
                        ),
                        atom_context(atom),
                    );
                }
            }
            EntityKind::Branched => {
                let (auth_seq_id, _, _, _) = atom.author_ids();
                if atom.label_seq_id().is_some()
                    || !self.model.branch_scheme.iter().any(|row| {
                        equal(&row.asym_id, atom.label_ids().2)
                            && equal(&row.entity_id, entity.id())
                            && equal(&row.component_id, component_id)
                            && equal(&row.auth_seq_id, auth_seq_id)
                    })
                {
                    self.error(
                        "PDBX_ATOM_BRANCH_IDENTITY",
                        format!(
                            "atom {:?} must use label_seq_id=. and resolve through the branch scheme",
                            atom.id()
                        ),
                        atom_context(atom),
                    );
                }
            }
        }
    }

    pub(super) fn connections(&mut self) {
        for connection in self.model.connections() {
            let (id, _) = connection.identity();
            let (first, second) = connection.endpoints();
            self.connection_endpoint(id, "first", first);
            self.connection_endpoint(id, "second", second);
        }
    }

    fn connection_endpoint(
        &mut self,
        connection_id: &str,
        side: &str,
        endpoint: &ConnectionEndpoint,
    ) {
        let count = self
            .model
            .atom_sites()
            .iter()
            .filter(|atom| {
                let (atom_id, component_id, asym_id, _) = atom.label_ids();
                equal(asym_id, endpoint.asym_id())
                    && equal(atom_id, endpoint.atom_id())
                    && endpoint
                        .component_id()
                        .is_none_or(|expected| equal(expected, component_id))
                    && endpoint
                        .sequence_id()
                        .is_none_or(|expected| atom.label_seq_id() == Some(expected))
            })
            .count();
        if count == 0 {
            self.error(
                "PDBX_CONNECTION_ENDPOINT",
                format!(
                    "connection {connection_id:?} {side} endpoint does not resolve to an atom site"
                ),
                vec![
                    format!("connection={connection_id}"),
                    format!("endpoint={side}"),
                ],
            );
        }
    }
}

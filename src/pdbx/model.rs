use crate::cif::CifDocument;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum EntityKind {
    Polymer,
    NonPolymer,
    Branched,
    Water,
}

impl EntityKind {
    pub(super) const fn as_str(self) -> &'static str {
        match self {
            Self::Polymer => "polymer",
            Self::NonPolymer => "non-polymer",
            Self::Branched => "branched",
            Self::Water => "water",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct PolymerMonomer {
    pub number: i64,
    pub component_id: String,
    pub heterogeneous: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct BranchedNode {
    pub number: i64,
    pub component_id: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct BranchedLink {
    pub first_node: i64,
    pub second_node: i64,
    pub first_atom_id: String,
    pub second_atom_id: String,
    pub order: Option<String>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Entity {
    pub id: String,
    pub kind: EntityKind,
    pub description: Option<String>,
    pub polymer_type: Option<String>,
    pub sequence: Vec<PolymerMonomer>,
    pub component_id: Option<String>,
    pub branch_nodes: Vec<BranchedNode>,
    pub branch_links: Vec<BranchedLink>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct AsymUnit {
    pub id: String,
    pub entity_id: String,
    pub details: Option<String>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum ComponentResolution {
    Embedded,
    MinimalRegistry,
    LocalCcd,
    Unresolved,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ComponentAtom {
    pub atom_id: String,
    pub element: String,
    pub formal_charge: Option<i32>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ComponentBond {
    pub first_atom_id: String,
    pub second_atom_id: String,
    pub order: String,
    pub aromatic: Option<bool>,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct ComponentDefinition {
    pub id: String,
    pub name: Option<String>,
    pub component_type: Option<String>,
    pub formula: Option<String>,
    pub formula_weight: Option<f64>,
    pub formal_charge: Option<i32>,
    pub atoms: Vec<ComponentAtom>,
    pub bonds: Vec<ComponentBond>,
    pub resolution: ComponentResolution,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct AtomSite {
    pub group_pdb: String,
    pub id: String,
    pub element: String,
    pub label_atom_id: String,
    pub label_alt_id: Option<String>,
    pub label_comp_id: String,
    pub label_asym_id: String,
    pub label_entity_id: String,
    pub label_seq_id: Option<i64>,
    pub insertion_code: Option<String>,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f64,
    pub b_iso: Option<f64>,
    pub formal_charge: Option<i32>,
    pub auth_seq_id: String,
    pub auth_comp_id: String,
    pub auth_asym_id: String,
    pub auth_atom_id: String,
    pub model_number: i64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ConnectionEndpoint {
    pub asym_id: String,
    pub component_id: Option<String>,
    pub sequence_id: Option<i64>,
    pub atom_id: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Connection {
    pub id: String,
    pub kind: String,
    pub first: ConnectionEndpoint,
    pub second: ConnectionEndpoint,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct PolySchemeRow {
    pub asym_id: String,
    pub entity_id: String,
    pub seq_id: i64,
    pub component_id: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct NonPolySchemeRow {
    pub asym_id: String,
    pub entity_id: String,
    pub component_id: String,
    pub auth_seq_id: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct BranchSchemeRow {
    pub asym_id: String,
    pub entity_id: String,
    pub number: i64,
    pub component_id: String,
    pub auth_seq_id: String,
}

/// An immutable, source-backed PDBx coordinate model.
#[derive(Clone, Debug, PartialEq)]
pub struct PdbxModel {
    pub(super) document: CifDocument,
    pub(super) entry_id: String,
    pub(super) audit_dictionary: Option<(String, String)>,
    pub(super) entities: Vec<Entity>,
    pub(super) asym_units: Vec<AsymUnit>,
    pub(super) components: Vec<ComponentDefinition>,
    pub(super) atom_sites: Vec<AtomSite>,
    pub(super) connections: Vec<Connection>,
    pub(super) poly_scheme: Vec<PolySchemeRow>,
    pub(super) nonpoly_scheme: Vec<NonPolySchemeRow>,
    pub(super) branch_scheme: Vec<BranchSchemeRow>,
}

impl PdbxModel {
    /// Return `_entry.id`.
    #[must_use]
    pub fn entry_id(&self) -> &str {
        &self.entry_id
    }

    /// Return the number of coordinate entities.
    #[must_use]
    pub fn entity_count(&self) -> usize {
        self.entities.len()
    }

    /// Return the number of asymmetric-unit instances.
    #[must_use]
    pub fn asym_unit_count(&self) -> usize {
        self.asym_units.len()
    }

    /// Return the number of component definitions and unresolved placeholders.
    #[must_use]
    pub fn component_count(&self) -> usize {
        self.components.len()
    }

    /// Return the number of coordinate atom sites.
    #[must_use]
    pub fn atom_site_count(&self) -> usize {
        self.atom_sites.len()
    }

    /// Return the number of explicit inter-site connections.
    #[must_use]
    pub fn connection_count(&self) -> usize {
        self.connections.len()
    }

    /// Iterate over canonical `_entity.type` values in source order.
    pub fn entity_kinds(&self) -> impl ExactSizeIterator<Item = &'static str> + '_ {
        self.entities.iter().map(|entity| entity.kind.as_str())
    }

    /// Iterate over resolved component identifiers in stable order.
    pub fn component_ids(&self) -> impl ExactSizeIterator<Item = &str> {
        self.components
            .iter()
            .map(|component| component.id.as_str())
    }

    pub(crate) const fn source_document(&self) -> &CifDocument {
        &self.document
    }

    pub(crate) fn coordinate_entity_ids(&self) -> impl Iterator<Item = &str> {
        self.entities.iter().map(|entity| entity.id.as_str())
    }

    pub(crate) fn asym_entity_ids(&self) -> impl Iterator<Item = (&str, &str)> {
        self.asym_units
            .iter()
            .map(|asym| (asym.id.as_str(), asym.entity_id.as_str()))
    }

    pub(crate) fn coordinate_model_numbers(&self) -> impl Iterator<Item = i64> + '_ {
        self.atom_sites.iter().map(|atom| atom.model_number)
    }

    pub(crate) fn polymer_sites(&self) -> impl Iterator<Item = (i64, &str, i64, &str)> {
        self.atom_sites.iter().filter_map(|atom| {
            atom.label_seq_id.map(|sequence_id| {
                (
                    atom.model_number,
                    atom.label_asym_id.as_str(),
                    sequence_id,
                    atom.label_comp_id.as_str(),
                )
            })
        })
    }
}

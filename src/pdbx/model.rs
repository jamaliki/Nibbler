use crate::cif::CifDocument;

/// The PDBx entity classification, kept distinct from asym-unit instances.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EntityKind {
    /// A covalently linked monomer sequence.
    Polymer,
    /// A ligand, ion, or other discrete non-polymer component.
    NonPolymer,
    /// A graph-shaped branched entity such as a glycan.
    Branched,
    /// A water entity, normally using component `HOH`.
    Water,
}

impl EntityKind {
    /// Return the canonical PDBx `_entity.type` value.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Polymer => "polymer",
            Self::NonPolymer => "non-polymer",
            Self::Branched => "branched",
            Self::Water => "water",
        }
    }
}

/// One position in a complete polymer sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PolymerMonomer {
    pub(super) number: i64,
    pub(super) component_id: String,
    pub(super) heterogeneous: bool,
}

impl PolymerMonomer {
    /// Return the one-based or dictionary-supplied sequence number.
    #[must_use]
    pub const fn number(&self) -> i64 {
        self.number
    }

    /// Return the chemical component identifier, including modifications such as MSE.
    #[must_use]
    pub fn component_id(&self) -> &str {
        &self.component_id
    }

    /// Return whether this position declares sequence microheterogeneity.
    #[must_use]
    pub const fn heterogeneous(&self) -> bool {
        self.heterogeneous
    }
}

/// One node in a branched-entity component graph.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BranchedNode {
    pub(super) number: i64,
    pub(super) component_id: String,
}

impl BranchedNode {
    /// Return the branch-local node number.
    #[must_use]
    pub const fn number(&self) -> i64 {
        self.number
    }

    /// Return the node's component identifier.
    #[must_use]
    pub fn component_id(&self) -> &str {
        &self.component_id
    }
}

/// One typed covalent edge in a branched entity.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BranchedLink {
    pub(super) first_node: i64,
    pub(super) second_node: i64,
    pub(super) first_atom_id: String,
    pub(super) second_atom_id: String,
    pub(super) order: Option<String>,
}

impl BranchedLink {
    /// Return the two branch-local node numbers.
    #[must_use]
    pub const fn nodes(&self) -> (i64, i64) {
        (self.first_node, self.second_node)
    }

    /// Return the atom identifiers forming the branch bond.
    #[must_use]
    pub fn atom_ids(&self) -> (&str, &str) {
        (&self.first_atom_id, &self.second_atom_id)
    }

    /// Return the dictionary bond-order value when supplied.
    #[must_use]
    pub fn order(&self) -> Option<&str> {
        self.order.as_deref()
    }
}

/// One unique chemical entity, separate from its coordinate instances.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Entity {
    pub(super) id: String,
    pub(super) kind: EntityKind,
    pub(super) description: Option<String>,
    pub(super) polymer_type: Option<String>,
    pub(super) sequence: Vec<PolymerMonomer>,
    pub(super) component_id: Option<String>,
    pub(super) branch_nodes: Vec<BranchedNode>,
    pub(super) branch_links: Vec<BranchedLink>,
}

impl Entity {
    /// Return the label-space entity identifier.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Return the semantic entity kind.
    #[must_use]
    pub const fn kind(&self) -> EntityKind {
        self.kind
    }

    /// Return the source description without using it for identity.
    #[must_use]
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }

    /// Return the PDBx polymer type for polymer entities.
    #[must_use]
    pub fn polymer_type(&self) -> Option<&str> {
        self.polymer_type.as_deref()
    }

    /// Return the complete polymer sequence in sequence-number order.
    #[must_use]
    pub fn sequence(&self) -> &[PolymerMonomer] {
        &self.sequence
    }

    /// Return the single component used by a non-polymer or water entity.
    #[must_use]
    pub fn component_id(&self) -> Option<&str> {
        self.component_id.as_deref()
    }

    /// Return the component nodes of a branched entity.
    #[must_use]
    pub fn branch_nodes(&self) -> &[BranchedNode] {
        &self.branch_nodes
    }

    /// Return the covalent edges of a branched entity.
    #[must_use]
    pub fn branch_links(&self) -> &[BranchedLink] {
        &self.branch_links
    }
}

/// One label-space instance of an entity.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AsymUnit {
    pub(super) id: String,
    pub(super) entity_id: String,
    pub(super) details: Option<String>,
}

impl AsymUnit {
    /// Return the label asym identifier.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Return the referenced label entity identifier.
    #[must_use]
    pub fn entity_id(&self) -> &str {
        &self.entity_id
    }

    /// Return optional source details.
    #[must_use]
    pub fn details(&self) -> Option<&str> {
        self.details.as_deref()
    }
}

/// Provenance of a component definition used by the semantic model.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ComponentResolution {
    /// Defined by `_chem_comp` and optional detail categories in the input.
    Embedded,
    /// Supplied by Nibbler's immutable minimal component registry.
    MinimalRegistry,
    /// Supplied by a caller-selected immutable local CCD cache.
    LocalCcd,
    /// Referenced by the model but not chemically defined.
    Unresolved,
}

/// One atom in a chemical component definition.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ComponentAtom {
    pub(super) atom_id: String,
    pub(super) element: String,
    pub(super) formal_charge: Option<i32>,
}

impl ComponentAtom {
    /// Return the component-local atom identifier.
    #[must_use]
    pub fn atom_id(&self) -> &str {
        &self.atom_id
    }

    /// Return the element symbol.
    #[must_use]
    pub fn element(&self) -> &str {
        &self.element
    }

    /// Return the formal charge when explicitly defined.
    #[must_use]
    pub const fn formal_charge(&self) -> Option<i32> {
        self.formal_charge
    }
}

/// One bond in a chemical component definition.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ComponentBond {
    pub(super) first_atom_id: String,
    pub(super) second_atom_id: String,
    pub(super) order: String,
    pub(super) aromatic: Option<bool>,
}

impl ComponentBond {
    /// Return the two component-local atom identifiers.
    #[must_use]
    pub fn atom_ids(&self) -> (&str, &str) {
        (&self.first_atom_id, &self.second_atom_id)
    }

    /// Return the dictionary bond-order code.
    #[must_use]
    pub fn order(&self) -> &str {
        &self.order
    }

    /// Return the aromatic flag when explicitly defined.
    #[must_use]
    pub const fn aromatic(&self) -> Option<bool> {
        self.aromatic
    }
}

/// One resolved or explicitly unresolved chemical component definition.
#[derive(Clone, Debug, PartialEq)]
pub struct ComponentDefinition {
    pub(super) id: String,
    pub(super) name: Option<String>,
    pub(super) component_type: Option<String>,
    pub(super) formula: Option<String>,
    pub(super) formula_weight: Option<f64>,
    pub(super) formal_charge: Option<i32>,
    pub(super) atoms: Vec<ComponentAtom>,
    pub(super) bonds: Vec<ComponentBond>,
    pub(super) resolution: ComponentResolution,
}

impl ComponentDefinition {
    /// Return the component identifier.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Return the component name when known.
    #[must_use]
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// Return the PDBx component type when known.
    #[must_use]
    pub fn component_type(&self) -> Option<&str> {
        self.component_type.as_deref()
    }

    /// Return the chemical formula when known.
    #[must_use]
    pub fn formula(&self) -> Option<&str> {
        self.formula.as_deref()
    }

    /// Return the formula weight in daltons when known.
    #[must_use]
    pub const fn formula_weight(&self) -> Option<f64> {
        self.formula_weight
    }

    /// Return the component formal charge when explicitly defined.
    #[must_use]
    pub const fn formal_charge(&self) -> Option<i32> {
        self.formal_charge
    }

    /// Return component atoms in source order.
    #[must_use]
    pub fn atoms(&self) -> &[ComponentAtom] {
        &self.atoms
    }

    /// Return component bonds in source order.
    #[must_use]
    pub fn bonds(&self) -> &[ComponentBond] {
        &self.bonds
    }

    /// Return how this definition was resolved.
    #[must_use]
    pub const fn resolution(&self) -> ComponentResolution {
        self.resolution
    }
}

/// One coordinate atom with label and author namespaces retained separately.
#[derive(Clone, Debug, PartialEq)]
pub struct AtomSite {
    pub(super) group_pdb: String,
    pub(super) id: String,
    pub(super) element: String,
    pub(super) label_atom_id: String,
    pub(super) label_alt_id: Option<String>,
    pub(super) label_comp_id: String,
    pub(super) label_asym_id: String,
    pub(super) label_entity_id: String,
    pub(super) label_seq_id: Option<i64>,
    pub(super) insertion_code: Option<String>,
    pub(super) x: f64,
    pub(super) y: f64,
    pub(super) z: f64,
    pub(super) occupancy: f64,
    pub(super) b_iso: Option<f64>,
    pub(super) formal_charge: Option<i32>,
    pub(super) auth_seq_id: String,
    pub(super) auth_comp_id: String,
    pub(super) auth_asym_id: String,
    pub(super) auth_atom_id: String,
    pub(super) model_number: i64,
}

impl AtomSite {
    /// Return the source atom-site identifier.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Return the `ATOM` or `HETATM` interoperability value.
    #[must_use]
    pub fn group_pdb(&self) -> &str {
        &self.group_pdb
    }

    /// Return the element symbol.
    #[must_use]
    pub fn element(&self) -> &str {
        &self.element
    }

    /// Return label-space atom, component, asym, and entity identifiers.
    #[must_use]
    pub fn label_ids(&self) -> (&str, &str, &str, &str) {
        (
            &self.label_atom_id,
            &self.label_comp_id,
            &self.label_asym_id,
            &self.label_entity_id,
        )
    }

    /// Return the optional label sequence identifier.
    #[must_use]
    pub const fn label_seq_id(&self) -> Option<i64> {
        self.label_seq_id
    }

    /// Return the label alternate identifier, with CIF missing states represented by `None`.
    #[must_use]
    pub fn label_alt_id(&self) -> Option<&str> {
        self.label_alt_id.as_deref()
    }

    /// Return the Cartesian coordinate in angstroms.
    #[must_use]
    pub const fn coordinates(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Return occupancy and optional isotropic displacement value.
    #[must_use]
    pub const fn occupancy_and_b_iso(&self) -> (f64, Option<f64>) {
        (self.occupancy, self.b_iso)
    }

    /// Return the coordinate formal charge when supplied.
    #[must_use]
    pub const fn formal_charge(&self) -> Option<i32> {
        self.formal_charge
    }

    /// Return author-space sequence, component, asym, and atom identifiers.
    #[must_use]
    pub fn author_ids(&self) -> (&str, &str, &str, &str) {
        (
            &self.auth_seq_id,
            &self.auth_comp_id,
            &self.auth_asym_id,
            &self.auth_atom_id,
        )
    }

    /// Return the positive PDB model number.
    #[must_use]
    pub const fn model_number(&self) -> i64 {
        self.model_number
    }
}

/// One label-space endpoint of an explicit inter-site connection.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConnectionEndpoint {
    pub(super) asym_id: String,
    pub(super) component_id: Option<String>,
    pub(super) sequence_id: Option<i64>,
    pub(super) atom_id: String,
}

impl ConnectionEndpoint {
    /// Return the label asym-unit identifier.
    #[must_use]
    pub fn asym_id(&self) -> &str {
        &self.asym_id
    }

    /// Return the label component identifier when supplied.
    #[must_use]
    pub fn component_id(&self) -> Option<&str> {
        self.component_id.as_deref()
    }

    /// Return the label sequence identifier when applicable.
    #[must_use]
    pub const fn sequence_id(&self) -> Option<i64> {
        self.sequence_id
    }

    /// Return the label atom identifier.
    #[must_use]
    pub fn atom_id(&self) -> &str {
        &self.atom_id
    }
}

/// One explicit inter-site connection from `_struct_conn`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Connection {
    pub(super) id: String,
    pub(super) kind: String,
    pub(super) first: ConnectionEndpoint,
    pub(super) second: ConnectionEndpoint,
}

impl Connection {
    /// Return the connection identifier and PDBx connection type.
    #[must_use]
    pub fn identity(&self) -> (&str, &str) {
        (&self.id, &self.kind)
    }

    /// Return the two label-space endpoints.
    #[must_use]
    pub const fn endpoints(&self) -> (&ConnectionEndpoint, &ConnectionEndpoint) {
        (&self.first, &self.second)
    }
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
    pub(super) block_code: String,
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
    /// Return the source data-block code.
    #[must_use]
    pub fn block_code(&self) -> &str {
        &self.block_code
    }

    /// Return `_entry.id` without conflating it with the block code.
    #[must_use]
    pub fn entry_id(&self) -> &str {
        &self.entry_id
    }

    /// Return entities in source order.
    #[must_use]
    pub fn entities(&self) -> &[Entity] {
        &self.entities
    }

    /// Return asym-unit instances in source order.
    #[must_use]
    pub fn asym_units(&self) -> &[AsymUnit] {
        &self.asym_units
    }

    /// Return component definitions in stable identifier order.
    #[must_use]
    pub fn components(&self) -> &[ComponentDefinition] {
        &self.components
    }

    /// Return coordinate atom sites in source order.
    #[must_use]
    pub fn atom_sites(&self) -> &[AtomSite] {
        &self.atom_sites
    }

    /// Return explicit inter-site connections in source order.
    #[must_use]
    pub fn connections(&self) -> &[Connection] {
        &self.connections
    }

    pub(crate) const fn source_document(&self) -> &CifDocument {
        &self.document
    }
}

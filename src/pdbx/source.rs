//! Construction of one immutable semantic model from generic CIF category rows.

use std::collections::{BTreeMap, BTreeSet};

use crate::cif::{BlockKind, CifDocument};

use super::category::{CategoryIndex, Row, case_key};
use super::component_source::{component_definitions, definitions_conflict};
use super::error::SemanticError;
use super::fields::{
    duplicate, optional_bool, optional_charge, optional_float, optional_integer, optional_text,
    required_float, required_integer, required_single, required_text, row_context,
};
use super::model::{
    AsymUnit, AtomSite, BranchSchemeRow, BranchedLink, BranchedNode, ComponentResolution,
    Connection, ConnectionEndpoint, Entity, EntityKind, NonPolySchemeRow, PdbxModel, PolySchemeRow,
    PolymerMonomer,
};
use super::registry::{ComponentRegistry, resolve_builtin, unresolved};

/// Build a source-backed semantic PDBx model without performing profile validation.
///
/// Unknown components remain explicit with [`ComponentResolution::Unresolved`].
///
/// # Errors
///
/// Returns [`SemanticError`] when required structural identifiers or typed atom-site values
/// are absent or malformed.
pub fn build_model(document: &CifDocument) -> Result<PdbxModel, SemanticError> {
    build_model_with_registry(document, None)
}

/// Build a source-backed PDBx model using an optional caller-selected local CCD registry.
///
/// # Errors
///
/// Returns [`SemanticError`] for structural failures or incompatible embedded and cached
/// component definitions.
pub fn build_model_with_registry(
    document: &CifDocument,
    registry: Option<&ComponentRegistry>,
) -> Result<PdbxModel, SemanticError> {
    let [block] = document.blocks() else {
        return Err(SemanticError::new(
            "PDBX_MODEL_BLOCK_COUNT",
            "a PDBx model requires exactly one data block",
            Vec::new(),
        ));
    };
    if block.kind() != BlockKind::Data {
        return Err(SemanticError::new(
            "PDBX_MODEL_DATA_BLOCK",
            "a PDBx model cannot be constructed from a global block",
            Vec::new(),
        ));
    }
    let block_code = block.code().unwrap_or_default().to_owned();
    let categories = CategoryIndex::new(block);
    let entry_id = required_single(categories.rows("entry"), "entry", "id")?;
    let audit_dictionary = audit_conform(&categories)?;

    let mut components = component_definitions(&categories, ComponentResolution::Embedded)?;

    let entities = entities(&categories)?;
    let asym_units = asym_units(&categories)?;
    let atom_sites = atom_sites(&categories)?;
    let poly_scheme = poly_scheme(&categories)?;
    let nonpoly_scheme = nonpoly_scheme(&categories)?;
    let branch_scheme = branch_scheme(&categories)?;
    let connections = connections(&categories)?;

    let referenced = referenced_components(&entities, &atom_sites);
    for component_id in referenced {
        let key = case_key(&component_id);
        if let (Some(embedded), Some(cached)) = (
            components.get(&key),
            registry.and_then(|value| value.get(&component_id)),
        ) {
            if definitions_conflict(embedded, cached) {
                return Err(SemanticError::new(
                    "PDBX_COMPONENT_CONFLICT",
                    format!(
                        "embedded and local CCD definitions for component {component_id:?} conflict"
                    ),
                    vec![format!("component={component_id}")],
                ));
            }
        }
        components.entry(key).or_insert_with(|| {
            registry
                .and_then(|value| value.get(&component_id))
                .cloned()
                .or_else(|| resolve_builtin(&component_id))
                .unwrap_or_else(|| unresolved(component_id))
        });
    }

    Ok(PdbxModel {
        document: document.clone(),
        block_code,
        entry_id,
        audit_dictionary,
        entities,
        asym_units,
        components: components.into_values().collect(),
        atom_sites,
        connections,
        poly_scheme,
        nonpoly_scheme,
        branch_scheme,
    })
}

/// Load component definitions from a local, already-parsed CCD CIF artifact.
///
/// # Errors
///
/// Returns [`SemanticError`] when component identifiers or typed chemical fields are invalid.
pub fn build_component_registry(
    document: &CifDocument,
) -> Result<ComponentRegistry, SemanticError> {
    let mut definitions = BTreeMap::new();
    for block in document
        .blocks()
        .iter()
        .filter(|block| block.kind() == BlockKind::Data)
    {
        let categories = CategoryIndex::new(block);
        for (key, definition) in component_definitions(&categories, ComponentResolution::LocalCcd)?
        {
            if definitions.insert(key, definition).is_some() {
                return Err(SemanticError::new(
                    "PDBX_COMPONENT_DUPLICATE",
                    "the local CCD cache defines a component more than once",
                    Vec::new(),
                ));
            }
        }
    }
    Ok(ComponentRegistry::new(definitions))
}

fn audit_conform(
    categories: &CategoryIndex<'_>,
) -> Result<Option<(String, String)>, SemanticError> {
    let rows = categories.rows("audit_conform");
    let Some(row) = rows.first() else {
        return Ok(None);
    };
    Ok(Some((
        required_text(row, "audit_conform", "dict_name", 0)?,
        required_text(row, "audit_conform", "dict_version", 0)?,
    )))
}

fn entities(categories: &CategoryIndex<'_>) -> Result<Vec<Entity>, SemanticError> {
    let polymer_rows = keyed_rows(categories, "entity_poly", "entity_id")?;
    let nonpoly_rows = keyed_rows(categories, "pdbx_entity_nonpoly", "entity_id")?;
    let mut sequences: BTreeMap<String, Vec<PolymerMonomer>> = BTreeMap::new();
    for (row_index, row) in categories.rows("entity_poly_seq").iter().enumerate() {
        let entity_id = required_text(row, "entity_poly_seq", "entity_id", row_index)?;
        sequences
            .entry(case_key(&entity_id))
            .or_default()
            .push(PolymerMonomer {
                number: required_integer(row, "entity_poly_seq", "num", row_index)?,
                component_id: required_text(row, "entity_poly_seq", "mon_id", row_index)?,
                heterogeneous: optional_bool(row, "hetero").unwrap_or(false),
            });
    }
    for sequence in sequences.values_mut() {
        sequence.sort_by_key(|monomer| monomer.number);
    }
    let mut branch_nodes: BTreeMap<String, Vec<BranchedNode>> = BTreeMap::new();
    for (row_index, row) in categories
        .rows("pdbx_entity_branch_list")
        .iter()
        .enumerate()
    {
        let entity_id = required_text(row, "pdbx_entity_branch_list", "entity_id", row_index)?;
        branch_nodes
            .entry(case_key(&entity_id))
            .or_default()
            .push(BranchedNode {
                number: required_integer(row, "pdbx_entity_branch_list", "num", row_index)?,
                component_id: required_text(row, "pdbx_entity_branch_list", "comp_id", row_index)?,
            });
    }
    let branch_links = branch_links(categories)?;

    let mut output = Vec::new();
    for (row_index, row) in categories.rows("entity").iter().enumerate() {
        let id = required_text(row, "entity", "id", row_index)?;
        let kind_text = required_text(row, "entity", "type", row_index)?;
        let kind = match case_key(&kind_text).as_str() {
            "polymer" => EntityKind::Polymer,
            "non-polymer" => EntityKind::NonPolymer,
            "branched" => EntityKind::Branched,
            "water" => EntityKind::Water,
            _ => {
                return Err(SemanticError::new(
                    "PDBX_ENTITY_TYPE",
                    format!("entity {id:?} has unsupported type {kind_text:?}"),
                    row_context("entity", "type", row_index),
                ));
            }
        };
        let key = case_key(&id);
        let polymer_type = polymer_rows
            .get(&key)
            .and_then(|row| optional_text(row, "type"));
        let component_id = nonpoly_rows
            .get(&key)
            .and_then(|row| optional_text(row, "comp_id"));
        output.push(Entity {
            id,
            kind,
            description: optional_text(row, "pdbx_description"),
            polymer_type,
            sequence: sequences.remove(&key).unwrap_or_default(),
            component_id,
            branch_nodes: branch_nodes.remove(&key).unwrap_or_default(),
            branch_links: branch_links.get(&key).cloned().unwrap_or_default(),
        });
    }
    Ok(output)
}

fn branch_links(
    categories: &CategoryIndex<'_>,
) -> Result<BTreeMap<String, Vec<BranchedLink>>, SemanticError> {
    let mut output: BTreeMap<String, Vec<BranchedLink>> = BTreeMap::new();
    for (row_index, row) in categories
        .rows("pdbx_entity_branch_link")
        .iter()
        .enumerate()
    {
        let entity_id = required_text(row, "pdbx_entity_branch_link", "entity_id", row_index)?;
        output
            .entry(case_key(&entity_id))
            .or_default()
            .push(BranchedLink {
                first_node: required_integer(
                    row,
                    "pdbx_entity_branch_link",
                    "entity_branch_list_num_1",
                    row_index,
                )?,
                second_node: required_integer(
                    row,
                    "pdbx_entity_branch_link",
                    "entity_branch_list_num_2",
                    row_index,
                )?,
                first_atom_id: required_text(
                    row,
                    "pdbx_entity_branch_link",
                    "atom_id_1",
                    row_index,
                )?,
                second_atom_id: required_text(
                    row,
                    "pdbx_entity_branch_link",
                    "atom_id_2",
                    row_index,
                )?,
                order: optional_text(row, "value_order"),
            });
    }
    Ok(output)
}

fn asym_units(categories: &CategoryIndex<'_>) -> Result<Vec<AsymUnit>, SemanticError> {
    categories
        .rows("struct_asym")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(AsymUnit {
                id: required_text(row, "struct_asym", "id", row_index)?,
                entity_id: required_text(row, "struct_asym", "entity_id", row_index)?,
                details: optional_text(row, "details"),
            })
        })
        .collect()
}

fn atom_sites(categories: &CategoryIndex<'_>) -> Result<Vec<AtomSite>, SemanticError> {
    categories
        .rows("atom_site")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(AtomSite {
                group_pdb: required_text(row, "atom_site", "group_pdb", row_index)?,
                id: required_text(row, "atom_site", "id", row_index)?,
                element: required_text(row, "atom_site", "type_symbol", row_index)?,
                label_atom_id: required_text(row, "atom_site", "label_atom_id", row_index)?,
                label_alt_id: optional_text(row, "label_alt_id"),
                label_comp_id: required_text(row, "atom_site", "label_comp_id", row_index)?,
                label_asym_id: required_text(row, "atom_site", "label_asym_id", row_index)?,
                label_entity_id: required_text(row, "atom_site", "label_entity_id", row_index)?,
                label_seq_id: optional_integer(row, "atom_site", "label_seq_id", row_index)?,
                insertion_code: optional_text(row, "pdbx_pdb_ins_code"),
                x: required_float(row, "atom_site", "cartn_x", row_index)?,
                y: required_float(row, "atom_site", "cartn_y", row_index)?,
                z: required_float(row, "atom_site", "cartn_z", row_index)?,
                occupancy: required_float(row, "atom_site", "occupancy", row_index)?,
                b_iso: optional_float(row, "atom_site", "b_iso_or_equiv", row_index)?,
                formal_charge: optional_charge(row, "atom_site", "pdbx_formal_charge", row_index)?,
                auth_seq_id: required_text(row, "atom_site", "auth_seq_id", row_index)?,
                auth_comp_id: required_text(row, "atom_site", "auth_comp_id", row_index)?,
                auth_asym_id: required_text(row, "atom_site", "auth_asym_id", row_index)?,
                auth_atom_id: required_text(row, "atom_site", "auth_atom_id", row_index)?,
                model_number: required_integer(row, "atom_site", "pdbx_pdb_model_num", row_index)?,
            })
        })
        .collect()
}

fn poly_scheme(categories: &CategoryIndex<'_>) -> Result<Vec<PolySchemeRow>, SemanticError> {
    categories
        .rows("pdbx_poly_seq_scheme")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(PolySchemeRow {
                asym_id: required_text(row, "pdbx_poly_seq_scheme", "asym_id", row_index)?,
                entity_id: required_text(row, "pdbx_poly_seq_scheme", "entity_id", row_index)?,
                seq_id: required_integer(row, "pdbx_poly_seq_scheme", "seq_id", row_index)?,
                component_id: required_text(row, "pdbx_poly_seq_scheme", "mon_id", row_index)?,
            })
        })
        .collect()
}

fn nonpoly_scheme(categories: &CategoryIndex<'_>) -> Result<Vec<NonPolySchemeRow>, SemanticError> {
    categories
        .rows("pdbx_nonpoly_scheme")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(NonPolySchemeRow {
                asym_id: required_text(row, "pdbx_nonpoly_scheme", "asym_id", row_index)?,
                entity_id: required_text(row, "pdbx_nonpoly_scheme", "entity_id", row_index)?,
                component_id: required_text(row, "pdbx_nonpoly_scheme", "mon_id", row_index)?,
                auth_seq_id: required_text(row, "pdbx_nonpoly_scheme", "auth_seq_num", row_index)?,
            })
        })
        .collect()
}

fn branch_scheme(categories: &CategoryIndex<'_>) -> Result<Vec<BranchSchemeRow>, SemanticError> {
    categories
        .rows("pdbx_branch_scheme")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(BranchSchemeRow {
                asym_id: required_text(row, "pdbx_branch_scheme", "asym_id", row_index)?,
                entity_id: required_text(row, "pdbx_branch_scheme", "entity_id", row_index)?,
                number: required_integer(row, "pdbx_branch_scheme", "num", row_index)?,
                component_id: required_text(row, "pdbx_branch_scheme", "mon_id", row_index)?,
                auth_seq_id: match optional_text(row, "auth_seq_num") {
                    Some(value) => value,
                    None => required_text(row, "pdbx_branch_scheme", "pdb_seq_num", row_index)?,
                },
            })
        })
        .collect()
}

fn connections(categories: &CategoryIndex<'_>) -> Result<Vec<Connection>, SemanticError> {
    categories
        .rows("struct_conn")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(Connection {
                id: required_text(row, "struct_conn", "id", row_index)?,
                kind: required_text(row, "struct_conn", "conn_type_id", row_index)?,
                first: connection_endpoint(row, row_index, "ptnr1")?,
                second: connection_endpoint(row, row_index, "ptnr2")?,
            })
        })
        .collect()
}

fn connection_endpoint(
    row: &Row<'_>,
    row_index: usize,
    prefix: &str,
) -> Result<ConnectionEndpoint, SemanticError> {
    Ok(ConnectionEndpoint {
        asym_id: required_text(
            row,
            "struct_conn",
            &format!("{prefix}_label_asym_id"),
            row_index,
        )?,
        component_id: optional_text(row, &format!("{prefix}_label_comp_id")),
        sequence_id: optional_integer(
            row,
            "struct_conn",
            &format!("{prefix}_label_seq_id"),
            row_index,
        )?,
        atom_id: required_text(
            row,
            "struct_conn",
            &format!("{prefix}_label_atom_id"),
            row_index,
        )?,
    })
}

fn referenced_components(entities: &[Entity], atoms: &[AtomSite]) -> BTreeSet<String> {
    let mut output = BTreeSet::new();
    for entity in entities {
        output.extend(
            entity
                .sequence
                .iter()
                .map(|monomer| monomer.component_id.clone()),
        );
        output.extend(entity.component_id.iter().cloned());
        output.extend(
            entity
                .branch_nodes
                .iter()
                .map(|node| node.component_id.clone()),
        );
    }
    output.extend(atoms.iter().map(|atom| atom.label_comp_id.clone()));
    output
}

fn keyed_rows<'a>(
    categories: &'a CategoryIndex<'a>,
    category: &str,
    key_item: &str,
) -> Result<BTreeMap<String, &'a Row<'a>>, SemanticError> {
    let mut output = BTreeMap::new();
    for (row_index, row) in categories.rows(category).iter().enumerate() {
        let key = required_text(row, category, key_item, row_index)?;
        if output.insert(case_key(&key), row).is_some() {
            return Err(duplicate(category, key_item, &key, row_index));
        }
    }
    Ok(output)
}

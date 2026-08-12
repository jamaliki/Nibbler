//! Construction of one immutable semantic model from generic CIF category rows.

use std::collections::{BTreeMap, BTreeSet};

use crate::cif::{BlockKind, CategoryView, CifDocument};

use super::component_source::{component_definitions, definitions_conflict};
use super::error::SemanticError;
use super::fields::{
    case_key, category_rows, duplicate, optional_bool, optional_charge, optional_float,
    optional_integer, optional_text, required_float, required_integer, required_single,
    required_text, row_context,
};
use super::model::{
    AsymUnit, AtomSite, BranchSchemeRow, BranchedLink, BranchedNode, ComponentResolution,
    Connection, ConnectionEndpoint, Entity, EntityKind, NonPolySchemeRow, PdbxModel, PolySchemeRow,
    PolymerMonomer,
};
use super::registry::{ComponentRegistry, resolve_builtin, unresolved};

/// Build a source-backed semantic PDBx model without performing profile validation.
///
/// # Errors
///
/// Returns [`SemanticError`] when required structural identifiers or typed values are invalid.
pub fn build_model(document: &CifDocument) -> Result<PdbxModel, SemanticError> {
    build_model_with_registry(document, None)
}

/// Build a source-backed PDBx model using an optional local CCD registry.
///
/// # Errors
///
/// Returns [`SemanticError`] for structural failures or incompatible component definitions.
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
    let categories = CategoryView::new(block.entries());
    let entry_id = required_single(categories, "entry", "id")?;
    let audit_dictionary = audit_conform(categories)?;
    let mut components = component_definitions(categories, ComponentResolution::Embedded)?;
    let entities = entities(categories)?;
    let asym_units = asym_units(categories)?;
    let atom_sites = atom_sites(categories)?;

    for component_id in referenced_components(&entities, &atom_sites) {
        let key = case_key(&component_id);
        if let (Some(embedded), Some(cached)) = (
            components.get(&key),
            registry.and_then(|value| value.get(&component_id)),
        ) && definitions_conflict(embedded, cached)
        {
            return Err(SemanticError::new(
                "PDBX_COMPONENT_CONFLICT",
                format!(
                    "embedded and local CCD definitions for component {component_id:?} conflict"
                ),
                vec![format!("component={component_id}")],
            ));
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
        entry_id,
        audit_dictionary,
        entities,
        asym_units,
        components: components.into_values().collect(),
        atom_sites,
        connections: connections(categories)?,
        poly_scheme: poly_scheme(categories)?,
        nonpoly_scheme: nonpoly_scheme(categories)?,
        branch_scheme: branch_scheme(categories)?,
    })
}

/// Load component definitions from an already-parsed local CCD artifact.
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
        for (key, definition) in component_definitions(
            CategoryView::new(block.entries()),
            ComponentResolution::LocalCcd,
        )? {
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

fn audit_conform(categories: CategoryView<'_>) -> Result<Option<(String, String)>, SemanticError> {
    let Some(row) =
        category_rows(categories, "audit_conform", ["dict_name", "dict_version"]).next()
    else {
        return Ok(None);
    };
    Ok(Some((
        required_text(&row, "audit_conform", "dict_name", 0)?,
        required_text(&row, "audit_conform", "dict_version", 1)?,
    )))
}

fn entities(categories: CategoryView<'_>) -> Result<Vec<Entity>, SemanticError> {
    let polymer_rows = keyed_optional_text(categories, "entity_poly", "entity_id", "type")?;
    let nonpoly_rows =
        keyed_optional_text(categories, "pdbx_entity_nonpoly", "entity_id", "comp_id")?;
    let mut sequences: BTreeMap<String, Vec<PolymerMonomer>> = BTreeMap::new();
    for row in category_rows(
        categories,
        "entity_poly_seq",
        ["entity_id", "num", "mon_id", "hetero"],
    ) {
        let entity_id = required_text(&row, "entity_poly_seq", "entity_id", 0)?;
        sequences
            .entry(case_key(&entity_id))
            .or_default()
            .push(PolymerMonomer {
                number: required_integer(&row, "entity_poly_seq", "num", 1)?,
                component_id: required_text(&row, "entity_poly_seq", "mon_id", 2)?,
                heterogeneous: optional_bool(&row, 3).unwrap_or(false),
            });
    }
    for sequence in sequences.values_mut() {
        sequence.sort_by_key(|monomer| monomer.number);
    }

    let mut branch_nodes: BTreeMap<String, Vec<BranchedNode>> = BTreeMap::new();
    for row in category_rows(
        categories,
        "pdbx_entity_branch_list",
        ["entity_id", "num", "comp_id"],
    ) {
        let entity_id = required_text(&row, "pdbx_entity_branch_list", "entity_id", 0)?;
        branch_nodes
            .entry(case_key(&entity_id))
            .or_default()
            .push(BranchedNode {
                number: required_integer(&row, "pdbx_entity_branch_list", "num", 1)?,
                component_id: required_text(&row, "pdbx_entity_branch_list", "comp_id", 2)?,
            });
    }
    let branch_links = branch_links(categories)?;

    category_rows(categories, "entity", ["id", "type", "pdbx_description"])
        .map(|row| {
            let id = required_text(&row, "entity", "id", 0)?;
            let kind_text = required_text(&row, "entity", "type", 1)?;
            let kind = match case_key(&kind_text).as_str() {
                "polymer" => EntityKind::Polymer,
                "non-polymer" => EntityKind::NonPolymer,
                "branched" => EntityKind::Branched,
                "water" => EntityKind::Water,
                _ => {
                    return Err(SemanticError::new(
                        "PDBX_ENTITY_TYPE",
                        format!("entity {id:?} has unsupported type {kind_text:?}"),
                        row_context("entity", "type", row.row_index()),
                    ));
                }
            };
            let key = case_key(&id);
            Ok(Entity {
                id,
                kind,
                description: optional_text(&row, 2),
                polymer_type: polymer_rows.get(&key).cloned().flatten(),
                sequence: sequences.remove(&key).unwrap_or_default(),
                component_id: nonpoly_rows.get(&key).cloned().flatten(),
                branch_nodes: branch_nodes.remove(&key).unwrap_or_default(),
                branch_links: branch_links.get(&key).cloned().unwrap_or_default(),
            })
        })
        .collect()
}

fn keyed_optional_text(
    categories: CategoryView<'_>,
    category: &str,
    key_item: &str,
    value_item: &str,
) -> Result<BTreeMap<String, Option<String>>, SemanticError> {
    let mut output = BTreeMap::new();
    for row in category_rows(categories, category, [key_item, value_item]) {
        let key = required_text(&row, category, key_item, 0)?;
        if output
            .insert(case_key(&key), optional_text(&row, 1))
            .is_some()
        {
            return Err(duplicate(category, key_item, &key, row.row_index()));
        }
    }
    Ok(output)
}

fn branch_links(
    categories: CategoryView<'_>,
) -> Result<BTreeMap<String, Vec<BranchedLink>>, SemanticError> {
    let category = "pdbx_entity_branch_link";
    let mut output: BTreeMap<String, Vec<BranchedLink>> = BTreeMap::new();
    for row in category_rows(
        categories,
        category,
        [
            "entity_id",
            "entity_branch_list_num_1",
            "entity_branch_list_num_2",
            "atom_id_1",
            "atom_id_2",
            "value_order",
        ],
    ) {
        let entity_id = required_text(&row, category, "entity_id", 0)?;
        output
            .entry(case_key(&entity_id))
            .or_default()
            .push(BranchedLink {
                first_node: required_integer(&row, category, "entity_branch_list_num_1", 1)?,
                second_node: required_integer(&row, category, "entity_branch_list_num_2", 2)?,
                first_atom_id: required_text(&row, category, "atom_id_1", 3)?,
                second_atom_id: required_text(&row, category, "atom_id_2", 4)?,
                order: optional_text(&row, 5),
            });
    }
    Ok(output)
}

fn asym_units(categories: CategoryView<'_>) -> Result<Vec<AsymUnit>, SemanticError> {
    category_rows(categories, "struct_asym", ["id", "entity_id", "details"])
        .map(|row| {
            Ok(AsymUnit {
                id: required_text(&row, "struct_asym", "id", 0)?,
                entity_id: required_text(&row, "struct_asym", "entity_id", 1)?,
                details: optional_text(&row, 2),
            })
        })
        .collect()
}

fn atom_sites(categories: CategoryView<'_>) -> Result<Vec<AtomSite>, SemanticError> {
    let category = "atom_site";
    category_rows(
        categories,
        category,
        [
            "group_pdb",
            "id",
            "type_symbol",
            "label_atom_id",
            "label_alt_id",
            "label_comp_id",
            "label_asym_id",
            "label_entity_id",
            "label_seq_id",
            "pdbx_pdb_ins_code",
            "cartn_x",
            "cartn_y",
            "cartn_z",
            "occupancy",
            "b_iso_or_equiv",
            "pdbx_formal_charge",
            "auth_seq_id",
            "auth_comp_id",
            "auth_asym_id",
            "auth_atom_id",
            "pdbx_pdb_model_num",
        ],
    )
    .map(|row| {
        Ok(AtomSite {
            group_pdb: required_text(&row, category, "group_pdb", 0)?,
            id: required_text(&row, category, "id", 1)?,
            element: required_text(&row, category, "type_symbol", 2)?,
            label_atom_id: required_text(&row, category, "label_atom_id", 3)?,
            label_alt_id: optional_text(&row, 4),
            label_comp_id: required_text(&row, category, "label_comp_id", 5)?,
            label_asym_id: required_text(&row, category, "label_asym_id", 6)?,
            label_entity_id: required_text(&row, category, "label_entity_id", 7)?,
            label_seq_id: optional_integer(&row, category, "label_seq_id", 8)?,
            insertion_code: optional_text(&row, 9),
            x: required_float(&row, category, "cartn_x", 10)?,
            y: required_float(&row, category, "cartn_y", 11)?,
            z: required_float(&row, category, "cartn_z", 12)?,
            occupancy: required_float(&row, category, "occupancy", 13)?,
            b_iso: optional_float(&row, category, "b_iso_or_equiv", 14)?,
            formal_charge: optional_charge(&row, category, "pdbx_formal_charge", 15)?,
            auth_seq_id: required_text(&row, category, "auth_seq_id", 16)?,
            auth_comp_id: required_text(&row, category, "auth_comp_id", 17)?,
            auth_asym_id: required_text(&row, category, "auth_asym_id", 18)?,
            auth_atom_id: required_text(&row, category, "auth_atom_id", 19)?,
            model_number: required_integer(&row, category, "pdbx_pdb_model_num", 20)?,
        })
    })
    .collect()
}

fn poly_scheme(categories: CategoryView<'_>) -> Result<Vec<PolySchemeRow>, SemanticError> {
    let category = "pdbx_poly_seq_scheme";
    category_rows(
        categories,
        category,
        ["asym_id", "entity_id", "seq_id", "mon_id"],
    )
    .map(|row| {
        Ok(PolySchemeRow {
            asym_id: required_text(&row, category, "asym_id", 0)?,
            entity_id: required_text(&row, category, "entity_id", 1)?,
            seq_id: required_integer(&row, category, "seq_id", 2)?,
            component_id: required_text(&row, category, "mon_id", 3)?,
        })
    })
    .collect()
}

fn nonpoly_scheme(categories: CategoryView<'_>) -> Result<Vec<NonPolySchemeRow>, SemanticError> {
    let category = "pdbx_nonpoly_scheme";
    category_rows(
        categories,
        category,
        ["asym_id", "entity_id", "mon_id", "auth_seq_num"],
    )
    .map(|row| {
        Ok(NonPolySchemeRow {
            asym_id: required_text(&row, category, "asym_id", 0)?,
            entity_id: required_text(&row, category, "entity_id", 1)?,
            component_id: required_text(&row, category, "mon_id", 2)?,
            auth_seq_id: required_text(&row, category, "auth_seq_num", 3)?,
        })
    })
    .collect()
}

fn branch_scheme(categories: CategoryView<'_>) -> Result<Vec<BranchSchemeRow>, SemanticError> {
    let category = "pdbx_branch_scheme";
    category_rows(
        categories,
        category,
        [
            "asym_id",
            "entity_id",
            "num",
            "mon_id",
            "auth_seq_num",
            "pdb_seq_num",
        ],
    )
    .map(|row| {
        Ok(BranchSchemeRow {
            asym_id: required_text(&row, category, "asym_id", 0)?,
            entity_id: required_text(&row, category, "entity_id", 1)?,
            number: required_integer(&row, category, "num", 2)?,
            component_id: required_text(&row, category, "mon_id", 3)?,
            auth_seq_id: match optional_text(&row, 4) {
                Some(value) => value,
                None => required_text(&row, category, "pdb_seq_num", 5)?,
            },
        })
    })
    .collect()
}

fn connections(categories: CategoryView<'_>) -> Result<Vec<Connection>, SemanticError> {
    let category = "struct_conn";
    category_rows(
        categories,
        category,
        [
            "id",
            "conn_type_id",
            "ptnr1_label_asym_id",
            "ptnr1_label_comp_id",
            "ptnr1_label_seq_id",
            "ptnr1_label_atom_id",
            "ptnr2_label_asym_id",
            "ptnr2_label_comp_id",
            "ptnr2_label_seq_id",
            "ptnr2_label_atom_id",
        ],
    )
    .map(|row| {
        Ok(Connection {
            id: required_text(&row, category, "id", 0)?,
            kind: required_text(&row, category, "conn_type_id", 1)?,
            first: ConnectionEndpoint {
                asym_id: required_text(&row, category, "ptnr1_label_asym_id", 2)?,
                component_id: optional_text(&row, 3),
                sequence_id: optional_integer(&row, category, "ptnr1_label_seq_id", 4)?,
                atom_id: required_text(&row, category, "ptnr1_label_atom_id", 5)?,
            },
            second: ConnectionEndpoint {
                asym_id: required_text(&row, category, "ptnr2_label_asym_id", 6)?,
                component_id: optional_text(&row, 7),
                sequence_id: optional_integer(&row, category, "ptnr2_label_seq_id", 8)?,
                atom_id: required_text(&row, category, "ptnr2_label_atom_id", 9)?,
            },
        })
    })
    .collect()
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

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::cif::parse;

    use super::*;

    const CHEMISTRY: &[u8] = include_bytes!("../../tests/fixtures/chemistry/ligand_ion_water.cif");
    const GLYCAN: &[u8] = include_bytes!("../../tests/fixtures/chemistry/branched_glycan.cif");

    #[test]
    fn retains_coordinate_and_chemistry_records_after_decoding() {
        let source = String::from_utf8(CHEMISTRY.to_vec())
            .expect("fixture is UTF-8")
            .replace("PG A ATP B 2 . ?", "PG A ATP B 2 . A");
        let document = parse(source.as_bytes()).expect("fixture parses");
        let model = build_model(&document).expect("fixture builds");

        assert_eq!(
            model.audit_dictionary,
            Some(("mmcif_pdbx.dic".to_owned(), "5.416".to_owned()))
        );
        assert_eq!(
            model.entities[0],
            Entity {
                id: "1".to_owned(),
                kind: EntityKind::Polymer,
                description: Some("ALA-MSE test peptide".to_owned()),
                polymer_type: Some("polypeptide(L)".to_owned()),
                sequence: vec![
                    PolymerMonomer {
                        number: 1,
                        component_id: "ALA".to_owned(),
                        heterogeneous: false,
                    },
                    PolymerMonomer {
                        number: 2,
                        component_id: "MSE".to_owned(),
                        heterogeneous: false,
                    },
                ],
                component_id: None,
                branch_nodes: Vec::new(),
                branch_links: Vec::new(),
            }
        );
        assert_eq!(model.entities[3].component_id.as_deref(), Some("HOH"));
        assert_eq!(
            model.asym_units[1],
            AsymUnit {
                id: "B".to_owned(),
                entity_id: "2".to_owned(),
                details: Some(
                    "ligand with alternate locations and negative author residue number".to_owned()
                ),
            }
        );
        assert_eq!(
            model.atom_sites[2],
            AtomSite {
                group_pdb: "HETATM".to_owned(),
                id: "3".to_owned(),
                element: "P".to_owned(),
                label_atom_id: "PG".to_owned(),
                label_alt_id: Some("A".to_owned()),
                label_comp_id: "ATP".to_owned(),
                label_asym_id: "B".to_owned(),
                label_entity_id: "2".to_owned(),
                label_seq_id: None,
                insertion_code: Some("A".to_owned()),
                x: 7.0,
                y: 1.0,
                z: 0.0,
                occupancy: 0.6,
                b_iso: Some(20.0),
                formal_charge: None,
                auth_seq_id: "-1".to_owned(),
                auth_comp_id: "ATP".to_owned(),
                auth_asym_id: "L".to_owned(),
                auth_atom_id: "PG".to_owned(),
                model_number: 1,
            }
        );
        assert_eq!(model.atom_sites[4].formal_charge, Some(2));
        assert_eq!(
            model.connections[0],
            Connection {
                id: "ZN1".to_owned(),
                kind: "metalc".to_owned(),
                first: ConnectionEndpoint {
                    asym_id: "B".to_owned(),
                    component_id: Some("ATP".to_owned()),
                    sequence_id: None,
                    atom_id: "PG".to_owned(),
                },
                second: ConnectionEndpoint {
                    asym_id: "C".to_owned(),
                    component_id: Some("ZN".to_owned()),
                    sequence_id: None,
                    atom_id: "ZN".to_owned(),
                },
            }
        );
        assert_eq!(
            model.poly_scheme[0],
            PolySchemeRow {
                asym_id: "A".to_owned(),
                entity_id: "1".to_owned(),
                seq_id: 1,
                component_id: "ALA".to_owned(),
            }
        );
        assert_eq!(
            model.nonpoly_scheme[0],
            NonPolySchemeRow {
                asym_id: "B".to_owned(),
                entity_id: "2".to_owned(),
                component_id: "ATP".to_owned(),
                auth_seq_id: "-1".to_owned(),
            }
        );

        let ccd = parse(
            br#"data_TST
_chem_comp.id TST
_chem_comp.name 'TEST COMPONENT'
_chem_comp.type non-polymer
_chem_comp.formula 'C N'
_chem_comp.formula_weight 26.02
_chem_comp.pdbx_formal_charge 1
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
TST C1 C 1
TST N1 N 0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
TST C1 N1 doub y
"#,
        )
        .expect("component fixture parses");
        let registry = build_component_registry(&ccd).expect("component fixture builds");
        let component = registry.get("TST").expect("component exists");
        assert_eq!(component.id, "TST");
        assert_eq!(component.name.as_deref(), Some("TEST COMPONENT"));
        assert_eq!(component.component_type.as_deref(), Some("non-polymer"));
        assert_eq!(component.formula.as_deref(), Some("C N"));
        assert_eq!(component.formula_weight, Some(26.02));
        assert_eq!(component.formal_charge, Some(1));
        assert_eq!(component.resolution, ComponentResolution::LocalCcd);
        assert_eq!(
            component.atoms,
            [
                super::super::model::ComponentAtom {
                    atom_id: "C1".to_owned(),
                    element: "C".to_owned(),
                    formal_charge: Some(1),
                },
                super::super::model::ComponentAtom {
                    atom_id: "N1".to_owned(),
                    element: "N".to_owned(),
                    formal_charge: Some(0),
                },
            ]
        );
        assert_eq!(
            component.bonds,
            [super::super::model::ComponentBond {
                first_atom_id: "C1".to_owned(),
                second_atom_id: "N1".to_owned(),
                order: "doub".to_owned(),
                aromatic: Some(true),
            }]
        );
    }

    #[test]
    fn retains_branched_entity_and_scheme_records_after_decoding() {
        let document = parse(GLYCAN).expect("fixture parses");
        let model = build_model(&document).expect("fixture builds");
        let entity = &model.entities[0];

        assert_eq!(entity.kind, EntityKind::Branched);
        assert_eq!(
            entity.branch_nodes,
            [
                BranchedNode {
                    number: 1,
                    component_id: "NAG".to_owned(),
                },
                BranchedNode {
                    number: 2,
                    component_id: "MAN".to_owned(),
                },
            ]
        );
        assert_eq!(
            entity.branch_links,
            [BranchedLink {
                first_node: 1,
                second_node: 2,
                first_atom_id: "O4".to_owned(),
                second_atom_id: "C1".to_owned(),
                order: Some("sing".to_owned()),
            }]
        );
        assert_eq!(
            model.branch_scheme[1],
            BranchSchemeRow {
                asym_id: "G".to_owned(),
                entity_id: "1".to_owned(),
                number: 2,
                component_id: "MAN".to_owned(),
                auth_seq_id: "2".to_owned(),
            }
        );
    }
}

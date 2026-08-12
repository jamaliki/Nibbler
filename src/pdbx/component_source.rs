//! Chemical-component extraction and cross-provenance conflict checks.

use std::collections::BTreeMap;

use crate::cif::CategoryView;

use super::error::SemanticError;
use super::fields::{
    case_key, category_rows, duplicate, optional_bool, optional_charge, optional_float,
    optional_text, required_text,
};
use super::model::{ComponentAtom, ComponentBond, ComponentDefinition, ComponentResolution};

pub(super) fn component_definitions(
    categories: CategoryView<'_>,
    resolution: ComponentResolution,
) -> Result<BTreeMap<String, ComponentDefinition>, SemanticError> {
    let mut output = BTreeMap::new();
    for row in category_rows(
        categories,
        "chem_comp",
        [
            "id",
            "name",
            "type",
            "formula",
            "formula_weight",
            "pdbx_formal_charge",
        ],
    ) {
        let id = required_text(&row, "chem_comp", "id", 0)?;
        let definition = ComponentDefinition {
            id: id.clone(),
            name: optional_text(&row, 1),
            component_type: optional_text(&row, 2),
            formula: optional_text(&row, 3),
            formula_weight: optional_float(&row, "chem_comp", "formula_weight", 4)?,
            formal_charge: optional_charge(&row, "chem_comp", "pdbx_formal_charge", 5)?,
            atoms: Vec::new(),
            bonds: Vec::new(),
            resolution,
        };
        if output.insert(case_key(&id), definition).is_some() {
            return Err(duplicate("chem_comp", "id", &id, row.row_index()));
        }
    }
    attach_component_atoms(categories, &mut output, resolution)?;
    attach_component_bonds(categories, &mut output, resolution)?;
    Ok(output)
}

pub(super) fn definitions_conflict(
    embedded: &ComponentDefinition,
    cached: &ComponentDefinition,
) -> bool {
    option_text_conflicts(embedded.name.as_deref(), cached.name.as_deref())
        || option_text_conflicts(
            embedded.component_type.as_deref(),
            cached.component_type.as_deref(),
        )
        || option_text_conflicts(embedded.formula.as_deref(), cached.formula.as_deref())
        || option_value_conflicts(embedded.formula_weight, cached.formula_weight)
        || option_value_conflicts(embedded.formal_charge, cached.formal_charge)
        || embedded.atoms.iter().any(|atom| {
            cached
                .atoms
                .iter()
                .find(|candidate| candidate.atom_id.eq_ignore_ascii_case(&atom.atom_id))
                .is_some_and(|candidate| {
                    !candidate.element.eq_ignore_ascii_case(&atom.element)
                        || option_value_conflicts(atom.formal_charge, candidate.formal_charge)
                })
        })
}

fn attach_component_atoms(
    categories: CategoryView<'_>,
    components: &mut BTreeMap<String, ComponentDefinition>,
    resolution: ComponentResolution,
) -> Result<(), SemanticError> {
    for row in category_rows(
        categories,
        "chem_comp_atom",
        ["comp_id", "atom_id", "type_symbol", "charge"],
    ) {
        let component_id = required_text(&row, "chem_comp_atom", "comp_id", 0)?;
        let definition = components
            .entry(case_key(&component_id))
            .or_insert_with(|| component_shell(component_id, resolution));
        definition.atoms.push(ComponentAtom {
            atom_id: required_text(&row, "chem_comp_atom", "atom_id", 1)?,
            element: required_text(&row, "chem_comp_atom", "type_symbol", 2)?,
            formal_charge: optional_charge(&row, "chem_comp_atom", "charge", 3)?,
        });
    }
    Ok(())
}

fn attach_component_bonds(
    categories: CategoryView<'_>,
    components: &mut BTreeMap<String, ComponentDefinition>,
    resolution: ComponentResolution,
) -> Result<(), SemanticError> {
    for row in category_rows(
        categories,
        "chem_comp_bond",
        [
            "comp_id",
            "atom_id_1",
            "atom_id_2",
            "value_order",
            "pdbx_aromatic_flag",
        ],
    ) {
        let component_id = required_text(&row, "chem_comp_bond", "comp_id", 0)?;
        let definition = components
            .entry(case_key(&component_id))
            .or_insert_with(|| component_shell(component_id, resolution));
        definition.bonds.push(ComponentBond {
            first_atom_id: required_text(&row, "chem_comp_bond", "atom_id_1", 1)?,
            second_atom_id: required_text(&row, "chem_comp_bond", "atom_id_2", 2)?,
            order: required_text(&row, "chem_comp_bond", "value_order", 3)?,
            aromatic: optional_bool(&row, 4),
        });
    }
    Ok(())
}

fn component_shell(component_id: String, resolution: ComponentResolution) -> ComponentDefinition {
    ComponentDefinition {
        id: component_id,
        name: None,
        component_type: None,
        formula: None,
        formula_weight: None,
        formal_charge: None,
        atoms: Vec::new(),
        bonds: Vec::new(),
        resolution,
    }
}

fn option_text_conflicts(left: Option<&str>, right: Option<&str>) -> bool {
    matches!((left, right), (Some(left), Some(right)) if !left.eq_ignore_ascii_case(right))
}

fn option_value_conflicts<T: PartialEq>(left: Option<T>, right: Option<T>) -> bool {
    matches!((left, right), (Some(left), Some(right)) if left != right)
}

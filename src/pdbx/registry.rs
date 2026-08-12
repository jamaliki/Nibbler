//! Immutable component registries used only for deterministic resolution.

use std::collections::BTreeMap;

use super::model::{ComponentAtom, ComponentDefinition, ComponentResolution};

/// Caller-selected component definitions loaded from a local CCD CIF artifact.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ComponentRegistry {
    pub(super) definitions: BTreeMap<String, ComponentDefinition>,
}

impl ComponentRegistry {
    pub(super) const fn new(definitions: BTreeMap<String, ComponentDefinition>) -> Self {
        Self { definitions }
    }

    pub(super) fn get(&self, component_id: &str) -> Option<&ComponentDefinition> {
        self.definitions.get(&component_id.to_ascii_lowercase())
    }

    /// Return the number of component definitions in the immutable registry.
    #[must_use]
    pub fn len(&self) -> usize {
        self.definitions.len()
    }

    /// Return whether the registry contains no definitions.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.definitions.is_empty()
    }
}

pub(super) fn resolve_builtin(component_id: &str) -> Option<ComponentDefinition> {
    let upper = component_id.to_ascii_uppercase();
    let record = builtin_record(&upper)?;
    let atoms = record.element.map_or_else(Vec::new, |element| {
        vec![ComponentAtom {
            atom_id: upper,
            element: element.to_owned(),
            formal_charge: record.charge,
        }]
    });
    Some(ComponentDefinition {
        id: component_id.to_owned(),
        name: record.name.map(str::to_owned),
        component_type: Some(record.component_type.to_owned()),
        formula: record.formula.map(str::to_owned),
        formula_weight: record.weight,
        formal_charge: record.charge,
        atoms,
        bonds: Vec::new(),
        resolution: ComponentResolution::MinimalRegistry,
    })
}

pub(super) fn unresolved(component_id: String) -> ComponentDefinition {
    ComponentDefinition {
        id: component_id,
        name: None,
        component_type: None,
        formula: None,
        formula_weight: None,
        formal_charge: None,
        atoms: Vec::new(),
        bonds: Vec::new(),
        resolution: ComponentResolution::Unresolved,
    }
}

struct BuiltinRecord {
    name: Option<&'static str>,
    component_type: &'static str,
    formula: Option<&'static str>,
    weight: Option<f64>,
    charge: Option<i32>,
    element: Option<&'static str>,
}

fn builtin_record(component_id: &str) -> Option<BuiltinRecord> {
    let ion = |name, formula, weight, charge, element| BuiltinRecord {
        name: Some(name),
        component_type: "non-polymer",
        formula: Some(formula),
        weight: Some(weight),
        charge: Some(charge),
        element: Some(element),
    };
    match component_id {
        "HOH" => Some(BuiltinRecord {
            name: Some("WATER"),
            component_type: "non-polymer",
            formula: Some("H2 O"),
            weight: Some(18.015),
            charge: None,
            element: None,
        }),
        "ZN" => Some(ion("ZINC ION", "Zn 2", 65.409, 2, "Zn")),
        "MG" => Some(ion("MAGNESIUM ION", "Mg 2", 24.305, 2, "Mg")),
        "CA" => Some(ion("CALCIUM ION", "Ca 2", 40.078, 2, "Ca")),
        "NA" => Some(ion("SODIUM ION", "Na 1", 22.990, 1, "Na")),
        "K" => Some(ion("POTASSIUM ION", "K 1", 39.098, 1, "K")),
        "CL" => Some(ion("CHLORIDE ION", "Cl 1", 35.453, -1, "Cl")),
        "ALA" | "ARG" | "ASN" | "ASP" | "CYS" | "GLN" | "GLU" | "GLY" | "HIS" | "ILE" | "LEU"
        | "LYS" | "MET" | "PHE" | "PRO" | "SER" | "THR" | "TRP" | "TYR" | "VAL" => {
            Some(BuiltinRecord {
                name: None,
                component_type: "L-peptide linking",
                formula: None,
                weight: None,
                charge: None,
                element: None,
            })
        }
        _ => None,
    }
}

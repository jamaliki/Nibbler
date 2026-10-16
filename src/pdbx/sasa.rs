//! Ligand surface areas read directly from a parsed atom_site table.

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};

use crate::cif::{BlockCategory, CategoryView, CifDocument, CifValueRef, numeric};
use crate::geometry::{sasa, validate_sasa_options};

use super::fields::category_rows;

mod model;
pub use model::model_ligand_sasa;

fn number(table: &BlockCategory<'_>, row: usize, column: Option<usize>) -> Option<f64> {
    match table.value(row, column?)? {
        CifValueRef::Integer(n, _) => Some(n as f64),
        CifValueRef::Float(n, _, _) => n.is_finite().then_some(n),
        CifValueRef::Text(s) => numeric::parse_float(s.as_str()),
        _ => None,
    }
}

fn model_number(
    table: &BlockCategory<'_>,
    row: usize,
    column: Option<usize>,
) -> Result<i64, String> {
    let Some(column) = column else {
        return Ok(1);
    };
    let value = match table.value(row, column) {
        Some(CifValueRef::Integer(n, _)) => Some(n),
        Some(CifValueRef::Text(s)) => numeric::parse_integer(s.as_str()),
        _ => None,
    };
    value.ok_or_else(|| format!("Invalid model number in atom_site row {row}"))
}

fn radius(element: &str, overrides: &HashMap<String, f64>) -> Result<f64, String> {
    let bytes = element.as_bytes();
    let symbol = [
        bytes.first().copied().unwrap_or(0).to_ascii_uppercase(),
        bytes.get(1).copied().unwrap_or(0).to_ascii_uppercase(),
    ];
    if !overrides.is_empty()
        && let Some(r) = overrides.get(&element.to_ascii_uppercase())
    {
        return Ok(*r);
    }
    let r = match symbol {
        [b'C', 0] => 1.70,
        [b'N', 0] => 1.55,
        [b'O', 0] => 1.52,
        [b'F', 0] => 1.47,
        [b'P', 0] | [b'S', 0] => 1.80,
        [b'C', b'L'] => 1.75,
        [b'B', b'R'] => 1.85,
        [b'I', 0] => 1.98,
        [b'S', b'E'] => 1.90,
        [b'H', b'E'] => 1.40,
        [b'N', b'A'] => 2.27,
        [b'M', b'G'] => 1.73,
        [b'K', 0] => 2.75,
        [b'C', b'A'] => 2.31,
        [b'N', b'I'] => 1.63,
        [b'C', b'U'] => 1.40,
        [b'Z', b'N'] => 1.39,
        [b'C', b'D'] => 1.58,
        [b'H', b'G'] => 1.55,
        _ => {
            return Err(format!(
                "No SASA radius for element {element:?}; supply radii"
            ));
        }
    };
    if bytes.len() > 2 {
        return Err(format!("Invalid element {element:?}"));
    }
    Ok(r)
}

struct Atoms {
    coordinates: Vec<[f64; 3]>,
    radii: Vec<f64>,
    rows: Vec<usize>,
    selected: Vec<usize>,
}

fn check_alternate<'a>(
    seen: &mut HashMap<[Cow<'a, str>; 3], Cow<'a, str>>,
    key: [Cow<'a, str>; 3],
    alternate: Cow<'a, str>,
) -> Result<(), String> {
    if key[0].is_empty() || key[1].is_empty() || seen.get(&key).is_some_and(|old| *old != alternate)
    {
        return Err("SASA requires one conformer; call mmcif.select_conformer first".into());
    }
    seen.entry(key).or_insert(alternate);
    Ok(())
}

/// Return zero-based atom_site rows and ligand SASA in square angstroms.
///
/// By default, non-polymer and branched entities are ligands. Water, H, D, T,
/// and zero-occupancy sites neither receive areas nor occlude other atoms.
/// Explicit `atom_rows` selects outputs without changing the occluding context.
/// The input must contain one conformer; biological assemblies must be explicit.
///
/// # Errors
///
/// Rejects invalid coordinates, occupancy, radii or options, alternate locations,
/// an absent model, and explicit rows outside the retained heavy-atom context.
/// Automatic ligand selection requires entity types and label_entity_id values.
pub fn ligand_sasa(
    document: &CifDocument,
    atom_rows: Option<&[usize]>,
    model: i64,
    probe_radius: f64,
    point_count: usize,
    radii: &HashMap<String, f64>,
) -> Result<(Vec<usize>, Vec<f64>), String> {
    validate_options(probe_radius, point_count, radii)?;
    let [block] = document.blocks() else {
        return Err("Ligand SASA requires one data block".into());
    };
    super::selected_document::check_tables(block.entries()).map_err(|e| e.to_string())?;
    let table = BlockCategory::find(block.entries(), "atom_site")
        .ok_or("Ligand SASA requires atom_site")?;
    let types: HashMap<_, _> =
        category_rows(CategoryView::new(block.entries()), "entity", ["id", "type"])
            .filter_map(|row| Some((row.text(0)?, row.text(1)?.to_ascii_lowercase())))
            .collect();
    let atoms = read_atoms(&table, &types, atom_rows, model, radii)?;
    let areas = sasa(
        &atoms.coordinates,
        &atoms.radii,
        &atoms.selected,
        probe_radius,
        point_count,
    )?;
    Ok((atoms.rows, areas))
}

fn read_atoms(
    table: &BlockCategory<'_>,
    types: &HashMap<String, String>,
    atom_rows: Option<&[usize]>,
    model: i64,
    radii: &HashMap<String, f64>,
) -> Result<Atoms, String> {
    let xyz = ["Cartn_x", "Cartn_y", "Cartn_z"].map(|name| table.column(name));
    let element = table
        .column("type_symbol")
        .ok_or("SASA requires type_symbol")?;
    let entity = table.column("label_entity_id");
    let component = table.column("label_comp_id");
    let occupancy = table.column("occupancy");
    let alternate = table.column("label_alt_id");
    let model_column = table.column("pdbx_PDB_model_num");
    let identity =
        ["label_asym_id", "auth_seq_id", "pdbx_PDB_ins_code"].map(|name| table.column(name));
    let sequence = table.column("label_seq_id");
    let mut alternates = HashMap::new();
    let requested: HashSet<_> = atom_rows.unwrap_or_default().iter().copied().collect();
    if atom_rows.is_some_and(|rows| rows.len() != requested.len()) {
        return Err("atom_rows must contain distinct source rows".into());
    }
    let mut result = Atoms {
        coordinates: Vec::new(),
        radii: Vec::new(),
        rows: Vec::new(),
        selected: Vec::new(),
    };
    let mut found = false;
    for row in 0..table.row_count() {
        if model_number(table, row, model_column)? != model {
            continue;
        }
        found = true;
        let occ = occupancy
            .map_or(Some(1.0), |c| number(table, row, Some(c)))
            .filter(|n| (0.0..=1.0).contains(n))
            .ok_or_else(|| format!("Invalid occupancy in atom_site row {row}"))?;
        if occ == 0.0 {
            continue;
        }
        let symbol = table
            .text(row, element)
            .ok_or_else(|| format!("Missing element in atom_site row {row}"))?;
        if ["H", "D", "T"]
            .iter()
            .any(|s| symbol.eq_ignore_ascii_case(s))
        {
            continue;
        }
        let kind = entity
            .and_then(|c| table.text(row, c))
            .and_then(|id| types.get(id.as_ref()));
        let water = kind.is_some_and(|t| t == "water")
            || component.and_then(|c| table.text(row, c)).is_some_and(|s| {
                ["HOH", "DOD", "WAT"]
                    .iter()
                    .any(|name| s.eq_ignore_ascii_case(name))
            });
        if water {
            continue;
        }
        if let Some(label) = alternate.and_then(|c| table.text(row, c)) {
            let mut key = identity.map(|c| c.and_then(|c| table.text(row, c)).unwrap_or_default());
            if key[1].is_empty() {
                key[1] = sequence
                    .and_then(|c| table.text(row, c))
                    .unwrap_or_default();
            }
            check_alternate(&mut alternates, key, label)?;
        }
        let selected = if atom_rows.is_some() {
            requested.contains(&row)
        } else {
            match kind.map(String::as_str) {
                Some("non-polymer" | "branched") => true,
                Some("polymer") => false,
                _ => return Err("Automatic ligand SASA requires entity types and label_entity_id; supply atom_rows".into()),
            }
        };
        if selected {
            result.selected.push(result.coordinates.len());
            result.rows.push(row);
        }
        let position = xyz.map(|c| number(table, row, c));
        let [Some(x), Some(y), Some(z)] = position else {
            return Err(format!("Invalid SASA coordinates in atom_site row {row}"));
        };
        result.coordinates.push([x, y, z]);
        result.radii.push(radius(&symbol, radii)?);
    }
    if !found {
        return Err(format!("Deposited model {model} has no atoms"));
    }
    if atom_rows.is_some_and(|rows| rows.len() != result.rows.len()) {
        return Err(
            "atom_rows must select occupied, non-water heavy atoms in the chosen model".into(),
        );
    }
    Ok(result)
}

fn validate_options(
    probe_radius: f64,
    point_count: usize,
    radii: &HashMap<String, f64>,
) -> Result<(), String> {
    validate_sasa_options(probe_radius, point_count)?;
    if radii
        .values()
        .any(|r| !r.is_finite() || *r <= 0.0 || *r > 100.0)
    {
        return Err("SASA radii must be finite and in (0, 100] A".into());
    }
    Ok(())
}

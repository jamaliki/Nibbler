use super::{Atoms, check_alternate, radius};
use crate::geometry::sasa;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};

/// Calculate ligand areas from an existing semantic model's decoded coordinates.
///
/// Uses the same selection, radii and occlusion rules as [`super::ligand_sasa`].
///
/// # Errors
///
/// Returns the same input and conformer errors as [`super::ligand_sasa`].
pub fn model_ligand_sasa(
    model: &crate::pdbx::PdbxModel,
    atom_rows: Option<&[usize]>,
    number: i64,
    probe_radius: f64,
    point_count: usize,
    radii: &HashMap<String, f64>,
) -> Result<(Vec<usize>, Vec<f64>), String> {
    use crate::pdbx::model::EntityKind;
    super::validate_options(probe_radius, point_count, radii)?;
    let kinds: HashMap<_, _> = model
        .entities
        .iter()
        .map(|e| (e.id.as_str(), e.kind))
        .collect();
    let requested: HashSet<_> = atom_rows.unwrap_or_default().iter().copied().collect();
    if atom_rows.is_some_and(|rows| rows.len() != requested.len()) {
        return Err("atom_rows must contain distinct source rows".into());
    }
    let mut atoms = Atoms {
        coordinates: Vec::new(),
        radii: Vec::new(),
        rows: Vec::new(),
        selected: Vec::new(),
    };
    let mut alternates = HashMap::new();
    let mut found = false;
    for (row, atom) in model
        .atom_sites
        .iter()
        .enumerate()
        .filter(|(_, a)| a.model_number == number)
    {
        found = true;
        if !atom.occupancy.is_finite() || !(0.0..=1.0).contains(&atom.occupancy) {
            return Err(format!("Invalid occupancy in atom_site row {row}"));
        }
        if atom.occupancy == 0.0
            || ["H", "D", "T"]
                .iter()
                .any(|s| atom.element.eq_ignore_ascii_case(s))
        {
            continue;
        }
        let kind = kinds.get(atom.label_entity_id.as_str()).copied();
        if kind == Some(EntityKind::Water)
            || ["HOH", "DOD", "WAT"]
                .iter()
                .any(|s| atom.label_comp_id.eq_ignore_ascii_case(s))
        {
            continue;
        }
        if let Some(label) = &atom.label_alt_id {
            check_alternate(
                &mut alternates,
                [
                    Cow::Borrowed(&atom.label_asym_id),
                    Cow::Borrowed(&atom.auth_seq_id),
                    Cow::Borrowed(atom.insertion_code.as_deref().unwrap_or_default()),
                ],
                Cow::Borrowed(label),
            )?;
        }
        let selected = if atom_rows.is_some() {
            requested.contains(&row)
        } else {
            match kind {
                Some(EntityKind::NonPolymer | EntityKind::Branched) => true,
                Some(EntityKind::Polymer) => false,
                _ => return Err("Automatic ligand SASA requires entity types and label_entity_id; supply atom_rows".into()),
            }
        };
        if selected {
            atoms.rows.push(row);
            atoms.selected.push(atoms.coordinates.len());
        }
        atoms.coordinates.push([atom.x, atom.y, atom.z]);
        atoms.radii.push(radius(&atom.element, radii)?);
    }
    if !found {
        return Err(format!("Deposited model {number} has no atoms"));
    }
    if atom_rows.is_some_and(|rows| rows.len() != atoms.rows.len()) {
        return Err(
            "atom_rows must select occupied, non-water heavy atoms in the chosen model".into(),
        );
    }
    let areas = sasa(
        &atoms.coordinates,
        &atoms.radii,
        &atoms.selected,
        probe_radius,
        point_count,
    )?;
    Ok((atoms.rows, areas))
}

//! Decode source sites and authored bond exclusions without reading CCD metadata.

use super::{
    ConformerAtom, ConformerPolicy, ConformerSelection, fail, select_with_links,
};
use crate::pdbx::{PdbxModel, SemanticError};
use crate::pdbx::{fields, model, source};
use std::collections::{BTreeMap, HashMap};

/// Select one deposited model without changing the source-backed model.
///
/// # Errors
///
/// Returns the same errors as [`select_atom_conformers`], or an absent model error.
pub fn select_conformer(
    model: &PdbxModel,
    number: Option<i64>,
    seed: u64,
) -> Result<ConformerSelection, SemanticError> {
    select_conformer_with_policy(model, number, seed, ConformerPolicy::Occupancy)
}

/// Select a deposited model with explicit alternative ordering.
///
/// # Errors
///
/// Returns the validation and search errors from [`select_conformer`].
pub fn select_conformer_with_policy(
    model: &PdbxModel,
    number: Option<i64>,
    seed: u64,
    policy: ConformerPolicy,
) -> Result<ConformerSelection, SemanticError> {
    if let Some(block) = model.document.blocks().first() {
        crate::pdbx::selected_document::check_tables(block.entries())?;
    }
    select_sites(
        &model.atom_sites,
        &model.connections,
        &model.document,
        number,
        seed,
        policy,
    )
}

/// Select directly from an atom_site document without requiring entity or CCD metadata.
///
/// # Errors
///
/// Rejects invalid coordinate fields and the same conflicts as [`select_conformer`].
pub fn select_conformer_document(
    document: &crate::cif::CifDocument,
    number: Option<i64>,
    seed: u64,
) -> Result<ConformerSelection, SemanticError> {
    select_conformer_document_with_policy(document, number, seed, ConformerPolicy::Occupancy)
}

/// Select source rows with an explicit alternative-ordering policy.
///
/// # Errors
///
/// Returns the validation and search errors from [`select_conformer_document`].
pub fn select_conformer_document_with_policy(
    document: &crate::cif::CifDocument,
    number: Option<i64>,
    seed: u64,
    policy: ConformerPolicy,
) -> Result<ConformerSelection, SemanticError> {
    let [block] = document.blocks() else {
        return Err(fail("Conformer selection requires one data block"));
    };
    let categories = crate::cif::CategoryView::new(block.entries());
    crate::pdbx::selected_document::check_tables(block.entries())?;
    let atoms = source::atom_sites(categories)?;
    let connections = source::connections(categories)?;
    select_sites(&atoms, &connections, document, number, seed, policy)
}

fn select_sites(
    atom_sites: &[model::AtomSite],
    connections: &[model::Connection],
    document: &crate::cif::CifDocument,
    number: Option<i64>,
    seed: u64,
    policy: ConformerPolicy,
) -> Result<ConformerSelection, SemanticError> {
    let number = number
        .or_else(|| atom_sites.iter().map(|a| a.model_number).min())
        .ok_or_else(|| fail("The model has no atom sites"))?;
    if !atom_sites
        .iter()
        .any(|a| a.model_number == number && a.label_alt_id.is_some())
    {
        let rows: Vec<_> = atom_sites
            .iter()
            .enumerate()
            .filter_map(|(i, a)| (a.model_number == number).then_some(i))
            .collect();
        if rows.is_empty() {
            return Err(fail(format!("Deposited model {number} has no atoms")));
        }
        return Ok(ConformerSelection {
            schema_version: policy.schema(),
            atom_rows: rows,
            assignments: BTreeMap::new(),
            search_nodes: 0,
            seed,
            model_number: Some(number),
        });
    }
    let mut residues = BTreeMap::new();
    let mut atoms = Vec::new();
    let mut rows = Vec::new();
    for (row, a) in atom_sites
        .iter()
        .enumerate()
        .filter(|(_, a)| a.model_number == number)
    {
        let key = (
            &a.label_asym_id,
            a.label_seq_id,
            &a.auth_seq_id,
            &a.insertion_code,
        );
        let count = residues.len();
        let residue = *residues.entry(key).or_insert(count);
        atoms.push(ConformerAtom {
            residue,
            chain: a.label_asym_id.to_string(),
            sequence: a.label_seq_id,
            atom_name: a.label_atom_id.as_str().to_owned(),
            alternate: a.label_alt_id.as_deref().map(str::to_owned),
            coordinates: [a.x, a.y, a.z],
            occupancy: a.occupancy,
            clash_eligible: !["H", "D", "T"]
                .iter()
                .any(|element| a.element.eq_ignore_ascii_case(element)),
        });
        rows.push(row);
    }
    if atoms.is_empty() {
        return Err(fail(format!("Deposited model {number} has no atoms")));
    }
    let mut index: HashMap<(&str, &str), Vec<usize>> = HashMap::new();
    for (i, &r) in rows.iter().enumerate() {
        let a = &atom_sites[r];
        index
            .entry((&a.label_asym_id, &a.label_atom_id))
            .or_default()
            .push(i);
    }
    let mut qualifiers = HashMap::new();
    if let Some(block) = document.blocks().first() {
        for row in fields::category_rows(
            crate::cif::CategoryView::new(block.entries()),
            "struct_conn",
            [
                "id",
                "pdbx_ptnr1_label_alt_id",
                "pdbx_ptnr2_label_alt_id",
                "ptnr1_auth_seq_id",
                "ptnr2_auth_seq_id",
                "ptnr1_symmetry",
                "ptnr2_symmetry",
                "pdbx_ptnr1_PDB_ins_code",
                "pdbx_ptnr2_PDB_ins_code",
            ],
        ) {
            if let Some(id) = fields::optional_text(&row, 0) {
                qualifiers.insert(
                    id,
                    [
                        fields::optional_text(&row, 1),
                        fields::optional_text(&row, 2),
                        fields::optional_text(&row, 3),
                        fields::optional_text(&row, 4),
                        fields::optional_text(&row, 5),
                        fields::optional_text(&row, 6),
                        fields::optional_text(&row, 7),
                        fields::optional_text(&row, 8),
                    ],
                );
            }
        }
    }
    let mut excluded = Vec::new();
    let mut linked = Vec::new();
    for c in connections
        .iter()
        .filter(|c| matches!(c.kind.as_str(), "covale" | "disulf" | "metalc"))
    {
        let q = qualifiers.get(&c.id);
        let image = |i: usize| q.and_then(|q| q[i].as_deref()).unwrap_or("1_555");
        if image(4) != image(5) {
            continue;
        }
        let matches = |endpoint: &model::ConnectionEndpoint, partner: usize| {
            let alt = q.and_then(|q| q[partner].as_deref());
            let author = q.and_then(|q| q[partner + 2].as_deref());
            let insertion = q.and_then(|q| q[partner + 6].as_deref());
            index
                .get(&(endpoint.asym_id.as_str(), endpoint.atom_id.as_str()))
                .into_iter()
                .flatten()
                .filter_map(|&i| {
                    let r = rows[i];
                    let a = &atom_sites[r];
                    (a.label_asym_id == endpoint.asym_id
                        && a.label_atom_id == endpoint.atom_id
                        && endpoint
                            .sequence_id
                            .is_none_or(|s| Some(s) == a.label_seq_id)
                        && endpoint
                            .component_id
                            .as_ref()
                            .is_none_or(|s| s == &a.label_comp_id)
                        && author.is_none_or(|s| s == a.auth_seq_id)
                        && insertion.is_none_or(|s| a.insertion_code.as_deref() == Some(s))
                        && alt.is_none_or(|s| a.label_alt_id.as_deref().is_none_or(|v| v == s)))
                    .then_some(i)
                })
                .collect::<Vec<_>>()
        };
        for a in matches(&c.first, 0) {
            for b in matches(&c.second, 1) {
                excluded.push([a, b]);
                if q.is_some_and(|values| values[0].is_some() && values[1].is_some())
                    && atoms[a].alternate.is_some() && atoms[b].alternate.is_some() {
                    linked.push([a, b]);
                }
            }
        }
    }
    let mut selected = select_with_links(&atoms, &excluded, &linked, seed, policy)?;
    let mut identities = HashMap::new();
    for &i in &selected.atom_rows {
        let component = atom_sites[rows[i]].label_comp_id.as_str();
        if identities
            .insert(atoms[i].residue, component)
            .is_some_and(|old| old != component)
        {
            return Err(fail(
                "Selected shared and alternate atoms have conflicting component identities",
            ));
        }
    }
    selected.model_number = Some(number);
    selected.atom_rows = selected
        .atom_rows
        .into_iter()
        .filter_map(|i| rows.get(i).copied())
        .collect();
    Ok(selected)
}

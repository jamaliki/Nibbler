//! Apply source-row diagnostics to only the affected training targets.

use super::{Schema, Tokens, record, source};
use crate::cif::CifDocument;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};

#[cfg(test)]
mod tests;

/// A stable heavy-atom key and the targets that lack source support.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct TargetUncertainty {
    pub parent_key: [String; 5],
    pub hydrogens_unknown: bool,
    pub chemical_unknown: bool,
    #[serde(default)]
    pub reasons: Vec<String>,
}

#[cfg(feature = "reduce3")]
pub(crate) fn from_output(document: &CifDocument, model: i64,
    output: &reduce3::pipeline::Output, required: Option<&[u32]>, allow_bonds: bool) -> Result<Vec<TargetUncertainty>, String> {
    let mut rows = BTreeMap::new();
    for site in &output.hydrogen_uncertainties {
        // An unobserved bridge coordinate does not change its known chemistry.
        if site.parent_sources.len() != 1 { continue; }
        use reduce3::pipeline::HydrogenUncertaintyReason as Reason;
        let reason = match site.reason {
            Reason::MissingHeavyNeighbor => "missing_heavy_neighbor",
            Reason::InsufficientRestraints => "insufficient_restraints",
            Reason::MissingObservedBridge => "missing_observed_bridge",
        };
        report(&mut rows, site.parent_sources[0], true, reason);
    }
    for site in output.bond_uncertainties.iter().filter(|_| allow_bonds) {
        use reduce3::interp::target_uncertainty::BondUncertaintyReason as Reason;
        let reason = match site.reason {
            Reason::AmbiguousEndpoint => "ambiguous_bond_endpoint",
            Reason::UninterpretedBond => "uninterpreted_bond",
            Reason::MultipleSulfurPartners => "multiple_sulfur_partners",
            Reason::UnsupportedResidueEnvironment => "unsupported_residue_environment",
            Reason::DegenerateHeavyBond => "degenerate_heavy_bond",
        };
        report(&mut rows, site.parent_source, site.hydrogens_unknown, reason);
    }
    if let Some(required) = required {
        let required = required.iter().copied().collect::<std::collections::BTreeSet<_>>();
        rows.retain(|row, _| required.contains(row));
    }
    if rows.is_empty() { return Ok(Vec::new()); }
    let source = source::read(document, model)?;
    let atoms = source.atoms.iter().map(|atom| (atom.source_row, atom)).collect::<HashMap<_, _>>();
    rows.into_iter().map(|(row, (hydrogens_unknown, reasons))| {
        let atom = atoms.get(&row).ok_or("An uncertain target lacks its original heavy source row")?;
        let mut parent_key = record::token_key(atom, true);
        // Protein target keys use label sequence rather than auth sequence.
        if atom.protein { parent_key[1] = atom.sequence.to_string(); }
        Ok(TargetUncertainty { parent_key, hydrogens_unknown, chemical_unknown: true,
            reasons: reasons.into_iter().map(str::to_owned).collect() })
    }).collect()
}

#[cfg(feature = "reduce3")]
fn report(rows: &mut BTreeMap<u32, (bool, std::collections::BTreeSet<&'static str>)>,
    row: u32, hydrogens_unknown: bool, reason: &'static str) {
    let entry = rows.entry(row).or_default();
    entry.0 |= hydrogens_unknown;
    entry.1.insert(reason);
}

/// Preserve observed coordinates and clear only uncertain target labels.
pub(crate) fn apply(tokens: &mut Tokens, schema: &Schema, values: &[TargetUncertainty]) -> Result<(), String> {
    if values.is_empty() { return Ok(()); }
    let lookup = tokens.keys.iter().enumerate().map(|(row, key)| (key.clone(), row)).collect::<HashMap<_, _>>();
    for value in values {
        let mut key = value.parent_key.clone();
        let parent_name = key[4].clone();
        let ligand_row = lookup.get(&key).copied();
        key[4] = "protein".into();
        let row = ligand_row.or_else(|| lookup.get(&key).copied())
            .ok_or("An uncertain target parent is absent from packed tokens")?;
        if tokens.ligand[row] != 0 {
            if value.chemical_unknown { tokens.charge_known[row] = 0; }
            if value.hydrogens_unknown { tokens.hydrogen_known[row] = 0; }
        } else {
            apply_protein(tokens, schema, row, &parent_name, value)?;
        }
    }
    Ok(())
}

fn apply_protein(tokens: &mut Tokens, schema: &Schema, row: usize, parent_name: &str,
    value: &TargetUncertainty) -> Result<(), String> {
    let aa = tokens.identity[row] as usize;
    if aa >= 20 { return Ok(()); }
    let residue = &schema.residues[aa];
    let parent = residue.names[..15].iter().position(|name| name == parent_name)
        .ok_or("An uncertain protein parent is absent from the slot schema")?;
    if value.chemical_unknown {
        tokens.state_known[row] = 0;
        let flag = match parent_name { "N" => Some(0), "C" | "OXT" => Some(1), "SG" => Some(2), _ => None };
        if let Some(flag) = flag { tokens.topology[row * 3 + flag] = -1; }
    }
    if value.hydrogens_unknown || value.chemical_unknown {
        for slot in 15..32 {
            if residue.parents[slot] == parent as i64 && tokens.masks[row * 32 + slot] == 0 {
                tokens.presence[row * 32 + slot] = 0;
            }
        }
    }
    Ok(())
}

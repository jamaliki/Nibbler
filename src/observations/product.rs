//! Apply occurrence-specific product targets before selecting complete tokens.

use super::{Schema, Tokens};
use serde::Deserialize;
use std::collections::HashMap;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ProductTarget {
    parent_key: [String; 5],
    sequence_id: Option<i64>,
    formal_charge: Option<i16>,
    hybridization: i64,
    hydrogens: Option<usize>,
    complete: bool,
    modified: bool,
}

pub(crate) fn apply(
    tokens: &mut Tokens,
    schema: &Schema,
    values: &[ProductTarget],
) -> Result<(), String> {
    let normalize = |mut key: [String; 5]| {
        if matches!(key[2].as_str(), "." | "?") {
            key[2].clear();
        }
        key
    };
    let lookup = tokens
        .keys
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, key)| (normalize(key), i))
        .collect::<HashMap<_, _>>();
    let prior = tokens.hydrogen_known.clone();
    for row in 0..tokens.identity.len() {
        if tokens.ligand[row] != 0 && tokens.identity[row] != 65 {
            tokens.hydrogen_known[row] = 0;
            tokens.charge_known[row] = 0;
        }
    }
    for value in values {
        if let Some(&row) = lookup.get(&normalize(value.parent_key.clone())) {
            let count = tokens.masks[row * 32 + 15..(row + 1) * 32]
                .iter()
                .filter(|&&present| present != 0)
                .count();
            tokens.charges[row] = value.formal_charge.map_or(0., f32::from);
            tokens.hybridization[row] = value.hybridization;
            let known = prior[row] != 0 && value.complete && value.hydrogens == Some(count);
            tokens.hydrogen_known[row] = u8::from(known);
            tokens.charge_known[row] = u8::from(known && value.formal_charge.is_some());
        } else if value.modified && value.parent_key[3] != "CYS" {
            let mut key = value.parent_key.clone();
            key[1] = value
                .sequence_id
                .map_or_else(|| key[1].clone(), |value| value.to_string());
            key[4] = "protein".into();
            if let Some(&row) = lookup.get(&normalize(key)) {
                tokens.state_known[row] = 0;
                let identity = tokens.identity[row] as usize;
                if let Some(residue) = schema.residues.get(identity) {
                    for slot in 15..32 {
                        if residue.parents[slot] >= 0 {
                            tokens.presence[row * 32 + slot] = tokens.masks[row * 32 + slot];
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

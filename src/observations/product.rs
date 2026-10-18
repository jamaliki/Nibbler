//! Apply occurrence-specific product targets before selecting complete tokens.

use super::record::same_residue as same_key_residue;
use super::{Schema, Tokens};
use serde::Deserialize;
use crate::hash::FastMap;
use crate::pdbx::ShortText;

/// One product target; native callers borrow the key text from their view.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ProductTarget<K = String> {
    parent_key: [K; 5],
    sequence_id: Option<i64>,
    formal_charge: Option<i16>,
    hybridization: i64,
    hydrogens: Option<usize>,
    complete: bool,
    modified: bool,
}

impl<K> ProductTarget<K> {
    #[cfg(feature = "reduce3")]
    pub(crate) fn new(parent_key: [K; 5], sequence_id: Option<i64>, formal_charge: Option<i32>,
        hybridization: i64, hydrogens: Option<usize>, complete: bool, modified: bool) -> Self {
        Self { parent_key, sequence_id, formal_charge: formal_charge.and_then(|value| i16::try_from(value).ok()),
            hybridization, hydrogens, complete, modified }
    }
}

/// The first four key fields, where "." or "?" insertion codes match an empty one.
fn residue<K: AsRef<str>>(key: &[K; 5]) -> [&str; 4] {
    let mut prefix = [&key[0], &key[1], &key[2], &key[3]].map(AsRef::as_ref);
    if matches!(prefix[2], "." | "?") {
        prefix[2] = "";
    }
    prefix
}

/// Whether two residue parts are equal.
///
/// Callers may borrow one residue's fields from one place, and then equal addresses
/// settle the test without reading the text.
fn same_residue(left: &[&str; 4], right: &[&str; 4]) -> bool {
    left.iter().zip(right).all(|(left, right)| std::ptr::eq(*left, *right) || left == right)
}

/// Token rows by key, found through the residue part of the key and then the name.
///
/// Tokens of one residue are mostly adjacent, and targets of one residue arrive together
/// in token order. So most tokens and targets reuse the previous residue without hashing,
/// and most names match at the row after the last match.
struct Lookup<'k> {
    keys: &'k [[ShortText; 5]],
    residues: FastMap<[&'k str; 4], usize>,
    /// Rows of residue `r` are `rows[starts[r]..starts[r + 1]]`, in token order.
    starts: Vec<usize>,
    rows: Vec<usize>,
    /// Whether the names of a residue's rows differ from each other.
    distinct: Vec<bool>,
    /// The residue of the last match and the position after it.
    hint: (usize, usize),
}

impl<'k> Lookup<'k> {
    fn new(keys: &'k [[ShortText; 5]]) -> Self {
        // Count the runs first, so that the map never grows.
        let same: Vec<bool> = (0..keys.len()).map(|row| row > 0 && same_key_residue(&keys[row], &keys[row - 1]))
            .collect();
        let runs = same.iter().filter(|&&same| !same).count();
        let mut residues = FastMap::with_capacity_and_hasher(runs, Default::default());
        let mut owner = Vec::with_capacity(keys.len());
        for (row, key) in keys.iter().enumerate() {
            let id = match owner.last() {
                Some(&id) if same[row] => id,
                _ => {
                    let next = residues.len();
                    *residues.entry(residue(key)).or_insert(next)
                }
            };
            owner.push(id);
        }
        let mut starts = vec![0; residues.len() + 1];
        for &id in &owner {
            starts[id + 1] += 1;
        }
        for id in 0..residues.len() {
            starts[id + 1] += starts[id];
        }
        let mut filled = starts.clone();
        let mut rows = vec![0; keys.len()];
        for (row, &id) in owner.iter().enumerate() {
            rows[filled[id]] = row;
            filled[id] += 1;
        }
        let names = |range: &[usize]| rows[range[0]..range[1]].iter().map(|&row| &keys[row][4]);
        // Copies of one residue type list the same names in the same order, so a residue
        // whose names equal those of a recent residue takes its answer and skips the test
        // of every pair. Short residues take the pair test, which costs little for them.
        let mut recent: Vec<(&[usize], bool)> = Vec::new();
        let distinct = starts.windows(2).map(|range| {
            let members = &rows[range[0]..range[1]];
            let long = members.len() >= 8;
            if long && let Some(&(_, value)) = recent.iter()
                .find(|(other, _)| other[1] - other[0] == members.len() && names(other).eq(names(range)))
            {
                return value;
            }
            let value = members.iter().enumerate()
                .all(|(index, &row)| members[index + 1..].iter().all(|&other| keys[other][4] != keys[row][4]));
            if long {
                if recent.len() == 16 {
                    recent.pop();
                }
                recent.push((range, value));
            }
            value
        }).collect();
        Self { keys, residues, starts, rows, distinct, hint: (usize::MAX, 0) }
    }

    /// The residue whose keys start with `prefix`.
    fn residue(&self, prefix: [&str; 4]) -> Option<usize> {
        let residues: &FastMap<[&str; 4], usize> = &self.residues;
        residues.get(&prefix).copied()
    }

    /// The last row of residue `id` whose name equals `name`.
    fn find(&mut self, id: usize, name: &str) -> Option<usize> {
        let hint = if self.hint.0 == id { self.hint.1 } else { 0 };
        let members = &self.rows[self.starts[id]..self.starts[id + 1]];
        let name = ShortText::new(name);
        let found = if self.distinct[id] {
            let (head, tail) = members.split_at(hint.min(members.len()));
            tail.iter().chain(head).position(|&row| self.keys[row][4] == name)
                .map(|position| (position + hint) % members.len())
        } else {
            members.iter().rposition(|&row| self.keys[row][4] == name)
        };
        let position = found?;
        self.hint = (id, position + 1);
        Some(members[position])
    }
}

pub(crate) fn apply<K: AsRef<str>>(
    tokens: &mut Tokens,
    schema: &Schema,
    values: &[ProductTarget<K>],
) -> Result<(), String> {
    let mut lookup = Lookup::new(&tokens.keys);
    // Targets of one residue arrive together, so most reuse the residue found before.
    let mut current: Option<([&str; 4], Option<usize>)> = None;
    let prior = tokens.hydrogen_known.clone();
    for row in 0..tokens.identity.len() {
        if tokens.ligand[row] != 0 && tokens.identity[row] != 65 {
            tokens.hydrogen_known[row] = 0;
            tokens.charge_known[row] = 0;
        }
    }
    for value in values {
        let prefix = residue(&value.parent_key);
        let id = match current {
            Some((previous, id)) if same_residue(&previous, &prefix) => id,
            _ => {
                let id = lookup.residue(prefix);
                current = Some((prefix, id));
                id
            }
        };
        if let Some(row) = id.and_then(|id| lookup.find(id, value.parent_key[4].as_ref())) {
            let count = tokens.masks[row * 32 + 15..(row + 1) * 32]
                .iter()
                .filter(|&&present| present != 0)
                .count();
            tokens.charges[row] = value.formal_charge.map_or(0., f32::from);
            tokens.hybridization[row] = value.hybridization;
            let known = prior[row] != 0 && value.complete && value.hydrogens == Some(count);
            tokens.hydrogen_known[row] = u8::from(known);
            tokens.charge_known[row] = u8::from(known && value.formal_charge.is_some());
        } else if value.modified && value.parent_key[3].as_ref() != "CYS" {
            let sequence = value.sequence_id.map(|value| value.to_string());
            let mut key = residue(&value.parent_key);
            if let Some(sequence) = &sequence {
                key[1] = sequence;
            }
            if let Some(row) = lookup.residue(key).and_then(|id| lookup.find(id, "protein")) {
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

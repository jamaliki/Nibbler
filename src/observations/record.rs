//! Owned dense arrays, complete-residue selection and identity alignment.

use super::source::Atom;
use super::{Schema, TokenTable, residue};
use crate::pdbx::ShortText;
use std::collections::{BTreeSet, HashMap};

pub(super) fn token_key(atom: &Atom, ligand: bool) -> [String; 5] {
    [
        atom.chain.as_str().to_owned(),
        if ligand {
            atom.auth_sequence.as_str().to_owned()
        } else {
            atom.sequence.to_string()
        },
        atom.insertion.as_str().to_owned(),
        atom.component.as_str().to_owned(),
        if ligand {
            atom.name.as_str().to_owned()
        } else {
            "protein".into()
        },
    ]
}

/// Whether two token keys share their residue part, the first four fields. A slice compare
/// calls a loop that is not inlined; four field compares inline to fixed-size compares.
#[inline]
pub(super) fn same_residue(left: &[ShortText; 5], right: &[ShortText; 5]) -> bool {
    left[0] == right[0] && left[1] == right[1] && left[2] == right[2] && left[3] == right[3]
}

/// The token key of `token_key`, held as inline text so that it needs no heap allocation.
fn short_key(atom: &Atom, ligand: bool) -> [ShortText; 5] {
    [
        atom.chain.clone(),
        if ligand {
            atom.auth_sequence.clone()
        } else {
            ShortText::new(&atom.sequence.to_string())
        },
        atom.insertion.clone(),
        atom.component.clone(),
        if ligand {
            atom.name.clone()
        } else {
            ShortText::new("protein")
        },
    ]
}

pub(crate) struct Tokens {
    pub positions: Vec<f32>,
    pub masks: Vec<u8>,
    pub presence: Vec<u8>,
    pub identity: Vec<i64>,
    pub sequence: Vec<i64>,
    pub polymer_segments: Vec<i64>,
    pub chains: Vec<i64>,
    pub ligand: Vec<u8>,
    pub charges: Vec<f32>,
    pub charge_known: Vec<u8>,
    pub hydrogen_known: Vec<u8>,
    pub hybridization: Vec<i64>,
    pub covalent: Vec<u8>,
    pub state: Vec<i64>,
    pub state_known: Vec<u8>,
    pub topology: Vec<i64>,
    pub attachments: Vec<(usize, usize)>,
    pub bridge_edges: Vec<(usize, usize, u8)>,
    pub crop_positions: Vec<f32>,
    pub crop_center_eligible: Vec<u8>,
    pub keys: Vec<[ShortText; 5]>,
}

impl Tokens {
    pub fn empty() -> Self {
        Self {
            positions: Vec::new(),
            masks: Vec::new(),
            presence: Vec::new(),
            identity: Vec::new(),
            sequence: Vec::new(),
            polymer_segments: Vec::new(),
            chains: Vec::new(),
            ligand: Vec::new(),
            charges: Vec::new(),
            charge_known: Vec::new(),
            hydrogen_known: Vec::new(),
            hybridization: Vec::new(),
            covalent: Vec::new(),
            state: Vec::new(),
            state_known: Vec::new(),
            topology: Vec::new(),
            attachments: Vec::new(),
            bridge_edges: Vec::new(),
            crop_positions: Vec::new(),
            crop_center_eligible: Vec::new(),
            keys: Vec::new(),
        }
    }

    /// Reserve room for `additional` more tokens.
    pub(super) fn reserve(&mut self, additional: usize) {
        self.positions.reserve(additional * 96);
        self.masks.reserve(additional * 32);
        self.presence.reserve(additional * 32);
        self.identity.reserve(additional);
        self.sequence.reserve(additional);
        self.polymer_segments.reserve(additional);
        self.chains.reserve(additional);
        self.ligand.reserve(additional);
        self.charges.reserve(additional);
        self.charge_known.reserve(additional);
        self.hydrogen_known.reserve(additional);
        self.hybridization.reserve(additional);
        self.covalent.reserve(additional);
        self.state.reserve(additional);
        self.state_known.reserve(additional);
        self.topology.reserve(additional * 3);
        self.keys.reserve(additional);
        self.crop_positions.reserve(additional * 3);
        self.crop_center_eligible.reserve(additional);
    }

    pub(super) fn push(
        &mut self,
        atom: &Atom,
        identity: usize,
        packed: &residue::Packed,
        ligand: bool,
        features: [i16; 4],
    ) {
        self.positions.extend_from_slice(packed.position.as_flattened());
        self.masks.extend(packed.mask.map(u8::from));
        self.presence.extend(packed.presence.map(u8::from));
        self.state.push(packed.state);
        self.state_known.push(u8::from(packed.state_known));
        self.topology.extend_from_slice(&packed.topology);
        self.crop_positions.extend_from_slice(&packed.position[0]);
        self.push_atom(atom, identity, ligand, features);
    }

    /// Push one ligand atom as its own token. The rows equal those `push` writes for
    /// `Packed::empty([-1; 3])` with the atom in slot 0, but skip building the 32-slot copy.
    pub(super) fn push_ligand(&mut self, atom: &Atom, identity: usize, features: [i16; 4]) {
        self.positions.extend_from_slice(&atom.position);
        self.positions.extend_from_slice(&[0.; 93]);
        self.masks.push(1);
        self.masks.extend_from_slice(&[0; 31]);
        self.presence.extend_from_slice(&[0; 32]);
        self.state.push(0);
        self.state_known.push(0);
        self.topology.extend_from_slice(&[-1; 3]);
        self.crop_positions.extend_from_slice(&atom.position);
        self.push_atom(atom, identity, true, features);
    }

    /// Push the fields that do not depend on the packed slots.
    fn push_atom(&mut self, atom: &Atom, identity: usize, ligand: bool, features: [i16; 4]) {
        self.identity.push(identity as i64);
        self.sequence.push(atom.sequence);
        self.polymer_segments.push(0);
        self.chains.push(0);
        self.ligand.push(u8::from(ligand));
        self.charges.push(f32::from(features[0]));
        self.charge_known.push(0);
        self.hydrogen_known.push(0);
        self.hybridization.push(i64::from(features[1]));
        self.covalent.push(0);
        self.keys.push(short_key(atom, ligand));
        self.crop_center_eligible.push(1);
    }

    pub fn select(&mut self, indices: &[usize]) -> Result<(), String> {
        fn selected<T: Clone>(values: &[T], rows: &[usize], width: usize) -> Vec<T> {
            rows.iter()
                .flat_map(|&row| values[row * width..(row + 1) * width].iter().cloned())
                .collect()
        }
        let mut inverse = vec![usize::MAX; self.identity.len()];
        for (index, &row) in indices.iter().enumerate() {
            if row >= inverse.len() || inverse[row] != usize::MAX {
                return Err("Token selection indices are absent or repeated".into());
            }
            inverse[row] = index;
        }
        for &(hydrogen, metal, _) in &self.bridge_edges {
            if inverse[hydrogen] != usize::MAX && inverse[metal] == usize::MAX {
                return Err("Token selection retains a bridge hydrogen without both metals".into());
            }
        }
        self.bridge_edges = self.bridge_edges.iter().filter_map(|&(hydrogen, metal, kind)| {
            (inverse[hydrogen] != usize::MAX).then_some((inverse[hydrogen], inverse[metal], kind))
        }).collect();
        self.attachments = self
            .attachments
            .iter()
            .filter_map(|&(a, b)| {
                (inverse[a] != usize::MAX && inverse[b] != usize::MAX)
                    .then_some((inverse[a], inverse[b]))
            })
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        self.positions = selected(&self.positions, indices, 96);
        self.masks = selected(&self.masks, indices, 32);
        self.presence = selected(&self.presence, indices, 32);
        self.topology = selected(&self.topology, indices, 3);
        self.identity = selected(&self.identity, indices, 1);
        self.sequence = selected(&self.sequence, indices, 1);
        self.polymer_segments = selected(&self.polymer_segments, indices, 1);
        self.chains = selected(&self.chains, indices, 1);
        self.ligand = selected(&self.ligand, indices, 1);
        self.charges = selected(&self.charges, indices, 1);
        self.charge_known = selected(&self.charge_known, indices, 1);
        self.hydrogen_known = selected(&self.hydrogen_known, indices, 1);
        self.hybridization = selected(&self.hybridization, indices, 1);
        self.state = selected(&self.state, indices, 1);
        self.state_known = selected(&self.state_known, indices, 1);
        self.keys = selected(&self.keys, indices, 1);
        self.crop_positions = selected(&self.crop_positions, indices, 3);
        self.crop_center_eligible = selected(&self.crop_center_eligible, indices, 1);
        self.covalent = vec![0; indices.len()];
        for &(a, _) in &self.attachments {
            self.covalent[a] = 1;
        }
        Ok(())
    }

    pub fn complete(&mut self, schema: &Schema) -> Result<(), String> {
        let indices: Vec<_> = (0..self.identity.len())
            .filter(|&row| {
                if self.ligand[row] != 0 {
                    return true;
                }
                let aa = self.identity[row] as usize;
                aa < 20
                    && (0..14).all(|slot| {
                        schema.residues[aa].names[slot].is_empty()
                            || self.masks[row * 32 + slot] != 0
                                && residue::finite(
                                    self.positions[row * 96 + slot * 3..][..3]
                                        .try_into()
                                        .unwrap_or([f32::NAN; 3]),
                                )
                    })
            })
            .collect();
        self.select(&indices)?;
        super::continuity::retained_proteins(self);
        Ok(())
    }

    pub fn protein_gaps(&mut self) {
        super::continuity::retained_proteins(self);
    }

    pub fn align(&mut self, expected: &str) -> Result<(), String> {
        let wanted: Vec<[String; 5]> =
            serde_json::from_str(expected).map_err(|error| error.to_string())?;
        let lookup: HashMap<_, _> = self
            .keys
            .iter()
            .enumerate()
            .map(|(row, key)| (key.each_ref().map(ShortText::as_str), row))
            .collect();
        let unique: BTreeSet<_> = wanted.iter().collect();
        if lookup.len() != self.keys.len()
            || unique.len() != wanted.len()
            || wanted.len() != self.keys.len()
            || wanted.iter().any(|key| !lookup.contains_key(&key.each_ref().map(String::as_str)))
        {
            return Err("PDB token identities differ from the heavy index".into());
        }
        let indices: Vec<_> = wanted.iter().map(|key| lookup[&key.each_ref().map(String::as_str)]).collect();
        self.select(&indices)?;
        renumber_ligands(self);
        Ok(())
    }

    pub fn validate(&self, schema: &Schema, table: &TokenTable) -> Result<(), String> {
        let mut bridges: HashMap<usize, BTreeSet<usize>> = HashMap::new();
        for &(hydrogen, metal, kind) in &self.bridge_edges {
            if hydrogen >= self.identity.len() || metal >= self.identity.len()
                || self.identity[hydrogen] != table.hydrogen() || self.ligand[hydrogen] == 0
                || self.ligand[metal] == 0 || hydrogen == metal || kind != 1 {
                return Err("Ligand bridge edges have invalid endpoints or type".into());
            }
            if !bridges.entry(hydrogen).or_default().insert(metal) {
                return Err("Ligand bridge edges contain a repeated metal endpoint".into());
            }
        }
        for (row, &identity) in self.identity.iter().enumerate() {
            if identity == table.hydrogen() && self.ligand[row] != 0
                && bridges.get(&row).is_none_or(|metals| metals.len() != 2) {
                return Err("Ligand bridge hydrogen requires exactly two metal edges".into());
            }
        }
        for row in 0..self.identity.len() {
            let aa = self.identity[row] as usize;
            let mask = &self.masks[row * 32..(row + 1) * 32];
            let known = &self.presence[row * 32..(row + 1) * 32];
            let topology = &self.topology[row * 3..(row + 1) * 3];
            // The test has no branch per slot, so the compiler can check a row in vector steps.
            let finite = mask.iter().zip(self.positions[row * 96..(row + 1) * 96].chunks_exact(3))
                .fold(true, |finite, (&observed, xyz)| {
                    finite & ((observed == 0) | xyz.iter().fold(true, |all, value| all & value.is_finite()))
                });
            if !finite {
                return Err("Protonated coordinates contain non-finite observed atoms".into());
            }
            if aa >= 20 {
                continue;
            }
            let allowed = |state: &super::schema::State| {
                state
                    .topology()
                    .iter()
                    .zip(topology)
                    .all(|(&expected, &value)| value < 0 || value == expected)
            };
            if known[15..].iter().any(|&value| value != 0)
                && !schema.residues[aa].states.iter().any(|state| {
                    allowed(state)
                        && state.present.iter().zip(mask).zip(known).all(
                            |((&present, &observed), &known)| {
                                known == 0 || present == (observed != 0)
                            },
                        )
                })
            {
                let detail = serde_json::json!({"row":row,"aatype":aa,"token_key":self.keys[row],
                    "topology":topology,"mask":mask,"atom_presence_known":known});
                return Err(format!(
                    "Atom-presence labels contradict every allowed chemical state: {detail}"
                ));
            }
            if self.state_known[row] != 0 {
                let state = &schema.residues[aa].states[self.state[row] as usize];
                if !allowed(state)
                    || state
                        .present
                        .iter()
                        .zip(mask)
                        .any(|(&present, &observed)| present != (observed != 0))
                {
                    return Err(
                        "Known chemical state disagrees with observed atoms or topology".into(),
                    );
                }
            }
        }
        Ok(())
    }
}

pub(super) fn renumber_ligands(tokens: &mut Tokens) {
    let mut counters = HashMap::new();
    for row in 0..tokens.identity.len() {
        if tokens.ligand[row] != 0 {
            let (next, segment) = counters.entry(tokens.chains[row])
                .or_insert((0, tokens.polymer_segments[row]));
            if *segment != tokens.polymer_segments[row] {
                *next += super::continuity::INDEX_SEPARATION - 1;
                *segment = tokens.polymer_segments[row];
            }
            tokens.sequence[row] = *next;
            *next += 1;
        }
    }
}

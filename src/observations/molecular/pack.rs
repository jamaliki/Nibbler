//! Typed Atom84 rows: protein residues, nucleotide residues and ligand atoms.
//!
//! Each step mirrors `pack_molecular` in kaveh-novo, including its claim order,
//! so source ownership and exclusion reasons match the Python packer exactly.

use super::arrays::{self, Array};
use super::json;
use super::snapshot::Snapshot;
use super::{Failure, MolecularSchema, SLOTS};
use crate::observations::{TokenTable, Tokens};
use crate::observations::record::same_residue;
use crate::observations::source::Source;
use crate::hash::{FastMap, FastSet};
use crate::pdbx::ShortText;
use std::borrow::Cow;
use std::collections::{BTreeMap, HashSet};
use std::fmt::Write;

const PROTEIN_BACKBONE: [&str; 4] = ["N", "CA", "C", "O"];
const CANONICAL_NUCLEOTIDES: [&str; 8] = ["DA", "DC", "DG", "DT", "A", "C", "G", "U"];

/// One typed token with the defaults of the Python row dictionary.
#[derive(Clone)]
pub(super) struct Row {
    pub pos: [[f32; 3]; SLOTS],
    pub mask: [bool; SLOTS],
    pub aatype: i64,
    pub res_idx: i64,
    pub chain_idx: i64,
    pub is_ligand: bool,
    pub chemical_state: i64,
    pub state_known: bool,
    pub topology: [i64; 3],
    pub sequence_target_known: bool,
    pub ligand_charge: f32,
    pub ligand_charge_known: bool,
    pub ligand_hydrogen_known: bool,
    pub ligand_hybridization: i64,
    pub covalently_bound_to_protein: bool,
    pub atom_presence_known: [bool; SLOTS],
    pub crop_positions: [f32; 3],
    pub crop_center_eligible: bool,
    pub token_kind: u8,
    pub nucleotide_identity: u8,
    pub token_atom_mask: [bool; SLOTS],
    pub nucleotide_atom_mask: [bool; SLOTS],
    pub nucleotide_atom_elements: [i16; SLOTS],
    pub ncaa_parent_token: i64,
    pub molecular_group: i64,
}

impl Default for Row {
    fn default() -> Self {
        Self {
            pos: [[0.; 3]; SLOTS],
            mask: [false; SLOTS],
            aatype: 20,
            res_idx: 0,
            chain_idx: 0,
            is_ligand: false,
            chemical_state: 0,
            state_known: false,
            topology: [-1; 3],
            sequence_target_known: true,
            ligand_charge: 0.,
            ligand_charge_known: false,
            ligand_hydrogen_known: false,
            ligand_hybridization: 0,
            covalently_bound_to_protein: false,
            atom_presence_known: [false; SLOTS],
            crop_positions: [0.; 3],
            crop_center_eligible: true,
            token_kind: 0,
            nucleotide_identity: 0,
            token_atom_mask: [true; SLOTS],
            nucleotide_atom_mask: [false; SLOTS],
            nucleotide_atom_elements: [0; SLOTS],
            ncaa_parent_token: -1,
            molecular_group: -1,
        }
    }
}

/// Audit counts in the key order of the Python packer.
pub(crate) struct Stats {
    pub tokens: usize,
    pub source_atoms: usize,
    pub unsupported_known_charge_targets: usize,
    pub unknown_sequence_tokens: usize,
    pub represented_atoms: usize,
    pub exclusions: Vec<(String, usize)>,
    pub nucleotide_residues: usize,
    pub unknown_source_atom_tokens: usize,
    pub ncaa_sidechain_tokens: usize,
    pub modified_nucleotide_atom_tokens: usize,
    pub fully_atomized_modified_nucleotides: usize,
}

/// Rows, source ownership and the edge targets that index them.
///
/// Keys and slot names borrow the legacy tokens, the source atoms and the schema.
pub(super) struct Packed<'a> {
    pub rows: Vec<Row>,
    pub keys: Vec<[&'a str; 5]>,
    pub owner: Vec<i64>,
    pub slot: Vec<i64>,
    pub reason: Vec<u8>,
    pub atom_sasa: Vec<[f32; SLOTS]>,
    pub atom_sasa_mask: Vec<[bool; SLOTS]>,
    pub attachment_pairs: Vec<[i64; 4]>,
    pub attachment_order_known: Vec<bool>,
    pub ligand_bonds: Vec<[i64; 2]>,
    pub ligand_bond_order: Vec<i8>,
    pub ligand_bond_known: Vec<bool>,
    pub attachment_known: Vec<bool>,
    pub protein_attachment_index: Vec<[i64; 2]>,
    pub legacy_token_index: Vec<i64>,
    /// Slot names per token, in the key order of the metadata JSON.
    pub names: Vec<(usize, Vec<&'a str>)>,
    pub bridge_index: Vec<[i64; 2]>,
    pub bridge_type: Vec<u8>,
    pub stats: Stats,
}

fn insertion(value: &str) -> &str {
    if matches!(value, "" | "." | "?") { "" } else { value }
}

/// Upper-case text, borrowed when it already is.
fn upper(text: &str) -> Cow<'_, str> {
    if text.bytes().any(|byte| !byte.is_ascii() || byte.is_ascii_lowercase()) {
        Cow::Owned(text.to_uppercase())
    } else {
        Cow::Borrowed(text)
    }
}

/// Whether the upper-cased element of `atom` equals `element`.
fn element_is(source: &Source, atom: usize, element: &str) -> bool {
    upper(&source.atoms[atom].element) == element
}

/// One residue identity as a single string; mmCIF values never contain NUL.
fn group_key(out: &mut String, chain: &str, sequence: &str, code: &str, component: &str) {
    out.clear();
    for (index, part) in [chain, sequence, insertion(code), component].into_iter().enumerate() {
        if index > 0 {
            out.push('\0');
        }
        out.push_str(part);
    }
}

/// Source residue groups in order of first appearance and their identity lookup.
struct Groups<'a> {
    order: Vec<usize>,
    /// Atoms of each group in source order, as the source holds them.
    atoms: &'a [Vec<usize>],
    /// Keys are inline text, so most need no heap allocation; `ShortText` hashes as `str`.
    lookup: FastMap<ShortText, usize>,
}

impl<'a> Groups<'a> {
    fn read(source: &'a Source, snapshot: &Snapshot) -> Result<Self, Failure> {
        let atoms = source.groups.as_slice();
        let mut seen = vec![false; atoms.len()];
        let mut order = Vec::new();
        // The snapshot copies each atom's group into a dense list, which is much faster to
        // read than the atoms.
        for &group in &snapshot.group {
            let group = group as usize;
            if !seen[group] {
                seen[group] = true;
                order.push(group);
            }
        }
        let mut lookup = FastMap::with_capacity_and_hasher(order.len() * 2, Default::default());
        let mut aliases = Vec::with_capacity(order.len());
        let mut number = String::new();
        let mut key = String::new();
        for &group in &order {
            let atom = &source.atoms[atoms[group][0]];
            let sequence = if snapshot.is_protein[atoms[group][0]] {
                number.clear();
                let _ = write!(number, "{}", snapshot.sequence[atoms[group][0]]);
                number.as_str()
            } else {
                atom.auth_sequence.as_str()
            };
            group_key(&mut key, &atom.chain, sequence, &atom.insertion, &atom.component);
            if lookup.get(key.as_str()).is_some_and(|&previous| previous != group) {
                return Err("Two source residue groups share a token identity".into());
            }
            lookup.insert(ShortText::new(&key), group);
            group_key(&mut key, &atom.chain, &atom.auth_sequence, &atom.insertion, &atom.component);
            aliases.push((ShortText::new(&key), group));
        }
        for (key, group) in aliases {
            lookup.entry(key).or_insert(group);
        }
        Ok(Self { order, atoms, lookup })
    }
}

/// The set bits of `bits`, lowest first.
fn slots(mut bits: u128) -> impl Iterator<Item = usize> {
    std::iter::from_fn(move || {
        (bits != 0).then(|| {
            let slot = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            slot
        })
    })
}

/// Source atom ownership with one owner per token slot.
struct Ownership {
    owner: Vec<i64>,
    slot: Vec<i64>,
    reason: Vec<u8>,
    target_known: Vec<bool>,
    claimed: Vec<u128>,
}

impl Ownership {
    fn new(snapshot: &Snapshot) -> Self {
        let count = snapshot.row.len();
        let first: Vec<bool> = snapshot.target_known.iter().zip(&snapshot.occupancy)
            .map(|(&known, &occupancy)| known && occupancy > 0.).collect();
        let target_known: Vec<bool> = (0..count).map(|atom| {
            let parent = snapshot.hydrogen_parent[atom];
            first[atom] && (parent < 0 || first[parent as usize])
        }).collect();
        let mut reason = vec![7_u8; count];
        for atom in 0..count {
            if !snapshot.mask[atom] {
                reason[atom] = 6;
            } else if snapshot.occupancy[atom] <= 0. {
                reason[atom] = 8;
            } else if !target_known[atom] {
                reason[atom] = 4;
            }
        }
        Self { owner: vec![-1; count], slot: vec![-1; count], reason, target_known, claimed: Vec::new() }
    }

    fn is_claimed(&self, token: usize, slot: usize) -> bool {
        self.claimed.get(token).is_some_and(|bits| bits >> slot & 1 == 1)
    }

    fn claim(&mut self, rows: &mut [Row], snapshot: &Snapshot, atom: usize, token: usize, slot: usize,
        add: bool) -> Result<(), Failure> {
        if !self.target_known[atom] {
            return Ok(());
        }
        let row = &mut rows[token];
        let position = snapshot.pos[atom];
        if add {
            row.pos[slot] = position;
            row.mask[slot] = true;
        }
        if !row.mask[slot] || row.pos[slot] != position {
            return Err("A source atom does not match its molecular target coordinate".into());
        }
        if self.is_claimed(token, slot) {
            return Err("A molecular atom slot already has a source owner".into());
        }
        self.owner[atom] = token as i64;
        self.slot[atom] = slot as i64;
        self.reason[atom] = 0;
        if self.claimed.len() <= token {
            self.claimed.resize(token + 1, 0);
        }
        self.claimed[token] |= 1 << slot;
        Ok(())
    }

    fn match_legacy(&mut self, rows: &mut [Row], snapshot: &Snapshot, source: &Source, names: &[String],
        atoms: &[usize], token: usize) -> Result<(), Failure> {
        // Slot sets are bit masks, visited from the lowest slot up.
        let candidates = (0..SLOTS).filter(|&slot| rows[token].mask[slot]).fold(0_u128, |bits, slot| bits | 1 << slot);
        for &atom in atoms {
            let position = snapshot.pos[atom];
            let name = &source.atoms[atom].name;
            if self.reason[atom] == 8 {
                for slot in slots(candidates) {
                    if rows[token].pos[slot] == position && names[slot] == *name {
                        rows[token].mask[slot] = false;
                    }
                }
            }
            if self.owner[atom] >= 0 || !self.target_known[atom] {
                continue;
            }
            let mut matches = 0_u128;
            for slot in slots(candidates) {
                if rows[token].pos[slot] == position && !self.is_claimed(token, slot) {
                    matches |= 1 << slot;
                }
            }
            // The name only chooses between slots that share the position.
            let chosen = if matches.count_ones() == 1 {
                matches
            } else {
                let named = slots(matches).filter(|&slot| names[slot] == *name)
                    .fold(0_u128, |bits, slot| bits | 1 << slot);
                if named.count_ones() == 1 { named } else { 0 }
            };
            if let Some(slot) = slots(chosen).next() {
                self.claim(rows, snapshot, atom, token, slot, false)?;
            }
        }
        Ok(())
    }
}

/// Expand one legacy Atom32 token by name, starting unsupplied protein slots at CA.
fn legacy_row(legacy: &Tokens, old: usize, molecular: &MolecularSchema, bridges: bool) -> Row {
    let ligand = legacy.ligand[old] != 0;
    let aatype = legacy.identity[old];
    let kind = if ligand { 21 } else { aatype.clamp(0, 20) as usize };
    let (index, defined) = (&molecular.index[kind], &molecular.defined[kind]);
    let values = &legacy.positions[old * 96..old * 96 + 96];
    let masks = &legacy.masks[old * 32..old * 32 + 32];
    let presence = &legacy.presence[old * 32..old * 32 + 32];
    let scale = if ligand { 0_f32 } else { 1. };
    let ca = [values[3] * scale, values[4] * scale, values[5] * scale];
    let mut row = Row { pos: [ca; SLOTS], ..Row::default() };
    for atom in 0..32 {
        let weight = if defined[atom] { 1_f32 } else { 0. };
        for axis in 0..3 {
            row.pos[index[atom]][axis] += (values[atom * 3 + axis] - ca[axis]) * weight;
        }
        row.mask[index[atom]] |= masks[atom] != 0 && defined[atom];
        row.atom_presence_known[index[atom]] |= presence[atom] != 0 && defined[atom];
    }
    for atom in 0..32 {
        if masks[atom] != 0 && defined[atom] {
            row.pos[index[atom]] = [values[atom * 3], values[atom * 3 + 1], values[atom * 3 + 2]];
        }
    }
    row.aatype = aatype;
    row.res_idx = legacy.sequence[old];
    row.chain_idx = legacy.chains[old];
    row.is_ligand = ligand;
    row.chemical_state = legacy.state[old];
    row.state_known = legacy.state_known[old] != 0;
    row.topology = [legacy.topology[old * 3], legacy.topology[old * 3 + 1], legacy.topology[old * 3 + 2]];
    row.ligand_charge = legacy.charges[old];
    row.ligand_charge_known = legacy.charge_known[old] != 0;
    row.ligand_hydrogen_known = legacy.hydrogen_known[old] != 0;
    row.ligand_hybridization = legacy.hybridization[old];
    row.covalently_bound_to_protein = legacy.covalent[old] != 0;
    if bridges {
        row.crop_positions = [legacy.crop_positions[old * 3], legacy.crop_positions[old * 3 + 1],
            legacy.crop_positions[old * 3 + 2]];
        row.crop_center_eligible = legacy.crop_center_eligible[old] != 0;
    }
    row.token_kind = if ligand { 3 } else { 0 };
    row.token_atom_mask = if ligand { row.mask } else { [true; SLOTS] };
    row
}

/// Canonical sidechain completeness of one legacy token.
fn complete_heavy(legacy: &Tokens, old: usize, molecular: &MolecularSchema) -> bool {
    let aatype = legacy.identity[old];
    let required = &molecular.required_heavy[aatype.clamp(0, 19) as usize];
    (0..20).contains(&aatype) && (0..14).all(|atom| {
        !required[atom] || (legacy.masks[old * 32 + atom] != 0
            && legacy.positions[(old * 32 + atom) * 3..(old * 32 + atom) * 3 + 3].iter().all(|value| value.is_finite()))
    })
}

fn polymer_kind(value: &str) -> u8 {
    let value = value.to_lowercase();
    if value.contains("deoxyribonucleotide") {
        1
    } else if value.contains("ribonucleotide") {
        2
    } else {
        0
    }
}

/// Heavy-atom neighbours from source bonds.
struct Neighbors {
    offsets: Vec<usize>,
    edges: Vec<usize>,
}

impl Neighbors {
    /// Neighbours of the atoms marked in `wanted`; other atoms get none.
    fn read(snapshot: &Snapshot, wanted: &[bool]) -> Self {
        let count = snapshot.row.len();
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        for &[first, second] in &snapshot.bonds {
            if first < 0 || second < 0 {
                continue;
            }
            let (first, second) = (first as usize, second as usize);
            if !snapshot.hydrogen[first] && !snapshot.hydrogen[second] {
                if wanted[first] {
                    pairs.push((first, second));
                }
                if wanted[second] {
                    pairs.push((second, first));
                }
            }
        }
        pairs.sort_unstable();
        pairs.dedup();
        let mut offsets = vec![0; count + 1];
        for &(first, _) in &pairs {
            offsets[first + 1] += 1;
        }
        for index in 0..count {
            offsets[index + 1] += offsets[index];
        }
        Self { offsets, edges: pairs.into_iter().map(|(_, second)| second).collect() }
    }

    /// Sorted neighbours of `atom`.
    fn of(&self, atom: usize) -> &[usize] {
        &self.edges[self.offsets[atom]..self.offsets[atom + 1]]
    }
}

/// Backbone role atoms, indexed like the schema's backbone names.
type Roles = Vec<Option<usize>>;

/// The most name lists each component keeps for reuse.
const LAYOUTS: usize = 16;

fn role(backbone: &[String], name: &str) -> usize {
    backbone.iter().position(|value| value == name).unwrap_or(usize::MAX)
}

/// A unique ribose or deoxyribose backbone, found from covalent identities alone.
fn modified_roles(atoms: &[usize], neighbors: &Neighbors, source: &Source,
    backbone: &[String]) -> Roles {
    let mut members: Vec<usize> = atoms.to_vec();
    members.sort_unstable();
    let inside = |atom: usize| members.binary_search(&atom).is_ok();
    let adjacent = |atom: usize, element: &str, excluded: &[usize]| -> Vec<usize> {
        neighbors.of(atom).iter().copied()
            .filter(|&other| inside(other) && !excluded.contains(&other) && element_is(source, other, element)).collect()
    };
    let index = |name: &str| role(backbone, name);
    let mut candidates: Vec<Roles> = Vec::new();
    for &oxygen in &members {
        if !element_is(source, oxygen, "O") {
            continue;
        }
        for c1 in adjacent(oxygen, "C", &[]) {
            for c2 in adjacent(c1, "C", &[oxygen]) {
                for c3 in adjacent(c2, "C", &[oxygen, c1]) {
                    for c4 in adjacent(c3, "C", &[oxygen, c1, c2]) {
                        let ring = [oxygen, c1, c2, c3, c4];
                        if neighbors.of(c4).binary_search(&oxygen).is_err() {
                            continue;
                        }
                        for c5 in adjacent(c4, "C", &ring) {
                            for o5 in adjacent(c5, "O", &ring) {
                                let phosphorus = adjacent(o5, "P", &ring);
                                let mut roles: Roles = vec![None; backbone.len()];
                                for (name, atom) in [("C4'", c4), ("C5'", c5), ("O5'", o5), ("O4'", oxygen),
                                    ("C1'", c1), ("C2'", c2), ("C3'", c3)] {
                                    if let Some(slot) = roles.get_mut(index(name)) {
                                        *slot = Some(atom);
                                    }
                                }
                                let mut branched = false;
                                for (name, parent) in [("O2'", c2), ("O3'", c3)] {
                                    let attached = adjacent(parent, "O", &ring);
                                    if attached.len() > 1 {
                                        branched = true;
                                        break;
                                    }
                                    if let (Some(&atom), Some(slot)) = (attached.first(), roles.get_mut(index(name))) {
                                        *slot = Some(atom);
                                    }
                                }
                                if branched || phosphorus.len() > 1 {
                                    continue;
                                }
                                if let Some(&atom) = phosphorus.first() {
                                    if let Some(slot) = roles.get_mut(index("P")) {
                                        *slot = Some(atom);
                                    }
                                    let mut oxygens = adjacent(atom, "O", &[o5]);
                                    if oxygens.len() > 3 {
                                        continue;
                                    }
                                    oxygens.sort_by(|&a, &b| source.atoms[a].name.cmp(&source.atoms[b].name));
                                    for (name, atom) in ["OP1", "OP2", "OP3"].into_iter().zip(oxygens) {
                                        if let Some(slot) = roles.get_mut(index(name)) {
                                            *slot = Some(atom);
                                        }
                                    }
                                }
                                if !candidates.contains(&roles) {
                                    candidates.push(roles);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    let o3 = index("O3'");
    let phosphorus = index("P");
    let score = |roles: &Roles| {
        let connected = roles.get(o3).copied().flatten().is_some_and(|atom| {
            neighbors.of(atom).iter().any(|&other| !inside(other) && element_is(source, other, "P"))
        });
        (u8::from(connected), u8::from(roles.get(phosphorus).copied().flatten().is_some()))
    };
    let Some(best) = candidates.iter().map(score).max() else {
        return vec![None; backbone.len()];
    };
    let selected: Vec<&Roles> = candidates.iter().filter(|roles| score(roles) == best).collect();
    if selected.len() == 1 { selected[0].clone() } else { vec![None; backbone.len()] }
}

/// Backbone roles of every nucleotide group.
fn backbone_roles(groups: &Groups, kinds: &[u8], source: &Source, snapshot: &Snapshot,
    backbone: &[String]) -> Vec<Roles> {
    // Only modified nucleotides search the bond graph, and only from their own atoms.
    let modified = |group: usize| kinds[group] != 0
        && !CANONICAL_NUCLEOTIDES.contains(&source.atoms[groups.atoms[group][0]].component.as_str());
    let mut wanted = vec![false; snapshot.row.len()];
    for &group in groups.order.iter().filter(|&&group| modified(group)) {
        for &atom in &groups.atoms[group] {
            wanted[atom] = true;
        }
    }
    let neighbors = Neighbors::read(snapshot, &wanted);
    // Inline text compares without a `memcmp` call.
    let names: Vec<ShortText> = backbone.iter().map(|name| ShortText::new(name)).collect();
    // Residues of one component mostly list the same names in the same order, so they
    // share the atom position of each role. Each component keeps a few name lists.
    let mut layouts: FastMap<&ShortText, Vec<(Vec<&ShortText>, Roles)>> = FastMap::default();
    let mut result = vec![Vec::new(); groups.atoms.len()];
    for &group in &groups.order {
        if kinds[group] == 0 {
            continue;
        }
        let atoms = &groups.atoms[group];
        result[group] = if !modified(group) {
            let known = layouts.entry(&source.atoms[atoms[0]].component).or_default();
            let found = known.iter().find(|(order, _)| order.len() == atoms.len()
                && atoms.iter().zip(order).all(|(&atom, name)| source.atoms[atom].name == **name));
            let positions = match found {
                Some((_, positions)) => positions,
                None => {
                    let mut positions = vec![None; backbone.len()];
                    for (position, &atom) in atoms.iter().enumerate() {
                        if let Some(role) = names.iter().position(|name| *name == source.atoms[atom].name) {
                            positions[role] = Some(position);
                        }
                    }
                    if known.len() >= LAYOUTS {
                        known.pop();
                    }
                    known.push((atoms.iter().map(|&atom| &source.atoms[atom].name).collect(), positions));
                    &known[known.len() - 1].1
                }
            };
            positions.iter().map(|position| position.map(|position| atoms[position])).collect()
        } else {
            modified_roles(atoms, &neighbors, source, backbone)
        };
    }
    result
}

/// Deposited nucleotide numbering with a 64-position offset after broken O3'-P links.
fn nucleotide_indices(groups: &Groups, kinds: &[u8], roles: &[Roles], source: &Source, snapshot: &Snapshot,
    backbone: &[String]) -> Vec<i64> {
    let threshold = f64::from(2.4_f32);
    let (o3, phosphorus) = (role(backbone, "O3'"), role(backbone, "P"));
    let position = |roles: &Roles, index: usize| -> Option<[f64; 3]> {
        let atom = roles.get(index).copied().flatten()?;
        if !snapshot.mask[atom] || snapshot.occupancy[atom] <= 0. {
            return None;
        }
        let value = snapshot.pos[atom];
        value.iter().all(|axis| axis.is_finite()).then(|| value.map(f64::from))
    };
    let mut result = vec![0_i64; groups.atoms.len()];
    // Only lookups by chain read this map, so its order never matters.
    let mut previous: FastMap<&str, (i64, Option<[f64; 3]>, i64)> = FastMap::default();
    for &group in &groups.order {
        if !matches!(kinds[group], 1 | 2) {
            continue;
        }
        let first = groups.atoms[group][0];
        let chain = source.atoms[first].chain.as_str();
        let sequence = snapshot.sequence[first];
        let mut offset = 0;
        if let Some(&(last, end, carried)) = previous.get(chain) {
            offset = carried;
            if sequence == last + 1 {
                let start = position(&roles[group], phosphorus);
                let broken = match (end, start) {
                    (Some(a), Some(b)) => {
                        let delta = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
                        (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]).sqrt() > threshold
                    }
                    _ => true,
                };
                if broken {
                    offset += 64;
                }
            }
        }
        result[group] = sequence + offset;
        previous.insert(chain, (sequence, position(&roles[group], o3), offset));
    }
    result
}

/// The fixed canonical nucleotide schema: backbone, sorted heavy atoms, sorted hydrogens.
fn nucleotide_names<'a>(molecular: &'a MolecularSchema, component: &str, observed: &[&'a str]) -> Vec<&'a str> {
    let mut names: FastSet<&str> = observed.iter().copied().collect();
    let base = component.chars().last().map(String::from).unwrap_or_default();
    if let Some(values) = molecular.bases.get(&base) {
        names.extend(values.iter().map(String::as_str));
    }
    names.extend(molecular.sugar_hydrogens.iter().map(String::as_str));
    names.insert(if component.starts_with('D') { "H2''" } else { "HO2'" });
    for name in &molecular.backbone {
        names.remove(name.as_str());
    }
    let hydrogen = |name: &&str| name.starts_with(['H', 'D', 'T']);
    let mut heavy: Vec<&str> = names.iter().copied().filter(|name| !hydrogen(name)).collect();
    let mut light: Vec<&str> = names.iter().copied().filter(hydrogen).collect();
    heavy.sort_unstable();
    light.sort_unstable();
    molecular.backbone.iter().map(String::as_str).chain(heavy).chain(light).collect()
}

/// The schema of a canonical nucleotide whose atom names all belong to it, with each name's first slot.
struct Canonical<'a> {
    names: Vec<&'a str>,
    slots: FastMap<&'a str, usize>,
    /// Recent atom name lists of this component with their slots, or `None` when some name
    /// is not in the schema.
    layouts: Vec<(Vec<&'a ShortText>, Option<Vec<usize>>)>,
}

struct Context<'a> {
    source: &'a Source,
    snapshot: &'a Snapshot,
    molecular: &'a MolecularSchema,
    table: &'a TokenTable,
    rows: Vec<Row>,
    /// Slot names by token.
    names: Vec<Option<Vec<&'a str>>>,
    ownership: Ownership,
    canonical: FastMap<&'a str, Canonical<'a>>,
    /// Atomic numbers by source element text.
    elements: Vec<(ShortText, i16)>,
}

/// Append the empty names of slots `start..end`, with a comma before each but the first slot.
fn push_empty(text: &mut String, empty: &str, start: usize, end: usize) {
    if start < end {
        text.push_str(&empty[usize::from(start == 0)..3 * (end - start)]);
    }
}

/// Append `value` in decimal.
fn push_decimal(text: &mut String, value: usize) {
    let mut digits = [0_u8; 20];
    let mut first = digits.len();
    let mut rest = value;
    loop {
        first -= 1;
        digits[first] = b'0' + (rest % 10) as u8;
        rest /= 10;
        if rest == 0 {
            break;
        }
    }
    text.extend(digits[first..].iter().map(|&digit| char::from(digit)));
}

/// The slot names of `token`, which may not be set yet.
fn names_of<'n, 'a>(names: &'n mut Vec<Option<Vec<&'a str>>>, token: usize) -> &'n mut Option<Vec<&'a str>> {
    if names.len() <= token {
        names.resize(token + 1, None);
    }
    &mut names[token]
}

impl<'a> Context<'a> {
    fn push(&mut self, row: Row) -> usize {
        self.rows.push(row);
        self.rows.len() - 1
    }

    fn claim(&mut self, atom: usize, token: usize, slot: usize, add: bool) -> Result<(), Failure> {
        self.ownership.claim(&mut self.rows, self.snapshot, atom, token, slot, add)
    }

    /// Return no backbone row when a component requires full atom-level packing.
    #[allow(clippy::too_many_arguments)]
    fn add_nucleotide(&mut self, atoms: &[usize], component: &'a str, kind: u8, chain: i64, residue: i64,
        roles: &Roles) -> Result<Option<usize>, Failure> {
        let molecular = self.molecular;
        let canonical = molecular.nucleotides.iter().position(|name| name == component);
        let selected: Vec<usize>;
        let atoms: &[usize] = if canonical.is_some() {
            atoms
        } else {
            if roles.iter().all(Option::is_none) {
                return Ok(None);
            }
            selected = roles.iter().filter_map(|atom| *atom).collect();
            &selected
        };
        let (schema, slots) = match canonical {
            Some(_) => match self.canonical_slots(atoms, component)? {
                Some(found) => found,
                None => self.observed_slots(atoms, component)?,
            },
            None => {
                self.unique_names(atoms)?;
                let mut slots: FastMap<usize, usize> = FastMap::default();
                for (position, atom) in roles.iter().enumerate() {
                    if let Some(atom) = atom {
                        slots.insert(*atom, position);
                    }
                }
                let slots = atoms.iter().map(|atom| slots.get(atom).copied().unwrap_or(0)).collect();
                (molecular.backbone.iter().map(String::as_str).collect(), slots)
            }
        };
        if schema.len() > SLOTS {
            return Err("A canonical nucleotide exceeds its coordinate schema".into());
        }
        let token = self.push(Row {
            aatype: 20,
            token_kind: kind,
            chain_idx: chain,
            molecular_group: chain,
            res_idx: residue,
            nucleotide_identity: canonical.map_or(8 + kind, |index| index as u8 + 1),
            token_atom_mask: [false; SLOTS],
            ..Row::default()
        });
        for (&atom, &slot) in atoms.iter().zip(&slots) {
            self.claim(atom, token, slot, true)?;
            self.rows[token].nucleotide_atom_elements[slot] = self.atomic_number(atom)?;
        }
        let row = &mut self.rows[token];
        row.nucleotide_atom_mask = row.mask;
        row.token_atom_mask = row.mask;
        let anchor = if row.mask[1] {
            1
        } else {
            row.mask.iter().position(|&value| value)
                .ok_or_else(|| Failure::Index("index 0 is out of bounds for axis 0 with size 0".into()))?
        };
        row.crop_positions = row.pos[anchor];
        row.crop_center_eligible = row.mask[1];
        *names_of(&mut self.names, token) = Some(schema);
        Ok(Some(token))
    }

    /// Reject repeated atom names in one nucleotide.
    fn unique_names(&self, atoms: &[usize]) -> Result<(), Failure> {
        let mut seen = FastSet::with_capacity_and_hasher(atoms.len(), Default::default());
        if !atoms.iter().all(|&atom| seen.insert(self.source.atoms[atom].name.as_str())) {
            return Err("A selected nucleotide has repeated atom names".into());
        }
        Ok(())
    }

    /// The schema and atom slots of a canonical nucleotide whose names all belong to its fixed schema.
    ///
    /// The schema joins the observed names with the fixed names of the component, so when every
    /// observed name is a fixed one, the schema is the same for every residue of that component.
    /// Returns `None` when some name is not a fixed one.
    #[allow(clippy::type_complexity)]
    fn canonical_slots(&mut self, atoms: &[usize], component: &'a str)
        -> Result<Option<(Vec<&'a str>, Vec<usize>)>, Failure> {
        let molecular = self.molecular;
        let fixed = self.canonical.entry(component).or_insert_with(|| {
            let names = nucleotide_names(molecular, component, &[]);
            let mut slots = FastMap::default();
            for (position, &name) in names.iter().enumerate().rev() {
                slots.insert(name, position);
            }
            Canonical { names, slots, layouts: Vec::new() }
        });
        if fixed.names.len() > 128 {
            return Ok(None);
        }
        let source = self.source;
        // Residues of one component usually list the same names in the same order, and
        // comparing the names costs less than hashing them.
        let known = fixed.layouts.iter().find(|(names, _)| names.len() == atoms.len()
            && names.iter().zip(atoms).all(|(name, &atom)| **name == source.atoms[atom].name));
        let slots = match known {
            Some((_, slots)) => slots.clone(),
            None => {
                let slots = atoms.iter().map(|&atom| fixed.slots.get(source.atoms[atom].name.as_str()).copied())
                    .collect::<Option<Vec<usize>>>();
                if fixed.layouts.len() == LAYOUTS {
                    fixed.layouts.pop();
                }
                fixed.layouts.push((atoms.iter().map(|&atom| &source.atoms[atom].name).collect(), slots.clone()));
                slots
            }
        };
        let Some(slots) = slots else {
            return Ok(None);
        };
        // Each name has one slot, so a repeated slot means a repeated name.
        let mut used = 0_u128;
        for &slot in &slots {
            if used >> slot & 1 == 1 {
                return Err("A selected nucleotide has repeated atom names".into());
            }
            used |= 1 << slot;
        }
        // The full output pads every name list to `SLOTS`, so room for that avoids a second allocation.
        let mut names = Vec::with_capacity(SLOTS.max(fixed.names.len()));
        names.extend_from_slice(&fixed.names);
        Ok(Some((names, slots)))
    }

    /// The schema and atom slots of a canonical nucleotide built from its observed names.
    #[allow(clippy::type_complexity)]
    fn observed_slots(&self, atoms: &[usize], component: &str) -> Result<(Vec<&'a str>, Vec<usize>), Failure> {
        let source = self.source;
        let names: Vec<&'a str> = atoms.iter().map(|&atom| source.atoms[atom].name.as_str()).collect();
        let schema = nucleotide_names(self.molecular, component, &names);
        self.unique_names(atoms)?;
        let mut named: FastMap<&str, usize> = FastMap::default();
        for (position, &name) in schema.iter().enumerate().rev() {
            named.insert(name, position);
        }
        let slots = names.iter().map(|name| named.get(name).copied().unwrap_or(0)).collect();
        Ok((schema, slots))
    }

    /// The atomic number of a nucleotide atom, with 1 for every hydrogen isotope.
    fn atomic_number(&mut self, atom: usize) -> Result<i16, Failure> {
        let element = &self.source.atoms[atom].element;
        if let Some(&(_, number)) = self.elements.iter().find(|(known, _)| known == element) {
            return Ok(number);
        }
        let upper = upper(element);
        let number = if matches!(upper.as_ref(), "H" | "D" | "T") {
            1
        } else {
            *self.molecular.atomic_numbers.get(upper.as_ref())
                .ok_or_else(|| Failure::Value(format!("Element '{upper}' not found")))?
        };
        self.elements.push((element.clone(), number));
        Ok(number)
    }

    fn add_ncaa_atoms(&mut self, atoms: &[usize], host: i64, chain: i64, group: usize) -> Result<(), Failure> {
        let mut ordinal = 0;
        for &atom in atoms {
            let element = upper(&self.source.atoms[atom].element);
            if self.ownership.owner[atom] >= 0 || matches!(element.as_ref(), "H" | "D" | "T")
                || !self.ownership.target_known[atom] {
                continue;
            }
            let aatype = self.table.element(&element);
            let token = self.push(Row {
                aatype,
                is_ligand: true,
                token_kind: 3,
                chain_idx: chain,
                res_idx: ordinal,
                molecular_group: group as i64,
                ncaa_parent_token: host,
                token_atom_mask: [false; SLOTS],
                crop_positions: self.snapshot.pos[atom],
                ligand_charge: self.snapshot.charge[atom],
                ligand_charge_known: self.snapshot.charge_known[atom],
                ligand_hybridization: self.snapshot.hybridization[atom],
                ..Row::default()
            });
            self.claim(atom, token, 0, true)?;
            *names_of(&mut self.names, token) = Some(vec![self.source.atoms[atom].name.as_str()]);
            ordinal += 1;
        }
        Ok(())
    }

    fn attach_hydrogens(&mut self, extension: &FastSet<usize>) -> Result<(), Failure> {
        let start = self.molecular.hydrogen_start;
        for atom in 0..self.ownership.owner.len() {
            if self.ownership.owner[atom] >= 0 || !matches!(self.ownership.reason[atom], 4 | 7)
                || !self.snapshot.hydrogen[atom] {
                continue;
            }
            let parent = self.snapshot.hydrogen_parent[atom];
            let host = if parent < 0 { -1 } else { self.ownership.owner[parent as usize] };
            if host < 0 || !extension.contains(&(host as usize)) {
                if self.ownership.reason[atom] == 7 {
                    self.ownership.reason[atom] = 4;
                }
                continue;
            }
            let token = host as usize;
            let row = &self.rows[token];
            let end = if row.is_ligand { start + 3 } else { SLOTS };
            let free: Vec<usize> = (start..end).filter(|&slot| !row.mask[slot]).collect();
            let name = self.source.atoms[atom].name.as_str();
            let preferred = free.iter().copied().find(|&slot| self.molecular.atom_names[slot] == name);
            let Some(slot) = preferred.or(free.first().copied()) else {
                return Err("A ligand atom has more than three attached hydrogen targets".into());
            };
            self.claim(atom, token, slot, true)?;
            let row = &mut self.rows[token];
            if matches!(row.token_kind, 1 | 2) && self.ownership.owner[atom] >= 0 {
                let slot = self.ownership.slot[atom] as usize;
                row.nucleotide_atom_mask[slot] = true;
                row.token_atom_mask[slot] = true;
                row.nucleotide_atom_elements[slot] = 1;
            }
        }
        Ok(())
    }
}

/// The residue part of a source atom's token key.
fn key_text(source: &Source, atom: usize) -> [&str; 4] {
    let atom = &source.atoms[atom];
    [&atom.chain, &atom.auth_sequence, &atom.insertion, &atom.component].map(|text| text.as_str())
}

/// Legacy tokens grouped by their source residue group, in token order.
struct ByGroup {
    offsets: Vec<usize>,
    tokens: Vec<usize>,
}

impl ByGroup {
    fn new(old_group: &[Option<usize>], groups: usize) -> Self {
        let mut offsets = vec![0_usize; groups + 1];
        for &group in old_group.iter().flatten() {
            offsets[group + 1] += 1;
        }
        for group in 0..groups {
            offsets[group + 1] += offsets[group];
        }
        let mut cursor = offsets.clone();
        let mut tokens = vec![0_usize; offsets[groups]];
        for (old, group) in old_group.iter().enumerate() {
            if let Some(group) = *group {
                tokens[cursor[group]] = old;
                cursor[group] += 1;
            }
        }
        Self { offsets, tokens }
    }

    fn of(&self, group: usize) -> &[usize] {
        &self.tokens[self.offsets[group]..self.offsets[group + 1]]
    }
}

/// Pack legacy tokens and source atoms into typed rows with one owner per slot.
///
/// With `index`, skip the slot names, SASA rows and ligand connectivity, which only
/// the full output reads; every check that can fail still runs.
pub(super) fn pack<'a>(legacy: &'a Tokens, source: &'a Source, snapshot: &'a Snapshot,
    molecular: &'a MolecularSchema, table: &'a TokenTable, index: bool) -> Result<Packed<'a>, Failure> {
    let groups = Groups::read(source, snapshot)?;
    let count = legacy.identity.len();
    let bridges = !legacy.bridge_edges.is_empty();
    let mut text = String::new();
    // Tokens of one residue are adjacent, so reuse the previous lookup while the residue key repeats.
    let mut previous: Option<(&[ShortText; 5], Option<usize>)> = None;
    let old_group: Vec<Option<usize>> = legacy.keys.iter().map(|key| {
        if let Some((last, group)) = previous
            && same_residue(last, key) {
            return group;
        }
        group_key(&mut text, &key[0], &key[1], &key[2], &key[3]);
        let group = groups.lookup.get(text.as_str()).copied();
        previous = Some((key, group));
        group
    }).collect();
    let old_by_group = ByGroup::new(&old_group, groups.atoms.len());
    let mut kinds = vec![0_u8; groups.atoms.len()];
    // Groups of one entity are mostly adjacent, so reuse the kind while the entity repeats.
    let mut entity: Option<(&ShortText, u8)> = None;
    for &group in &groups.order {
        let atom = &source.atoms[groups.atoms[group][0]];
        kinds[group] = match entity {
            Some((last, kind)) if *last == atom.entity => kind,
            _ => {
                let kind = polymer_kind(source.entity_types.get(atom.entity.as_str()).map_or("", String::as_str));
                entity = Some((&atom.entity, kind));
                kind
            }
        };
    }
    let roles = backbone_roles(&groups, &kinds, source, snapshot, &molecular.backbone);
    let residues = nucleotide_indices(&groups, &kinds, &roles, source, snapshot, &molecular.backbone);
    let mut keep: Vec<bool> = (0..count).map(|old| {
        let unknown = legacy.identity[old] == 20 && (0..4).all(|atom| legacy.masks[old * 32 + atom] != 0);
        legacy.ligand[old] != 0 || complete_heavy(legacy, old, molecular) || unknown
    }).collect();
    let mut vacant: Vec<Vec<usize>> = vec![Vec::new(); groups.atoms.len()];
    for &group in &groups.order {
        vacant[group] = groups.atoms[group].iter().copied().filter(|&atom| {
            let occupancy = snapshot.occupancy[atom];
            occupancy.is_nan() || occupancy <= 0.
        }).collect();
    }
    for (old, group) in old_group.iter().enumerate() {
        let Some(group) = group else { continue };
        for &atom in &vacant[*group] {
            let name = source.atoms[atom].name.as_str();
            if (legacy.ligand[old] != 0 && name == legacy.keys[old][4])
                || (legacy.ligand[old] == 0 && !snapshot.hydrogen[atom]
                    && (legacy.identity[old] < 20 || PROTEIN_BACKBONE.contains(&name))) {
                keep[old] = false;
            }
        }
    }
    let mut context = Context {
        source, snapshot, molecular, table, rows: Vec::with_capacity(count + groups.order.len()),
        names: Vec::new(), ownership: Ownership::new(snapshot), canonical: FastMap::default(), elements: Vec::new(),
    };
    let mut keys: Vec<[&'a str; 5]> = Vec::with_capacity(count);
    let mut old_to_new = vec![-1_i64; count];
    let mut hosts: Vec<(usize, usize)> = Vec::new();
    let mut host_index: FastMap<usize, usize> = FastMap::default();
    let mut set_host = |hosts: &mut Vec<(usize, usize)>, group: usize, token: usize| {
        match host_index.get(&group) {
            Some(&index) => hosts[index].1 = token,
            None => {
                host_index.insert(group, hosts.len());
                hosts.push((group, token));
            }
        }
    };
    for (old, key) in legacy.keys.iter().enumerate() {
        let group = old_group[old];
        if group.is_some_and(|group| kinds[group] != 0) {
            continue;
        }
        let ligand = legacy.ligand[old] != 0;
        if !keep[old] {
            if let (Some(group), false) = (group, ligand) {
                let replacement = if legacy.identity[old] == 20 { 3 } else { 2 };
                for &atom in &groups.atoms[group] {
                    if context.ownership.reason[atom] != 8 {
                        context.ownership.reason[atom] = replacement;
                    }
                }
            }
            continue;
        }
        let mut row = legacy_row(legacy, old, molecular, bridges);
        row.sequence_target_known = ligand || !molecular.unknown_components.contains(key[3].as_str());
        if let Some(group) = group {
            row.molecular_group = group as i64;
        }
        let token = context.push(row);
        old_to_new[old] = token as i64;
        keys.push(key.each_ref().map(ShortText::as_str));
        if let Some(group) = group {
            if !ligand {
                set_host(&mut hosts, group, token);
            }
            context.ownership.match_legacy(&mut context.rows, snapshot, source, &molecular.atom_names,
                &groups.atoms[group], token)?;
        }
    }
    let mut chains: FastMap<&str, i64> = FastMap::default();
    for (index, key) in keys.iter().enumerate() {
        chains.insert(key[0], context.rows[index].chain_idx);
    }
    let mut atomized: Vec<usize> = Vec::new();
    for &group in &groups.order {
        if kinds[group] == 0 {
            continue;
        }
        let atoms = &groups.atoms[group];
        if !atoms.iter().any(|&atom| context.ownership.target_known[atom]) {
            continue;
        }
        let key = key_text(source, atoms[0]);
        let chain = match chains.get(key[0]) {
            Some(&chain) => chain,
            None => {
                let chain = chains.values().copied().max().unwrap_or(-1) + 1;
                chains.insert(key[0], chain);
                chain
            }
        };
        let canonical = molecular.nucleotides.iter().any(|name| name == key[3]);
        match context.add_nucleotide(atoms, key[3], kinds[group], chain, residues[group], &roles[group])? {
            Some(token) => {
                keys.push([key[0], key[1], key[2], key[3], "nucleotide"]);
                for &old in old_by_group.of(group) {
                    old_to_new[old] = token as i64;
                }
                if !canonical {
                    set_host(&mut hosts, group, token);
                }
            }
            None if !canonical => atomized.push(group),
            None => {}
        }
    }
    atomized.sort_unstable();
    atomized.dedup();
    let mut cohort = context.rows.iter().map(|row| row.chain_idx).max().unwrap_or(-1) + 1;
    let mut extension: FastSet<usize> = FastSet::default();
    let entries: Vec<(usize, i64)> = hosts.iter().map(|&(group, token)| (group, token as i64))
        .chain(atomized.iter().map(|&group| (group, -1))).collect();
    for (group, host) in entries {
        if host >= 0 && context.rows[host as usize].aatype != 20 {
            continue;
        }
        let before = context.rows.len();
        let first = groups.atoms[group][0];
        let key = if host >= 0 {
            let key = keys[host as usize];
            [key[0], key[1], key[2], key[3]]
        } else {
            key_text(source, first)
        };
        context.add_ncaa_atoms(&groups.atoms[group], host, cohort, group)?;
        cohort += 1;
        for token in before..context.rows.len() {
            let name = context.names.get(token).and_then(Option::as_ref).and_then(|names| names.first()).copied()
                .unwrap_or_default();
            keys.push([key[0], key[1], key[2], key[3], name]);
        }
        extension.extend(before..context.rows.len());
        if host >= 0 {
            extension.insert(host as usize);
        }
    }
    context.attach_hydrogens(&extension)?;
    let mut by_name: FastMap<&str, i64> = FastMap::default();
    for &group in &groups.order {
        if kinds[group] == 0 {
            continue;
        }
        let atoms = &groups.atoms[group];
        let owner = &context.ownership.owner;
        let first = atoms.first().map_or(-1, |&atom| owner[atom]);
        if atoms.iter().all(|&atom| owner[atom] == first) {
            // With one owner for every atom, a token maps to that owner when some atom has its
            // name. Tokens follow atom order, so the search starts after the last match.
            let mut cursor = 0;
            for &old in old_by_group.of(group) {
                let name = &legacy.keys[old][4];
                let found = (cursor..atoms.len()).chain(0..cursor)
                    .find(|&index| source.atoms[atoms[index]].name == *name);
                old_to_new[old] = match found {
                    Some(index) => {
                        cursor = index + 1;
                        first
                    }
                    None => -1,
                };
            }
            continue;
        }
        by_name.clear();
        for &atom in &groups.atoms[group] {
            by_name.insert(source.atoms[atom].name.as_str(), context.ownership.owner[atom]);
        }
        for &old in old_by_group.of(group) {
            old_to_new[old] = by_name.get(legacy.keys[old][4].as_str()).copied().unwrap_or(-1);
        }
    }
    if !index {
        for atom in 0..context.ownership.owner.len() {
            let token = context.ownership.owner[atom];
            if token < 0 {
                continue;
            }
            let token = token as usize;
            let row = &context.rows[token];
            if !(row.is_ligand || row.aatype == 20) {
                continue;
            }
            let kind = row.token_kind;
            let names = names_of(&mut context.names, token).get_or_insert_with(|| {
                if kind == 0 { molecular.atom_names.iter().map(String::as_str).collect() } else { vec![""; SLOTS] }
            });
            if names.len() < SLOTS {
                names.resize(SLOTS, "");
            }
            // Keep an equal name, so that equal lists keep the same `&str` values for `metadata`.
            let name = &source.atoms[atom].name;
            let slot = &mut names[context.ownership.slot[atom] as usize];
            if *slot != *name {
                *slot = name.as_str();
            }
        }
    }
    for &group in &groups.order {
        let unassigned: Vec<usize> = groups.atoms[group].iter().copied()
            .filter(|&atom| context.ownership.reason[atom] == 7).collect();
        if unassigned.is_empty() {
            continue;
        }
        let heavy: Vec<usize> = unassigned.iter().copied().filter(|&atom| !snapshot.hydrogen[atom]).collect();
        for &atom in &unassigned {
            if snapshot.hydrogen[atom] {
                context.ownership.reason[atom] = 4;
            }
        }
        if let (Some(&index), false) = (host_index.get(&group), heavy.is_empty()) {
            let host = hosts[index].1;
            let names: Vec<&str> = heavy.iter().map(|&atom| source.atoms[atom].name.as_str()).collect();
            return Err(Failure::Value(format!("A protein source heavy atom has no molecular target slot: {}; {}",
                json::repr_sequence(keys[host], false), json::repr_sequence(names, false))));
        }
        for atom in heavy {
            context.ownership.reason[atom] = 1;
        }
    }
    let Context { mut rows, names, ownership, .. } = context;
    let mut unusual = 0;
    for row in &mut rows {
        if row.ligand_charge_known && row.ligand_charge.abs() > 6. {
            row.ligand_charge_known = false;
            unusual += 1;
        }
        if row.is_ligand {
            row.token_atom_mask = row.mask;
        }
    }
    let Ownership { owner, slot, reason, .. } = ownership;
    let (atom_sasa, atom_sasa_mask) = ownership_payload(snapshot, &owner, &slot, &reason, &rows, index)?;
    let kind: Vec<u8> = rows.iter().map(|row| row.token_kind).collect();
    let (attachment_pairs, attachment_order_known) = attachment_targets(snapshot, &owner, &slot, &kind)?;
    let (ligand_bonds, ligand_bond_order, ligand_bond_known, attachment_known) = if index {
        Default::default()
    } else {
        connectivity_targets(snapshot, &owner, &slot, &kind)
    };
    for atom in 0..owner.len() {
        if owner[atom] < 0 || slot[atom] != 0 {
            continue;
        }
        let token = owner[atom] as usize;
        if rows[token].ncaa_parent_token >= 0 {
            let observed = rows[token].mask[37..40].iter().filter(|&&value| value).count() as i64;
            let hydrogens = snapshot.hydrogen_count[atom];
            rows[token].ligand_hydrogen_known = snapshot.hydrogen_known[atom] && (0..=3).contains(&hydrogens)
                && observed == hydrogens;
        }
    }
    let mut protein: Vec<[i64; 2]> = attachment_pairs.iter().filter(|pair| kind[pair[1] as usize] == 0)
        .map(|pair| [pair[0], pair[1]]).collect();
    for &(first, second) in &legacy.attachments {
        let mapped = [old_to_new[first], old_to_new[second]];
        if mapped.iter().all(|&token| token >= 0) && rows[mapped[0] as usize].is_ligand
            && kind[mapped[1] as usize] == 0 {
            protein.push(mapped);
        }
    }
    protein.sort_unstable();
    protein.dedup();
    for row in &mut rows {
        row.covalently_bound_to_protein = false;
    }
    for pair in &protein {
        rows[pair[0] as usize].covalently_bound_to_protein = true;
    }
    let mut bridge_index = Vec::with_capacity(legacy.bridge_edges.len());
    for &(first, second, _) in &legacy.bridge_edges {
        let mapped = [old_to_new[first], old_to_new[second]];
        if mapped.iter().any(|&token| token < 0) {
            return Err("Molecular packing removed a metal-hydride bridge endpoint".into());
        }
        bridge_index.push(mapped);
    }
    let mut counts = [0_usize; 9];
    for &value in &reason {
        if let Some(count) = counts.get_mut(value as usize) {
            *count += 1;
        }
    }
    let nucleotide_groups: HashSet<i64> = groups.order.iter().filter(|&&group| kinds[group] != 0)
        .map(|&group| group as i64).collect();
    let parent_of = |row: &Row| (row.ncaa_parent_token >= 0).then(|| &rows[row.ncaa_parent_token as usize]);
    let stats = Stats {
        tokens: keys.len(),
        source_atoms: owner.len(),
        unsupported_known_charge_targets: unusual,
        unknown_sequence_tokens: rows.iter().filter(|row| !row.sequence_target_known).count(),
        represented_atoms: owner.iter().filter(|&&token| token >= 0).count(),
        exclusions: counts.iter().enumerate().skip(1).filter(|&(_, &count)| count > 0)
            .map(|(index, &count)| (molecular.exclusions[index].clone(), count)).collect(),
        nucleotide_residues: rows.iter().filter(|row| matches!(row.token_kind, 1 | 2)).count(),
        unknown_source_atom_tokens: rows.iter().filter(|row| parent_of(row)
            .is_some_and(|parent| !parent.sequence_target_known)).count(),
        ncaa_sidechain_tokens: rows.iter().filter(|row| parent_of(row)
            .is_some_and(|parent| parent.sequence_target_known && parent.token_kind == 0)).count(),
        modified_nucleotide_atom_tokens: rows.iter()
            .filter(|row| row.is_ligand && nucleotide_groups.contains(&row.molecular_group)).count(),
        fully_atomized_modified_nucleotides: atomized.len(),
    };
    Ok(Packed {
        rows, keys, owner, slot, reason, atom_sasa, atom_sasa_mask, attachment_pairs, attachment_order_known,
        ligand_bonds, ligand_bond_order, ligand_bond_known, attachment_known,
        protein_attachment_index: protein, legacy_token_index: old_to_new,
        names: names.into_iter().enumerate().filter_map(|(token, names)| names.map(|names| (token, names))).collect(),
        bridge_index,
        bridge_type: legacy.bridge_edges.iter().map(|&(_, _, kind)| kind).collect(), stats,
    })
}

type Areas = (Vec<[f32; SLOTS]>, Vec<[bool; SLOTS]>);

/// Require exact observed geometry for each represented source atom; without `index`, place source SASA.
fn ownership_payload(snapshot: &Snapshot, owner: &[i64], slot: &[i64], reason: &[u8],
    rows: &[Row], index: bool) -> Result<Areas, Failure> {
    let mut seen = vec![0_u128; rows.len()];
    let mut duplicate = false;
    for atom in 0..owner.len() {
        let represented = owner[atom] >= 0;
        let valid = if represented {
            let (token, index) = (owner[atom] as usize, slot[atom]);
            (0..SLOTS as i64).contains(&index) && token < rows.len() && reason[atom] == 0
                && rows[token].mask[index as usize] && rows[token].pos[index as usize] == snapshot.pos[atom]
        } else {
            slot[atom] == -1 && reason[atom] != 0
        };
        if !valid {
            return Err("Molecular targets changed a represented source atom or lost its exclusion reason".into());
        }
        if represented {
            let bit = 1_u128 << slot[atom];
            duplicate |= seen[owner[atom] as usize] & bit != 0;
            seen[owner[atom] as usize] |= bit;
        }
    }
    if duplicate {
        return Err("Two source atoms cannot own the same modeled atom slot".into());
    }
    if snapshot.sasa.iter().zip(&snapshot.sasa_known).any(|(&value, &known)| known && !(value.is_finite() && value >= 0.)) {
        return Err("Known source SASA values must be finite, nonnegative and explicitly masked".into());
    }
    if index {
        return Ok((Vec::new(), Vec::new()));
    }
    let mut area = vec![[0_f32; SLOTS]; rows.len()];
    let mut mask = vec![[false; SLOTS]; rows.len()];
    for atom in 0..owner.len() {
        if owner[atom] >= 0 {
            area[owner[atom] as usize][slot[atom] as usize] = snapshot.sasa[atom];
            mask[owner[atom] as usize][slot[atom] as usize] = snapshot.sasa_known[atom];
        }
    }
    Ok((area, mask))
}

/// Ligand-token to polymer-token attachment targets with their slot and order.
fn attachment_targets(snapshot: &Snapshot, owner: &[i64], slot: &[i64],
    kind: &[u8]) -> Result<(Vec<[i64; 4]>, Vec<bool>), Failure> {
    let mut rows: BTreeMap<(i64, i64, i64), (i64, bool)> = BTreeMap::new();
    for (edge, &[first, second]) in snapshot.bonds.iter().enumerate() {
        if first < 0 || second < 0 || owner[first as usize] < 0 || owner[second as usize] < 0 {
            continue;
        }
        let (mut first, mut second) = (first as usize, second as usize);
        let (mut a, mut b) = (owner[first], owner[second]);
        if kind[a as usize] <= 2 && kind[b as usize] == 3 {
            std::mem::swap(&mut first, &mut second);
            std::mem::swap(&mut a, &mut b);
        }
        let _ = first;
        if kind[a as usize] != 3 || kind[b as usize] > 2 {
            continue;
        }
        let key = (a, b, slot[second]);
        let order = i64::from(snapshot.bond_order[edge]);
        let known = snapshot.bond_order_known[edge];
        let previous = rows.get(&key).copied();
        if let Some((value, true)) = previous
            && known && value != order {
            return Err("Source bonds disagree on a protein attachment order".into());
        }
        if previous.is_none() || known {
            rows.insert(key, (order, known));
        }
    }
    Ok(rows.into_iter().map(|((a, b, slot), (order, known))| ([a, b, slot, order], known)).unzip())
}

type Connectivity = (Vec<[i64; 2]>, Vec<i8>, Vec<bool>, Vec<bool>);

/// Bonds inside one ligand component and whether each atom token's attachments are known.
fn connectivity_targets(snapshot: &Snapshot, owner: &[i64], slot: &[i64], kind: &[u8]) -> Connectivity {
    let mut pairs = Vec::new();
    let mut orders = Vec::new();
    let mut known = Vec::new();
    for (edge, &[first, second]) in snapshot.bonds.iter().enumerate() {
        if first < 0 || second < 0 {
            continue;
        }
        let (first, second) = (first as usize, second as usize);
        let (a, b) = (owner[first], owner[second]);
        if a < 0 || b < 0 {
            continue;
        }
        if a != b && kind[a as usize] == 3 && kind[b as usize] == 3 && snapshot.group[first] == snapshot.group[second]
            && slot[first] == 0 && slot[second] == 0 {
            pairs.push([a, b]);
            orders.push(snapshot.bond_order[edge]);
            known.push(snapshot.bond_order_known[edge]);
        }
    }
    let mut attachment = vec![false; kind.len()];
    for atom in 0..owner.len() {
        if owner[atom] >= 0 && slot[atom] == 0 {
            let token = owner[atom] as usize;
            attachment[token] = snapshot.connectivity_known[atom] && kind[token] == 3;
        }
    }
    for &[first, second] in &snapshot.bonds {
        if first < 0 || second < 0 {
            continue;
        }
        let (a, b) = (owner[first as usize], owner[second as usize]);
        if a >= 0 && b < 0 {
            attachment[a as usize] = false;
        }
        if b >= 0 && a < 0 {
            attachment[b as usize] = false;
        }
    }
    if snapshot.connections.iter().flatten().any(|&atom| atom < 0) {
        attachment.fill(false);
    }
    (pairs, orders, known, attachment)
}

/// Row fields in the order of the Python row defaults.
const ROW_FIELDS: [&str; 25] = ["pos", "mask", "aatype", "res_idx", "chain_idx", "is_ligand", "chemical_state",
    "state_known", "topology", "sequence_target_known", "ligand_charge", "ligand_charge_known",
    "ligand_hydrogen_known", "ligand_hybridization", "covalently_bound_to_protein", "atom_presence_known",
    "crop_positions", "crop_center_eligible", "token_kind", "nucleotide_identity", "token_atom_mask",
    "nucleotide_atom_mask", "nucleotide_atom_elements", "ncaa_parent_token", "molecular_group"];

/// Row fields that the heavy index stores or hashes into its crop layout.
const INDEX_FIELDS: [&str; 8] = ["mask", "aatype", "res_idx", "chain_idx", "is_ligand", "crop_positions",
    "crop_center_eligible", "token_kind"];

impl Packed<'_> {
    /// One row field as a NumPy buffer.
    fn row_array(&self, name: &'static str) -> Array {
        let rows = &self.rows;
        let count = rows.len();
        let int = |f: fn(&Row) -> i64| arrays::int64(&rows.iter().map(f).collect::<Vec<_>>(), vec![count]);
        let flag = |f: fn(&Row) -> bool| arrays::boolean(&rows.iter().map(f).collect::<Vec<_>>(), vec![count]);
        // A flattening iterator gives no exact length, so reserve the lists first.
        let slots = |f: fn(&Row) -> &[bool; SLOTS]| {
            let mut values = Vec::with_capacity(count * SLOTS);
            for row in rows {
                values.extend_from_slice(f(row));
            }
            arrays::boolean(&values, vec![count, SLOTS])
        };
        match name {
            "pos" => {
                let mut values = Vec::with_capacity(count * SLOTS);
                for row in rows {
                    values.extend_from_slice(&row.pos);
                }
                arrays::points(&values, vec![count, SLOTS, 3])
            }
            "mask" => slots(|row| &row.mask),
            "aatype" => int(|row| row.aatype),
            "res_idx" => int(|row| row.res_idx),
            "chain_idx" => int(|row| row.chain_idx),
            "is_ligand" => flag(|row| row.is_ligand),
            "chemical_state" => int(|row| row.chemical_state),
            "state_known" => flag(|row| row.state_known),
            "topology" => arrays::int64(&rows.iter().flat_map(|row| row.topology).collect::<Vec<_>>(), vec![count, 3]),
            "sequence_target_known" => flag(|row| row.sequence_target_known),
            "ligand_charge" => arrays::float32(&rows.iter().map(|row| row.ligand_charge).collect::<Vec<_>>(),
                vec![count]),
            "ligand_charge_known" => flag(|row| row.ligand_charge_known),
            "ligand_hydrogen_known" => flag(|row| row.ligand_hydrogen_known),
            "ligand_hybridization" => int(|row| row.ligand_hybridization),
            "covalently_bound_to_protein" => flag(|row| row.covalently_bound_to_protein),
            "atom_presence_known" => slots(|row| &row.atom_presence_known),
            "crop_positions" => arrays::points(&rows.iter().map(|row| row.crop_positions).collect::<Vec<_>>(),
                vec![count, 3]),
            "crop_center_eligible" => flag(|row| row.crop_center_eligible),
            "token_kind" => arrays::uint8(rows.iter().map(|row| row.token_kind).collect(), vec![count]),
            "nucleotide_identity" => arrays::uint8(rows.iter().map(|row| row.nucleotide_identity).collect(),
                vec![count]),
            "token_atom_mask" => slots(|row| &row.token_atom_mask),
            "nucleotide_atom_mask" => slots(|row| &row.nucleotide_atom_mask),
            "nucleotide_atom_elements" => {
                let mut values = Vec::with_capacity(count * SLOTS);
                for row in rows {
                    values.extend_from_slice(&row.nucleotide_atom_elements);
                }
                arrays::int16(&values, vec![count, SLOTS])
            }
            "ncaa_parent_token" => int(|row| row.ncaa_parent_token),
            _ => int(|row| row.molecular_group),
        }
    }

    /// Row fields in the order of the Python row defaults.
    pub(super) fn row_arrays(&self) -> Vec<(&'static str, Array)> {
        ROW_FIELDS.iter().map(|&name| (name, self.row_array(name))).collect()
    }

    /// The row fields of the heavy index.
    pub(super) fn index_arrays(&self) -> Vec<(&'static str, Array)> {
        INDEX_FIELDS.iter().map(|&name| (name, self.row_array(name))).collect()
    }

    fn metadata(&self, molecular: &MolecularSchema) -> String {
        let mut text = String::with_capacity(molecular.schema_json.len() + self.names.len() * 512);
        text.push_str("{\"schema\":");
        text.push_str(&molecular.schema_json);
        text.push_str(",\"token_slot_names\":{");
        // Most slots have no name, so write each run of empty names as one copy.
        let longest = self.names.iter().map(|(_, names)| names.len()).max().unwrap_or(0);
        let empty = ",\"\"".repeat(longest);
        // Tokens of one nucleotide type share a list of the same `&str` values, so write each
        // such list once and copy its text.
        let mut written: Vec<(&[&str], String)> = Vec::new();
        for (index, (token, names)) in self.names.iter().enumerate() {
            if index > 0 {
                text.push(',');
            }
            text.push('"');
            push_decimal(&mut text, *token);
            text.push_str("\":");
            let same = |list: &[&str]| list.len() == names.len()
                && list.iter().zip(names).all(|(left, right)| std::ptr::eq(*left, *right));
            if let Some((_, list)) = written.iter().find(|(list, _)| same(list)) {
                text.push_str(list);
                continue;
            }
            let first = text.len();
            text.push('[');
            let mut start = 0;
            for (slot, name) in names.iter().enumerate() {
                if name.is_empty() {
                    continue;
                }
                push_empty(&mut text, &empty, start, slot);
                if slot > 0 {
                    text.push(',');
                }
                json::string(&mut text, name);
                start = slot + 1;
            }
            push_empty(&mut text, &empty, start, names.len());
            text.push(']');
            // A list with one name almost always holds a ligand atom's own name, which no other
            // token shares.
            if written.len() < 256 && names.iter().filter(|name| !name.is_empty()).nth(1).is_some() {
                written.push((names, text[first..].to_owned()));
            }
        }
        text.push_str("}}");
        text
    }

    /// Ownership, attachment and bridge targets after the source snapshot.
    pub(super) fn target_arrays(&self, molecular: &MolecularSchema) -> Vec<(&'static str, Array)> {
        let atoms = self.owner.len();
        let tokens = self.rows.len();
        let pairs: Vec<i64> = self.attachment_pairs.iter().flatten().copied().collect();
        let areas: Vec<f32> = self.atom_sasa.iter().flatten().copied().collect();
        let area_mask: Vec<bool> = self.atom_sasa_mask.iter().flatten().copied().collect();
        vec![
            ("source_atom_token", arrays::int64(&self.owner, vec![atoms])),
            ("source_atom_slot", arrays::int64(&self.slot, vec![atoms])),
            ("source_atom_in_scope", arrays::boolean(&self.owner.iter().map(|&token| token >= 0)
                .collect::<Vec<_>>(), vec![atoms])),
            ("source_exclusion_reason", arrays::uint8(self.reason.clone(), vec![atoms])),
            ("atom_sasa", arrays::float32(&areas, vec![tokens, SLOTS])),
            ("atom_sasa_mask", arrays::boolean(&area_mask, vec![tokens, SLOTS])),
            ("attachment_pairs", arrays::int64(&pairs, vec![self.attachment_pairs.len(), 4])),
            ("attachment_order_known", arrays::boolean(&self.attachment_order_known,
                vec![self.attachment_order_known.len()])),
            ("ligand_bonds", arrays::pairs(&self.ligand_bonds)),
            ("ligand_bond_order", arrays::int8(&self.ligand_bond_order, vec![self.ligand_bond_order.len()])),
            ("ligand_bond_known", arrays::boolean(&self.ligand_bond_known, vec![self.ligand_bond_known.len()])),
            ("attachment_known", arrays::boolean(&self.attachment_known, vec![tokens])),
            ("protein_attachment_index", arrays::columns(&self.protein_attachment_index)),
            ("protein_attachment_offsets", arrays::int64(&[0, self.protein_attachment_index.len() as i64], vec![2])),
            ("legacy_token_index", arrays::int64(&self.legacy_token_index, vec![self.legacy_token_index.len()])),
            ("sample_offsets", arrays::int64(&[0, tokens as i64], vec![2])),
            ("molecular_schema", arrays::text(molecular.schema_sha256.clone())),
            ("molecular_metadata", arrays::text(self.metadata(molecular))),
            ("ligand_bridge_index", arrays::columns(&self.bridge_index)),
            ("ligand_bridge_type", arrays::uint8(self.bridge_type.clone(), vec![self.bridge_type.len()])),
            ("ligand_bridge_offsets", arrays::int64(&[0, self.bridge_index.len() as i64], vec![2])),
        ]
    }
}

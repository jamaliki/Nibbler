//! Protein residue tokens and heavy ligand atom tokens.

use super::record::{Tokens, renumber_ligands, same_residue};
use super::source::{Atom, Source};
use super::{
    HydrogenParents, LigandFeatures, ResidueParents, Schema, TokenTable, residue, topology,
};
use crate::hash::{FastMap, FastSet};
use crate::pdbx::ShortText;
use std::collections::{BTreeSet, HashMap};

fn hydrogen_overflow(
    source: &Source,
    features: &LigandFeatures,
    parent: usize,
    assigned: &[usize],
) -> String {
    let atom = &source.atoms[parent];
    let children: Vec<_> = assigned
        .iter()
        .map(|&row| {
            let child = &source.atoms[row];
            serde_json::json!({"name":child.name,"element":child.element,
            "position":child.position,"distance":residue::distance(atom.position,child.position)})
        })
        .collect();
    let ccd = features
        .get(atom.component.as_str())
        .and_then(|component| component.get(atom.name.as_str()));
    let partners: Vec<_> = source
        .bonds
        .iter()
        .filter_map(|&(a, b)| {
            let other = if a == parent {
                b
            } else if b == parent {
                a
            } else {
                return None;
            };
            let other = &source.atoms[other];
            Some(
                serde_json::json!({"chain":other.chain,"component":other.component,
            "atom_name":other.name,"auth_sequence":other.auth_sequence,"position":other.position}),
            )
        })
        .collect();
    let detail = serde_json::json!({"parent_key":[atom.chain.clone(),atom.auth_sequence.clone(),
        atom.insertion.clone(),atom.component.clone(),atom.name.clone()],"parent_element":atom.element,
        "parent_position":atom.position,"capacity":3,"hydrogens":children,"covalent_partners":partners,
        "ccd_features":ccd.map(|values| serde_json::json!({"formal_charge":values[0],
            "hybridization":values[1],"total_hydrogen":values[2],"formal_charge_known":values[3] != 0}))});
    format!("Ligand parent has more than three hydrogens; unsupported Atom32 ligand: {detail}")
}

fn parent_error(atom: &Atom, parent: &str, reason: &str) -> String {
    let detail = serde_json::json!({"hydrogen_key":[atom.chain.clone(),atom.auth_sequence.clone(),
        atom.insertion.clone(),atom.component.clone(),atom.name.clone()],
        "position":atom.position,"parent":parent,"reason":reason});
    format!("Invalid ligand hydrogen parent: {detail}")
}

fn declared_parents(
    source: &Source,
    bridges: &BTreeSet<usize>,
) -> Result<FastMap<usize, Option<usize>>, String> {
    let mut result = FastMap::default();
    for &(child, parent) in &source.hydrogen_bonds {
        if bridges.contains(&child) {
            continue;
        }
        if let Some(previous) = result.insert(child, parent)
            && previous != parent
        {
            return Err(parent_error(
                &source.atoms[child],
                "struct_conn",
                "conflicting covalent parents",
            ));
        }
    }
    Ok(result)
}

#[cfg(feature = "reduce3")]
fn add_bridges(
    source: &Source,
    sites: &[super::bridges::Site],
    excluded: &[String],
    hydrogen: usize,
    tokens: &mut Tokens,
    mapping: &mut [usize],
) -> Result<HashMap<usize, usize>, String> {
    let mut counts = HashMap::new();
    for site in sites {
        let first = &source.atoms[site.metal_rows[0]];
        if excluded.iter().any(|name| *name == first.component) {
            continue;
        }
        let token = tokens.identity.len();
        let mut packed = residue::Packed::empty([-1; 3]);
        if let Some(row) = site.hydrogen_row {
            packed.position[0] = source.atoms[row].position;
            packed.mask[0] = true;
            mapping[row] = token;
        }
        tokens.push(
            first,
            hydrogen,
            &packed,
            true,
            [site.formal_charge, 0, 0, 1],
        );
        // A hydride token is never a crop centre.
        tokens.crop_center_eligible[token] = 0;
        tokens.keys[token][4] = ShortText::new(&site.hydrogen_name);
        tokens.charge_known[token] = 1;
        tokens.hydrogen_known[token] = 1;
        if site.hydrogen_row.is_none() {
            let [first, second] = site.metal_rows.map(|row| source.atoms[row].position);
            for axis in 0..3 {
                tokens.crop_positions[token * 3 + axis] = first[axis] * 0.5 + second[axis] * 0.5;
            }
        }
        for edge in site.edges(token, mapping)? {
            tokens
                .bridge_edges
                .push((edge.hydrogen, edge.metal, edge.kind as u8));
        }
        for row in site.metal_rows {
            *counts.entry(row).or_insert(0) += 1;
        }
    }
    Ok(counts)
}

fn nucleotide_reference(component: &str, name: &str) -> Option<String> {
    #[cfg(feature = "reduce3")]
    {
        let class = reduce3::resclass::get_class_ext(component, true);
        (class.is_rna_dna() || class == reduce3::resclass::ResClass::Ccp4MonLibRnaDna)
            .then(|| reduce3::names::na_reference_name(name))
    }
    #[cfg(not(feature = "reduce3"))]
    {
        let _ = (component, name);
        None
    }
}

fn ccd_parent<'a>(
    parents: &'a HashMap<String, String>,
    atom: &Atom,
) -> Result<Option<&'a str>, String> {
    if let Some(parent) = parents.get(atom.name.as_str()) {
        return Ok(Some(parent));
    }
    let name = if matches!(atom.element.as_str(), "D" | "T") {
        atom.name.replacen(atom.element.as_str(), "H", 1)
    } else {
        atom.name.as_str().to_owned()
    };
    if let Some(parent) = parents.get(&name) {
        return Ok(Some(parent));
    }
    let Some(reference) = nucleotide_reference(&atom.component, &name) else {
        return Ok(None);
    };
    let mut found: Option<&str> = None;
    for (child, parent) in parents {
        if nucleotide_reference(&atom.component, child).as_ref() == Some(&reference) {
            if let Some(previous) = found
                && nucleotide_reference(&atom.component, previous)
                    != nucleotide_reference(&atom.component, parent)
            {
                return Err(parent_error(
                    atom,
                    parent,
                    "nucleotide aliases have conflicting parents",
                ));
            }
            found = Some(parent);
        }
    }
    Ok(found)
}

fn dictionary_parent(
    atom: &Atom,
    name: &str,
    text: &ShortText,
    heavy: &[(usize, &ShortText)],
) -> Result<usize, String> {
    if let Some(&(row, _)) = heavy.iter().find(|(_, known)| **known == *text) {
        return Ok(row);
    }
    if let Some(reference) = nucleotide_reference(&atom.component, name) {
        let mut matching = heavy.iter().filter(|(_, known)| {
            nucleotide_reference(&atom.component, known).as_ref() == Some(&reference)
        });
        if let Some(&(row, _)) = matching.next() {
            if matching.next().is_some() {
                return Err(parent_error(
                    atom,
                    name,
                    "Dictionary parent alias is ambiguous",
                ));
            }
            return Ok(row);
        }
    }
    Err(parent_error(atom, name, "Dictionary parent is absent"))
}

/// The heavy atoms of one group with their names, in row order, read once per group.
struct HeavyNames<'a> {
    group: usize,
    atoms: Vec<(usize, &'a ShortText)>,
}

impl<'a> HeavyNames<'a> {
    fn of(
        &mut self,
        source: &'a Source,
        group: usize,
        rows: &[usize],
        hydrogen: &[bool],
    ) -> &[(usize, &'a ShortText)] {
        if self.group != group {
            self.group = group;
            self.atoms.clear();
            self.atoms.extend(
                rows.iter()
                    .filter(|&&row| !hydrogen[row])
                    .map(|&row| (row, &source.atoms[row].name)),
            );
        }
        &self.atoms
    }
}

/// Hydrogen names mapped to their dictionary parent names, for one component.
type ParentNames = HashMap<String, String>;

/// Hydrogen name and element, and the parent name as `str` and as inline text.
type ParentEntry<'a> = (&'a ShortText, &'a ShortText, Option<(&'a str, ShortText)>);

/// The dictionary parent names of one component's hydrogens, with the names already
/// looked up. Copies of a residue list their hydrogens in the same order, so a search that
/// starts after the last match usually succeeds at once and skips hashing the name.
struct ParentTable<'a> {
    map: &'a ParentNames,
    component: &'a ShortText,
    seen: Vec<ParentEntry<'a>>,
    cursor: usize,
}

impl<'a> ParentTable<'a> {
    /// Names past this count go straight to the hash map, which bounds the search.
    const LIMIT: usize = 256;
    /// The number of dictionaries to keep tables for.
    const TABLES: usize = 64;

    /// The result of `ccd_parent` for `atom`, which must belong to `component`.
    fn parent(&mut self, atom: &'a Atom) -> Result<Option<(&'a str, ShortText)>, String> {
        let count = self.seen.len();
        if let Some(index) = (self.cursor..count)
            .chain(0..self.cursor)
            .find(|&index| *self.seen[index].0 == atom.name && *self.seen[index].1 == atom.element)
        {
            self.cursor = index + 1;
            return Ok(self.seen[index].2.clone());
        }
        let parent = ccd_parent(self.map, atom)?.map(|name| (name, ShortText::new(name)));
        if count < Self::LIMIT {
            self.seen.push((&atom.name, &atom.element, parent.clone()));
            self.cursor = count + 1;
        }
        Ok(parent)
    }
}

/// The heavy parent of hydrogen `row` in group `group`, whose atoms are `rows`.
///
/// `dictionary` gives the dictionary parent name, as `ccd_parent` does. `heavy` caches the
/// heavy atoms of the group with their positions once the nearest-atom rule needs them.
#[allow(clippy::too_many_arguments)]
fn ligand_parent<'a>(
    source: &'a Source,
    dictionary: impl FnOnce() -> Result<Option<(&'a str, ShortText)>, String>,
    declared: &FastMap<usize, Option<usize>>,
    row: usize,
    group: usize,
    rows: &[usize],
    hydrogen: &[bool],
    heavy: &mut Option<Vec<(usize, [f32; 3])>>,
    names: &mut HeavyNames<'a>,
) -> Result<Option<usize>, String> {
    let atom = &source.atoms[row];
    let parent = if let Some(&parent) = declared.get(&row) {
        parent.ok_or_else(|| parent_error(atom, "struct_conn", "bonded parent is absent"))?
    } else if let Some((name, text)) = dictionary()? {
        dictionary_parent(atom, name, &text, names.of(source, group, rows, hydrogen))?
    } else {
        let heavy = heavy.get_or_insert_with(|| {
            rows.iter()
                .copied()
                .filter(|&i| !hydrogen[i])
                .map(|i| (i, source.atoms[i].position))
                .collect()
        });
        return Ok(residue::nearest(atom.position, heavy));
    };
    let target = &source.atoms[parent];
    let distance = residue::distance(atom.position, target.position);
    if !(residue::finite(atom.position)
        && residue::finite(target.position)
        && distance > 0.5
        && distance < 1.6
        && !target.protein)
    {
        return Err(parent_error(
            atom,
            &target.name,
            "bonded parent has invalid geometry or is a protein atom",
        ));
    }
    Ok(Some(parent))
}

fn add_proteins(
    source: &Source,
    schema: &Schema,
    complete: bool,
    protein_only: bool,
    policy: &str,
    tokens: &mut Tokens,
    mapping: &mut [usize],
) -> Result<(), String> {
    let topology = if protein_only {
        topology::protein(source, policy)?
    } else {
        topology::pdb(source)
    };
    let mut groups: Vec<_> = source
        .groups
        .iter()
        .enumerate()
        .filter_map(|(group, rows)| {
            rows.iter()
                .copied()
                .find(|&row| {
                    source.atoms[row].protein
                        && source.atoms[row].name == "N"
                        && residue::finite(source.atoms[row].position)
                })
                .map(|center| (center, group))
        })
        .collect();
    if protein_only {
        groups = source
            .groups
            .iter()
            .enumerate()
            .filter_map(|(group, rows)| {
                rows.iter()
                    .copied()
                    .find(|&row| source.atoms[row].protein)
                    .map(|first| (first, group))
            })
            .collect();
    }
    groups.sort_unstable();
    tokens.reserve(groups.len());
    // Slot names as `ShortText`, so that name tests compare fixed-size arrays.
    let names: Vec<Vec<ShortText>> = schema
        .residues
        .iter()
        .map(|residue| {
            residue
                .names
                .iter()
                .map(|name| ShortText::new(name))
                .collect()
        })
        .collect();
    let mut atoms: Vec<&Atom> = Vec::new();
    for (_, group) in groups {
        let rows = &source.groups[group];
        let first = &source.atoms[rows[0]];
        let aa = schema.identity(&first.component);
        if protein_only && aa >= 20 {
            continue;
        }
        atoms.clear();
        atoms.extend(rows.iter().map(|&row| &source.atoms[row]));
        let mut packed = if aa < 20 {
            residue::pack(
                &schema.residues[aa],
                &names[aa],
                &atoms,
                topology[group],
                (!protein_only && topology[group][2] == 1).then_some(true),
                complete,
                true,
            )
        } else {
            residue::Packed::empty([-1; 3])
        };
        if aa >= 20 {
            for atom in &atoms {
                if let Some(slot) = ["N", "CA", "C", "O"]
                    .iter()
                    .position(|&name| name == atom.name)
                    && residue::finite(atom.position)
                {
                    packed.position[slot] = atom.position;
                    packed.mask[slot] = true;
                }
            }
        }
        if protein_only && !packed.mask[1] {
            continue;
        }
        if !protein_only && first.sequence < 1 {
            return Err("Protein token lacks a positive label_seq_id".into());
        }
        let token = tokens.identity.len();
        for &row in rows {
            let atom = &source.atoms[row];
            if !(atom.element == "H" || atom.element == "D") && residue::finite(atom.position) {
                mapping[row] = token;
            }
        }
        tokens.push(first, aa, &packed, false, [0; 4]);
    }
    Ok(())
}

/// The dictionary features of one ligand component, with the atom names already looked up.
/// Copies of a residue list their atoms in the same order, so a search that starts after
/// the last match usually succeeds at once and skips hashing the name.
struct FeatureTable<'a> {
    component: &'a ShortText,
    names: Option<&'a HashMap<String, [i16; 4]>>,
    seen: Vec<(&'a ShortText, [i16; 4])>,
    cursor: usize,
}

impl<'a> FeatureTable<'a> {
    /// Names past this count go straight to the hash map, which bounds the search.
    const LIMIT: usize = 256;
    /// The number of components to keep tables for.
    const TABLES: usize = 64;

    fn feature(&mut self, name: &'a ShortText) -> [i16; 4] {
        let count = self.seen.len();
        if let Some(index) = (self.cursor..count)
            .chain(0..self.cursor)
            .find(|&index| *self.seen[index].0 == *name)
        {
            self.cursor = index + 1;
            return self.seen[index].1;
        }
        let feature = self
            .names
            .and_then(|names| names.get(name.as_str()))
            .copied()
            .unwrap_or([0; 4]);
        if count < Self::LIMIT {
            self.seen.push((name, feature));
            self.cursor = count + 1;
        }
        feature
    }
}

#[allow(clippy::too_many_arguments)]
fn add_ligands(
    source: &Source,
    features: &LigandFeatures,
    parents: &HydrogenParents,
    residue_parents: &ResidueParents,
    #[cfg(feature = "reduce3")] bridge_sites: &[super::bridges::Site],
    token_table: &TokenTable,
    excluded: &[String],
    tokens: &mut Tokens,
    mapping: &mut [usize],
    complete: bool,
) -> Result<(), String> {
    #[cfg(feature = "reduce3")]
    let bridge_rows = bridge_sites
        .iter()
        .filter_map(|site| site.hydrogen_row)
        .collect();
    #[cfg(not(feature = "reduce3"))]
    let bridge_rows = BTreeSet::new();
    let declared = declared_parents(source, &bridge_rows)?;
    // Atoms of one residue share a component, so test each run of one component once.
    let mut previous: Option<(&ShortText, bool)> = None;
    // Per row: whether the atom is a hydrogen, and whether it is a ligand hydrogen that
    // needs a parent. Per group: whether it holds such a hydrogen.
    let mut hydrogen: Vec<bool> = Vec::with_capacity(source.atoms.len());
    let mut orphan: Vec<bool> = Vec::with_capacity(source.atoms.len());
    let mut pending = vec![false; source.groups.len()];
    // The loop below pushes each heavy ligand atom as it reads it, so it reads every atom
    // once. Hydrogens make up about half of a protonated structure, so half the atoms
    // bounds the new tokens of most inputs; more tokens only grow the arrays.
    tokens.reserve(source.atoms.len() / 2);
    // Rows of the heavy ligand atoms, in order, with their dictionary features.
    let mut ligand_rows: Vec<(usize, [i16; 4])> = Vec::with_capacity(source.atoms.len() / 2);
    let mut identities: Vec<(&ShortText, usize)> = Vec::new();
    let mut tables: Vec<FeatureTable> = Vec::new();
    let mut current = 0;
    for (row, atom) in source.atoms.iter().enumerate() {
        let skipped = match previous {
            Some((component, flag)) if *component == atom.component => flag,
            _ => {
                let flag = excluded.iter().any(|name| *name == atom.component);
                previous = Some((&atom.component, flag));
                flag
            }
        };
        let is_hydrogen = residue::hydrogen(&atom.element);
        hydrogen.push(is_hydrogen);
        let needs_parent = !atom.protein && is_hydrogen && !skipped && !bridge_rows.contains(&row);
        orphan.push(needs_parent);
        pending[atom.group] |= needs_parent;
        if atom.protein || is_hydrogen || !residue::finite(atom.position) || skipped {
            continue;
        }
        mapping[row] = tokens.identity.len();
        let identity = match identities
            .iter()
            .find(|(element, _)| **element == atom.element)
        {
            Some(&(_, identity)) => identity,
            None => {
                // The table checks that every token lies above the residue tokens.
                let identity = token_table.element(atom.element.as_str()) as usize;
                identities.push((&atom.element, identity));
                identity
            }
        };
        if tables
            .get(current)
            .is_none_or(|table| *table.component != atom.component)
        {
            current = match tables
                .iter()
                .position(|table| *table.component == atom.component)
            {
                Some(index) => index,
                None => {
                    let names = features.get(atom.component.as_str());
                    let table = FeatureTable {
                        component: &atom.component,
                        names,
                        seen: Vec::new(),
                        cursor: 0,
                    };
                    // Past `TABLES` components, reuse the last place so the search stays short.
                    if tables.len() == FeatureTable::TABLES {
                        tables.pop();
                    }
                    tables.push(table);
                    tables.len() - 1
                }
            };
        }
        let feature = tables[current].feature(&atom.name);
        ligand_rows.push((row, feature));
        tokens.push_ligand(atom, identity, feature);
        if complete {
            *tokens
                .hydrogen_known
                .last_mut()
                .ok_or("Missing ligand token")? = 1;
        }
    }
    #[cfg(feature = "reduce3")]
    let bridge_counts = add_bridges(
        source,
        bridge_sites,
        excluded,
        token_table.hydrogen() as usize,
        tokens,
        mapping,
    )?;
    #[cfg(not(feature = "reduce3"))]
    let bridge_counts: HashMap<usize, usize> = HashMap::new();
    // Hydrogen rows given a parent, in order, and each one's parent.
    let mut assigned: Vec<usize> = Vec::new();
    let mut parent_of = vec![usize::MAX; source.atoms.len()];
    let mut names = HeavyNames {
        group: usize::MAX,
        atoms: Vec::new(),
    };
    let mut tables: Vec<ParentTable> = Vec::new();
    for (group, rows) in source.groups.iter().enumerate() {
        if !pending[group] {
            continue;
        }
        let mut heavy = None;
        // Atoms of one group share a component, so this map rarely changes within it.
        let mut dictionary: Option<(&ShortText, Option<&ParentNames>, usize)> = None;
        for &row in rows.iter().filter(|&&row| orphan[row]) {
            let atom = &source.atoms[row];
            let (map, table) = match dictionary {
                Some((component, map, table)) if *component == atom.component => (map, table),
                _ => {
                    let map = residue_parents
                        .get(&atom.group)
                        .map(AsRef::as_ref)
                        .or_else(|| parents.get(atom.component.as_str()));
                    let table = map.map_or(usize::MAX, |map| {
                        let known = tables.iter().position(|table| {
                            std::ptr::eq(table.map, map) && *table.component == atom.component
                        });
                        known.unwrap_or_else(|| {
                            // Past `TABLES` dictionaries, reuse the last place so the search stays short.
                            if tables.len() == ParentTable::TABLES {
                                tables.pop();
                            }
                            tables.push(ParentTable {
                                map,
                                component: &atom.component,
                                seen: Vec::new(),
                                cursor: 0,
                            });
                            tables.len() - 1
                        })
                    });
                    dictionary = Some((&atom.component, map, table));
                    (map, table)
                }
            };
            let lookup = || match tables.get_mut(table) {
                Some(table) => table.parent(atom),
                None => Ok(None),
            };
            let parent = ligand_parent(
                source, lookup, &declared, row, group, rows, &hydrogen, &mut heavy, &mut names,
            )?;
            if let Some(parent) = parent
                && mapping[parent] == usize::MAX
                && (declared.contains_key(&row)
                    || map
                        .map(|map| ccd_parent(map, atom))
                        .transpose()?
                        .flatten()
                        .is_some())
            {
                return Err(parent_error(
                    &source.atoms[row],
                    &source.atoms[parent].name,
                    "bonded parent has no ligand token",
                ));
            }
            let parent = parent.filter(|&parent| mapping[parent] != usize::MAX);
            let Some(parent) = parent else {
                for &other in rows {
                    if mapping[other] != usize::MAX && !hydrogen[other] {
                        tokens.hydrogen_known[mapping[other]] = 0;
                    }
                }
                continue;
            };
            let token = mapping[parent];
            parent_of[row] = parent;
            assigned.push(row);
            let count = tokens.masks[token * 32 + 15..token * 32 + 18]
                .iter()
                .filter(|&&value| value != 0)
                .count();
            if count >= 3 {
                let rows: Vec<usize> = assigned
                    .iter()
                    .copied()
                    .filter(|&other| parent_of[other] == parent)
                    .collect();
                return Err(hydrogen_overflow(source, features, parent, &rows));
            }
            let slot = 15 + count;
            tokens.masks[token * 32 + slot] = 1;
            tokens.positions[token * 96 + slot * 3..][..3]
                .copy_from_slice(&source.atoms[row].position);
        }
    }
    // Only heavy ligand atoms map to ligand tokens; bridge tokens map hydrogens.
    for &(row, feature) in &ligand_rows {
        let token = mapping[row];
        if token == usize::MAX || tokens.ligand[token] == 0 {
            continue;
        }
        let count = tokens.masks[token * 32 + 15..(token + 1) * 32]
            .iter()
            .filter(|&&value| value != 0)
            .count()
            + bridge_counts.get(&row).copied().unwrap_or(0);
        tokens.charge_known[token] = u8::from(
            tokens.hydrogen_known[token] != 0 && feature[3] != 0 && count == feature[2] as usize,
        );
    }
    Ok(())
}

/// Whether no two keys are equal.
///
/// Keys of one residue are adjacent and share their first four fields. When no such
/// residue part starts two separate runs, two equal keys must lie in one run, so it is
/// enough to compare the names within each run. Copies of one residue type list the same
/// names, so a run whose names equal those of a recent run takes its answer.
fn distinct_keys(keys: &[[ShortText; 5]]) -> bool {
    // Find the runs first, so that the set never grows.
    let mut starts: Vec<usize> = (0..keys.len())
        .filter(|&row| row == 0 || !same_residue(&keys[row], &keys[row - 1]))
        .collect();
    starts.push(keys.len());
    let mut residues: FastSet<&[ShortText]> =
        FastSet::with_capacity_and_hasher(starts.len(), Default::default());
    // Recent runs with distinct names.
    let mut recent: Vec<&[[ShortText; 5]]> = Vec::new();
    for range in starts.windows(2) {
        let run = &keys[range[0]..range[1]];
        if !residues.insert(&run[0][..4]) {
            let mut seen = FastSet::with_capacity_and_hasher(keys.len(), Default::default());
            return keys.iter().all(|key| seen.insert(key));
        }
        if run.len() < 2 {
            continue;
        }
        let same = |other: &&[[ShortText; 5]]| {
            other.len() == run.len()
                && other
                    .iter()
                    .zip(run)
                    .all(|(left, right)| left[4] == right[4])
        };
        if run.len() >= 8 && recent.iter().any(same) {
            continue;
        }
        let distinct = if run.len() <= 64 {
            run.iter()
                .enumerate()
                .all(|(index, key)| run[index + 1..].iter().all(|other| other[4] != key[4]))
        } else {
            let mut names = FastSet::with_capacity_and_hasher(run.len(), Default::default());
            run.iter().all(|key| names.insert(&key[4]))
        };
        if !distinct {
            return false;
        }
        if run.len() >= 8 {
            if recent.len() == 16 {
                recent.pop();
            }
            recent.push(run);
        }
    }
    true
}

#[allow(clippy::too_many_arguments)]
pub(super) fn pack(
    source: &Source,
    schema: &Schema,
    features: &LigandFeatures,
    parents: &HydrogenParents,
    residue_parents: &ResidueParents,
    #[cfg(feature = "reduce3")] bridge_sites: &[super::bridges::Site],
    token_table: &TokenTable,
    excluded: &[String],
    options: &super::PackOptions<'_>,
) -> Result<Tokens, String> {
    let super::PackOptions {
        complete,
        topology_policy: policy,
        protein_only,
        ..
    } = *options;
    if protein_only && !source.has_protein_annotations {
        return Err("Protein export requires entity_poly protein annotations".into());
    }
    let mut result = Tokens::empty();
    let mut mapping = vec![usize::MAX; source.atoms.len()];
    add_proteins(
        source,
        schema,
        complete,
        protein_only,
        policy,
        &mut result,
        &mut mapping,
    )?;
    if !protein_only {
        add_ligands(
            source,
            features,
            parents,
            residue_parents,
            #[cfg(feature = "reduce3")]
            bridge_sites,
            token_table,
            excluded,
            &mut result,
            &mut mapping,
            complete,
        )?;
    }
    // Protein exports number chains by first appearance, others by sorted name.
    let mut chains: FastMap<&str, i64> = FastMap::default();
    if protein_only {
        for key in &result.keys {
            let next = chains.len() as i64;
            chains.entry(key[0].as_str()).or_insert(next);
        }
    } else {
        // Chains repeat over long runs of tokens, so skip each repeat before sorting.
        let mut names: Vec<&ShortText> = Vec::new();
        for key in &result.keys {
            if names.last().is_none_or(|last| **last != key[0]) {
                names.push(&key[0]);
            }
        }
        names.sort_unstable();
        names.dedup();
        chains.extend(
            names
                .into_iter()
                .enumerate()
                .map(|(index, name)| (name.as_str(), index as i64)),
        );
    }
    let mut previous: Option<(&ShortText, i64)> = None;
    for (row, key) in result.keys.iter().enumerate() {
        let chain = match previous {
            Some((last, chain)) if *last == key[0] => chain,
            _ => *chains.get(key[0].as_str()).ok_or("Missing token chain")?,
        };
        previous = Some((&key[0], chain));
        result.chains[row] = chain;
    }
    if !protein_only {
        let mut pairs = BTreeSet::new();
        for &(a, b) in &source.bonds {
            let (left, right) = (mapping[a], mapping[b]);
            if left == usize::MAX
                || right == usize::MAX
                || source.atoms[a].protein == source.atoms[b].protein
            {
                continue;
            }
            pairs.insert(if source.atoms[a].protein {
                (right, left)
            } else {
                (left, right)
            });
        }
        result.attachments = pairs.into_iter().collect();
        for &(ligand, _) in &result.attachments {
            result.covalent[ligand] = 1;
        }
    }
    super::continuity::apply(source, &mapping, &mut result);
    renumber_ligands(&mut result);
    if !distinct_keys(&result.keys) {
        return Err("PDB tokens have duplicate mmCIF identities".into());
    }
    if protein_only && result.identity.is_empty() {
        return Err("No canonical protein residues with CA coordinates".into());
    }
    Ok(result)
}

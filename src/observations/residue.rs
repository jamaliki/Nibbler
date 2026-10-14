//! Pack observed atoms without inferring missing coordinates.

use super::schema::Residue;
use super::source::Atom;

pub(super) struct Packed {
    pub position: [[f32; 3]; 32],
    pub mask: [bool; 32],
    pub presence: [bool; 32],
    pub topology: [i64; 3],
    pub state: i64,
    pub state_known: bool,
}

impl Packed {
    pub fn empty(topology: [i64; 3]) -> Self {
        Self {
            position: [[0.; 3]; 32],
            mask: [false; 32],
            presence: [false; 32],
            topology,
            state: 0,
            state_known: false,
        }
    }
}

pub(super) fn hydrogen(element: &str) -> bool {
    matches!(element, "H" | "D" | "T")
}

fn protein_hydrogen(element: &str) -> bool {
    element.eq_ignore_ascii_case("H")
        || element.eq_ignore_ascii_case("D")
        || element.eq_ignore_ascii_case("T")
}

pub(super) fn finite(xyz: [f32; 3]) -> bool {
    xyz.iter().all(|value| value.is_finite())
}

pub(super) fn distance(first: [f32; 3], second: [f32; 3]) -> f32 {
    let [x, y, z] = std::array::from_fn(|axis| first[axis] - second[axis]);
    ((x * x + y * y) + z * z).sqrt()
}

pub(super) fn nearest(position: [f32; 3], heavy: &[(usize, [f32; 3])]) -> Option<usize> {
    if !finite(position) {
        return None;
    }
    let mut first = (usize::MAX, f32::INFINITY);
    let mut second = f32::INFINITY;
    for &(index, xyz) in heavy {
        let value = distance(xyz, position);
        if !value.is_finite() {
            return None;
        }
        if value < first.1 {
            second = first.1;
            first = (index, value);
        } else if value < second {
            second = value;
        }
    }
    (first.1 > 0.5 && first.1 < 1.6 && second - first.1 >= 0.2).then_some(first.0)
}

fn protein_parent(schema: &Residue, atom: &Atom, heavy: &[(usize, [f32; 3])]) -> Option<usize> {
    let named_parent = |name: &str| {
        schema
            .names
            .iter()
            .zip(&schema.parents)
            .find_map(|(slot, &parent)| (slot == name && parent >= 0).then_some(parent as usize))
    };
    let parent = named_parent(&atom.name)
        .or_else(|| match atom.name.as_str() {
            "H1" | "D1" | "T1" => named_parent("H"),
            _ => None,
        })
        .or_else(|| {
            if atom.element.eq_ignore_ascii_case("D") || atom.element.eq_ignore_ascii_case("T") {
                named_parent(
                    &atom
                        .name
                        .replacen(&atom.element.to_ascii_uppercase(), "H", 1),
                )
            } else {
                None
            }
        });
    if let Some(parent) = parent {
        return heavy.iter().find_map(|&(slot, xyz)| {
            let value = distance(atom.position, xyz);
            (slot == parent && value > 0.5 && value < 1.6).then_some(slot)
        });
    }
    nearest(atom.position, heavy)
}

pub(super) fn pack(
    schema: &Residue,
    atoms: &[&Atom],
    topology: [i64; 3],
    linked: Option<bool>,
    complete: bool,
    constrain_termini: bool,
) -> Packed {
    let mut result = Packed::empty(topology);
    let mut invalid = false;
    for atom in atoms.iter().filter(|atom| !protein_hydrogen(&atom.element)) {
        let slot = schema.names[..15]
            .iter()
            .position(|name| !name.is_empty() && *name == atom.name);
        let Some(slot) = slot else {
            invalid = true;
            continue;
        };
        if result.mask[slot] || !finite(atom.position) {
            invalid = true;
            continue;
        }
        result.position[slot] = atom.position;
        result.mask[slot] = true;
    }
    let heavy: Vec<_> = (0..15)
        .filter(|&slot| result.mask[slot])
        .map(|slot| (slot, result.position[slot]))
        .collect();
    let mut groups: [Vec<&Atom>; 15] = std::array::from_fn(|_| Vec::new());
    for atom in atoms.iter().filter(|atom| protein_hydrogen(&atom.element)) {
        let Some(parent) = protein_parent(schema, atom, &heavy) else {
            invalid = true;
            continue;
        };
        if !schema.parents.contains(&(parent as i64)) {
            invalid = true;
            continue;
        }
        groups[parent].push(atom);
    }
    for (parent, group) in groups.iter_mut().enumerate() {
        let slots: Vec<_> = schema
            .parents
            .iter()
            .enumerate()
            .filter_map(|(slot, &value)| (value == parent as i64).then_some(slot))
            .collect();
        if group.len() > slots.len() {
            invalid = true;
            continue;
        }
        group.sort_by(|first, second| first.name.cmp(&second.name));
        for (&slot, atom) in slots.iter().zip(group.iter()) {
            result.position[slot] = atom.position;
            result.mask[slot] = true;
        }
    }
    result.presence = result.mask;
    if !complete || invalid {
        return result;
    }
    let supported = schema.states.iter().any(|state| {
        (15..32)
            .filter(|&slot| schema.parents[slot] >= 4 && schema.parents[slot] < 14)
            .all(|slot| state.present[slot] == result.mask[slot])
    });
    if supported && (4..14).all(|slot| schema.names[slot].is_empty() || result.mask[slot]) {
        for slot in 15..32 {
            if schema.parents[slot] >= 4 && schema.parents[slot] < 14 {
                result.presence[slot] = true;
            }
        }
    }
    let matches: Vec<_> = schema
        .states
        .iter()
        .enumerate()
        .filter(|(_, state)| {
            linked.is_none_or(|value| state.linked == value)
                && state.present.iter().zip(result.mask).all(|(&a, b)| a == b)
                && (!constrain_termini || state.topology()[..2] == topology[..2])
        })
        .collect();
    if matches.len() == 1 {
        let (index, state) = matches[0];
        result.state = index as i64;
        result.state_known = true;
        result.topology = state.topology();
        result.presence.fill(true);
    }
    result
}

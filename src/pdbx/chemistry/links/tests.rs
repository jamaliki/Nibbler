//! Tests for link inference from distances and deposited connections.

use super::super::{BondSource, PreparedChemistry, prepare_product};
use crate::cif::parse;
use crate::pdbx::build_model;

/// A Gly-Gly dipeptide, as the header of an mmCIF block that lacks only `atom_site`.
const HEADER: &str = "data_link
_entry.id LINK
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'GG peptide'
loop_
_entity_poly.entity_id
_entity_poly.type
1 polypeptide(L)
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1 GLY n
1 2 GLY n
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_chem_comp.id
_chem_comp.type
GLY 'L-peptide linking'
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_leaving_atom_flag
GLY N N 0 N
GLY CA C 0 N
GLY C C 0 N
GLY O O 0 N
GLY OXT O 0 Y
GLY H H 0 N
GLY H2 H 0 Y
GLY HA2 H 0 N
GLY HA3 H 0 N
GLY HXT H 0 Y
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
GLY N CA SING
GLY N H SING
GLY N H2 SING
GLY CA C SING
GLY CA HA2 SING
GLY CA HA3 SING
GLY C O DOUB
GLY C OXT SING
GLY OXT HXT SING
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
";

/// The first residue of every pair below, as `(N, CA, C)`.
const RESIDUE_1: [[f64; 3]; 3] = [[-0.523, 1.363, 0.], [0., 0., 0.], [1.52, 0., 0.]];

// Second residues, as `(N, CA, C)`, with the C(1) to N(2) and CA-CA distances in A. The
// same pairs as the chain-break tests of `observations::continuity`.
const BONDED: [[f64; 3]; 3] = [[2.103, 1.195, 0.], [3.556, 1.34, 0.], [4.443, 0.106, 0.]]; // C-N 1.33, CA-CA 3.80
const BONDED_STRETCHED: [[f64; 3]; 3] =
    [[2.357, 1.072, 0.], [3.795, 1.326, 0.], [4.773, 0.162, 0.]]; // C-N 1.36, CA-CA 4.02
const BROKEN_CLOSE: [[f64; 3]; 3] = [[2.835, 2.696, 0.], [3.528, 1.411, 0.], [2.73, 0.118, 0.]]; // C-N 3.00, CA-CA 3.80
const BROKEN_FAR: [[f64; 3]; 3] = [[3.054, 3.146, 0.], [4.48, 2.832, 0.], [4.939, 1.383, 0.]]; // C-N 3.50, CA-CA 5.30
const LINK_2_39: [[f64; 3]; 3] = [[2.568, 2.148, 0.], [3.702, 1.228, 0.], [3.462, -0.273, 0.]]; // C-N 2.39, CA-CA 3.90
const LINK_2_41: [[f64; 3]; 3] = [[2.576, 2.166, 0.], [3.7, 1.233, 0.], [3.444, -0.265, 0.]]; // C-N 2.41, CA-CA 3.90

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}

const NONE: &[usize] = &[];
const N: usize = 0;
const C: usize = 2;

/// The product of the dipeptide, with backbone atoms of either residue left unobserved
/// where `drop_left` and `drop_right` list their `(N, CA, C)` slots.
fn dipeptide(
    second: [[f64; 3]; 3],
    drop_left: &[usize],
    drop_right: &[usize],
) -> PreparedChemistry {
    let mut text = HEADER.to_owned();
    let mut id = 0;
    for (number, atoms, drop) in [(1, RESIDUE_1, drop_left), (2, second, drop_right)] {
        for (slot, name) in ["N", "CA", "C"].into_iter().enumerate() {
            if !drop.contains(&slot) {
                let [x, y, z] = atoms[slot];
                id += 1;
                let element = &name[..1];
                text += &format!(
                    "ATOM {id} {element} {name} . GLY A 1 {number} ? {x} {y} {z} 1.00 10.00 ? {number} GLY A {name} 1\n"
                );
            }
        }
    }
    let model =
        build_model(&parse(text.as_bytes()).expect("fixture parses")).expect("fixture builds");
    prepare_product(&model, None, 1, None, false).expect("fixture gives a product")
}

/// Whether the product joins C of the first residue to N of the second by a sequence link.
fn linked(product: &PreparedChemistry) -> bool {
    let site = |number: &str, name: &str| {
        product
            .atoms()
            .iter()
            .position(|atom| atom.auth_sequence_id == number && atom.atom_id == name)
            .expect("backbone atom")
    };
    let (c, n) = (site("1", "C"), site("2", "N"));
    product.bonds().iter().any(|bond| {
        bond.active
            && bond.source == BondSource::PolymerSequence
            && bond.atom_indices.contains(&c)
            && bond.atom_indices.contains(&n)
    })
}

/// Whether a diagnostic says that the pair was left unlinked.
fn excluded(product: &PreparedChemistry) -> bool {
    product
        .diagnostics()
        .iter()
        .any(|text| text.starts_with("Excluded distant"))
}

#[test]
fn fixtures_have_the_distances_they_are_named_for() {
    for (second, link, ca) in [
        (BONDED, 1.33, 3.80),
        (BONDED_STRETCHED, 1.36, 4.02),
        (BROKEN_CLOSE, 3.00, 3.80),
        (BROKEN_FAR, 3.50, 5.30),
        (LINK_2_39, 2.39, 3.90),
        (LINK_2_41, 2.41, 3.90),
    ] {
        assert!((distance(RESIDUE_1[2], second[0]) - link).abs() < 0.005);
        assert!((distance(RESIDUE_1[1], second[1]) - ca).abs() < 0.005);
    }
}

#[test]
fn a_peptide_bond_is_linked_whatever_its_ca_distance() {
    // 1SFV ILE13 to PRO14 has C-N 1.36 A and CA-CA 4.02 A: bonded, so linked.
    for second in [BONDED, BONDED_STRETCHED] {
        let product = dipeptide(second, NONE, NONE);
        assert!(linked(&product));
        assert!(!excluded(&product), "{:?}", product.diagnostics());
    }
}

#[test]
fn a_linked_pair_loses_its_terminal_atoms() {
    let atoms = |product: &PreparedChemistry, number: &str, name: &str| {
        product
            .atoms()
            .iter()
            .find(|atom| atom.auth_sequence_id == number && atom.atom_id == name)
            .map(|atom| atom.active)
    };
    let stretched = dipeptide(BONDED_STRETCHED, NONE, NONE);
    let broken = dipeptide(BROKEN_CLOSE, NONE, NONE);
    // The carboxylate oxygen and one amino hydrogen leave when the link forms.
    assert_eq!(atoms(&stretched, "1", "OXT"), Some(false));
    assert_eq!(atoms(&stretched, "2", "H2"), Some(false));
    assert_eq!(atoms(&broken, "1", "OXT"), Some(true));
    assert_eq!(atoms(&broken, "2", "H2"), Some(true));
}

#[test]
fn a_missing_peptide_bond_is_not_linked_whether_ca_atoms_are_near_or_far() {
    // CA-CA below 4.0 A is still no link when C and N are 3.0 A apart.
    for second in [BROKEN_CLOSE, BROKEN_FAR] {
        let product = dipeptide(second, NONE, NONE);
        assert!(!linked(&product));
        assert!(excluded(&product));
    }
}

#[test]
fn the_peptide_link_limit_is_2_4_angstroms_between_c_and_n() {
    // Both pairs have the same CA-CA, 3.9 A, so only C-N tells them apart.
    assert!(linked(&dipeptide(LINK_2_39, NONE, NONE)));
    assert!(!linked(&dipeptide(LINK_2_41, NONE, NONE)));
}

#[test]
fn an_unobserved_c_or_n_falls_back_to_the_ca_distance() {
    for (left, right) in [(NONE, &[N][..]), (&[C][..], NONE), (&[C][..], &[N][..])] {
        assert!(linked(&dipeptide(BONDED, left, right)));
        // Without the bond atoms nothing shows the 4.02 A pair is bonded.
        assert!(!linked(&dipeptide(BONDED_STRETCHED, left, right)));
        assert!(!linked(&dipeptide(BROKEN_FAR, left, right)));
    }
}

#[test]
fn only_the_linking_atoms_matter_when_other_backbone_atoms_are_missing() {
    // The left residue has no N and the right no C, but its C and N are both observed.
    assert!(linked(&dipeptide(BONDED_STRETCHED, &[N], &[C])));
    assert!(!linked(&dipeptide(BROKEN_CLOSE, &[N], &[C])));
}

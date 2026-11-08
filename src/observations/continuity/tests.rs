//! Tests for sequence gaps kept and added at polymer breaks.

use super::*;
use crate::observations::{residue::Packed, source};

fn input(kind: &str, rows: &str) -> Source {
    let raw = format!(
        "data_gap\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n1 '{kind}'\n\
        loop_\n_atom_site.label_entity_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n\
        _atom_site.auth_seq_id\n_atom_site.label_comp_id\n_atom_site.label_atom_id\n\
        _atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
        _atom_site.pdbx_PDB_model_num\n{rows}"
    );
    source::read(
        &crate::cif::parse(raw.as_bytes()).expect("Fixture must parse"),
        1,
    )
    .expect("Fixture must read")
}

fn packed(source: &Source, protein: bool) -> Tokens {
    let mut tokens = Tokens::empty();
    let mut chains = Vec::new();
    let mut mapping = vec![usize::MAX; source.atoms.len()];
    for rows in &source.groups {
        for &row in rows {
            mapping[row] = tokens.identity.len();
            let mut packed = Packed::empty([-1; 3]);
            if protein {
                // One token per residue, with the backbone atoms the fixture lists.
                for &atom in rows {
                    let slot = ["N", "CA", "C"]
                        .iter()
                        .position(|name| source.atoms[atom].name == *name);
                    if let Some(slot) = slot {
                        packed.mask[slot] = true;
                        packed.position[slot] = source.atoms[atom].position;
                    }
                }
            }
            tokens.push(&source.atoms[row], 0, &packed, !protein, [0; 4]);
            let chain = &source.atoms[row].chain;
            if !chains.contains(chain) {
                chains.push(chain.clone());
            }
            let token = tokens.identity.len() - 1;
            tokens.chains[token] = chains
                .iter()
                .position(|name| name == chain)
                .expect("Fixture chain must exist") as i64;
            if protein {
                break;
            }
        }
    }
    apply(source, &mapping, &mut tokens);
    super::super::record::renumber_ligands(&mut tokens);
    retained_proteins(&mut tokens);
    tokens
}

#[test]
fn protein_gaps_accumulate_without_changing_source_keys_or_existing_gaps() {
    let source = input(
        "polypeptide(L)",
        concat!(
            "1 A 1 101 ALA CA C 0 0 0 1\n",
            "1 A 2 102 ALA CA C 4 0 0 1\n",
            "1 A 3 103 ALA CA C 8.01 0 0 1\n",
            "1 A 4 104 ALA CA C 11.8 0 0 1\n",
            "1 A 8 108 ALA CA C 30 0 0 1\n",
            "1 A 9 109 ALA CA C 40 0 0 1\n",
            "1 B 1 101 ALA CA C 90 0 0 1\n",
            "1 B 2 102 ALA CA C 93.8 0 0 1\n"
        ),
    );
    let mut tokens = packed(&source, true);
    assert_eq!(tokens.sequence, [1, 2, 67, 68, 72, 137, 1, 2]);
    assert_eq!(tokens.keys[2][1], "3");
    tokens.select(&[1, 3]).expect("Crop must retain index gaps");
    assert_eq!(tokens.sequence, [2, 68]);
}

#[test]
fn nucleotide_gaps_survive_alignment_for_rna_dna_and_prime_aliases() {
    for kind in [
        "polyribonucleotide",
        "polydeoxyribonucleotide",
        "polydeoxyribonucleotide/polyribonucleotide hybrid",
    ] {
        for oxygen in ["O3'", "O3*"] {
            let source = input(
                kind,
                &format!(
                    "1 A 1 101 A P P 0 0 0 1\n1 A 1 101 A \"{oxygen}\" O 2 0 0 1\n\
                 1 A 2 102 A P P 3.6 0 0 1\n1 A 2 102 A \"{oxygen}\" O 6 0 0 1\n\
                 1 A 3 103 A P P 9 0 0 1\n1 A 3 103 A \"{oxygen}\" O 11 0 0 1\n\
                 1 A 5 105 A P P 12.6 0 0 1\n"
                ),
            );
            let mut tokens = packed(&source, false);
            assert_eq!(tokens.sequence, [0, 1, 2, 3, 68, 69, 70]);
            let expected = serde_json::to_string(&tokens.keys).expect("Keys must encode");
            tokens
                .align(&expected)
                .expect("Alignment must preserve polymer breaks");
            assert_eq!(tokens.sequence, [0, 1, 2, 3, 68, 69, 70]);
            assert_eq!(tokens.keys[4][1], "103");
        }
    }
}

#[test]
fn missing_nucleotide_link_atoms_break_adjacency_but_nonpolymers_do_not() {
    let rows = "1 A 1 1 A P P 0 0 0 1\n1 A 2 2 A P P 1.6 0 0 1\n";
    assert_eq!(
        packed(&input("polyribonucleotide", rows), false).sequence,
        [0, 65]
    );
    assert_eq!(packed(&input("other", rows), false).sequence, [0, 1]);
}

#[test]
fn a_larger_deposited_gap_stays_larger_than_the_clipping_range() {
    let source = input(
        "polypeptide(L)",
        "1 A 1 1 ALA CA C 0 0 0 1\n1 A 100 100 ALA CA C 90 0 0 1\n",
    );
    assert_eq!(packed(&source, true).sequence, [1, 100]);
}

#[test]
fn existing_small_gaps_stay_exactly_the_same_size() {
    let source = input(
        "polypeptide(L)",
        "1 A 61 245 ALA CA C 0 0 0 1\n1 A 76 260 ALA CA C 28 0 0 1\n",
    );
    assert_eq!(packed(&source, true).sequence, [61, 76]);
}

#[test]
fn completeness_keeps_existing_gaps_and_only_widens_new_consecutive_breaks() {
    let source = input(
        "polypeptide(L)",
        concat!(
            "1 A 1 1 ALA CA C 0 0 0 1\n1 A 2 2 ALA CA C 3.8 0 0 1\n",
            "1 A 3 3 ALA CA C 7.6 0 0 1\n1 A 4 4 ALA CA C 11.4 0 0 1\n",
            "1 A 5 5 ALA CA C 20 0 0 1\n"
        ),
    );
    let mut tokens = packed(&source, true);
    for row in 0..5 {
        tokens.masks[row * 32 + 1] = 1;
        tokens.masks[row * 32 + 4] = u8::from(row != 1);
        tokens.positions[row * 96 + 3..][..3].copy_from_slice(&source.atoms[row].position);
    }
    let mut names = vec![String::new(); 32];
    names[1] = "CA".into();
    names[4] = "CB".into();
    let schema = crate::observations::Schema {
        slots: 32,
        hydrogen_start: 15,
        residues: vec![crate::observations::schema::Residue {
            residue: "ALA".into(),
            names,
            parents: vec![-1; 32],
            states: Vec::new(),
        }],
    };
    tokens
        .complete(&schema)
        .expect("Complete-residue selection must preserve source identities");
    assert_eq!(tokens.sequence, [1, 3, 4, 69]);
    assert_eq!(tokens.keys[1][1], "3");
    tokens
        .complete(&schema)
        .expect("Repeated completeness checks must preserve index gaps");
    assert_eq!(tokens.sequence, [1, 3, 4, 69]);
}

/// The first residue of every pair below, as `(N, CA, C)`.
const RESIDUE_1: [[f32; 3]; 3] = [[-0.523, 1.363, 0.], [0., 0., 0.], [1.52, 0., 0.]];

// Second residues, as `(N, CA, C)`, with the C(1) to N(2) and CA-CA distances in A.
const BONDED: [[f32; 3]; 3] = [[2.103, 1.195, 0.], [3.556, 1.34, 0.], [4.443, 0.106, 0.]]; // C-N 1.33, CA-CA 3.80
const BONDED_STRETCHED: [[f32; 3]; 3] =
    [[2.357, 1.072, 0.], [3.795, 1.326, 0.], [4.773, 0.162, 0.]]; // C-N 1.36, CA-CA 4.02
const BROKEN_CLOSE: [[f32; 3]; 3] = [[2.835, 2.696, 0.], [3.528, 1.411, 0.], [2.73, 0.118, 0.]]; // C-N 3.00, CA-CA 3.80
const BROKEN_FAR: [[f32; 3]; 3] = [[3.054, 3.146, 0.], [4.48, 2.832, 0.], [4.939, 1.383, 0.]]; // C-N 3.50, CA-CA 5.30
const LINK_2_39: [[f32; 3]; 3] = [[2.568, 2.148, 0.], [3.702, 1.228, 0.], [3.462, -0.273, 0.]]; // C-N 2.39, CA-CA 3.90
const LINK_2_41: [[f32; 3]; 3] = [[2.576, 2.166, 0.], [3.7, 1.233, 0.], [3.444, -0.265, 0.]]; // C-N 2.41, CA-CA 3.90

/// Backbone rows of one residue; a `None` atom is left unobserved.
fn backbone(chain: &str, number: i64, atoms: [Option<[f32; 3]>; 3]) -> String {
    let mut rows = String::new();
    for (name, position) in ["N", "CA", "C"].iter().zip(atoms) {
        if let Some([x, y, z]) = position {
            let element = &name[..1];
            rows += &format!("1 {chain} {number} {number} ALA {name} {element} {x} {y} {z} 1\n");
        }
    }
    rows
}

/// Gaps of a two-residue chain, with either residue's atoms optionally removed.
fn pair_gap(second: [[f32; 3]; 3], drop_left: &[usize], drop_right: &[usize]) -> i64 {
    let keep = |atoms: [[f32; 3]; 3], drop: &[usize]| {
        std::array::from_fn(|slot| (!drop.contains(&slot)).then_some(atoms[slot]))
    };
    let rows =
        backbone("A", 1, keep(RESIDUE_1, drop_left)) + &backbone("A", 2, keep(second, drop_right));
    let sequence = packed(&input("polypeptide(L)", &rows), true).sequence;
    sequence[1] - sequence[0]
}

const NONE: &[usize] = &[];
const N: usize = 0;
const C: usize = 2;

#[test]
fn a_peptide_bond_keeps_residues_consecutive_whatever_their_ca_distance() {
    // 1SFV ILE13 to PRO14 has C-N 1.36 A and CA-CA 4.02 A: bonded, so no break.
    assert_eq!(pair_gap(BONDED, NONE, NONE), 1);
    assert_eq!(pair_gap(BONDED_STRETCHED, NONE, NONE), 1);
}

#[test]
fn a_missing_peptide_bond_breaks_whether_ca_atoms_are_near_or_far() {
    // CA-CA below the old 4.0 A limit is still a break when C and N are 3.0 A apart.
    assert_eq!(pair_gap(BROKEN_CLOSE, NONE, NONE), INDEX_SEPARATION);
    assert_eq!(pair_gap(BROKEN_FAR, NONE, NONE), INDEX_SEPARATION);
}

#[test]
fn the_peptide_link_limit_is_2_4_angstroms_between_c_and_n() {
    // Both pairs have the same CA-CA, 3.9 A, so only C-N tells them apart.
    assert_eq!(pair_gap(LINK_2_39, NONE, NONE), 1);
    assert_eq!(pair_gap(LINK_2_41, NONE, NONE), INDEX_SEPARATION);
}

#[test]
fn an_unobserved_c_or_n_falls_back_to_the_ca_distance() {
    for (left, right) in [(NONE, &[N][..]), (&[C][..], NONE), (&[C][..], &[N][..])] {
        assert_eq!(pair_gap(BONDED, left, right), 1);
        // Without the bond atoms nothing shows the 4.02 A pair is bonded.
        assert_eq!(pair_gap(BONDED_STRETCHED, left, right), INDEX_SEPARATION);
        assert_eq!(pair_gap(BROKEN_FAR, left, right), INDEX_SEPARATION);
    }
}

#[test]
fn only_the_linking_atoms_matter_when_other_backbone_atoms_are_missing() {
    // The left residue has no N and the right no C, but its C and N are both observed.
    assert_eq!(pair_gap(BONDED_STRETCHED, &[N], &[C]), 1);
    assert_eq!(pair_gap(BROKEN_CLOSE, &[N], &[C]), INDEX_SEPARATION);
}

#[test]
fn a_nucleotide_link_limit_is_still_2_4_angstroms_between_o3_and_p() {
    for (distance, expected) in [(2.39, [0, 1]), (2.41, [0, 65])] {
        let rows = format!("1 A 1 1 A \"O3'\" O 0 0 0 1\n1 A 2 2 A P P {distance} 0 0 1\n");
        assert_eq!(
            packed(&input("polyribonucleotide", &rows), false).sequence,
            expected
        );
    }
}

fn pair_topology(second: [[f32; 3]; 3]) -> [Vec<[i64; 3]>; 2] {
    let rows = backbone("A", 1, RESIDUE_1.map(Some)) + &backbone("A", 2, second.map(Some));
    let mut source = input("polypeptide(L)", &rows);
    source.bonds.push((2, 3));
    [
        crate::observations::topology::pdb(&source),
        crate::observations::topology::protein(&source, "observed_connectivity")
            .expect("Topology must resolve"),
    ]
}

#[test]
fn distant_backbone_atoms_do_not_establish_internal_peptide_topology() {
    // C and N are 3.0 A apart, so the declared bond between them is not accepted.
    for values in pair_topology(BROKEN_CLOSE) {
        assert_eq!(values[0][1], -1);
        assert_eq!(values[1][0], -1);
    }
}

#[test]
fn a_bonded_pair_with_a_stretched_ca_distance_is_internal_peptide_topology() {
    for values in pair_topology(BONDED_STRETCHED) {
        assert_eq!(values[0][1], 0);
        assert_eq!(values[1][0], 0);
    }
}

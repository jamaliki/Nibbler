//! One physical site stored twice as separate residues is resolved to one conformer.
//!
//! Alternate locations are normally letters inside a residue. Some entries store
//! the competing conformers as different residues, ligands or chains instead, each
//! carrying one letter. Those atoms are not a variable group, so selection used to
//! keep every one of them and the view held the site twice.

#![allow(clippy::expect_used)]

use _core::cif::{CifDocument, parse};
use _core::pdbx::{
    ConformerAtom, ConformerPolicy, select_atom_conformers, select_atom_conformers_with_policy,
    select_conformer_document, select_conformer_document_with_policy,
};
use std::collections::BTreeSet;

const TAGS: [&str; 20] = [
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_PDB_ins_code",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "occupancy",
    "B_iso_or_equiv",
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "pdbx_PDB_model_num",
];

/// 9KVL, neutron: MPD A501 (alternate A, 0.70) and MRD A502 (alternate B, 0.30)
/// are two enantiomers of one ligand bound at one site. Rows 0..22 and 22..44.
const MPD_AND_MRD_9KVL: &str = r##"
HETATM 5109 C C1 A MPD B 2 . ? -21.243 -14.923 17.830 0.70 89.84 501 MPD A C1 1
HETATM 5110 C C2 A MPD B 2 . ? -20.464 -14.534 16.580 0.70 68.70 501 MPD A C2 1
HETATM 5111 O O2 A MPD B 2 . ? -20.025 -15.755 15.932 0.70 55.82 501 MPD A O2 1
HETATM 5112 C CM A MPD B 2 . ? -21.380 -13.798 15.609 0.70 83.80 501 MPD A CM 1
HETATM 5113 C C3 A MPD B 2 . ? -19.280 -13.630 16.934 0.70 68.24 501 MPD A C3 1
HETATM 5114 C C4 A MPD B 2 . ? -18.005 -14.356 17.364 0.70 65.19 501 MPD A C4 1
HETATM 5115 O O4 A MPD B 2 . ? -17.557 -15.174 16.306 0.70 57.83 501 MPD A O4 1
HETATM 5116 C C5 A MPD B 2 . ? -16.905 -13.356 17.701 0.70 76.76 501 MPD A C5 1
HETATM 5117 D D11 A MPD B 2 . ? -20.612 -14.793 18.708 0.70 107.81 501 MPD A D11 1
HETATM 5118 D D12 A MPD B 2 . ? -22.125 -14.290 17.928 0.70 107.81 501 MPD A D12 1
HETATM 5119 D D13 A MPD B 2 . ? -21.555 -15.964 17.751 0.70 107.81 501 MPD A D13 1
HETATM 5120 D DO2 A MPD B 2 . ? -19.084 -15.690 15.674 0.70 66.98 501 MPD A DO2 1
HETATM 5121 D DM1 A MPD B 2 . ? -21.208 -12.725 15.687 0.70 100.56 501 MPD A DM1 1
HETATM 5122 D DM2 A MPD B 2 . ? -21.158 -14.120 14.593 0.70 100.56 501 MPD A DM2 1
HETATM 5123 D DM3 A MPD B 2 . ? -22.422 -14.020 15.837 0.70 100.56 501 MPD A DM3 1
HETATM 5124 D D31 A MPD B 2 . ? -19.587 -12.956 17.733 0.70 81.89 501 MPD A D31 1
HETATM 5125 D D32 A MPD B 2 . ? -19.045 -13.015 16.066 0.70 81.89 501 MPD A D32 1
HETATM 5126 D D4 A MPD B 2 . ? -18.209 -14.959 18.251 0.70 78.23 501 MPD A D4 1
HETATM 5127 D DO4 A MPD B 2 . ? -16.908 -14.684 15.759 0.70 69.40 501 MPD A DO4 1
HETATM 5128 D D51 A MPD B 2 . ? -16.539 -12.895 16.785 0.70 92.11 501 MPD A D51 1
HETATM 5129 D D52 A MPD B 2 . ? -16.081 -13.880 18.183 0.70 92.11 501 MPD A D52 1
HETATM 5130 D D53 A MPD B 2 . ? -17.288 -12.591 18.373 0.70 92.11 501 MPD A D53 1
HETATM 5131 C C1 B MRD C 3 . ? -16.851 -13.163 17.582 0.30 89.88 502 MRD A C1 1
HETATM 5132 C C2 B MRD C 3 . ? -17.865 -14.287 17.408 0.30 69.32 502 MRD A C2 1
HETATM 5133 O O2 B MRD C 3 . ? -17.478 -14.964 16.194 0.30 57.22 502 MRD A O2 1
HETATM 5134 C CM B MRD C 3 . ? -17.765 -15.277 18.560 0.30 84.23 502 MRD A CM 1
HETATM 5135 C C3 B MRD C 3 . ? -19.273 -13.703 17.270 0.30 68.48 502 MRD A C3 1
HETATM 5136 C C4 B MRD C 3 . ? -20.321 -14.707 16.780 0.30 64.96 502 MRD A C4 1
HETATM 5137 O O4 B MRD C 3 . ? -19.731 -15.724 16.001 0.30 57.64 502 MRD A O4 1
HETATM 5138 C C5 B MRD C 3 . ? -21.430 -14.029 15.980 0.30 77.35 502 MRD A C5 1
HETATM 5139 D D1C1 B MRD C 3 . ? -17.173 -12.498 18.381 0.30 107.85 502 MRD A D1C1 1
HETATM 5140 D D1C2 B MRD C 3 . ? -15.875 -13.588 17.817 0.30 107.85 502 MRD A D1C2 1
HETATM 5141 D D1C3 B MRD C 3 . ? -16.780 -12.593 16.656 0.30 107.85 502 MRD A D1C3 1
HETATM 5142 D D2 B MRD C 3 . ? -17.166 -14.281 15.570 0.30 68.67 502 MRD A D2 1
HETATM 5143 D DMC1 B MRD C 3 . ? -18.748 -15.695 18.773 0.30 101.08 502 MRD A DMC1 1
HETATM 5144 D DMC2 B MRD C 3 . ? -17.075 -16.075 18.289 0.30 101.08 502 MRD A DMC2 1
HETATM 5145 D DMC3 B MRD C 3 . ? -17.394 -14.759 19.443 0.30 101.08 502 MRD A DMC3 1
HETATM 5146 D D3C1 B MRD C 3 . ? -19.236 -12.864 16.576 0.30 82.17 502 MRD A D3C1 1
HETATM 5147 D D3C2 B MRD C 3 . ? -19.589 -13.311 18.235 0.30 82.17 502 MRD A D3C2 1
HETATM 5148 D D4 B MRD C 3 . ? -20.782 -15.160 17.658 0.30 77.95 502 MRD A D4 1
HETATM 5149 D DA B MRD C 3 . ? -19.121 -15.323 15.349 0.30 69.17 502 MRD A DA 1
HETATM 5150 D D5C1 B MRD C 3 . ? -21.900 -14.758 15.321 0.30 92.81 502 MRD A D5C1 1
HETATM 5151 D D5C2 B MRD C 3 . ? -21.015 -13.223 15.376 0.30 92.81 502 MRD A D5C2 1
HETATM 5152 D D5C3 B MRD C 3 . ? -22.175 -13.622 16.664 0.30 92.81 502 MRD A D5C3 1
"##;

/// 6S48: DG C1 (alternate A, 0.50) and DG H1 (alternate B, 0.50) are one 5'
/// nucleotide stored as two chains, 0.08 A apart. Rows 0..19 and 19..38.
const DG_C1_AND_H1_6S48: &str = r##"
ATOM 3992 O "O5'" A DG C 2 1 ? 11.235 36.033 41.774 0.50 60.60 1 DG C "O5'" 1
ATOM 3993 C "C5'" A DG C 2 1 ? 10.190 35.467 40.960 0.50 59.92 1 DG C "C5'" 1
ATOM 3994 C "C4'" A DG C 2 1 ? 9.660 36.503 39.996 0.50 58.47 1 DG C "C4'" 1
ATOM 3995 O "O4'" A DG C 2 1 ? 8.478 37.112 40.565 0.50 57.65 1 DG C "O4'" 1
ATOM 3996 C "C3'" A DG C 2 1 ? 10.627 37.647 39.687 0.50 58.28 1 DG C "C3'" 1
ATOM 3997 O "O3'" A DG C 2 1 ? 11.243 37.485 38.404 0.50 58.32 1 DG C "O3'" 1
ATOM 3998 C "C2'" A DG C 2 1 ? 9.744 38.880 39.653 0.50 56.72 1 DG C "C2'" 1
ATOM 3999 C "C1'" A DG C 2 1 ? 8.590 38.526 40.573 0.50 54.58 1 DG C "C1'" 1
ATOM 4000 N N9 A DG C 2 1 ? 8.756 38.957 41.957 0.50 50.15 1 DG C N9 1
ATOM 4001 C C8 A DG C 2 1 ? 9.151 38.195 43.031 0.50 48.02 1 DG C C8 1
ATOM 4002 N N7 A DG C 2 1 ? 9.192 38.869 44.148 0.50 46.79 1 DG C N7 1
ATOM 4003 C C5 A DG C 2 1 ? 8.793 40.150 43.792 0.50 45.39 1 DG C C5 1
ATOM 4004 C C6 A DG C 2 1 ? 8.639 41.318 44.584 0.50 43.93 1 DG C C6 1
ATOM 4005 O O6 A DG C 2 1 ? 8.830 41.457 45.799 0.50 42.38 1 DG C O6 1
ATOM 4006 N N1 A DG C 2 1 ? 8.220 42.402 43.818 0.50 43.93 1 DG C N1 1
ATOM 4007 C C2 A DG C 2 1 ? 7.980 42.368 42.467 0.50 44.61 1 DG C C2 1
ATOM 4008 N N2 A DG C 2 1 ? 7.581 43.520 41.909 0.50 44.01 1 DG C N2 1
ATOM 4009 N N3 A DG C 2 1 ? 8.116 41.285 41.719 0.50 46.61 1 DG C N3 1
ATOM 4010 C C4 A DG C 2 1 ? 8.521 40.220 42.443 0.50 46.67 1 DG C C4 1
ATOM 4661 O "O5'" B DG H 3 1 ? 11.138 35.940 41.829 0.50 64.51 1 DG H "O5'" 1
ATOM 4662 C "C5'" B DG H 3 1 ? 10.136 35.415 40.938 0.50 63.51 1 DG H "C5'" 1
ATOM 4663 C "C4'" B DG H 3 1 ? 9.662 36.493 39.990 0.50 62.00 1 DG H "C4'" 1
ATOM 4664 O "O4'" B DG H 3 1 ? 8.495 37.133 40.555 0.50 60.75 1 DG H "O4'" 1
ATOM 4665 C "C3'" B DG H 3 1 ? 10.677 37.605 39.723 0.50 61.89 1 DG H "C3'" 1
ATOM 4666 O "O3'" B DG H 3 1 ? 11.314 37.439 38.451 0.50 62.60 1 DG H "O3'" 1
ATOM 4667 C "C2'" B DG H 3 1 ? 9.839 38.870 39.690 0.50 59.99 1 DG H "C2'" 1
ATOM 4668 C "C1'" B DG H 3 1 ? 8.647 38.544 40.571 0.50 57.33 1 DG H "C1'" 1
ATOM 4669 N N9 B DG H 3 1 ? 8.783 38.968 41.960 0.50 52.25 1 DG H N9 1
ATOM 4670 C C8 B DG H 3 1 ? 9.165 38.204 43.037 0.50 49.89 1 DG H C8 1
ATOM 4671 N N7 B DG H 3 1 ? 9.184 38.873 44.158 0.50 48.54 1 DG H N7 1
ATOM 4672 C C5 B DG H 3 1 ? 8.786 40.154 43.800 0.50 47.04 1 DG H C5 1
ATOM 4673 C C6 B DG H 3 1 ? 8.614 41.318 44.592 0.50 45.43 1 DG H C6 1
ATOM 4674 O O6 B DG H 3 1 ? 8.784 41.454 45.811 0.50 43.71 1 DG H O6 1
ATOM 4675 N N1 B DG H 3 1 ? 8.204 42.403 43.824 0.50 45.40 1 DG H N1 1
ATOM 4676 C C2 B DG H 3 1 ? 7.986 42.373 42.469 0.50 46.13 1 DG H C2 1
ATOM 4677 N N2 B DG H 3 1 ? 7.592 43.526 41.908 0.50 45.48 1 DG H N2 1
ATOM 4678 N N3 B DG H 3 1 ? 8.140 41.294 41.719 0.50 48.22 1 DG H N3 1
ATOM 4679 C C4 B DG H 3 1 ? 8.536 40.228 42.446 0.50 48.44 1 DG H C4 1
"##;

fn document(rows: &str) -> CifDocument {
    let tags: String = TAGS
        .iter()
        .map(|tag| format!("_atom_site.{tag}\n"))
        .collect();
    let text = format!("data_site\nloop_\n{tags}{}\n", rows.trim());
    parse(text.as_bytes()).expect("source")
}

fn rows_of(selection: &[usize]) -> BTreeSet<usize> {
    selection.iter().copied().collect()
}

#[test]
fn ligand_conformers_stored_as_two_residues_keep_the_higher_occupancy_one() {
    let doc = document(MPD_AND_MRD_9KVL);
    for seed in 0..8 {
        let selection = select_conformer_document(&doc, None, seed).expect("selection");
        assert_eq!(
            selection.atom_rows,
            (0..22).collect::<Vec<_>>(),
            "seed {seed}"
        );
    }
}

#[test]
fn random_policy_picks_one_whole_ligand_and_replays() {
    let doc = document(MPD_AND_MRD_9KVL);
    let mut choices = BTreeSet::new();
    for seed in 0..32 {
        let choose = || {
            select_conformer_document_with_policy(
                &doc,
                None,
                seed,
                ConformerPolicy::RandomClashAware,
            )
            .expect("selection")
            .atom_rows
        };
        let rows = choose();
        assert_eq!(rows, choose());
        assert!(
            rows == (0..22).collect::<Vec<_>>() || rows == (22..44).collect::<Vec<_>>(),
            "seed {seed} kept {rows:?}"
        );
        choices.insert(rows[0]);
    }
    assert_eq!(choices.len(), 2, "the 0.30 ligand is a valid random choice");
}

#[test]
fn nucleotide_conformers_stored_as_two_chains_keep_one_chain() {
    let doc = document(DG_C1_AND_H1_6S48);
    let mut choices = BTreeSet::new();
    for seed in 0..16 {
        let selection = select_conformer_document(&doc, None, seed).expect("selection");
        assert!(
            selection.atom_rows == (0..19).collect::<Vec<_>>()
                || selection.atom_rows == (19..38).collect::<Vec<_>>(),
            "seed {seed} kept {:?}",
            selection.atom_rows
        );
        choices.insert(selection.atom_rows[0]);
    }
    assert_eq!(choices.len(), 2, "equal occupancy is broken by the seed");
}

#[test]
fn a_lone_alternate_residue_is_kept() {
    let rows: Vec<_> = DG_C1_AND_H1_6S48.trim().lines().take(19).collect();
    let doc = document(&rows.join("\n"));
    let selection = select_conformer_document(&doc, None, 0).expect("selection");
    assert_eq!(selection.atom_rows, (0..19).collect::<Vec<_>>());
}

fn atom(residue: usize, chain: &str, alt: Option<&str>, x: f64, occupancy: f64) -> ConformerAtom {
    ConformerAtom {
        residue,
        chain: chain.into(),
        sequence: None,
        atom_name: "C1".into(),
        alternate: alt.map(str::to_owned),
        coordinates: [x, 0., 0.],
        occupancy,
        clash_eligible: true,
    }
}

#[test]
fn single_label_residues_that_clash_are_alternatives_of_one_site() {
    let atoms = [
        atom(0, "A", Some("A"), 0., 0.6),
        atom(0, "A", Some("A"), 5., 0.6),
        atom(1, "B", Some("B"), 0.5, 0.4),
        atom(1, "B", Some("B"), 5.5, 0.4),
    ];
    let selection = select_atom_conformers(&atoms, &[], 0).expect("selection");
    assert_eq!(selection.atom_rows, [0, 1]);
    let flipped = [
        atom(0, "A", Some("A"), 0., 0.4),
        atom(0, "A", Some("A"), 5., 0.4),
        atom(1, "B", Some("B"), 0.5, 0.6),
        atom(1, "B", Some("B"), 5.5, 0.6),
    ];
    let selection = select_atom_conformers(&flipped, &[], 0).expect("selection");
    assert_eq!(selection.atom_rows, [2, 3]);
}

#[test]
fn a_clash_chain_across_residues_resolves_to_one_label() {
    // A(0) touches B(1), B(1) touches A(2); the label with more occupancy wins.
    let atoms = [
        atom(0, "A", Some("A"), 0., 0.5),
        atom(1, "B", Some("B"), 1., 0.6),
        atom(2, "C", Some("A"), 2., 0.5),
    ];
    let selection = select_atom_conformers(&atoms, &[], 0).expect("selection");
    assert_eq!(selection.atom_rows, [0, 2]);
}

#[test]
fn different_letters_that_do_not_touch_are_all_kept() {
    let atoms = [
        atom(0, "A", Some("A"), 0., 0.5),
        atom(1, "B", Some("B"), 10., 0.5),
    ];
    let selection = select_atom_conformers(&atoms, &[], 0).expect("selection");
    assert_eq!(selection.atom_rows, [0, 1]);
}

#[test]
fn the_same_letter_on_neighbouring_residues_is_not_a_conflict() {
    let atoms = [
        atom(0, "A", Some("A"), 0., 0.5),
        atom(1, "B", Some("A"), 1., 0.5),
    ];
    let selection = select_atom_conformers(&atoms, &[], 0).expect("selection");
    assert_eq!(selection.atom_rows, [0, 1]);
}

#[test]
fn bonded_atoms_with_different_letters_are_not_alternatives() {
    let atoms = [
        atom(0, "A", Some("A"), 0., 0.5),
        atom(1, "A", Some("B"), 1.4, 0.5),
    ];
    let bonded = select_atom_conformers(&atoms, &[[0, 1]], 0).expect("selection");
    assert_eq!(bonded.atom_rows, [0, 1]);
    let mut peptide = atoms.clone();
    peptide[0].atom_name = "C".into();
    peptide[0].sequence = Some(4);
    peptide[1].atom_name = "N".into();
    peptide[1].sequence = Some(5);
    let backbone = select_atom_conformers(&peptide, &[], 0).expect("selection");
    assert_eq!(backbone.atom_rows, [0, 1]);
}

#[test]
fn hydrogens_and_zero_occupancy_atoms_do_not_make_alternatives() {
    let mut atoms = [
        atom(0, "A", Some("A"), 0., 0.5),
        atom(1, "B", Some("B"), 0.5, 0.5),
        atom(2, "C", Some("C"), 0.5, 0.),
    ];
    atoms[1].clash_eligible = false;
    let selection = select_atom_conformers(&atoms, &[], 0).expect("selection");
    assert_eq!(selection.atom_rows, [0, 1, 2]);
}

#[test]
fn random_policy_resolves_cross_residue_alternatives_and_replays() {
    let atoms = [
        atom(0, "A", Some("A"), 0., 0.9),
        atom(1, "B", Some("B"), 0.5, 0.1),
        atom(2, "C", Some("C"), 0.5, 0.),
    ];
    let mut choices = BTreeSet::new();
    for seed in 0..32 {
        let choose = || {
            select_atom_conformers_with_policy(&atoms, &[], seed, ConformerPolicy::RandomClashAware)
                .expect("selection")
                .atom_rows
        };
        let rows = choose();
        assert_eq!(rows, choose());
        assert_eq!(rows.len(), 2, "seed {seed}");
        assert!(rows.contains(&2), "an unoccupied atom never competes");
        choices.insert(rows_of(&rows));
    }
    assert_eq!(choices.len(), 2);
}

#[test]
fn single_label_alternatives_still_block_a_clashing_variable_residue_variant() {
    // Residue 0 is a variable group; residue 1 is a lone shared atom with the
    // letter B. The existing contact rule is unchanged by the cluster pass.
    let atoms = [
        atom(0, "A", Some("A"), 0., 0.8),
        atom(0, "A", Some("B"), 10., 0.2),
        atom(1, "B", Some("B"), 0.5, 0.4),
    ];
    let selection = select_atom_conformers(&atoms, &[], 0).expect("selection");
    assert_eq!(selection.atom_rows, [1, 2]);
}

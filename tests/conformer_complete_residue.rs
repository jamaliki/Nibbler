//! A residue never takes an incomplete conformer while another one is complete.
//!
//! The search weighs conformers by occupancy, or at random under the random policy.
//! Neither looks at which atoms a conformer holds, so a seeded choice could keep a
//! conformer that lacks atoms another one has, or one that holds only the exchangeable
//! deuterium of a neutron structure. The canonical residue then arrives with heavy atoms
//! missing, and the packer drops it whole.

#![allow(clippy::expect_used)]

use _core::cif::{CifDocument, parse};
use _core::pdbx::{
    ConformerAtom, ConformerPolicy, select_atom_conformers, select_atom_conformers_with_policy,
    select_conformer_document_with_policy,
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

/// 1CJX GLU 313: conformer A (0.50) lacks OE1, which only conformer B (0.50) holds.
const GLU_1CJX: &str = r##"
ATOM 2453 N N . GLU A 1 313 ? -29.749 29.414 63.760 1.00 21.06 313 GLU A N 1
ATOM 2454 C CA . GLU A 1 313 ? -30.397 29.516 62.461 1.00 22.00 313 GLU A CA 1
ATOM 2455 C C . GLU A 1 313 ? -29.284 29.851 61.457 1.00 20.97 313 GLU A C 1
ATOM 2456 O O . GLU A 1 313 ? -28.187 30.179 61.910 1.00 20.97 313 GLU A O 1
ATOM 2457 C CB . GLU A 1 313 ? -31.411 30.646 62.401 1.00 23.03 313 GLU A CB 1
ATOM 2458 C CG A GLU A 1 313 ? -32.844 30.408 62.806 0.50 30.62 313 GLU A CG 1
ATOM 2459 C CG B GLU A 1 313 ? -32.568 30.458 63.353 0.50 27.11 313 GLU A CG 1
ATOM 2460 C CD A GLU A 1 313 ? -33.130 28.975 60.927 0.50 43.84 313 GLU A CD 1
ATOM 2461 C CD B GLU A 1 313 ? -33.680 31.464 63.167 0.50 29.56 313 GLU A CD 1
ATOM 2462 O OE1 B GLU A 1 313 ? -33.512 32.462 62.439 0.50 31.03 313 GLU A OE1 1
ATOM 2463 O OE2 A GLU A 1 313 ? -34.603 28.813 62.535 0.50 39.92 313 GLU A OE2 1
ATOM 2464 O OE2 B GLU A 1 313 ? -34.727 31.200 63.802 0.50 34.53 313 GLU A OE2 1
"##;

/// 5MOO LEU 67, neutron: heavy atoms are labelled A (0.45) and C (0.55); the exchangeable
/// amide deuterium is labelled B (0.88), a conformer with a single deuterium and no heavy
/// atom of its own.
const LEU_5MOO: &str = r##"
ATOM 803 N N . LEU A 1 50 ? -1.665 -25.250 2.895 1.00 18.73 67 LEU A N 1
ATOM 804 C CA A LEU A 1 50 ? -1.541 -26.547 2.249 0.45 19.61 67 LEU A CA 1
ATOM 805 C CA C LEU A 1 50 ? -1.540 -26.551 2.255 0.55 18.62 67 LEU A CA 1
ATOM 806 C C . LEU A 1 50 ? -2.837 -27.317 2.456 1.00 18.89 67 LEU A C 1
ATOM 807 O O . LEU A 1 50 ? -3.604 -27.045 3.379 1.00 21.24 67 LEU A O 1
ATOM 808 C CB A LEU A 1 50 ? -0.372 -27.350 2.830 0.45 22.85 67 LEU A CB 1
ATOM 809 C CB C LEU A 1 50 ? -0.399 -27.361 2.884 0.55 19.10 67 LEU A CB 1
ATOM 810 C CG A LEU A 1 50 ? 0.997 -26.662 2.817 0.45 27.11 67 LEU A CG 1
ATOM 811 C CG C LEU A 1 50 ? 0.956 -26.667 2.976 0.55 21.07 67 LEU A CG 1
ATOM 812 C CD1 A LEU A 1 50 ? 1.337 -26.058 4.181 0.45 30.09 67 LEU A CD1 1
ATOM 813 C CD1 C LEU A 1 50 ? 1.968 -27.524 3.738 0.55 25.00 67 LEU A CD1 1
ATOM 814 C CD2 A LEU A 1 50 ? 2.074 -27.630 2.383 0.45 33.41 67 LEU A CD2 1
ATOM 815 C CD2 C LEU A 1 50 ? 1.454 -26.326 1.602 0.55 31.64 67 LEU A CD2 1
ATOM 816 H H A LEU A 1 50 ? -1.659 -25.251 3.921 0.12 18.86 67 LEU A H 1
ATOM 817 D D B LEU A 1 50 ? -1.659 -25.251 3.921 0.88 18.86 67 LEU A D 1
ATOM 818 H HA . LEU A 1 50 ? -1.344 -26.421 1.184 1.00 19.12 67 LEU A HA 1
ATOM 819 H HB2 A LEU A 1 50 ? -0.609 -27.598 3.860 0.45 23.87 67 LEU A HB2 1
ATOM 820 H HB2 C LEU A 1 50 ? -0.696 -27.623 3.895 0.55 20.77 67 LEU A HB2 1
ATOM 821 H HB3 A LEU A 1 50 ? -0.271 -28.278 2.261 0.45 23.87 67 LEU A HB3 1
ATOM 822 H HB3 C LEU A 1 50 ? -0.254 -28.285 2.319 0.55 20.77 67 LEU A HB3 1
ATOM 823 H HG A LEU A 1 50 ? 0.968 -25.846 2.096 0.45 28.18 67 LEU A HG 1
ATOM 824 H HD11 A LEU A 1 50 ? 1.254 -26.840 4.936 0.45 32.44 67 LEU A HD11 1
ATOM 825 H HD11 C LEU A 1 50 ? 2.327 -28.327 3.095 0.55 27.41 67 LEU A HD11 1
ATOM 826 H HD12 A LEU A 1 50 ? 2.356 -25.671 4.156 0.45 32.44 67 LEU A HD12 1
ATOM 827 H HD13 A LEU A 1 50 ? 0.628 -25.261 4.400 0.45 32.44 67 LEU A HD13 1
ATOM 828 H HD13 C LEU A 1 50 ? 1.492 -27.932 4.628 0.55 27.41 67 LEU A HD13 1
ATOM 829 H HD21 A LEU A 1 50 ? 1.974 -28.550 2.958 0.45 36.09 67 LEU A HD21 1
ATOM 830 H HD21 C LEU A 1 50 ? 1.521 -27.244 1.015 0.55 35.41 67 LEU A HD21 1
ATOM 831 H HD22 A LEU A 1 50 ? 1.949 -27.834 1.320 0.45 36.09 67 LEU A HD22 1
ATOM 832 H HD22 C LEU A 1 50 ? 0.755 -25.633 1.135 0.55 35.41 67 LEU A HD22 1
ATOM 833 H HD23 A LEU A 1 50 ? 3.052 -27.180 2.565 0.45 36.09 67 LEU A HD23 1
ATOM 834 H HD23 C LEU A 1 50 ? 2.439 -25.870 1.686 0.55 35.41 67 LEU A HD23 1
"##;

fn document(rows: &str) -> CifDocument {
    let tags: String = TAGS
        .iter()
        .map(|tag| format!("_atom_site.{tag}\n"))
        .collect();
    let text = format!("data_site\nloop_\n{tags}{}\n", rows.trim());
    parse(text.as_bytes()).expect("source")
}

fn select(doc: &CifDocument, seed: u64, policy: ConformerPolicy) -> Vec<usize> {
    select_conformer_document_with_policy(doc, None, seed, policy)
        .expect("selection")
        .atom_rows
}

#[test]
fn a_residue_keeps_the_conformer_that_holds_every_heavy_atom() {
    let doc = document(GLU_1CJX);
    // Rows 0-4 are shared; B holds rows 6, 8, 9 and 11, A rows 5, 7 and 10.
    let complete = vec![0, 1, 2, 3, 4, 6, 8, 9, 11];
    for seed in 0..32 {
        assert_eq!(
            select(&doc, seed, ConformerPolicy::RandomClashAware),
            complete,
            "random, seed {seed}"
        );
        assert_eq!(
            select(&doc, seed, ConformerPolicy::Occupancy),
            complete,
            "occupancy, seed {seed}"
        );
    }
}

#[test]
fn a_conformer_without_heavy_atoms_is_never_the_residue() {
    let doc = document(LEU_5MOO);
    let mut letters = BTreeSet::new();
    for seed in 0..32 {
        for policy in [
            ConformerPolicy::RandomClashAware,
            ConformerPolicy::Occupancy,
        ] {
            let rows = select(&doc, seed, policy);
            // Rows 1 and 2 are the alpha carbon of A and C; row 14 is the lone deuterium of B.
            assert!(
                !rows.contains(&14),
                "seed {seed} kept the deuterium conformer"
            );
            assert!(rows.contains(&1) != rows.contains(&2));
            letters.insert(rows.contains(&1));
        }
    }
    assert_eq!(letters.len(), 2, "both complete conformers stay possible");
}

fn named(alt: Option<&str>, name: &str, x: f64, occupancy: f64) -> ConformerAtom {
    ConformerAtom {
        residue: 0,
        chain: "A".into(),
        sequence: Some(1),
        atom_name: name.into(),
        alternate: alt.map(str::to_owned),
        coordinates: [x, 0., 0.],
        occupancy,
        clash_eligible: true,
    }
}

#[test]
fn completeness_outranks_occupancy_only_between_unequal_conformers() {
    let mut residue = vec![
        named(Some("A"), "N", 0., 0.7),
        named(Some("A"), "CA", 3., 0.7),
        named(Some("B"), "N", 0., 0.4),
        named(Some("B"), "CA", 3., 0.4),
        named(Some("B"), "CB", 6., 0.4),
    ];
    let b_only = select_atom_conformers(&residue, &[], 0).expect("selection");
    assert_eq!(b_only.atom_rows, [2, 3, 4]);
    // When both conformers are complete, occupancy decides as before.
    residue.insert(2, named(Some("A"), "CB", 6., 0.7));
    let both = select_atom_conformers(&residue, &[], 0).expect("selection");
    assert_eq!(both.atom_rows, [0, 1, 2]);
    // When neither is complete, occupancy decides as before.
    let neither = [
        named(Some("A"), "N", 0., 0.6),
        named(Some("A"), "CA", 3., 0.6),
        named(Some("B"), "N", 0., 0.4),
        named(Some("B"), "CB", 6., 0.4),
    ];
    let partial = select_atom_conformers(&neither, &[], 0).expect("selection");
    assert_eq!(partial.atom_rows, [0, 1]);
}

#[test]
fn hydrogens_do_not_make_a_conformer_complete_and_a_clash_still_blocks_it() {
    let mut atoms = vec![
        named(Some("A"), "N", 0., 0.4),
        named(Some("A"), "CA", 3., 0.4),
        named(Some("B"), "N", 0., 0.6),
        named(Some("B"), "CA", 3., 0.6),
        named(Some("B"), "HA", 3.5, 0.6),
        named(Some("A"), "CB", 6., 0.4),
    ];
    atoms[4].clash_eligible = false;
    // B has no CB: the hydrogen does not stand in for it, so A is the complete one.
    let selection = select_atom_conformers(&atoms, &[], 0).expect("selection");
    assert_eq!(selection.atom_rows, [0, 1, 5]);
    // A shared atom of another residue sits on A's CB, which blocks the complete conformer.
    let mut other = named(None, "O", 6.2, 1.);
    other.residue = 1;
    other.chain = "B".into();
    atoms.push(other);
    let blocked = select_atom_conformers(&atoms, &[], 0).expect("selection");
    assert_eq!(blocked.atom_rows, [2, 3, 4, 6]);
    for seed in 0..8 {
        let random = select_atom_conformers_with_policy(
            &atoms,
            &[],
            seed,
            ConformerPolicy::RandomClashAware,
        )
        .expect("selection");
        assert_eq!(random.atom_rows, [2, 3, 4, 6], "seed {seed}");
    }
}

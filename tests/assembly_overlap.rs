//! Whole-ligand overlap checks keep source atom mappings and distinct molecules.

#![allow(clippy::expect_used)]

use _core::pdbx::{deduplicate_assembly_copies, deduplicate_assembly_copies_with_occupancy};

#[test]
fn permuted_symmetric_atoms_remove_one_complete_copy() {
    let positions = [[-2., 0., 0.], [2., 0., 0.], [2., 0., 0.], [-2., 0., 0.]];
    let result = deduplicate_assembly_copies(
        &positions,
        &[0, 1, 0, 1],
        &[0; 4],
        &[0, 0, 1, 1],
        &[true; 4],
        &[true; 4],
    )
    .expect("deduplication");
    assert_eq!(result.retained_rows, [0, 1]);
    assert_eq!(result.duplicate_pairs, [[2, 0], [3, 1]]);
}

#[test]
fn separate_source_units_and_protected_copies_stay_separate() {
    let positions = [[0., 0., 0.], [0., 0., 0.]];
    let result = deduplicate_assembly_copies(
        &positions,
        &[0, 1],
        &[0, 1],
        &[0, 1],
        &[true; 2],
        &[true; 2],
    )
    .expect("distinct units");
    assert_eq!(result.retained_rows, [0, 1]);
    let protected = deduplicate_assembly_copies(
        &positions,
        &[0, 0],
        &[0, 0],
        &[0, 1],
        &[false; 2],
        &[true; 2],
    )
    .expect("protected units");
    assert_eq!(protected.retained_rows, [0, 1]);
}

#[test]
fn distant_copies_and_zero_occupancy_sites_do_not_trigger_removal() {
    let positions = [[0., 0., 0.], [5., 0., 0.]];
    let result = deduplicate_assembly_copies(
        &positions,
        &[0, 0],
        &[0; 2],
        &[0, 1],
        &[true; 2],
        &[true; 2],
    )
    .expect("distant");
    assert!(result.duplicate_pairs.is_empty());
    let unresolved = deduplicate_assembly_copies(
        &[[0.; 3]; 2],
        &[0, 0],
        &[0; 2],
        &[0, 1],
        &[true; 2],
        &[false; 2],
    )
    .expect("unresolved");
    assert!(unresolved.duplicate_pairs.is_empty());
}

/// Copy 0 holds source atoms 0..3 at x = 0, 5 and 10. Copy 1 holds the same source
/// atoms with the first at `shift` and the others far away, so only the first pair can
/// touch: one of three atoms, below the half the whole-unit rule needs.
fn two_copies(shift: f64, occupancy: [f64; 6]) -> (Vec<[f64; 3]>, Vec<f64>) {
    let positions = vec![
        [0., 0., 0.],
        [5., 0., 0.],
        [10., 0., 0.],
        [shift, 0., 0.],
        [25., 0., 0.],
        [35., 0., 0.],
    ];
    (positions, occupancy.to_vec())
}

fn dedupe_with_occupancy(positions: &[[f64; 3]], occupancy: &[f64]) -> Vec<usize> {
    deduplicate_assembly_copies_with_occupancy(
        positions,
        &[0, 1, 2, 0, 1, 2],
        &[0; 6],
        &[0, 0, 0, 1, 1, 1],
        &[true; 6],
        occupancy,
    )
    .expect("deduplication")
    .retained_rows
}

#[test]
fn a_wholly_partial_copy_that_touches_an_earlier_copy_is_an_alternative() {
    let (positions, occupancy) = two_copies(0.5, [0.5; 6]);
    assert_eq!(dedupe_with_occupancy(&positions, &occupancy), [0, 1, 2]);
    let whole = deduplicate_assembly_copies(
        &positions,
        &[0, 1, 2, 0, 1, 2],
        &[0; 6],
        &[0, 0, 0, 1, 1, 1],
        &[true; 6],
        &[true; 6],
    )
    .expect("whole-unit rule");
    assert_eq!(whole.retained_rows.len(), 6);
}

#[test]
fn full_occupancy_or_distance_keeps_both_copies() {
    let (positions, full) = two_copies(0.5, [1.; 6]);
    assert_eq!(dedupe_with_occupancy(&positions, &full).len(), 6);
    let mixed = [0.5, 0.5, 1., 0.5, 0.5, 1.];
    assert_eq!(dedupe_with_occupancy(&positions, &mixed).len(), 6);
    let (apart, partial) = two_copies(2., [0.5; 6]);
    assert_eq!(dedupe_with_occupancy(&apart, &partial).len(), 6);
}

#[test]
fn unoccupied_atoms_do_not_make_an_alternative() {
    let (positions, _) = two_copies(0.5, [0.; 6]);
    let occupancy = [0.5, 0.5, 0.5, 0., 0.5, 0.5];
    assert_eq!(dedupe_with_occupancy(&positions, &occupancy).len(), 6);
}

#[test]
fn the_alternative_rule_maps_removed_atoms_to_the_earlier_copy() {
    let (positions, occupancy) = two_copies(0.5, [0.5; 6]);
    let result = deduplicate_assembly_copies_with_occupancy(
        &positions,
        &[0, 1, 2, 0, 1, 2],
        &[0; 6],
        &[0, 0, 0, 1, 1, 1],
        &[true; 6],
        &occupancy,
    )
    .expect("deduplication");
    assert_eq!(result.duplicate_pairs, [[3, 0], [4, 1], [5, 2]]);
}

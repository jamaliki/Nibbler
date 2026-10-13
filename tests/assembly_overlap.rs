//! Whole-ligand overlap checks keep source atom mappings and distinct molecules.

#![allow(clippy::expect_used)]

use _core::pdbx::deduplicate_assembly_copies;

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

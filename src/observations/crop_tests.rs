#![allow(clippy::unwrap_used)]

use super::{permutation, plan};

fn fixture(count: usize, ligand: &[u8], chains: &[i64]) -> super::CropPlans {
    let positions: Vec<u8> = (0..count)
        .flat_map(|index| [index as f32, 0., 0.])
        .flat_map(|value| value.to_ne_bytes())
        .collect();
    let chains: Vec<u8> = chains
        .iter()
        .flat_map(|value| value.to_ne_bytes())
        .collect();
    plan(&positions, &chains, ligand, &vec![1; count], 50, None).unwrap()
}

#[test]
fn stable_ties_and_short_chain_filters() {
    let result = fixture(60, &vec![0; 60], &vec![0; 60]);
    assert_eq!(result.centers.len(), 60);
    assert_eq!(result.ties[25], 0);
    assert_eq!(result.radii[25], 625.);
    assert_eq!(result.bits, vec![1; 60]);
    assert!((result.mass - 1.).abs() < 1e-14);
    let result = fixture(
        60,
        &vec![0; 60],
        &(0..60).map(|index| (index / 30) as i64).collect::<Vec<_>>(),
    );
    assert!(result.centers.is_empty());
}

#[test]
fn empty_and_protein_free_views_have_no_mass() {
    assert_eq!(plan(&[], &[], &[], &[], 50, None).unwrap().mass, 0.);
    assert_eq!(fixture(60, &vec![1; 60], &vec![0; 60]).mass, 0.);
}

#[test]
fn rejects_invalid_buffers() {
    assert!(plan(&[], &[], &[0], &[1], 50, None).is_err());
    assert!(plan(&[0; 12], &(-1_i64).to_ne_bytes(), &[0], &[1], 50, None).is_err());
    assert!(plan(&[0; 12], &(1_i64).to_ne_bytes(), &[0], &[1], 50, None).is_err());
    assert!(plan(&[0; 12], &(0_i64).to_ne_bytes(), &[2], &[1], 50, None).is_err());
    assert!(plan(&[0; 12], &(0_i64).to_ne_bytes(), &[0], &[1], 49, None).is_err());
    assert!(plan(&[0; 12], &(0_i64).to_ne_bytes(), &[0], &[1], 50, Some(&[])).is_err());
    assert!(plan(&[0; 12], &(0_i64).to_ne_bytes(), &[0], &[1], 50, Some(&[2])).is_err());
}

#[test]
fn masked_centers_stay_in_the_neighbor_budget_and_cohort() {
    let positions: Vec<u8> = (0..53).flat_map(|_| [0_f32; 3]).flat_map(f32::to_ne_bytes).collect();
    let chains: Vec<u8> = (0..53).map(|row| if row < 50 { 0_i64 } else { 1 }).flat_map(i64::to_ne_bytes).collect();
    let ligand: Vec<_> = (0..53).map(|row| u8::from(row >= 50)).collect();
    let mut eligible = vec![1; 53];
    eligible[52] = 0;
    let result = plan(&positions, &chains, &ligand, &vec![1; 53], 53, Some(&eligible)).unwrap();
    assert_eq!(result.centers, (0..52).collect::<Vec<_>>());
    assert_eq!(result.bits, vec![3; 52]);
    assert!((result.mass - 1.).abs() < 1e-14);
    assert!((result.probabilities[50] - (0.01 / 52. + 0.99 / 2.)).abs() < 1e-14);
    let limited = plan(&positions, &chains, &ligand, &vec![1; 53], 52, Some(&eligible)).unwrap();
    assert_eq!(limited.bits, vec![1; 52]);
    assert_eq!(limited.ties, vec![51; 52]);
    let none = plan(&positions, &chains, &ligand, &vec![1; 53], 53, Some(&vec![0; 53])).unwrap();
    assert_eq!(none.mass, 0.);
    assert!(none.centers.is_empty());
}

#[test]
fn ligand_permutation_keeps_protein_and_chain_coordinates() {
    let bytes = |values: &[f32]| {
        values
            .iter()
            .flat_map(|value| value.to_ne_bytes())
            .collect::<Vec<_>>()
    };
    let left = bytes(&[0., 0., 0., 1., 0., 0., 2., 0., 0.]);
    let right = bytes(&[0., 0., 0., 2., 0., 0., 1., 0., 0.]);
    let chains: Vec<u8> = [0_i64, 1, 1]
        .into_iter()
        .flat_map(|value| value.to_ne_bytes())
        .collect();
    assert_eq!(
        permutation(&left, &right, &chains, &[0, 1, 1]).unwrap(),
        vec![0, 2, 1]
    );
    assert!(permutation(&left, &right, &chains, &[0, 0, 1]).is_err());
    let different = bytes(&[0., 0., 0., 3., 0., 0., 1., 0., 0.]);
    assert!(permutation(&left, &different, &chains, &[0, 1, 1]).is_err());
}

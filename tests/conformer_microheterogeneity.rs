//! A residue with two chemical identities keeps one of them, whatever the seed.
//!
//! 2BWX CYS 249 is a microheterogeneous site: the backbone N, CA, C and O is shared
//! and labelled CYS (0.70), the sulfur side chain is CYS A (0.40) or CYS B (0.30), and
//! CSO C (0.30) is a second residue that labels its own backbone on top of the shared
//! one, with the extra sulfinic oxygen OD. The CSO label is the only one that holds every
//! heavy atom any label holds, so a search that ranks conformers by completeness takes it
//! on every seed, keeps the CYS backbone beside it and fails on two component identities.

#![allow(clippy::expect_used)]

use _core::cif::{CifDocument, parse};
use _core::pdbx::{ConformerPolicy, select_conformer_document_with_policy};
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

/// 2BWX chain A residue 249, model 1, as deposited.
const CYS_CSO_2BWX: &str = r##"
ATOM 2038 N N . CYS A 1 249 ? 68.746 44.381 71.143 0.70 21.04 249 CYS A N 1
ATOM 2039 C CA . CYS A 1 249 ? 68.957 43.022 71.606 0.70 21.28 249 CYS A CA 1
ATOM 2040 C C . CYS A 1 249 ? 70.359 42.507 71.362 0.70 19.80 249 CYS A C 1
ATOM 2041 O O . CYS A 1 249 ? 71.055 42.917 70.439 0.70 19.80 249 CYS A O 1
ATOM 2042 C CB A CYS A 1 249 ? 67.945 42.064 70.987 0.40 24.99 249 CYS A CB 1
ATOM 2043 C CB B CYS A 1 249 ? 67.928 42.101 70.948 0.30 23.34 249 CYS A CB 1
ATOM 2044 S SG A CYS A 1 249 ? 66.261 42.472 71.389 0.40 27.94 249 CYS A SG 1
ATOM 2045 S SG B CYS A 1 249 ? 67.977 40.404 71.507 0.30 26.46 249 CYS A SG 1
HETATM 2046 N N C CSO A 1 249 ? 68.746 44.381 71.143 0.30 21.04 249 CSO A N 1
HETATM 2047 C CA C CSO A 1 249 ? 68.957 43.022 71.606 0.30 21.28 249 CSO A CA 1
HETATM 2048 C CB C CSO A 1 249 ? 67.945 42.064 70.987 0.30 24.99 249 CSO A CB 1
HETATM 2049 S SG C CSO A 1 249 ? 66.261 42.472 71.389 0.30 27.94 249 CSO A SG 1
HETATM 2050 C C C CSO A 1 249 ? 70.359 42.507 71.362 0.30 19.80 249 CSO A C 1
HETATM 2051 O O C CSO A 1 249 ? 71.055 42.917 70.439 0.30 19.80 249 CSO A O 1
HETATM 2052 O OD C CSO A 1 249 ? 66.275 42.201 72.870 0.30 23.67 249 CSO A OD 1
"##;

/// The component of each row of [`CYS_CSO_2BWX`], in order.
const COMPONENTS: [&str; 15] = [
    "CYS", "CYS", "CYS", "CYS", "CYS", "CYS", "CYS", "CYS", "CSO", "CSO", "CSO", "CSO", "CSO",
    "CSO", "CSO",
];

fn document() -> CifDocument {
    let tags: String = TAGS
        .iter()
        .map(|tag| format!("_atom_site.{tag}\n"))
        .collect();
    let text = format!("data_site\nloop_\n{tags}{}\n", CYS_CSO_2BWX.trim());
    parse(text.as_bytes()).expect("source")
}

#[test]
fn a_microheterogeneous_residue_keeps_one_component_on_every_seed() {
    let doc = document();
    for policy in [
        ConformerPolicy::RandomClashAware,
        ConformerPolicy::Occupancy,
    ] {
        for seed in 0..32 {
            let rows = select_conformer_document_with_policy(&doc, None, seed, policy)
                .map_err(|error| format!("{policy:?}, seed {seed}: {error}"))
                .expect("selection")
                .atom_rows;
            let components: BTreeSet<_> = rows.iter().map(|&row| COMPONENTS[row]).collect();
            assert_eq!(components.len(), 1, "{policy:?}, seed {seed}: {rows:?}");
            // The kept component is whole: CYS holds six heavy atoms, CSO seven.
            let expected = if components.contains("CYS") { 6 } else { 7 };
            assert_eq!(rows.len(), expected, "{policy:?}, seed {seed}: {rows:?}");
        }
    }
}

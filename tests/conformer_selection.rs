//! Sparse conformer constraints preserve occupancy preferences and source identity.

#![allow(clippy::expect_used)]

use _core::pdbx::{ConformerAtom, select_atom_conformers};

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
fn independent_occupancy_choices_cannot_create_a_clash() {
    let atoms = [
        atom(0, "A", Some("A"), 0., 0.8),
        atom(0, "A", Some("B"), 10., 0.2),
        atom(1, "A", Some("A"), 0.5, 0.8),
        atom(1, "A", Some("B"), 20., 0.2),
    ];
    let selection = select_atom_conformers(&atoms, &[], 0).expect("assignment");
    assert_eq!(selection.atom_rows, [0, 3]);
}

#[test]
fn shared_sites_block_an_alternate_but_stay_selected() {
    let atoms = [
        atom(0, "A", Some("A"), 0., 0.8),
        atom(0, "A", Some("B"), 10., 0.2),
        atom(1, "B", None, 0.5, 1.),
    ];
    let selection = select_atom_conformers(&atoms, &[], 0).expect("assignment");
    assert_eq!(selection.atom_rows, [1, 2]);
    let bonded = select_atom_conformers(&atoms, &[[0, 2]], 0).expect("bonded contact");
    assert_eq!(bonded.atom_rows, [0, 2]);
}

#[test]
fn unsatisfiable_assignment_fails_instead_of_substituting_a_parent() {
    let atoms = [
        atom(0, "A", Some("A"), 0., 0.8),
        atom(0, "A", Some("B"), 1., 0.2),
        atom(1, "B", None, 0.5, 1.),
    ];
    let error = select_atom_conformers(&atoms, &[], 0).expect_err("clash");
    assert_eq!(error.code(), "PDBX_CONFORMER_CONFLICT");
}

#[test]
fn proximal_cross_chain_letters_agree_on_equal_occupancy() {
    let atoms = [
        atom(0, "A", Some("A"), 0., 0.5),
        atom(0, "A", Some("B"), 10., 0.5),
        atom(1, "B", Some("A"), 3., 0.5),
        atom(1, "B", Some("B"), 13., 0.5),
    ];
    let selection = select_atom_conformers(&atoms, &[], 0).expect("assignment");
    assert_eq!(selection.assignments[&0], selection.assignments[&1]);
    assert_eq!(selection.atom_rows, [0, 2]);
}

#[test]
fn invalid_inputs_fail_and_fixed_only_inputs_skip_search() {
    let atoms = [atom(0, "A", None, 0., 1.)];
    assert!(select_atom_conformers(&atoms, &[[0, 2]], 0).is_err());
    let selection = select_atom_conformers(&atoms, &[], 0).expect("fixed");
    assert_eq!(selection.atom_rows, [0]);
    assert_eq!(selection.search_nodes, 0);
}

#[test]
fn hydrogens_affect_occupancy_but_do_not_make_heavy_atom_clashes() {
    let mut atoms = [
        atom(0, "A", Some("A"), 0., 0.8),
        atom(0, "A", Some("B"), 10., 0.2),
        atom(1, "B", None, 0.5, 1.),
    ];
    atoms[0].clash_eligible = false;
    atoms[1].clash_eligible = false;
    assert_eq!(
        select_atom_conformers(&atoms, &[], 0)
            .expect("selection")
            .atom_rows,
        [0, 2]
    );
}

#[test]
fn symmetry_image_bond_does_not_exempt_an_asu_clash() {
    let tags = [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "auth_seq_id",
        "auth_comp_id",
        "auth_asym_id",
        "auth_atom_id",
        "pdbx_PDB_model_num",
    ]
    .map(|name| format!("_atom_site.{name}"))
    .join("\n");
    let text = format!(
        "data_images\nloop_\n{tags}\n\
        ATOM 1 C C1 A LIG A 1 . 0 0 0 .8 1 LIG A C1 1\n\
        ATOM 2 C C1 B LIG A 1 . 10 0 0 .2 1 LIG A C1 1\n\
        ATOM 3 C C1 . LIG B 2 . .5 0 0 1 1 LIG B C1 1\n\
        loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
        _struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_atom_id\n\
        _struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_atom_id\n\
        _struct_conn.ptnr1_symmetry\n_struct_conn.ptnr2_symmetry\n\
        link covale A C1 B C1 2_555 1_555\n"
    );
    let doc = _core::cif::parse(text.as_bytes()).expect("source");
    let result = _core::pdbx::select_conformer_document(&doc, None, 0).expect("selection");
    assert_eq!(result.atom_rows, [1, 2]);
    let same_image =
        _core::cif::parse(text.replace("2_555 1_555", "2_555 2_555").as_bytes()).expect("source");
    let result = _core::pdbx::select_conformer_document(&same_image, None, 0).expect("selection");
    assert_eq!(result.atom_rows, [0, 2]);
}

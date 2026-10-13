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

#[test]
fn scalar_links_and_anisotrop_follow_the_selected_source_rows() {
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
    let source = format!(
        "data_scalar\n_entry.id scalar\nloop_\n{tags}\n\
        ATOM 1 C C1 A LIG A 1 . 0 0 0 .8 1 LIG A C1 1\n\
        ATOM 2 C C1 B LIG A 1 . 10 0 0 .2 1 LIG A C1 1\n\
        ATOM 3 C C1 . LIG B 2 . .5 0 0 1 1 LIG B C1 1\n\
        _atom_site_anisotrop.id 1\n_atom_site_anisotrop.U[1][1] .01\n\
        _struct_conn.id link\n_struct_conn.conn_type_id covale\n\
        _struct_conn.ptnr1_label_asym_id A\n_struct_conn.ptnr1_label_atom_id C1\n\
        _struct_conn.pdbx_ptnr1_label_alt_id A\n\
        _struct_conn.ptnr2_label_asym_id B\n_struct_conn.ptnr2_label_atom_id C1\n"
    );
    let doc = _core::cif::parse(source.as_bytes()).expect("source");
    let selected = _core::pdbx::selected_atom_document(&doc, &[1, 2]).expect("scalar view");
    let text = _core::cif::write_canonical(&selected).expect("selected text");
    assert!(!text.contains("_struct_conn.id"));
    assert!(!text.contains("_atom_site_anisotrop.id"));
    assert!(text.contains("_entry.id scalar"));
    let retained = source
        .replace("_atom_site_anisotrop.id 1", "_atom_site_anisotrop.id 2")
        .replace(
            "_struct_conn.pdbx_ptnr1_label_alt_id A",
            "_struct_conn.pdbx_ptnr1_label_alt_id B",
        );
    let doc = _core::cif::parse(retained.as_bytes()).expect("retained source");
    let selected =
        _core::pdbx::selected_atom_document(&doc, &[1, 2]).expect("retained scalar view");
    let text = _core::cif::write_canonical(&selected).expect("retained text");
    assert!(text.contains("_struct_conn.id link"));
    assert!(text.contains("_atom_site_anisotrop.id 2"));
    let selected = _core::pdbx::select_conformer_document(&doc, None, 0).expect("scalar conformer");
    assert_eq!(selected.atom_rows, [1, 2]);
}

#[test]
fn a_scalar_atom_site_is_a_single_row_table() {
    let source = b"data_scalar\n_atom_site.group_PDB HETATM\n_atom_site.id 7\n\
        _atom_site.type_symbol C\n_atom_site.label_atom_id C1\n_atom_site.label_alt_id .\n\
        _atom_site.label_comp_id LIG\n_atom_site.label_asym_id A\n_atom_site.label_entity_id 1\n\
        _atom_site.label_seq_id .\n_atom_site.Cartn_x 0\n_atom_site.Cartn_y 0\n_atom_site.Cartn_z 0\n\
        _atom_site.occupancy 1\n_atom_site.auth_seq_id 1\n_atom_site.auth_comp_id LIG\n\
        _atom_site.auth_asym_id A\n_atom_site.auth_atom_id C1\n_atom_site.pdbx_PDB_model_num 1\n";
    let doc = _core::cif::parse(source).expect("scalar atom");
    assert_eq!(
        _core::pdbx::select_conformer_document(&doc, None, 0)
            .expect("scalar selection")
            .atom_rows,
        [0]
    );
    assert_eq!(
        _core::pdbx::selected_atom_document(&doc, &[0]).expect("all rows"),
        doc
    );
    let empty = _core::pdbx::selected_atom_document(&doc, &[]).expect("empty selection");
    assert!(
        !_core::cif::write_canonical(&empty)
            .expect("empty text")
            .contains("_atom_site.id")
    );
}

//! Linked-product rules preserve source identities and fail on unknown chemistry.

#![allow(clippy::expect_used)]

use _core::cif::parse;
use _core::pdbx::{PreparedChemistry, build_model, prepare_chemistry};

fn fixture(atoms: &str, bonds: &str, observed: &str, links: &str) -> String {
    format!(
        "data_product\n_entry.id product\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_leaving_atom_flag\n{atoms}\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n{bonds}\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\n{observed}\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.pdbx_value_order\n{links}\n"
    )
}

fn prepare(text: &str) -> PreparedChemistry {
    let document = parse(text.as_bytes()).expect("valid fixture");
    let model = build_model(&document).expect("model");
    let original = model.clone();
    let product = prepare_chemistry(&model, None, 1, None).expect("valid product");
    assert_eq!(model, original, "source stays unchanged");
    product
}

fn atom(graph: &PreparedChemistry, comp: &str, name: &str) -> usize {
    graph
        .atoms()
        .iter()
        .position(|a| a.component_id == comp && a.atom_id == name)
        .expect("atom")
}

fn order(graph: &PreparedChemistry, a: usize, b: usize) -> u8 {
    graph
        .bonds()
        .iter()
        .find(|edge| {
            edge.active && edge.atom_indices.contains(&a) && edge.atom_indices.contains(&b)
        })
        .expect("bond")
        .order
}

fn sulfur_addition() -> String {
    fixture(
        "L C C 0 N\nL N N 0 N\nL C2 C 0 N\nL HC H 0 N\nL HN H 0 N\nL H1 H 0 N\nL H2 H 0 N\nL H3 H 0 N\nD S S 0 N\nD C C 0 N\nD HS H 0 N\nD H1 H 0 N\nD H2 H 0 N\nD H3 H 0 N",
        "L C N DOUB\nL C C2 SING\nL C HC SING\nL N HN SING\nL C2 H1 SING\nL C2 H2 SING\nL C2 H3 SING\nD S C SING\nD S HS SING\nD C H1 SING\nD C H2 SING\nD C H3 SING",
        "HETATM 1 C C . L A 1 . 0 0 0 1 1 L A C 1\nHETATM 2 N N . L A 1 . 1.3 0 0 1 1 L A N 1\nHETATM 3 S S . D B 2 . -1.8 0 0 1 1 D B S 1",
        "link covale A L C B D S SING",
    )
}

#[test]
fn addition_with_carbon_hydrogen_displaces_it_first() {
    let g = prepare(&sulfur_addition());
    let c = atom(&g, "L", "C");
    assert_eq!(g.atoms()[c].hydrogen_count, Some(0));
    assert!(!g.atoms()[atom(&g, "L", "HC")].active);
    assert_eq!(order(&g, c, atom(&g, "L", "N")), 2);
}

#[test]
fn pi_addition_transfers_donor_hydrogen_without_inventing_coordinates() {
    let text = sulfur_addition()
        .replace("L C HC SING", "L C C3 SING")
        .replace("L HC H 0 N", "L C3 C 0 N");
    let g = prepare(&text);
    let c = atom(&g, "L", "C");
    let n = atom(&g, "L", "N");
    assert_eq!(order(&g, c, n), 1);
    assert_eq!(g.atoms()[n].hydrogen_count, Some(2));
    assert_eq!(g.atoms()[atom(&g, "D", "S")].hydrogen_count, Some(0));
    let transferred = &g.atoms()[atom(&g, "L", "H@N")];
    assert!(transferred.source_row.is_none() && transferred.coordinates.is_none());
    assert!(g.edits().iter().any(|e| e.reason == "pi_bond_addition"));
}

#[test]
fn coordination_does_not_remove_hydrogens_or_change_charge() {
    let g = prepare(&sulfur_addition().replace("link covale", "link metalc"));
    let s = atom(&g, "D", "S");
    assert_eq!(g.atoms()[s].formal_charge, Some(0));
    assert_eq!(g.atoms()[s].hydrogen_count, Some(1));
    assert!(g.edits().is_empty());
    assert_eq!(order(&g, s, atom(&g, "L", "C")), 0);
}

#[test]
fn unsupported_explicit_order_fails() {
    let text = sulfur_addition().replace("S SING\n", "S QUAD\n");
    let document = parse(text.as_bytes()).expect("fixture");
    let model = build_model(&document).expect("model");
    assert_eq!(
        prepare_chemistry(&model, None, 1, None)
            .expect_err("unsupported order")
            .code(),
        "PDBX_LINK_PRODUCT"
    );
}

#[test]
fn model_number_is_a_deposited_identifier() {
    let text = sulfur_addition()
        .replace(" L A C 1", " L A C 7")
        .replace(" L A N 1", " L A N 7")
        .replace(" D B S 1", " D B S 7");
    let document = parse(text.as_bytes()).expect("fixture");
    let model = build_model(&document).expect("model");
    assert!(prepare_chemistry(&model, None, 1, None).is_err());
    assert!(prepare_chemistry(&model, None, 7, None).is_ok());
}

#[test]
fn source_heavy_coordinates_survive_product_repairs() {
    let g = prepare(&sulfur_addition());
    assert_eq!(
        g.atoms()[atom(&g, "L", "C")].coordinates,
        Some([0., 0., 0.])
    );
    assert_eq!(
        g.atoms()[atom(&g, "D", "S")].coordinates,
        Some([-1.8, 0., 0.])
    );
}

#[test]
fn observed_deuterium_uses_hydrogen_template_without_losing_isotope() {
    let text = sulfur_addition().replace("link covale", "link metalc").replace(
        "HETATM 3 S S . D B 2 . -1.8 0 0 1 1 D B S 1",
        "HETATM 3 S S . D B 2 . -1.8 0 0 1 1 D B S 1\nHETATM 4 D HS . D B 2 . -2.8 0 0 1 1 D B HS 1");
    let g = prepare(&text);
    assert_eq!(g.atoms()[atom(&g, "D", "HS")].element, "D");
}

#[test]
fn malformed_leaving_flag_does_not_become_unknown() {
    let text = sulfur_addition().replace("L C C 0 N", "L C C 0 broken");
    let doc = parse(text.as_bytes()).expect("valid syntax");
    assert_eq!(
        build_model(&doc).expect_err("invalid flag").code(),
        "PDBX_ITEM_TYPE"
    );
}

#[test]
fn ccd_alias_keeps_the_observed_source_row() {
    let text = sulfur_addition().replace("link covale", "link metalc")
        .replace("_chem_comp_atom.pdbx_leaving_atom_flag", "_chem_comp_atom.alt_atom_id\n_chem_comp_atom.pdbx_leaving_atom_flag")
        .replace(" 0 N", " 0 . N")
        .replace("L C C 0 . N", "L C C 0 CX N")
        .replace("HETATM 1 C C .", "HETATM 1 C CX .")
        .replace("link metalc A L C", "link metalc A L CX");
    let graph = prepare(&text);
    let carbon = &graph.atoms()[atom(&graph, "L", "C")];
    assert_eq!(carbon.source_row, Some(0));
    assert_eq!(carbon.coordinates, Some([0., 0., 0.]));
}

#[test]
fn canonical_names_take_precedence_over_overlapping_aliases() {
    let text = sulfur_addition().replace("link covale", "link metalc")
        .replace("_chem_comp_atom.pdbx_leaving_atom_flag", "_chem_comp_atom.alt_atom_id\n_chem_comp_atom.pdbx_leaving_atom_flag")
        .replace(" 0 N", " 0 . N")
        .replace("L C C 0 . N", "L C C 0 C2 N");
    let graph = prepare(&text);
    assert_eq!(graph.atoms()[atom(&graph, "L", "C")].source_row, Some(0));
    assert_eq!(graph.atoms()[atom(&graph, "L", "C2")].source_row, None);
}

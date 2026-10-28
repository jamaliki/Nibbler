use super::read;
use crate::cif::parse;

const HEADER: &str = "data_test\nloop_\n_atom_site.label_atom_id\n_atom_site.type_symbol\n_atom_site.label_asym_id\n_atom_site.label_comp_id\n_atom_site.label_entity_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.label_alt_id\n_atom_site.group_PDB\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.pdbx_PDB_model_num\n";

#[test]
fn keeps_selected_atom_order_and_numeric_missing_kinds() {
    let document = parse(format!("{HEADER}N N A ALA 1 8 . . ATOM 1 1.25(3) 2 3 0.2 2\nCA C A ALA 1 8 . A ATOM 1 4 5 6 0.8 2\nN N A ALA 1 8 . A ATOM 1 9 9 9 1.0 2\nO O B HOH 2 9 ? . HETATM . ? 0 0 . 2\n").as_bytes()).expect("fixture");
    let source = read(&document, 2).expect("selected source");
    assert_eq!(source.atoms.len(), 3);
    assert_eq!(source.atoms[0].position, [1.25, 2.0, 3.0]);
    assert_eq!(source.atoms[1].name, "CA");
    assert_eq!(source.atoms[2].sequence, -1);
    assert_eq!(source.atoms[2].insertion, "?");
    assert!(source.atoms[2].position[0].is_nan());
    assert_eq!(source.groups, [vec![0, 1], vec![2]]);
    assert!(source.atoms[0].protein && !source.atoms[2].protein);
    assert!(!source.has_protein_annotations);
}

#[test]
fn keeps_occupancy_choice_and_entity_classification() {
    let document = parse(format!("{HEADER}N N A ALA 1 1 . . HETATM 1 0 0 0 0.2 1\nN N A ALA 1 1 . . HETATM 1 1 0 0 0.8 1\nCA C A ALA 1 1 . . HETATM 1 2 0 0 1 1\n#\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n1 'polypeptide(D)'\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n1 1\n1 7\n").as_bytes()).expect("fixture");
    let source = read(&document, 1).expect("source");
    assert_eq!(source.atoms[0].position[0], 1.0);
    assert!(source.atoms.iter().all(|atom| atom.protein));
    assert!(source.has_protein_annotations);
    assert_eq!(source.lengths["1"], 7);
}

#[test]
fn rejects_unselected_alternates_models_and_assemblies() {
    let input = format!(
        "{HEADER}N N A ALA 1 1 . A ATOM 1 0 0 0 1 1\nCA C A ALA 1 1 . B ATOM 1 1 0 0 1 1\n"
    );
    assert!(read(&parse(input.as_bytes()).expect("fixture"), 1).is_err());
    let input = input.replace(" . B ATOM", " . A ATOM");
    assert!(read(&parse(input.as_bytes()).expect("fixture"), 2).is_err());
    let input = format!("{input}_pdbx_struct_assembly_gen.assembly_id 1\n");
    assert!(read(&parse(input.as_bytes()).expect("fixture"), 1).is_err());
}

#[test]
fn rejects_malformed_present_coordinates_without_inventing_missingness() {
    let input = format!("{HEADER}N N A ALA 1 1 . . ATOM 1 wrong 0 0 1 1\n");
    let error = read(&parse(input.as_bytes()).expect("fixture"), 1)
        .err()
        .expect("error");
    assert!(error.contains("Invalid coordinate number"));
}

#[test]
fn keeps_only_identity_symmetry_covalent_links() {
    let input = format!(
        "{HEADER}SG S A CYS 1 1 . . ATOM 1 0 0 0 1 1\nSG S A CYS 1 2 . A ATOM 2 2 0 0 1 1\nC1 C B LIG 2 3 . . HETATM . 3 0 0 1 1\n#\nloop_\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr1_symmetry\n_struct_conn.ptnr2_symmetry\n_struct_conn.pdbx_ptnr2_label_alt_id\ndisulf A 1 CYS SG A 2 CYS SG 1_555 . A\ncovale A 1 CYS SG B . LIG C1 . ? .\nhydrog A 2 CYS SG B . LIG C1 . . .\ncovale A 2 CYS SG B . LIG C1 2_555 . .\n"
    );
    let source = read(&parse(input.as_bytes()).expect("fixture"), 1).expect("source");
    assert_eq!(source.bonds, [(0, 1), (0, 2)]);
}

/// Every value of a source, as text.
#[cfg(feature = "reduce3")]
fn describe(source: &super::Source) -> String {
    use std::collections::BTreeMap;
    let atoms: Vec<_> = source
        .atoms
        .iter()
        .map(|atom| {
            format!(
                "{} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {} {:?} {} {} {}",
                atom.source_row,
                atom.name,
                atom.element,
                atom.chain,
                atom.component,
                atom.entity,
                atom.auth_sequence,
                atom.auth_chain,
                atom.auth_component,
                atom.insertion,
                atom.alternate,
                atom.sequence,
                atom.position.map(f32::to_bits),
                atom.occupancy.to_bits(),
                atom.protein,
                atom.group
            )
        })
        .collect();
    let connections: Vec<_> = source
        .connections
        .iter()
        .map(|link| {
            format!(
                "{} {:?} {} {} {:?}",
                link.row, link.atoms, link.kind, link.order, link.symmetries
            )
        })
        .collect();
    format!(
        "{atoms:?} {:?} {:?} {:?} {connections:?} {} {:?} {:?}",
        source.groups,
        source.bonds,
        source.hydrogen_bonds,
        source.has_protein_annotations,
        source.lengths.iter().collect::<BTreeMap<_, _>>(),
        source.entity_types.iter().collect::<BTreeMap<_, _>>()
    )
}

#[cfg(feature = "reduce3")]
const COPIED: &str = "data_copy\n_entry.id COPY\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n1 'polypeptide(L)'\n\
loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n1 1 GLY\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n\
_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n\
_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n\
_atom_site.B_iso_or_equiv\n_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n\
_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\n\
ATOM 1 N N . GLY A 1 1 ? 0.0 0.0 0.0 1.00 10.0 ? 1 GLY A N 1\n\
ATOM 2 C CA . GLY A 1 1 ? 1.4 0.0 0.0 1.00 10.0 ? 1 GLY A CA 1\n\
ROWS\
HETATM 9 ZN ZN . ZN B 2 . ? ZINC 1.00 20.0 2 2 ZN B ZN 1\n\
loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_comp_id\n\
_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_asym_id\n\
_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\n\
metalc1 metalc A GLY 1 N B ZN . ZN\n";

#[cfg(feature = "reduce3")]
#[test]
fn copied_rows_give_the_source_their_own_model_gives() {
    let original = parse(
        COPIED
            .replace("ROWS", "")
            .replace("ZINC", "5.0 5.0 5.0")
            .as_bytes(),
    )
    .expect("fixture");
    let written = COPIED
        .replace(
            "ROWS",
            "ATOM 3 H HA2 . GLY A 1 1 ? 1.8 1.0 0.0 1.00 12.0 ? 1 GLY A HA2 1\n",
        )
        .replace("ZINC", "6.0 5.0 5.0");
    let document = parse(written.as_bytes()).expect("fixture");
    let model = crate::pdbx::build_model(&original).expect("model");
    let base = super::read_with_sites(&original, 1, model.loop_sites()).expect("source");
    let rows = [0, 1, 1, 2];
    let hydrogen = [1, 2, 3, 4, 10, 11, 12, 13, 14, 15, 19]
        .iter()
        .fold(0_u64, |mask, bit| mask | 1 << bit);
    let made = [
        1 << 1,
        1 << 1,
        hydrogen,
        1 << 1 | 1 << 10 | 1 << 11 | 1 << 12,
    ];
    let copied =
        super::Copied::new(&document, &model, &rows, &made, 1, &base).expect("copied rows");
    let own = crate::pdbx::build_model(&document).expect("model");
    let expected = super::read_with_sites(&document, 1, own.loop_sites()).expect("source");
    assert_eq!(describe(&copied.source), describe(&expected));
    assert_eq!(
        (
            copied.source.atoms.len(),
            copied.source.connections[0].atoms
        ),
        (4, [Some(0), Some(3)])
    );
    assert_eq!(copied.heavy(), [0, 1, 3]);
    assert_eq!(copied.position(3), [6.0, 5.0, 5.0]);
    assert!(copied.copies_name(1, 1) && !copied.copies_name(2, 1));
    // A row that writes a residue field leaves the document to its own model.
    let mut changed = made;
    changed[2] |= 1 << 6;
    assert!(super::Copied::new(&document, &model, &rows, &changed, 1, &base).is_none());
    // So does a source site that no row copies.
    assert!(super::Copied::new(&document, &model, &[0, 1, 1, 1], &made, 1, &base).is_none());
}

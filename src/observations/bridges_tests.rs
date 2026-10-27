use super::*;
use reduce3::monlib::{CompAtom, CompBond};

#[path = "bridge_token_tests.rs"]
mod token_tests;

/// kaveh-data's 66-token table: 43 elements from 21, other elements 64 and hydrogen 65.
fn token_table() -> crate::observations::TokenTable {
    let elements = "C O N P S MG ZN CA FE F CL MN K NI CU BR CO CD HG AS B I PT MO BE AL BA RU V CS XE SE W AU GD YB LI IR RB Y PB OS AG";
    crate::observations::TokenTable::new(elements.split_whitespace().map(str::to_owned).zip(21..), 64, 65)
        .expect("The 66-token table must be valid")
}

fn library() -> MonLib {
    let root = MonLib::locate(None).expect("Bridge source tests require REDUCE3_CHEM_DATA");
    MonLib::load(&root).expect("Bridge source tests require a readable dictionary")
}

fn source(component: &str, observed: bool) -> Source {
    let mut rows = match component {
        "NWN" => concat!(
            "FE FE 19.249 5.075 -13.055\nNI NI 17.788 5.315 -15.164\n",
            "C C1 21.026 5.409 -13.514\nN N1 22.114 5.592 -13.838\n",
            "C C2 19.620 3.270 -12.628\nN N2 19.817 2.145 -12.312\n",
            "C C3 19.499 5.427 -11.383\nO O3 19.663 5.636 -10.236\n",
        ),
        "EJ2" => concat!(
            "FE FE -4.600 -5.078 26.969\nNI NI -5.053 -2.663 26.213\n",
            "C C1 -3.495 -5.027 28.482\nN N1 -2.796 -4.998 29.405\n",
            "C C2 -3.227 -5.893 25.921\nN N2 -2.353 -6.357 25.319\n",
            "C C3 -5.225 -6.686 27.290\nO O3 -5.613 -7.733 27.462\n",
        ),
        _ => panic!("Unknown test component"),
    }.lines().map(|line| format!("{line} L 601 {component} 1\n")).collect::<String>();
    if observed { rows += &format!("H H -3.806 -3.506 26.700 L 601 {component} 1\n"); }
    let header = concat!(
        "data_bridge\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n",
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
        "_atom_site.label_asym_id\n_atom_site.auth_seq_id\n_atom_site.label_comp_id\n",
        "_atom_site.pdbx_PDB_model_num\n",
    );
    let document = cif::parse(format!("{header}{rows}").as_bytes()).expect("Source fixture must parse");
    super::super::source::read(&document, 1).expect("Source fixture must read")
}

#[test]
fn selected_nwn_and_ej2_dictionaries_keep_both_negative_hydride_bonds() {
    let library = library();
    for component in ["NWN", "EJ2"] {
        let input = source(component, false);
        let names: Vec<_> = input.atoms.iter().map(|atom| atom.name.as_str().to_owned()).collect();
        let comp = reduce3::hplace::hydrogen_dictionary(&library, component, &names, false).expect("Source dictionary must exist");
        let definitions = dictionary(&library, &comp).expect("Source hydride definition must read");
        assert_eq!(definitions.len(), 1);
        let definition = definitions.get("H").expect("Source dictionary must retain its hydride");
        assert_eq!(definition.metals, ["FE", "NI"]);
        assert_eq!(definition.metal_elements, ["FE", "NI"]);
        assert_eq!(definition.formal_charge, -1);
        assert_eq!(definition.kind, EdgeKind::MetalHydride);
    }
}

#[test]
fn absent_hydride_coordinates_keep_both_original_metal_dependencies() {
    let library = library();
    for component in ["NWN", "EJ2"] {
        let input = source(component, false);
        let original: Vec<_> = input.atoms.iter().map(|atom| atom.position.map(f32::to_bits)).collect();
        let sites = for_source(&library, &input).expect("Source bridge must resolve");
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].hydrogen_name, "H");
        assert_eq!(sites[0].hydrogen_row, None);
        assert_eq!(sites[0].metal_rows, [0, 1]);
        assert_eq!(sites[0].formal_charge, -1);
        assert_eq!(input.atoms.iter().map(|atom| atom.position.map(f32::to_bits)).collect::<Vec<_>>(), original);
        let mapping = [8, 9, 10, 11, 12, 13, 14, 15];
        assert_eq!(sites[0].edges(16, &mapping).expect("Both metal tokens must map"), [
            Edge { hydrogen: 16, metal: 8, kind: EdgeKind::MetalHydride },
            Edge { hydrogen: 16, metal: 9, kind: EdgeKind::MetalHydride },
        ]);
    }
}

#[test]
fn observed_hydride_keeps_its_exact_row_and_coordinate() {
    let library = library();
    let input = source("EJ2", true);
    let original = input.atoms[8].position.map(f32::to_bits);
    let sites = for_source(&library, &input).expect("Observed bridge must resolve");
    assert_eq!(sites[0].hydrogen_row, Some(8));
    assert_eq!(input.atoms[8].position.map(f32::to_bits), original);
    assert_eq!(sites[0].metal_rows, [0, 1]);
}

#[test]
fn observed_bridge_isotopes_keep_original_rows_coordinates_and_canonical_key() {
    let library = library();
    for (element, name) in [("H", "H"), ("D", "H"), ("T", "H"), ("D", "D"), ("T", "T")] {
        let mut input = source("EJ2", true);
        input.atoms[8].name = name.into();
        input.atoms[8].element = element.into();
        let bits = input.atoms[8].position.map(f32::to_bits);
        let sites = for_source(&library, &input).expect("Observed isotope bridge must resolve");
        assert_eq!(sites[0].hydrogen_row, Some(8));
        assert_eq!(sites[0].hydrogen_name, "H");
        assert_eq!(sites[0].metal_rows, [0, 1]);
        let schema = crate::observations::Schema { slots: 32, hydrogen_start: 15, residues: Vec::new() };
        let parents = crate::observations::parents::for_source(&library, &input).expect("Ordinary registry must resolve");
        let options = crate::observations::PackOptions { model: 1, complete: true,
            topology_policy: "observed_connectivity", protein_only: false, chem_data: None };
        let tokens = crate::observations::tokens::pack(&input, &schema, &HashMap::new(), &HashMap::new(),
            &parents, &sites, &token_table(), &[], &options).expect("Isotope must stay atom-level");
        assert_eq!(tokens.keys[8][4], "H");
        assert_eq!(tokens.identity[8], 65);
        assert_eq!(tokens.positions[8 * 96..][..3].iter().copied().map(f32::to_bits).collect::<Vec<_>>(), bits);
        assert_eq!(tokens.bridge_edges, [(8, 0, 1), (8, 1, 1)]);
        assert_eq!(tokens.masks[0..32][15..18], [0; 3]);
        assert_eq!(tokens.masks[32..64][15..18], [0; 3]);
        assert_eq!(input.atoms[8].name, name);
        assert_eq!(input.atoms[8].element, element);
    }
}

#[test]
fn ambiguous_or_mistyped_bridge_isotope_aliases_fail() {
    let library = library();
    for (element, name) in [("C", "D"), ("H", "D"), ("T", "D"), ("C", "T"), ("D", "T")] {
        let mut input = source("EJ2", true);
        input.atoms[8].name = name.into();
        input.atoms[8].element = element.into();
        assert!(for_source(&library, &input).is_err(), "{element} {name}");
    }
    for name in ["D", "T"] {
        let mut input = source("EJ2", true);
        input.atoms.push(super::super::source::Atom { source_row: 9, name: name.into(), element: name.into(),
            chain: "L".into(), component: "EJ2".into(), entity: ".".into(), auth_sequence: "601".into(),
            auth_chain: "L".into(), auth_component: "EJ2".into(),
            insertion: ".".into(), alternate: ".".into(), sequence: -1,
            position: input.atoms[8].position, occupancy: 1.0, protein: false, group: 0 });
        input.groups[0].push(9);
        assert!(for_source(&library, &input).is_err(), "H+{name}");
    }
}

#[test]
fn ambiguous_mistyped_or_nonfinite_observed_hydride_fails() {
    let library = library();
    for invalid in 0..3 {
        let mut input = source("EJ2", true);
        match invalid {
            0 => { input.groups[0].push(8); }
            1 => input.atoms[8].element = "C".into(),
            2 => input.atoms[8].position[0] = f32::NAN,
            _ => unreachable!(),
        }
        assert!(for_source(&library, &input).is_err(), "case {invalid}");
    }
}

#[test]
fn missing_ambiguous_or_mistyped_metal_endpoints_fail() {
    let library = library();
    for invalid in 0..4 {
        let mut input = source("EJ2", false);
        match invalid {
            0 => { input.groups[0].remove(1); }
            1 => { input.groups[0].push(1); }
            2 => input.atoms[1].element = "C".into(),
            3 => input.atoms[1].position[0] = f32::NAN,
            _ => unreachable!(),
        }
        assert!(for_source(&library, &input).is_err(), "case {invalid}");
    }
}

#[test]
fn emitting_only_one_metal_edge_is_an_error() {
    let library = library();
    let input = source("EJ2", false);
    let sites = for_source(&library, &input).expect("Bridge must resolve");
    for mapping in [[0, usize::MAX], [0, 0], [2, 1]] {
        assert!(sites[0].edges(2, &mapping).is_err());
    }
}

#[test]
fn ordinary_single_parent_and_nonmetal_conflicts_are_not_bridges() {
    let atom = |id: &str, element: &str| CompAtom { id: id.into(), type_symbol: element.into(), ..Default::default() };
    let bond = |first: &str, second: &str| CompBond { a1: first.into(), a2: second.into(), type_: "single".into(), ..Default::default() };
    let mut comp = Comp { atoms: vec![atom("H", "H"), atom("FE", "FE"), atom("NI", "NI")],
        bonds: vec![bond("H", "FE")], ..Default::default() };
    assert!(candidates(&comp).is_empty());
    comp.bonds.push(bond("H", "NI"));
    assert_eq!(candidates(&comp).len(), 1);
    comp.atoms[2].type_symbol = "C".into();
    assert!(candidates(&comp).is_empty());
    comp.atoms[2].type_symbol = "NI".into();
    comp.bonds[1].type_ = "double".into();
    assert!(candidates(&comp).is_empty());
    comp.bonds[1].type_ = "single".into();
    comp.atoms.push(atom("C", "C"));
    comp.bonds.push(bond("H", "C"));
    assert!(candidates(&comp).is_empty());
}

#[test]
fn both_metal_flags_and_known_negative_charge_are_required() {
    let library = library();
    let comp = library.ccd_comp("EJ2", false).expect("EJ2 must exist");
    let candidate = candidates(&comp).remove("H").expect("EJ2 must have two metal neighbors");
    let bytes = std::fs::read(&comp.source).expect("EJ2 source must read");
    let document = cif::parse(&bytes).expect("EJ2 source must parse");
    assert!(source_definition(&comp, "H", &candidate, &document));
    for (charge, second_order, second_flag) in [("0", "SING", "Y"), ("?", "SING", "Y"),
        ("-1", "DOUB", "Y"), ("-1", "SING", "N")] {
        let changed = format!(concat!(
            "data_test\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n",
            "_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\nEJ2 H H {}\n",
            "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n",
            "_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_metal_coordination_flag\n",
            "EJ2 FE H SING Y\nEJ2 NI H {} {}\n",
        ), charge, second_order, second_flag);
        let source = cif::parse(changed.as_bytes()).expect("Changed fixture must parse");
        assert!(!source_definition(&comp, "H", &candidate, &source));
    }
}

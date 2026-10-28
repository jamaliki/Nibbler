use super::{library, source, token_table};
use crate::observations::source::Source;
use crate::observations::{self, HydrogenParents, LigandFeatures, PackOptions, Schema, TokenTable};
use std::collections::HashMap;

fn schema() -> Schema {
    Schema {
        slots: 32,
        hydrogen_start: 15,
        residues: Vec::new(),
    }
}

/// The 66-token table with hydrogen moved to `hydrogen`, which no element token uses.
fn hydrogen_at(hydrogen: i64) -> TokenTable {
    let elements = "C O N P S MG ZN CA FE F CL MN K NI CU BR CO CD HG AS B I PT MO BE AL BA RU V CS XE SE W AU GD YB LI IR RB Y PB OS AG";
    TokenTable::new(
        elements.split_whitespace().map(str::to_owned).zip(21..),
        64,
        hydrogen,
    )
    .expect("The moved hydrogen token must give a valid table")
}

fn packed(source: &Source) -> observations::Tokens {
    let library = library();
    let sites = super::super::for_source(&library, source).expect("Bridge must resolve");
    let parents =
        observations::parents::for_source(&library, source).expect("Unique parents must resolve");
    let features: LigandFeatures = HashMap::from([(
        source.atoms[0].component.as_str().to_owned(),
        HashMap::from([("FE".into(), [0, 0, 1, 1]), ("NI".into(), [0, 0, 1, 1])]),
    )]);
    observations::tokens::pack(
        source,
        &schema(),
        &features,
        &HydrogenParents::new(),
        &parents,
        &sites,
        &token_table(),
        &[],
        &PackOptions {
            model: 1,
            complete: true,
            topology_policy: "observed_connectivity",
            protein_only: false,
            chem_data: None,
        },
    )
    .expect("Bridge tokens must pack")
}

#[test]
fn product_targets_preserve_the_checked_hydride_charge_and_both_bonds() {
    for component in ["NWN", "EJ2"] {
        for observed in [false, true] {
            let mut tokens = packed(&source(component, observed));
            let positions = tokens.positions.clone();
            observations::product::apply::<&str>(&mut tokens, &schema(), &token_table(), &[])
                .expect("Empty product targets must preserve checked bridge chemistry");
            assert_eq!(tokens.charges[8], -1.);
            assert_eq!(tokens.charge_known[8], 1);
            assert_eq!(tokens.hydrogen_known[8], 1);
            assert_eq!(tokens.bridge_edges, [(8, 0, 1), (8, 1, 1)]);
            assert_eq!(tokens.positions, positions);
            assert_eq!(tokens.charge_known[0..2], [0, 0]);
            tokens
                .validate(&schema(), &token_table())
                .expect("Product targets must retain valid hydride labels");
        }
    }
}

#[test]
fn absent_hydride_has_unknown_coordinates_and_known_charge_with_both_bonds() {
    for component in ["NWN", "EJ2"] {
        let input = source(component, false);
        let tokens = packed(&input);
        assert_eq!(tokens.identity.len(), 9);
        assert_eq!(tokens.identity[0..2], [29, 34]);
        assert_eq!(tokens.identity[8], 65);
        assert_eq!(tokens.keys[8][4], "H");
        assert_eq!(tokens.positions[8 * 96..], vec![0.; 96]);
        assert_eq!(tokens.masks[8 * 32..], vec![0; 32]);
        assert_eq!(tokens.presence[8 * 32..], vec![0; 32]);
        assert_eq!(tokens.charges[8], -1.);
        assert_eq!(tokens.charge_known[8], 1);
        assert_eq!(tokens.hydrogen_known[8], 1);
        assert_eq!(tokens.charge_known[0..2], [1, 1]);
        assert_eq!(tokens.bridge_edges, [(8, 0, 1), (8, 1, 1)]);
        assert_eq!(tokens.crop_center_eligible, [1, 1, 1, 1, 1, 1, 1, 1, 0]);
        for axis in 0..3 {
            assert_eq!(
                tokens.crop_positions[8 * 3 + axis],
                input.atoms[0].position[axis] * 0.5 + input.atoms[1].position[axis] * 0.5
            );
        }
        for (row, atom) in input.atoms.iter().enumerate() {
            assert_eq!(
                tokens.positions[row * 96..][..3]
                    .iter()
                    .copied()
                    .map(f32::to_bits)
                    .collect::<Vec<_>>(),
                atom.position.map(f32::to_bits)
            );
            assert_eq!(tokens.masks[row * 32 + 15..row * 32 + 18], [0; 3]);
        }
        tokens
            .validate(&schema(), &token_table())
            .expect("Bridge labels must validate");
    }
}

#[test]
fn observed_hydride_stays_atom_level_at_its_original_position() {
    let input = source("EJ2", true);
    let tokens = packed(&input);
    assert_eq!(tokens.positions[8 * 96..][..3], input.atoms[8].position);
    assert_eq!(tokens.crop_positions[8 * 3..], input.atoms[8].position);
    assert_eq!(tokens.masks[8 * 32], 1);
    assert_eq!(tokens.masks[8 * 32 + 1..], vec![0; 31]);
    assert_eq!(tokens.bridge_edges, [(8, 0, 1), (8, 1, 1)]);
}

#[test]
fn selected_registry_skips_only_proven_bridge_names() {
    let library = library();
    for component in ["NWN", "EJ2"] {
        let comp = library
            .ccd_comp(component, false)
            .expect("CCD component must exist");
        let parents =
            observations::parents::dictionary(&library, &comp).expect("Registry must resolve");
        assert!(!parents.contains_key("H"));
        if component == "NWN" {
            assert_eq!(parents.get("H1").map(String::as_str), Some("C3"));
            assert_eq!(parents.get("H2").map(String::as_str), Some("C3"));
            assert_eq!(parents.get("H3").map(String::as_str), Some("O3"));
        }
        let mut changed = (*comp).clone();
        changed.applied_mods.push("test-carbon-parent".into());
        changed
            .atoms
            .iter_mut()
            .find(|atom| atom.id == "NI")
            .expect("NI must exist")
            .type_symbol = "C".into();
        assert!(observations::parents::dictionary(&library, &changed).is_err());
    }
}

#[test]
fn observed_bridge_declarations_keep_both_metals_and_reject_other_parents() {
    let library = library();
    let mut input = source("EJ2", true);
    input.hydrogen_bonds = vec![(8, Some(0)), (8, Some(1)), (8, Some(0))];
    let sites = super::super::for_source(&library, &input)
        .expect("Repeated correct declarations must resolve");
    assert_eq!(sites[0].metal_rows, [0, 1]);
    assert_eq!(packed(&input).bridge_edges, [(8, 0, 1), (8, 1, 1)]);
    for parent in [None, Some(2)] {
        input.hydrogen_bonds.push((8, parent));
        assert!(super::super::for_source(&library, &input).is_err());
        input.hydrogen_bonds.pop();
    }
}

#[test]
fn key_alignment_and_selection_keep_both_edges_or_fail_without_mutation() {
    let input = source("EJ2", false);
    let mut tokens = packed(&input);
    let original = tokens.positions.clone();
    assert!(tokens.select(&[8, 0]).is_err());
    assert_eq!(tokens.positions, original);
    let mut wanted = tokens.keys.clone();
    wanted.reverse();
    tokens
        .align(&serde_json::to_string(&wanted).expect("Keys must encode"))
        .expect("Keys must align");
    assert_eq!(tokens.keys, wanted);
    assert_eq!(tokens.bridge_edges, [(0, 8, 1), (0, 7, 1)]);
    assert_eq!(tokens.crop_center_eligible[0], 0);
    tokens
        .complete(&schema())
        .expect("Complete filter must retain ligand cohort");
    tokens
        .select(&[0, 8, 7])
        .expect("Both metal dependencies must survive");
    assert_eq!(tokens.bridge_edges, [(0, 1, 1), (0, 2, 1)]);
    tokens
        .validate(&schema(), &token_table())
        .expect("Remapped labels must validate");
    tokens
        .select(&[1, 2])
        .expect("Dropping H must drop both edges");
    assert!(tokens.bridge_edges.is_empty());
}

#[test]
fn ordinary_single_parent_hydrogens_keep_their_slots_whatever_the_hydrogen_token() {
    let raw = b"data_test\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.label_asym_id\n_atom_site.auth_seq_id\n_atom_site.label_comp_id\n_atom_site.pdbx_PDB_model_num\nC C1 0 0 0 L 1 LIG 1\nH H1 0.9 0 0 L 1 LIG 1\n";
    let document = crate::cif::parse(raw).expect("Ordinary fixture must parse");
    let input = observations::source::read(&document, 1).expect("Ordinary fixture must read");
    let parents = HashMap::from([("LIG".into(), HashMap::from([("H1".into(), "C1".into())]))]);
    let features = HashMap::from([("LIG".into(), HashMap::from([("C1".into(), [0, 0, 1, 1])]))]);
    let options = PackOptions {
        model: 1,
        complete: true,
        topology_policy: "observed_connectivity",
        protein_only: false,
        chem_data: None,
    };
    let pack = |hydrogen| {
        observations::tokens::pack(
            &input,
            &schema(),
            &features,
            &parents,
            &HashMap::new(),
            &[],
            &hydrogen_at(hydrogen),
            &[],
            &options,
        )
        .expect("Ordinary H must pack")
    };
    let (before, after) = (pack(65), pack(90));
    assert_eq!(before.positions, after.positions);
    assert_eq!(before.masks, after.masks);
    assert_eq!(before.identity, after.identity);
    assert_eq!(before.charge_known, after.charge_known);
    assert_eq!(after.positions[45..48], [0.9, 0., 0.]);
    assert!(after.bridge_edges.is_empty());
}

#[test]
fn virtual_h_key_needs_both_original_metal_rows_and_takes_the_table_hydrogen_token() {
    let raw = b"data_test\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.label_asym_id\n_atom_site.auth_seq_id\n_atom_site.label_comp_id\n_atom_site.pdbx_PDB_model_num\nFE FE 0 0 0 L 601 EJ2 1\nNI NI 2 0 0 L 601 EJ2 1\nC C1 0 0 2 L 601 EJ2 1\nN N1 0 0 3 L 601 EJ2 1\nC C2 0 2 0 L 601 EJ2 1\nN N2 0 3 0 L 601 EJ2 1\nC C3 0 -2 0 L 601 EJ2 1\nO O3 0 -3 0 L 601 EJ2 1\n";
    let document = crate::cif::parse(raw).expect("Bridge fixture must parse");
    let keys = r#"[["L","601",".","EJ2","H"]]"#;
    let library = library();
    assert_eq!(
        observations::required_parent_sources(&document, keys, 1, Some(&library.root))
            .expect("Both parents must resolve"),
        [0, 1]
    );
    let source = observations::source::read(&document, 1).expect("Bridge fixture must read");
    let keys = observations::ParentKeys::parse(keys, &[]).expect("Bridge keys must parse");
    let retained = observations::RetainedParents::new(&source, &keys, Some(&library.root))
        .expect("Bridge parents must resolve");
    assert_eq!(retained.ligand_sources(), [0, 1]);
    assert_eq!(
        retained
            .bridge_sources(Some(&library.root))
            .expect("Bridge declaration must resolve"),
        [("H".into(), [0, 1])]
    );
    for chem_data in [Some(library.root.as_path()), None] {
        let options = PackOptions {
            model: 1,
            complete: false,
            topology_policy: "observed_connectivity",
            protein_only: false,
            chem_data,
        };
        let tokens = observations::pack(
            &document,
            &schema(),
            &HashMap::new(),
            &HashMap::new(),
            &hydrogen_at(90),
            &[],
            &options,
        )
        .expect("A table with any hydrogen token must pack the hydride");
        assert_eq!(tokens.identity, [29, 34, 21, 23, 21, 23, 21, 22, 90]);
        assert_eq!(tokens.keys[8][4], "H");
        tokens
            .validate(&schema(), &hydrogen_at(90))
            .expect("The packing table must validate the hydride");
        let error = tokens
            .validate(&schema(), &token_table())
            .expect_err("Another hydrogen token must not validate");
        assert!(error.contains("invalid endpoints"), "{error}");
    }
}

#[test]
fn compact_key_scan_matches_serde_or_defers_to_it() {
    let serde = |text: &str| serde_json::from_str::<Vec<[String; 5]>>(text).ok();
    for text in [
        "[]",
        r#"[["A","1","","ALA","protein"]]"#,
        r#"[["A","1","","ALA","protein"],["B","2","x","ATP","O5'"]]"#,
    ] {
        let scanned = observations::compact_keys(text).expect("Compact keys must scan");
        let scanned: Vec<[String; 5]> = scanned
            .into_iter()
            .map(|key| key.map(String::from))
            .collect();
        assert_eq!(Some(scanned), serde(text));
    }
    for text in [
        " []",
        r#"[["A","1","","ALA"]]"#,
        r#"[["A","1","","ALA","protein","x"]]"#,
        r#"[["A\"","1","","ALA","protein"]]"#,
        r#"[["A","1","","ALA","protein"],]"#,
        r#"[["A","1","","ALA","protein"] ]"#,
        "[[\"A\t\",\"1\",\"\",\"ALA\",\"protein\"]]",
        r#"[["A","1","","ALA","protein"]"#,
    ] {
        assert!(observations::compact_keys(text).is_none(), "{text}");
    }
}

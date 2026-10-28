use super::*;
use crate::observations::{
    residue::Packed,
    schema::{Residue, State},
    source::Atom,
};

fn atom() -> Atom {
    Atom {
        source_row: 8,
        occupancy: 1.0,
        name: "N".into(),
        element: "N".into(),
        chain: "A".into(),
        component: "ALA".into(),
        entity: "1".into(),
        auth_sequence: "93".into(),
        auth_chain: "A".into(),
        auth_component: "ALA".into(),
        insertion: "?".into(),
        alternate: ".".into(),
        sequence: 7,
        position: [1., 2., 3.],
        protein: true,
        group: 0,
    }
}

fn schema() -> Schema {
    let mut names = vec![String::new(); 32];
    for (slot, name) in [(0, "N"), (4, "CB"), (15, "H"), (16, "H2"), (17, "HB1")] {
        names[slot] = name.into();
    }
    let mut parents = vec![-1; 32];
    parents[15] = 0;
    parents[16] = 0;
    parents[17] = 4;
    Schema {
        slots: 32,
        hydrogen_start: 15,
        residues: vec![Residue {
            residue: "ALA".into(),
            names,
            parents,
            states: vec![State {
                present: vec![true; 32],
                n_terminal: true,
                c_terminal: false,
                linked: false,
            }],
        }],
    }
}

fn tokens() -> Tokens {
    let mut packed = Packed::empty([1, 0, 0]);
    packed.position[0] = [1., 2., 3.];
    packed.mask[0] = true;
    packed.mask[15] = true;
    packed.presence.fill(true);
    packed.state_known = true;
    let mut tokens = Tokens::empty();
    tokens.push(&atom(), 0, &packed, false, [0; 4]);
    let mut other = atom();
    other.sequence = 8;
    tokens.push(&other, 0, &packed, false, [0; 4]);
    tokens
}

fn uncertainty(hydrogens_unknown: bool) -> TargetUncertainty {
    TargetUncertainty {
        parent_key: ["A".into(), "7".into(), "?".into(), "ALA".into(), "N".into()],
        hydrogens_unknown,
        chemical_unknown: true,
        reasons: Vec::new(),
    }
}

#[test]
fn partial_protein_keeps_observed_h_and_other_parent_targets() {
    let mut tokens = tokens();
    let original = tokens.positions.clone();
    apply(&mut tokens, &schema(), &[uncertainty(true)]).unwrap();
    assert_eq!(tokens.positions, original);
    assert_eq!(tokens.masks[15], 1);
    assert_eq!(tokens.presence[15..18], [1, 0, 1]);
    assert_eq!(tokens.topology, [-1, 0, 0, 1, 0, 0]);
    assert_eq!(tokens.state_known, [0, 1]);
    assert_eq!(tokens.presence[32..], [1; 32]);
}

#[test]
fn chemical_only_report_keeps_coordinate_targets() {
    let mut tokens = tokens();
    let original = tokens.masks.clone();
    apply(&mut tokens, &schema(), &[uncertainty(false)]).unwrap();
    assert_eq!(tokens.masks, original);
    assert_eq!(tokens.presence[15], 1);
    assert_eq!(tokens.state_known, [0, 1]);
}

#[test]
fn ligand_partial_report_clears_only_that_parents_label_flags() {
    let mut tokens = tokens();
    tokens.ligand = vec![1, 1];
    tokens.keys[0][4] = "C1".into();
    tokens.charge_known.fill(1);
    tokens.hydrogen_known.fill(1);
    let mut value = uncertainty(true);
    value.parent_key[4] = "C1".into();
    let original = tokens.masks.clone();
    apply(&mut tokens, &schema(), &[value]).unwrap();
    assert_eq!(tokens.charge_known, [0, 1]);
    assert_eq!(tokens.hydrogen_known, [0, 1]);
    assert_eq!(tokens.masks, original);
}

#[test]
fn absent_parent_and_malformed_status_do_not_silently_mask() {
    let mut value = uncertainty(true);
    value.parent_key[0] = "missing".into();
    assert!(apply(&mut tokens(), &schema(), &[value]).is_err());
    assert!(
        serde_json::from_str::<TargetUncertainty>(
            r#"{"parent_key":[],"hydrogens_unknown":true,"chemical_unknown":true}"#
        )
        .is_err()
    );
    let mut json = serde_json::to_value(uncertainty(true)).unwrap();
    json["unknown_field"] = serde_json::Value::Bool(true);
    assert!(serde_json::from_value::<TargetUncertainty>(json).is_err());
}

use super::*;
use crate::observations::{residue::Packed, source};

fn input(kind: &str, rows: &str) -> Source {
    let raw = format!("data_gap\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n1 '{kind}'\n\
        loop_\n_atom_site.label_entity_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n\
        _atom_site.auth_seq_id\n_atom_site.label_comp_id\n_atom_site.label_atom_id\n\
        _atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
        _atom_site.pdbx_PDB_model_num\n{rows}");
    source::read(&crate::cif::parse(raw.as_bytes()).expect("Fixture must parse"), 1)
        .expect("Fixture must read")
}

fn packed(source: &Source, protein: bool) -> Tokens {
    let mut tokens = Tokens::empty();
    let mut chains = Vec::new();
    let mut mapping = vec![usize::MAX; source.atoms.len()];
    for rows in &source.groups {
        for &row in rows {
            mapping[row] = tokens.identity.len();
            let mut packed = Packed::empty([-1; 3]);
            if protein {
                packed.mask[1] = true;
                packed.position[1] = source.atoms[row].position;
            }
            tokens.push(&source.atoms[row], 0, &packed, !protein, [0; 4]);
            let chain = &source.atoms[row].chain;
            if !chains.contains(chain) { chains.push(chain.clone()); }
            let token = tokens.identity.len() - 1;
            tokens.chains[token] = chains.iter().position(|name| name == chain)
                .expect("Fixture chain must exist") as i64;
            if protein { break; }
        }
    }
    apply(source, &mapping, &mut tokens);
    super::super::record::renumber_ligands(&mut tokens);
    retained_proteins(&mut tokens);
    tokens
}

#[test]
fn protein_gaps_accumulate_without_changing_source_keys_or_existing_gaps() {
    let source = input("polypeptide(L)", concat!(
        "1 A 1 101 ALA CA C 0 0 0 1\n",
        "1 A 2 102 ALA CA C 4 0 0 1\n",
        "1 A 3 103 ALA CA C 8.01 0 0 1\n",
        "1 A 4 104 ALA CA C 11.8 0 0 1\n",
        "1 A 8 108 ALA CA C 30 0 0 1\n",
        "1 A 9 109 ALA CA C 40 0 0 1\n",
        "1 B 1 101 ALA CA C 90 0 0 1\n",
        "1 B 2 102 ALA CA C 93.8 0 0 1\n"));
    let mut tokens = packed(&source, true);
    assert_eq!(tokens.sequence, [1, 2, 67, 68, 72, 137, 1, 2]);
    assert_eq!(tokens.keys[2][1], "3");
    tokens.select(&[1, 3]).expect("Crop must retain index gaps");
    assert_eq!(tokens.sequence, [2, 68]);
}

#[test]
fn nucleotide_gaps_survive_alignment_for_rna_dna_and_prime_aliases() {
    for kind in ["polyribonucleotide", "polydeoxyribonucleotide", "polydeoxyribonucleotide/polyribonucleotide hybrid"] {
        for oxygen in ["O3'", "O3*"] {
            let source = input(kind, &format!(
                "1 A 1 101 A P P 0 0 0 1\n1 A 1 101 A \"{oxygen}\" O 2 0 0 1\n\
                 1 A 2 102 A P P 3.6 0 0 1\n1 A 2 102 A \"{oxygen}\" O 6 0 0 1\n\
                 1 A 3 103 A P P 9 0 0 1\n1 A 3 103 A \"{oxygen}\" O 11 0 0 1\n\
                 1 A 5 105 A P P 12.6 0 0 1\n"));
            let mut tokens = packed(&source, false);
            assert_eq!(tokens.sequence, [0, 1, 2, 3, 68, 69, 70]);
            let expected = serde_json::to_string(&tokens.keys).expect("Keys must encode");
            tokens.align(&expected).expect("Alignment must preserve polymer breaks");
            assert_eq!(tokens.sequence, [0, 1, 2, 3, 68, 69, 70]);
            assert_eq!(tokens.keys[4][1], "103");
        }
    }
}

#[test]
fn missing_nucleotide_link_atoms_break_adjacency_but_nonpolymers_do_not() {
    let rows = "1 A 1 1 A P P 0 0 0 1\n1 A 2 2 A P P 1.6 0 0 1\n";
    assert_eq!(packed(&input("polyribonucleotide", rows), false).sequence, [0, 65]);
    assert_eq!(packed(&input("other", rows), false).sequence, [0, 1]);
}

#[test]
fn a_larger_deposited_gap_stays_larger_than_the_clipping_range() {
    let source = input("polypeptide(L)",
        "1 A 1 1 ALA CA C 0 0 0 1\n1 A 100 100 ALA CA C 90 0 0 1\n");
    assert_eq!(packed(&source, true).sequence, [1, 100]);
}

#[test]
fn existing_small_gaps_stay_exactly_the_same_size() {
    let source = input("polypeptide(L)",
        "1 A 61 245 ALA CA C 0 0 0 1\n1 A 76 260 ALA CA C 28 0 0 1\n");
    assert_eq!(packed(&source, true).sequence, [61, 76]);
}

#[test]
fn completeness_keeps_existing_gaps_and_only_widens_new_consecutive_breaks() {
    let source = input("polypeptide(L)", concat!(
        "1 A 1 1 ALA CA C 0 0 0 1\n1 A 2 2 ALA CA C 3.8 0 0 1\n",
        "1 A 3 3 ALA CA C 7.6 0 0 1\n1 A 4 4 ALA CA C 11.4 0 0 1\n",
        "1 A 5 5 ALA CA C 20 0 0 1\n"));
    let mut tokens = packed(&source, true);
    for row in 0..5 {
        tokens.masks[row * 32 + 1] = 1;
        tokens.masks[row * 32 + 4] = u8::from(row != 1);
        tokens.positions[row * 96 + 3..][..3].copy_from_slice(&source.atoms[row].position);
    }
    let mut names = vec![String::new(); 32];
    names[1] = "CA".into();
    names[4] = "CB".into();
    let schema = crate::observations::Schema { slots: 32, hydrogen_start: 15,
        residues: vec![crate::observations::schema::Residue { residue: "ALA".into(), names,
            parents: vec![-1; 32], states: Vec::new() }] };
    tokens.complete(&schema).expect("Complete-residue selection must preserve source identities");
    assert_eq!(tokens.sequence, [1, 3, 4, 69]);
    assert_eq!(tokens.keys[1][1], "3");
    tokens.complete(&schema).expect("Repeated completeness checks must preserve index gaps");
    assert_eq!(tokens.sequence, [1, 3, 4, 69]);
}

#[test]
fn distant_ca_atoms_do_not_establish_internal_peptide_topology() {
    let mut source = input("polypeptide(L)", concat!(
        "1 A 1 101 ALA N N -1 0 0 1\n1 A 1 101 ALA CA C 0 0 0 1\n",
        "1 A 1 101 ALA C C 1.5 0 0 1\n1 A 2 102 ALA N N 2.8 0 0 1\n",
        "1 A 2 102 ALA CA C 5 0 0 1\n1 A 2 102 ALA C C 6.5 0 0 1\n"));
    source.bonds.push((2, 3));
    let pdb = crate::observations::topology::pdb(&source);
    let protein = crate::observations::topology::protein(&source, "observed_connectivity")
        .expect("Topology must resolve without false links");
    for values in [pdb, protein] {
        assert_eq!(values[0][1], -1);
        assert_eq!(values[1][0], -1);
    }
}

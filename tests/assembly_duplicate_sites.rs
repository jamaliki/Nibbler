//! Partially occupied copies of one unit that overlap are one site stored twice.
//!
//! 9CNT places inositol hexakisphosphate (IHP E, occupancy 0.17) beside a five-fold
//! axis. Operator 5 puts a copy on the identity copy: 17 of its 36 atoms lie within
//! 1 A of an atom of the first, below the half that the whole-ligand rule needs. The
//! two copies are alternatives of one site and only one belongs in the assembly.

#![allow(clippy::expect_used)]

use _core::cif::parse;
use _core::pdbx::{assembly_document_with_mapping, build_model};

const OPERATORS: &str = "\
1 'identity operation' 1_555 x,y,z 1.00000000 0.00000000 0.00000000 0.00000 0.00000000 1.00000000 0.00000000 0.00000 0.00000000 0.00000000 1.00000000 0.00000
2 'point symmetry operation' ? ? 0.30900000 0.95110000 0.00000000 -58.74760 -0.95110000 0.30900000 0.00000000 370.91670 0.00000000 0.00000000 1.00000000 0.00000
3 'point symmetry operation' ? ? -0.80900000 0.58780000 0.00000000 275.86090 -0.58780000 -0.80900000 0.00000000 541.40900 0.00000000 0.00000000 1.00000000 0.00000
4 'point symmetry operation' ? ? -0.80900000 -0.58780000 0.00000000 541.40870 0.58780000 -0.80900000 0.00000000 275.86260 0.00000000 0.00000000 1.00000000 0.00000
5 'point symmetry operation' ? ? 0.30900000 -0.95110000 0.00000000 370.91820 0.95110000 0.30900000 0.00000000 -58.74690 0.00000000 0.00000000 1.00000000 0.00000
";

/// 9CNT IHP E, auth A 301: 36 atoms, all at occupancy 0.17.
const IHP_9CNT: &str = r##"
HETATM 6945 C C1 . IHP E 2 . ? 226.242 224.975 141.911 0.17 20.00 301 IHP A C1 1
HETATM 6946 C C2 . IHP E 2 . ? 224.802 225.390 141.610 0.17 20.00 301 IHP A C2 1
HETATM 6947 C C3 . IHP E 2 . ? 224.538 226.855 141.912 0.17 20.00 301 IHP A C3 1
HETATM 6948 C C4 . IHP E 2 . ? 225.543 227.679 141.141 0.17 20.00 301 IHP A C4 1
HETATM 6949 C C5 . IHP E 2 . ? 226.829 227.456 141.896 0.17 20.00 301 IHP A C5 1
HETATM 6950 C C6 . IHP E 2 . ? 227.304 226.036 141.620 0.17 20.00 301 IHP A C6 1
HETATM 6951 O O11 . IHP E 2 . ? 226.522 223.867 141.108 0.17 20.00 301 IHP A O11 1
HETATM 6952 P P1 . IHP E 2 . ? 226.874 222.418 141.805 0.17 20.00 301 IHP A P1 1
HETATM 6953 O O21 . IHP E 2 . ? 228.278 222.414 142.336 0.17 20.00 301 IHP A O21 1
HETATM 6954 O O31 . IHP E 2 . ? 225.921 222.118 142.932 0.17 20.00 301 IHP A O31 1
HETATM 6955 O O41 . IHP E 2 . ? 226.753 221.341 140.764 0.17 20.00 301 IHP A O41 1
HETATM 6956 O O12 . IHP E 2 . ? 223.971 224.672 142.461 0.17 20.00 301 IHP A O12 1
HETATM 6957 P P2 . IHP E 2 . ? 222.745 223.791 141.831 0.17 20.00 301 IHP A P2 1
HETATM 6958 O O22 . IHP E 2 . ? 223.301 222.702 140.955 0.17 20.00 301 IHP A O22 1
HETATM 6959 O O32 . IHP E 2 . ? 221.925 223.204 142.947 0.17 20.00 301 IHP A O32 1
HETATM 6960 O O42 . IHP E 2 . ? 221.882 224.689 140.993 0.17 20.00 301 IHP A O42 1
HETATM 6961 O O13 . IHP E 2 . ? 223.234 227.199 141.548 0.17 20.00 301 IHP A O13 1
HETATM 6962 P P3 . IHP E 2 . ? 222.571 228.573 142.164 0.17 20.00 301 IHP A P3 1
HETATM 6963 O O23 . IHP E 2 . ? 223.166 229.788 141.513 0.17 20.00 301 IHP A O23 1
HETATM 6964 O O33 . IHP E 2 . ? 221.082 228.564 141.946 0.17 20.00 301 IHP A O33 1
HETATM 6965 O O43 . IHP E 2 . ? 222.862 228.669 143.640 0.17 20.00 301 IHP A O43 1
HETATM 6966 O O14 . IHP E 2 . ? 225.677 227.252 139.821 0.17 20.00 301 IHP A O14 1
HETATM 6967 P P4 . IHP E 2 . ? 225.107 228.245 138.638 0.17 20.00 301 IHP A P4 1
HETATM 6968 O O24 . IHP E 2 . ? 223.636 228.489 138.877 0.17 20.00 301 IHP A O24 1
HETATM 6969 O O34 . IHP E 2 . ? 225.869 229.551 138.677 0.17 20.00 301 IHP A O34 1
HETATM 6970 O O44 . IHP E 2 . ? 225.302 227.589 137.293 0.17 20.00 301 IHP A O44 1
HETATM 6971 O O15 . IHP E 2 . ? 227.749 228.410 141.464 0.17 20.00 301 IHP A O15 1
HETATM 6972 P P5 . IHP E 2 . ? 227.782 229.852 142.256 0.17 20.00 301 IHP A P5 1
HETATM 6973 O O25 . IHP E 2 . ? 226.472 230.580 142.144 0.17 20.00 301 IHP A O25 1
HETATM 6974 O O35 . IHP E 2 . ? 228.100 229.618 143.712 0.17 20.00 301 IHP A O35 1
HETATM 6975 O O45 . IHP E 2 . ? 228.848 230.716 141.643 0.17 20.00 301 IHP A O45 1
HETATM 6976 O O16 . IHP E 2 . ? 228.340 225.766 142.498 0.17 20.00 301 IHP A O16 1
HETATM 6977 P P6 . IHP E 2 . ? 229.883 225.700 141.956 0.17 20.00 301 IHP A P6 1
HETATM 6978 O O26 . IHP E 2 . ? 229.999 224.575 140.966 0.17 20.00 301 IHP A O26 1
HETATM 6979 O O36 . IHP E 2 . ? 230.779 225.426 143.133 0.17 20.00 301 IHP A O36 1
HETATM 6980 O O46 . IHP E 2 . ? 230.283 226.991 141.301 0.17 20.00 301 IHP A O46 1
"##;

const ATOM_TAGS: [&str; 20] = [
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

fn source(occupancy: &str) -> String {
    let tags: String = ATOM_TAGS
        .iter()
        .map(|tag| format!("_atom_site.{tag}\n"))
        .collect();
    let rows = IHP_9CNT.trim().replace(" 0.17 ", &format!(" {occupancy} "));
    format!(
        "data_9CNT\n_entry.id 9CNT\nloop_\n_entity.id\n_entity.type\n2 non-polymer\n\
        loop_\n_struct_asym.id\n_struct_asym.entity_id\nE 2\n\
        loop_\n_pdbx_struct_assembly.id\n_pdbx_struct_assembly.details\n1 'complete point assembly'\n\
        loop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n\
        _pdbx_struct_assembly_gen.asym_id_list\n1 '(1-5)' E\n\
        loop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.type\n_pdbx_struct_oper_list.name\n\
        _pdbx_struct_oper_list.symmetry_operation\n\
        _pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n\
        _pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n\
        _pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n\
        _pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n\
        _pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n\
        _pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n{OPERATORS}\
        loop_\n{tags}{rows}\n"
    )
}

fn chains(text: &str) -> (Vec<String>, Vec<[usize; 3]>) {
    let document = parse(text.as_bytes()).expect("source");
    let model = build_model(&document).expect("model");
    let (expanded, mapping) = assembly_document_with_mapping(&model, None).expect("assembly");
    let written = _core::cif::write_canonical(&expanded).expect("expanded text");
    let mut names: Vec<String> = written
        .lines()
        .filter(|line| line.starts_with("HETATM "))
        .filter_map(|line| line.split_whitespace().nth(6))
        .map(str::to_owned)
        .collect();
    names.dedup();
    (names, mapping)
}

#[test]
fn a_partially_occupied_ligand_copy_on_its_own_site_is_dropped() {
    let (names, mapping) = chains(&source("0.17"));
    assert_eq!(names, ["E"]);
    let fifth: Vec<_> = mapping
        .iter()
        .filter(|&&[_, removed, _]| removed == 5)
        .collect();
    assert_eq!(fifth.len(), 36);
    assert!(fifth.iter().all(|&&[_, _, kept]| kept == 1));
}

#[test]
fn the_same_ligand_at_full_occupancy_keeps_that_copy() {
    // The whole-ligand rule already drops copies 2 to 4, which overlap by half or
    // more. Copy 5 overlaps by less, so only occupancy marks it as an alternative.
    let (names, mapping) = chains(&source("1.00"));
    assert_eq!(names, ["E", "E-5"]);
    assert!(mapping.iter().all(|&[_, removed, _]| removed != 5));
}

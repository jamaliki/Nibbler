//! Source-backed PDBx semantic-model and strict-writer tests for Phase 4.

#![allow(clippy::expect_used)]

use _core::cif::{parse, write_canonical};
use _core::pdbx::{
    ComponentResolution, EntityKind, ProfileSeverity, build_component_registry, build_model,
    build_model_with_registry, canonical_document, validate_document, validate_model,
};

const CHEMISTRY: &[u8] = include_bytes!("fixtures/chemistry/ligand_ion_water.cif");
const GLYCAN: &[u8] = include_bytes!("fixtures/chemistry/branched_glycan.cif");

#[test]
fn chemistry_fixture_builds_without_losing_semantic_identity() {
    let document = parse(CHEMISTRY).expect("fixture must parse");
    let model = build_model(&document).expect("fixture must build");

    assert_eq!(model.entry_id(), "NIBBLER_LIGAND_ION_WATER");
    assert_eq!(model.entities().len(), 4);
    assert_eq!(model.asym_units().len(), 4);
    assert_eq!(model.atom_sites().len(), 6);
    assert_eq!(model.connections().len(), 1);
    assert_eq!(
        model
            .entities()
            .iter()
            .map(|entity| entity.kind())
            .collect::<Vec<_>>(),
        vec![
            EntityKind::Polymer,
            EntityKind::NonPolymer,
            EntityKind::NonPolymer,
            EntityKind::Water,
        ]
    );

    let polymer = &model.entities()[0];
    assert_eq!(polymer.sequence()[1].component_id(), "MSE");
    let ligand = &model.atom_sites()[2];
    assert_eq!(ligand.label_seq_id(), None);
    assert_eq!(ligand.author_ids().0, "-1");
    assert_eq!(model.atom_sites()[4].formal_charge(), Some(2));
    assert_eq!(model.entities()[3].component_id(), Some("HOH"));
    let connection = &model.connections()[0];
    assert_eq!(connection.identity(), ("ZN1", "metalc"));
    assert_eq!(connection.endpoints().0.asym_id(), "B");
    assert_eq!(connection.endpoints().1.asym_id(), "C");
}

#[test]
fn chemistry_fixture_passes_dictionary_and_semantic_profile() {
    let document = parse(CHEMISTRY).expect("fixture must parse");
    let report = validate_document(&document).expect("dictionary must load");

    assert_eq!(report.dictionary_version(), "5.416");
    assert!(report.coverage().contains(&"chemical-component-resolution"));
    assert!(report.is_valid(), "unexpected: {:?}", report.diagnostics());
}

#[test]
fn semantic_diagnostic_limit_reports_truncation() {
    let duplicate_atom = "ATOM 1 C CA . ALA A 1 1 ? 0.000 0.000 0.000 1.00 10.00 ? 10 ALA X CA 1\n";
    let source = String::from_utf8(CHEMISTRY.to_vec())
        .expect("fixture is UTF-8")
        .replace(
            "\n_struct_conn_type.id metalc",
            &format!(
                "\n{}\n_struct_conn_type.id metalc",
                duplicate_atom.repeat(10_001)
            ),
        );
    let document = parse(source.as_bytes()).expect("large duplicate fixture must parse");
    let model = build_model(&document).expect("duplicate atom identifiers remain representable");
    let report = validate_model(&model);

    assert_eq!(report.diagnostics().len(), 10_001);
    assert!(
        report.diagnostics()[..10_000]
            .iter()
            .all(|diagnostic| diagnostic.code() == "PDBX_ATOM_ID_DUPLICATE")
    );
    let truncated = &report.diagnostics()[10_000];
    assert_eq!(truncated.code(), "PDBX_DIAGNOSTICS_TRUNCATED");
    assert_eq!(truncated.severity(), ProfileSeverity::Warning);
}

#[test]
fn glycan_fixture_maps_atom_sites_through_the_branch_scheme() {
    let document = parse(GLYCAN).expect("glycan fixture must parse");
    let model = build_model(&document).expect("glycan fixture must build");
    let entity = &model.entities()[0];

    assert_eq!(entity.kind(), EntityKind::Branched);
    assert_eq!(entity.branch_nodes().len(), 2);
    assert_eq!(entity.branch_links().len(), 1);
    assert_eq!(model.atom_sites()[0].label_seq_id(), None);
    assert!(
        validate_document(&document)
            .expect("dictionary must load")
            .is_valid()
    );
}

#[test]
fn strict_profile_rejects_an_unresolved_component() {
    let source = String::from_utf8(CHEMISTRY.to_vec())
        .expect("fixture is UTF-8")
        .replace(
            "ATP non-polymer . \"ADENOSINE-5'-TRIPHOSPHATE\" 'C10 H16 N5 O13 P3' 507.181\n",
            "",
        );
    let document = parse(source.as_bytes()).expect("modified fixture must parse");
    let model = build_model(&document).expect("unresolved chemistry remains representable");
    let atp = model
        .components()
        .iter()
        .find(|component| component.id() == "ATP")
        .expect("referenced ATP placeholder");

    assert_eq!(atp.resolution(), ComponentResolution::Unresolved);
    assert!(
        validate_model(&model)
            .diagnostics()
            .iter()
            .any(|diagnostic| diagnostic.code() == "PDBX_COMPONENT_UNRESOLVED")
    );
}

#[test]
fn caller_selected_local_ccd_resolves_without_network_or_guessing() {
    let source = String::from_utf8(CHEMISTRY.to_vec())
        .expect("fixture is UTF-8")
        .replace(
            "ATP non-polymer . \"ADENOSINE-5'-TRIPHOSPHATE\" 'C10 H16 N5 O13 P3' 507.181\n",
            "",
        );
    let document = parse(source.as_bytes()).expect("modified fixture must parse");
    let ccd = parse(
        br#"data_ATP
_chem_comp.id ATP
_chem_comp.type non-polymer
_chem_comp.name "ADENOSINE-5'-TRIPHOSPHATE"
_chem_comp.formula 'C10 H16 N5 O13 P3'
_chem_comp.formula_weight 507.181
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
ATP PG P
"#,
    )
    .expect("local CCD fixture must parse");
    let registry = build_component_registry(&ccd).expect("local CCD must load");
    let model = build_model_with_registry(&document, Some(&registry)).expect("ATP must resolve");
    let atp = model
        .components()
        .iter()
        .find(|component| component.id() == "ATP")
        .expect("ATP definition");

    assert_eq!(registry.len(), 1);
    assert_eq!(atp.resolution(), ComponentResolution::LocalCcd);
    assert!(validate_model(&model).is_valid());
    let emitted = canonical_document(&model);
    assert!(
        validate_document(&emitted)
            .expect("dictionary must load")
            .is_valid()
    );
}

#[test]
fn local_ccd_conflicts_do_not_hide_behind_source_priority() {
    let document = parse(CHEMISTRY).expect("fixture must parse");
    let ccd = parse(b"data_ATP\n_chem_comp.id ATP\n_chem_comp.name 'INCOMPATIBLE NAME'\n")
        .expect("local CCD fixture must parse");
    let registry = build_component_registry(&ccd).expect("local CCD must load");
    let error = build_model_with_registry(&document, Some(&registry))
        .expect_err("incompatible definitions must fail");

    assert_eq!(error.code(), "PDBX_COMPONENT_CONFLICT");
}

#[test]
fn profile_catches_atom_entity_mismatch_without_using_group_pdb() {
    let source = String::from_utf8(CHEMISTRY.to_vec())
        .expect("fixture is UTF-8")
        .replace("HETATM 5 Zn ZN . ZN C 3 . ?", "ATOM 5 Zn ZN . ZN C 2 . ?");
    let document = parse(source.as_bytes()).expect("modified fixture must parse");
    let model = build_model(&document).expect("structural model still builds");
    let report = validate_model(&model);
    let codes = report
        .diagnostics()
        .iter()
        .map(|diagnostic| diagnostic.code())
        .collect::<Vec<_>>();

    assert!(codes.contains(&"PDBX_ATOM_ENTITY"));
}

#[test]
fn canonical_profile_output_uses_shared_serializer_and_reparses() {
    let document = parse(CHEMISTRY).expect("fixture must parse");
    let model = build_model(&document).expect("fixture must build");
    let ordered = canonical_document(&model);
    let output = write_canonical(&ordered).expect("profile document must serialize");
    let reparsed = parse(output.as_bytes()).expect("profile output must reparse");

    assert_eq!(reparsed, ordered);
    assert!(
        output.find("_entity.id").expect("entity category")
            < output.find("_chem_comp.id").expect("component category")
    );
    assert!(
        output.find("_atom_site.group_PDB").expect("group column")
            < output
                .find("_atom_site.Cartn_x")
                .expect("coordinate column")
    );
}

#[test]
fn connection_endpoints_must_resolve_to_emitted_atom_sites() {
    let source = String::from_utf8(CHEMISTRY.to_vec())
        .expect("fixture is UTF-8")
        .replace(
            "ZN1 metalc ATP B . PG ZN C . ZN",
            "ZN1 metalc ATP B . PG ZN C . ABSENT",
        );
    let document = parse(source.as_bytes()).expect("fixture with connection must parse");
    let model = build_model(&document).expect("connection must build");

    assert!(
        validate_model(&model)
            .diagnostics()
            .iter()
            .any(|diagnostic| diagnostic.code() == "PDBX_CONNECTION_ENDPOINT")
    );
}

#[test]
fn branched_entity_retains_graph_and_requires_instance_scheme() {
    let source = br#"data_glycan
_entry.id GLYCAN
loop_
_audit_conform.dict_name
_audit_conform.dict_version
mmcif_pdbx.dic 5.416
loop_
_entity.id
_entity.type
1 branched
loop_
_pdbx_entity_branch_list.entity_id
_pdbx_entity_branch_list.num
_pdbx_entity_branch_list.comp_id
1 1 NAG
1 2 MAN
loop_
_pdbx_entity_branch_link.entity_id
_pdbx_entity_branch_link.entity_branch_list_num_1
_pdbx_entity_branch_link.entity_branch_list_num_2
_pdbx_entity_branch_link.atom_id_1
_pdbx_entity_branch_link.atom_id_2
_pdbx_entity_branch_link.value_order
1 1 2 C1 O4 sing
loop_
_chem_comp.id
_chem_comp.type
NAG saccharide
MAN saccharide
_struct_asym.id G
_struct_asym.entity_id 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1 C C1 . NAG G 1 1 ? 0 0 0 1 10 ? 1 NAG G C1 1
HETATM 2 O O4 . MAN G 1 2 ? 1 0 0 1 10 ? 2 MAN G O4 1
"#;
    let document = parse(source).expect("glycan fixture must parse");
    let model = build_model(&document).expect("glycan fixture must build");
    let entity = &model.entities()[0];

    assert_eq!(entity.kind(), EntityKind::Branched);
    assert_eq!(entity.branch_nodes().len(), 2);
    assert_eq!(entity.branch_links().len(), 1);
    assert!(
        validate_model(&model)
            .diagnostics()
            .iter()
            .any(|diagnostic| diagnostic.code() == "PDBX_BRANCH_SCHEME")
    );
}

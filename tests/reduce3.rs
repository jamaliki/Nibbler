//! Reduce3 on parsed documents: the in-memory path must give exactly what the
//! `reduce3` program reads from and writes to mmCIF files.

#![cfg(feature = "reduce3")]
#![allow(clippy::expect_used)]

use std::path::Path;

use _core::cif::{CifDocument, CifEntry, decode_binary, encode_binary, parse, write_canonical};
use _core::reduce::{
    BlockSource, OptParams, Params, ReduceErrorCode, load_monomer_library, run, structure_document,
};

/// The first five residues of PDB entry 1CRN (crambin), heavy atoms only.
const CRAMBIN: &str = include_str!("reduce3/crambin_1-5.cif");

const CHEMISTRY: &str = include_str!("fixtures/chemistry/ligand_ion_water.cif");
const GLYCAN: &str = include_str!("fixtures/chemistry/branched_glycan.cif");
const PREDICTION: &str = include_str!("fixtures/modelcif/prediction_with_qa.cif");

fn models() -> [&'static str; 4] {
    [CRAMBIN, CHEMISTRY, GLYCAN, PREDICTION]
}

/// The model as the `reduce3` program reads and writes it.
fn file_round_trip(text: &str) -> String {
    let structure = reduce3::mmcif::read_mmcif(text).expect("reduce3 must read the model");
    reduce3::mmcif::write_mmcif(&structure)
}

fn entries(document: &CifDocument) -> &[CifEntry] {
    document.blocks().first().expect("one block").entries()
}

#[test]
fn parsed_blocks_build_the_same_model_as_files() {
    for text in models() {
        let document = parse(text.as_bytes()).expect("fixture must parse");
        let block = document.blocks().first().expect("one block");
        let structure = reduce3::mmcif::structure_from_cif(&BlockSource::new(block))
            .expect("the block must be readable");
        assert_eq!(
            reduce3::mmcif::write_mmcif(&structure),
            file_round_trip(text)
        );
    }
}

#[test]
fn binary_cif_documents_build_the_same_model() {
    let document = parse(CRAMBIN.as_bytes()).expect("fixture must parse");
    let decoded = decode_binary(&encode_binary(&document).expect("encodable")).expect("decodable");
    let block = decoded.blocks().first().expect("one block");
    let structure = reduce3::mmcif::structure_from_cif(&BlockSource::new(block))
        .expect("the block must be readable");
    assert_eq!(
        reduce3::mmcif::write_mmcif(&structure),
        file_round_trip(CRAMBIN)
    );
}

#[test]
fn output_documents_equal_a_parse_of_the_written_file() {
    for text in models() {
        let structure = reduce3::mmcif::read_mmcif(text).expect("reduce3 must read the model");
        let document = structure_document(&structure, "check").expect("representable");
        let written = parse(reduce3::mmcif::write_mmcif(&structure).as_bytes())
            .expect("reduce3 output must parse");
        assert_eq!(entries(&document), entries(&written));
        assert_eq!(
            document.blocks().first().and_then(|b| b.code()),
            Some("check")
        );

        let canonical = write_canonical(&document).expect("canonical output");
        let reparsed = parse(canonical.as_bytes()).expect("canonical output must parse");
        assert_eq!(entries(&reparsed), entries(&document));
        encode_binary(&document).expect("BinaryCIF output");
    }
}

#[test]
fn a_named_library_directory_must_exist() {
    let error = load_monomer_library(Some(Path::new("/nonexistent/chem_data")))
        .err()
        .expect("not a chem_data directory");
    assert_eq!(error.code(), ReduceErrorCode::MonomerLibraryNotFound);
}

/// Drop the timing lines of a report.
fn without_timings(report: &str) -> String {
    report
        .lines()
        .filter(|line| !line.contains("Time to"))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn hydrogens_are_added_in_memory() {
    let Ok(monomers) = load_monomer_library(None) else {
        eprintln!("skipped: chem_data not found (set REDUCE3_CHEM_DATA)");
        return;
    };
    let document = parse(CRAMBIN.as_bytes()).expect("fixture must parse");
    for (compat, flips) in [(false, false), (false, true), (true, false), (true, true)] {
        let params = Params {
            compat,
            opt: OptParams {
                add_flip_movers: flips,
                ..OptParams::default()
            },
            ..Params::default()
        };
        let reduction = run(&document, &monomers, &params).expect("Reduce3 must succeed");

        let expected = reduce3::pipeline::run(
            reduce3::mmcif::read_mmcif(CRAMBIN).expect("reduce3 must read the model"),
            &monomers,
            &params,
        )
        .expect("Reduce3 must succeed");
        let written = parse(reduce3::mmcif::write_mmcif(&expected.structure).as_bytes())
            .expect("reduce3 output must parse");
        assert_eq!(entries(reduction.document()), entries(&written));
        assert_eq!(
            without_timings(reduction.report()),
            without_timings(&expected.description)
        );
        assert_eq!(
            reduction.document().blocks().first().and_then(|b| b.code()),
            Some("crn5")
        );
        let atoms = entries(reduction.document())
            .iter()
            .find_map(|entry| match entry {
                CifEntry::Loop(cif_loop) if cif_loop.tags()[0].starts_with("_atom_site.") => {
                    Some(cif_loop.row_count())
                }
                _ => None,
            });
        assert_eq!(atoms, Some(66), "33 heavy atoms and 33 hydrogens");
    }

    let no_model = parse(b"data_empty\n_entry.id empty\n").expect("must parse");
    let error = run(&no_model, &monomers, &Params::default()).expect_err("no model");
    assert_eq!(error.code(), ReduceErrorCode::InvalidModel);
}

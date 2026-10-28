//! Hydrogen removal filters atom sites and their tensors and copies the rest.

use _core::cif::{parse, write_canonical};
use _core::pdbx::hydrogen_free_document;

const SOURCE: &[u8] = b"data_x\n_entry.id x\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n\
ATOM 1 N N\nATOM 2 H H\nATOM 3 C CA\nHETATM 4 d D1\nHETATM 5 T T1\nATOM 6 ? X\n\
loop_\n_atom_site_anisotrop.id\n_atom_site_anisotrop.type_symbol\n_atom_site_anisotrop.U[1][1]\n\
1 N 0.1\n2 H 0.2\n3 C 0.3\n9 C 0.4\n\
loop_\n_struct_conn.id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_atom_id\n\
c1 H CA\n";

fn canonical(raw: &[u8]) -> Result<String, Box<dyn std::error::Error>> {
    Ok(write_canonical(&parse(raw)?)?)
}

#[test]
fn removes_hydrogen_isotopes_and_orphan_tensors() -> Result<(), Box<dyn std::error::Error>> {
    let result = write_canonical(&hydrogen_free_document(&parse(SOURCE)?)?)?;
    let expected = canonical(
        b"data_x\n_entry.id x\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n\
ATOM 1 N N\nATOM 3 C CA\nATOM 6 ? X\n\
loop_\n_atom_site_anisotrop.id\n_atom_site_anisotrop.type_symbol\n_atom_site_anisotrop.U[1][1]\n\
1 N 0.1\n3 C 0.3\n\
loop_\n_struct_conn.id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_atom_id\n\
c1 H CA\n",
    )?;
    assert_eq!(result, expected);
    Ok(())
}

#[test]
fn keeps_single_rows_and_drops_emptied_tables() -> Result<(), Box<dyn std::error::Error>> {
    let raw = b"data_x\n\
loop_\n_atom_site.id\n_atom_site.type_symbol\n1 C\n2 H\n\
loop_\n_atom_site_anisotrop.id\n_atom_site_anisotrop.U[1][1]\n2 0.2\n";
    let result = write_canonical(&hydrogen_free_document(&parse(raw)?)?)?;
    assert_eq!(
        result,
        canonical(b"data_x\nloop_\n_atom_site.id\n_atom_site.type_symbol\n1 C\n")?
    );

    let scalar = b"data_x\n_atom_site.id 1\n_atom_site.type_symbol C\n\
_atom_site_anisotrop.id 7\n_atom_site_anisotrop.U[1][1] 0.2\n";
    let result = write_canonical(&hydrogen_free_document(&parse(scalar)?)?)?;
    assert_eq!(
        result,
        canonical(b"data_x\n_atom_site.id 1\n_atom_site.type_symbol C\n")?
    );
    Ok(())
}

#[test]
fn unchanged_documents_share_their_blocks() -> Result<(), Box<dyn std::error::Error>> {
    let raw = b"data_x\nloop_\n_atom_site.id\n_atom_site.type_symbol\n1 C\n2 O\n\
loop_\n_atom_site_anisotrop.id\n_atom_site_anisotrop.U[1][1]\n2 0.2\n";
    let source = parse(raw)?;
    let result = hydrogen_free_document(&source)?;
    assert!(std::ptr::eq(result.blocks(), source.blocks()));
    Ok(())
}

#[test]
fn rejects_documents_it_cannot_filter() -> Result<(), Box<dyn std::error::Error>> {
    for raw in [
        b"data_a\nloop_\n_atom_site.id\n_atom_site.type_symbol\n1 C\ndata_b\n_entry.id b\n"
            .as_slice(),
        b"data_a\n_entry.id a\n".as_slice(),
        b"data_a\nloop_\n_atom_site.id\n_atom_site.label_atom_id\n1 CA\n".as_slice(),
    ] {
        let result = hydrogen_free_document(&parse(raw)?);
        assert!(matches!(result, Err(ref error) if error.code() == "PDBX_HYDROGEN_REMOVAL"));
    }
    Ok(())
}

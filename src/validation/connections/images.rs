use reduce3::cell::{SymOp, crystal_operators};
use reduce3::cifsource::{CifSource, CifTable};
use reduce3::model::Structure;

use crate::reduce::BlockSource;

pub(super) fn relative(
    source: &BlockSource<'_>,
    connections: &impl CifTable,
    row: usize,
    structure: &Structure,
) -> Result<SymOp, String> {
    let first = endpoint(source, connections, row, 1, structure)?;
    let second = endpoint(source, connections, row, 2, structure)?;
    Ok(first
        .inverse()
        .ok_or("Connection symmetry operator is singular")?
        .compose(&second))
}

fn endpoint(
    source: &BlockSource<'_>,
    connections: &impl CifTable,
    row: usize,
    side: usize,
    structure: &Structure,
) -> Result<SymOp, String> {
    let code_value = connections
        .column(&format!("ptnr{side}_symmetry"))
        .map(|column| connections.cell(row, column))
        .unwrap_or(".".into());
    let value = code_value.trim();
    if value == "." {
        return SymOp::parse("x,y,z").ok_or_else(|| "Identity symmetry operator is invalid".into());
    }
    let (id, code) = value
        .split_once('_')
        .ok_or_else(|| format!("Invalid connection symmetry code {value}"))?;
    if id.is_empty()
        || !id.bytes().all(|c| c.is_ascii_digit())
        || code.len() != 3
        || !code.bytes().all(|c| c.is_ascii_digit())
    {
        return Err(format!("Invalid connection symmetry code {value}"));
    }
    let mut image = declared(source, id)?
        .or_else(|| {
            if id == "1" && code == "555" && crystal_operators(structure).is_none() {
                return SymOp::parse("x,y,z");
            }
            let index = id.parse::<usize>().ok()?.checked_sub(1)?;
            crystal_operators(structure)?.1.get(index).copied()
        })
        .ok_or_else(|| format!("Connection symmetry operator {id} is undefined"))?;
    for (translation, digit) in image.t.iter_mut().zip(code.bytes()) {
        *translation += (digit as i32 - b'5' as i32) as f64;
    }
    if !image.is_identity() && crystal_operators(structure).is_none() {
        return Err("Connection image needs a valid cell and space group".into());
    }
    if image.inverse().is_none() {
        return Err("Connection symmetry operator is singular".into());
    }
    Ok(image)
}

fn declared(source: &BlockSource<'_>, id: &str) -> Result<Option<SymOp>, String> {
    for (category, field) in [
        ("space_group_symop", "operation_xyz"),
        ("symmetry_equiv", "pos_as_xyz"),
    ] {
        let Some(table) = source.table(category) else {
            continue;
        };
        let (Some(identity), Some(operation)) = (table.column("id"), table.column(field)) else {
            return Err(format!("Incomplete {category} operator table"));
        };
        let rows: Vec<_> = (0..table.row_count())
            .filter(|&row| table.cell(row, identity).trim() == id)
            .collect();
        if rows.len() != 1 {
            return Err(format!(
                "Connection symmetry operator {id} is missing or ambiguous"
            ));
        }
        let image = SymOp::parse(&table.cell(rows[0], operation))
            .ok_or("Invalid declared symmetry operation")?;
        return Ok(Some(image));
    }
    Ok(None)
}

pub(super) fn same(first: SymOp, second: SymOp) -> bool {
    first
        .r
        .into_iter()
        .flatten()
        .chain(first.t)
        .zip(second.r.into_iter().flatten().chain(second.t))
        .all(|(a, b)| (a - b).abs() < 1e-10)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cif::parse;
    use reduce3::model::CrystalSymmetry;

    #[test]
    fn relative_images_preserve_cell_shift_and_declared_ids() -> Result<(), String> {
        let document = parse(b"data_test\nloop_\n_space_group_symop.id\n_space_group_symop.operation_xyz\n7 x,y,z\n9 -x,-y,z\nloop_\n_struct_conn.ptnr1_symmetry\n_struct_conn.ptnr2_symmetry\n7_655 9_555\n").map_err(|e|e.to_string())?;
        let source = BlockSource::new(&document.blocks()[0]);
        let structure = Structure {
            crystal: Some(CrystalSymmetry {
                cell: [20., 30., 40., 90., 90., 90.],
                space_group: "P 1".into(),
            }),
            ..Default::default()
        };
        let image = relative(
            &source,
            &source.table("struct_conn").ok_or("missing connections")?,
            0,
            &structure,
        )?;
        assert!(same(
            image,
            SymOp::parse("-x-1,-y,z").ok_or("invalid expected image")?
        ));
        Ok(())
    }

    #[test]
    fn unknown_and_malformed_images_cannot_become_identity() -> Result<(), String> {
        for code in ["?", "0_555", "bogus", "1_55", "1_55X", "2_555"] {
            let text = format!(
                "data_test\n_struct_conn.ptnr1_symmetry 1_555\n_struct_conn.ptnr2_symmetry {code}\n"
            );
            let document = parse(text.as_bytes()).map_err(|e| e.to_string())?;
            let source = BlockSource::new(&document.blocks()[0]);
            assert!(
                relative(
                    &source,
                    &source.table("struct_conn").ok_or("missing connections")?,
                    0,
                    &Structure::default()
                )
                .is_err()
            );
        }
        Ok(())
    }

    #[test]
    fn dummy_cells_keep_identity_connections() -> Result<(), String> {
        let document = parse(
            b"data_test\n_struct_conn.ptnr1_symmetry 1_555\n_struct_conn.ptnr2_symmetry 1_555\n",
        )
        .map_err(|error| error.to_string())?;
        let source = BlockSource::new(&document.blocks()[0]);
        let structure = Structure {
            crystal: Some(CrystalSymmetry {
                cell: [1., 1., 1., 90., 90., 90.],
                space_group: "P 1".into(),
            }),
            ..Default::default()
        };
        let image = relative(
            &source,
            &source.table("struct_conn").ok_or("missing connections")?,
            0,
            &structure,
        )?;
        assert!(image.is_identity());
        Ok(())
    }
}

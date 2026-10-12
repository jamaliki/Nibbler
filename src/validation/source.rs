use reduce3::cifsource::{CifSource, CifTable};

use crate::cif::{CifBlock, CifDocument, CifEntry};
use crate::reduce::BlockSource;

pub(super) fn validate_charges(block: &CifBlock) -> Result<(), String> {
    let source = BlockSource::new(block);
    let Some(table) = source.table("atom_site") else {
        return Ok(());
    };
    let Some(column) = table.column("pdbx_formal_charge") else {
        return Ok(());
    };
    for row in 0..table.row_count() {
        let charge = table.cell(row, column).parse::<i32>();
        if table.missing(row, column).is_none() && (charge.is_err() || charge == Ok(i32::MIN)) {
            return Err(format!("Atom row {row} has an invalid formal charge"));
        }
    }
    Ok(())
}

// Reduce3's preparation hierarchy omits these interpretation inputs. Check the
// original document before validation so they cannot disappear from coverage.
pub(super) fn unsupported(block: &CifBlock) -> Vec<String> {
    let mut reasons = Vec::new();
    if block
        .entries()
        .iter()
        .any(|entry| matches!(entry, CifEntry::Frame(_)))
    {
        reasons.push("Save-frame restraint inputs need reference interpretation".into());
    }
    if let Some(table) = BlockSource::new(block).table("struct_conn") {
        for name in ["ptnr1_symmetry", "ptnr2_symmetry"] {
            let Some(column) = table.column(name) else {
                continue;
            };
            for row in 0..table.row_count() {
                let value = table.cell(row, column);
                if !matches!(value.trim(), "1_555" | ".") {
                    reasons.push(format!(
                        "Connection row {row} has unsupported {name}={value}"
                    ));
                }
            }
        }
    }
    reasons
}

pub(super) fn has_restraints(document: &CifDocument) -> bool {
    document.blocks().iter().any(|data_block| {
        let declared = |tag: &str| {
            let tag = tag.to_ascii_lowercase();
            tag.starts_with("_lib_")
                || tag.starts_with("_chem_link")
                || tag.starts_with("_chem_mod")
        };
        if data_block.entries().iter().any(|entry| match entry {
            CifEntry::Item(item) => declared(item.tag()),
            CifEntry::Loop(table) => table.tags().iter().any(|tag| declared(tag)),
            CifEntry::Frame(_) => false,
        }) {
            return true;
        }
        let source = BlockSource::new(data_block);
        [
            ("lib_atom", "type"),
            ("chem_comp_bond", "value_dist"),
            ("chem_comp_bond", "value_dist_esd"),
            ("chem_comp_atom", "type_energy"),
            ("chem_comp_angle", "value_angle"),
            ("chem_comp_chir", "volume_sign"),
            ("chem_comp_plane_atom", "dist_esd"),
            ("chem_comp_tor", "value_angle"),
            ("chem_link_bond", "value_dist"),
            ("chem_link_angle", "value_angle"),
            ("chem_link_chir", "volume_sign"),
            ("chem_link_plane", "dist_esd"),
            ("chem_link_tor", "value_angle"),
            ("chem_mod_atom", "function"),
            ("chem_mod_bond", "function"),
            ("chem_mod_angle", "function"),
            ("chem_mod_chir", "function"),
            ("chem_mod_plane_atom", "function"),
            ("chem_mod_tor", "function"),
        ]
        .iter()
        .any(|(name, target)| {
            source
                .table(name)
                .is_some_and(|t| t.column(target).is_some())
        })
    })
}

#[cfg(test)]
mod tests {
    use super::{has_restraints, unsupported, validate_charges};
    use crate::cif::parse;

    #[test]
    fn connection_symmetry_cannot_disappear_into_the_preparation_model()
    -> Result<(), Box<dyn std::error::Error>> {
        let input = b"data_test\nloop_\n_struct_conn.ptnr1_symmetry\n_struct_conn.ptnr2_symmetry\n1_555 2_555\n";
        let document = parse(input)?;
        let reasons = unsupported(&document.blocks()[0]);
        assert_eq!(reasons.len(), 1);
        assert!(reasons[0].contains("ptnr2_symmetry=2_555"));
        Ok(())
    }

    #[test]
    fn unknown_energy_and_link_categories_cannot_disappear()
    -> Result<(), Box<dyn std::error::Error>> {
        for source in [
            b"data_test\nloop_\n_lib_bond.type\nC\n".as_slice(),
            b"data_test\nloop_\n_chem_link_plane_atom.dist_esd\n.02\n",
        ] {
            assert!(has_restraints(&parse(source)?));
        }
        Ok(())
    }

    #[test]
    fn an_invalid_formal_charge_cannot_become_neutral() -> Result<(), Box<dyn std::error::Error>> {
        let document = parse(b"data_test\nloop_\n_atom_site.pdbx_formal_charge\nbogus\n")?;
        assert!(validate_charges(&document.blocks()[0]).is_err());
        let document = parse(b"data_test\nloop_\n_atom_site.pdbx_formal_charge\n?\n.\n1\n-1\n")?;
        assert!(validate_charges(&document.blocks()[0]).is_ok());
        Ok(())
    }

    #[test]
    fn ccd_bond_order_metadata_is_not_a_custom_restraint() -> Result<(), Box<dyn std::error::Error>>
    {
        let document = parse(b"data_test\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nALA N CA SING\n")?;
        assert!(!has_restraints(&document));
        let document = parse(
            b"data_test\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.value_dist\nALA 1.5\n",
        )?;
        assert!(has_restraints(&document));
        Ok(())
    }
}

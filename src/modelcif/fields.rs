//! The error codes of ModelCIF item decoding at the generic-CIF trust boundary.

use crate::pdbx::fields::Decoder;

/// ModelCIF items: `MODELCIF_ITEM_REQUIRED` and `MODELCIF_ITEM_TYPE`.
pub(super) const MODELCIF: Decoder = Decoder::new("MODELCIF_ITEM_REQUIRED", "MODELCIF_ITEM_TYPE");

#[cfg(test)]
mod tests {
    use crate::cif::{CategoryView, parse};
    use crate::pdbx::fields::{PDBX, category_rows};

    use super::MODELCIF;

    #[test]
    fn each_profile_reports_its_own_codes_with_shared_messages() -> Result<(), String> {
        let document = parse(b"data_demo\n_ma_demo.count x\n_ma_demo.name ?\n")
            .map_err(|error| error.message().to_owned())?;
        let view = CategoryView::new(document.blocks()[0].entries());
        let row = category_rows(view, "ma_demo", ["count", "name"])
            .next()
            .ok_or("no row")?;
        for (decoder, prefix) in [(PDBX, "PDBX"), (MODELCIF, "MODELCIF")] {
            let typed = decoder
                .required_integer(&row, "ma_demo", "count", 0)
                .err()
                .ok_or("integer decoded")?;
            assert_eq!(typed.code(), format!("{prefix}_ITEM_TYPE"));
            assert_eq!(
                typed.message(),
                "value \"x\" for _ma_demo.count is not a integer"
            );
            assert_eq!(
                typed.context(),
                ["category=ma_demo", "item=_ma_demo.count", "row=1"]
            );
            let missing = decoder
                .required_text(&row, "ma_demo", "name", 1)
                .err()
                .ok_or("text decoded")?;
            assert_eq!(missing.code(), format!("{prefix}_ITEM_REQUIRED"));
            assert_eq!(
                missing.message(),
                "required item _ma_demo.name is absent or missing"
            );
            let float = decoder
                .optional_float(&row, "ma_demo", "count", 0)
                .err()
                .ok_or("float decoded")?;
            assert_eq!(
                float.message(),
                "value \"x\" for _ma_demo.count is not a finite float"
            );
        }
        Ok(())
    }
}

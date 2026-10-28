//! Ligand token identities owned by the caller.
//!
//! The packer derives no ligand token. The caller names the token of each element it
//! keeps, one token for every other element, and the token of hydrogen, which only a
//! bridging metal hydride carries as a token of its own. Residue tokens take 0-20.

use std::collections::{BTreeMap, HashMap};

/// The lowest token a ligand atom may carry; 0-20 are residues.
const FIRST_LIGAND_TOKEN: i64 = 21;

#[derive(Clone, Debug)]
pub(crate) struct TokenTable {
    elements: HashMap<String, i64>,
    other: i64,
    hydrogen: i64,
}

impl TokenTable {
    /// Check the table: upper-case element symbols other than hydrogen, and distinct tokens
    /// above the residue range for every element, for other elements and for hydrogen.
    pub(crate) fn new(
        elements: impl IntoIterator<Item = (String, i64)>,
        other: i64,
        hydrogen: i64,
    ) -> Result<Self, String> {
        let mut table = HashMap::new();
        for (element, token) in elements {
            if element.is_empty() || !element.bytes().all(|byte| byte.is_ascii_uppercase()) {
                return Err(format!(
                    "Ligand token table element {element:?} is not an upper-case symbol"
                ));
            }
            if matches!(element.as_str(), "H" | "D" | "T") {
                return Err(
                    "Hydrogen takes the table's hydrogen token, not an element entry".into(),
                );
            }
            if table.insert(element.clone(), token).is_some() {
                return Err(format!("Ligand token table lists {element} twice"));
            }
        }
        let mut owners = BTreeMap::new();
        let named = table
            .iter()
            .map(|(element, &token)| (element.as_str(), token));
        for (owner, token) in named.chain([("other elements", other), ("hydrogen", hydrogen)]) {
            if token < FIRST_LIGAND_TOKEN {
                return Err(format!(
                    "Ligand token {token} of {owner} is not above the residue tokens 0-20"
                ));
            }
            if let Some(previous) = owners.insert(token, owner) {
                let (first, second) = if previous < owner {
                    (previous, owner)
                } else {
                    (owner, previous)
                };
                return Err(format!(
                    "Ligand token {token} is given to both {first} and {second}"
                ));
            }
        }
        Ok(Self {
            elements: table,
            other,
            hydrogen,
        })
    }

    /// The token of a heavy ligand atom with this element symbol, or the token of other
    /// elements. The symbol is matched as written.
    pub(crate) fn element(&self, element: &str) -> i64 {
        self.elements.get(element).copied().unwrap_or(self.other)
    }

    /// The token of a bridging hydride hydrogen.
    pub(crate) fn hydrogen(&self) -> i64 {
        self.hydrogen
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use super::TokenTable;

    fn table(elements: &[(&str, i64)], other: i64, hydrogen: i64) -> Result<TokenTable, String> {
        TokenTable::new(
            elements
                .iter()
                .map(|&(element, token)| (element.to_owned(), token)),
            other,
            hydrogen,
        )
    }

    #[test]
    fn named_elements_take_their_tokens_and_others_share_one() {
        let table = table(&[("C", 21), ("SI", 66)], 64, 65).expect("Valid table");
        assert_eq!(table.element("C"), 21);
        assert_eq!(table.element("SI"), 66);
        assert_eq!(table.element("ZR"), 64);
        assert_eq!(table.element("Si"), 64);
        assert_eq!(table.hydrogen(), 65);
    }

    #[test]
    fn tables_with_shared_or_residue_tokens_are_refused() {
        for (elements, other, hydrogen, message) in [
            (&[("C", 21), ("O", 21)][..], 64, 65, "given to both C and O"),
            (
                &[("C", 64)][..],
                64,
                65,
                "given to both C and other elements",
            ),
            (&[("C", 65)][..], 64, 65, "given to both C and hydrogen"),
            (
                &[("C", 21)][..],
                65,
                65,
                "given to both hydrogen and other elements",
            ),
            (&[("C", 20)][..], 64, 65, "Ligand token 20 of C"),
            (&[("C", 21)][..], 0, 65, "Ligand token 0 of other elements"),
            (&[("C", 21)][..], 64, -1, "Ligand token -1 of hydrogen"),
            (&[("H", 66)][..], 64, 65, "Hydrogen takes"),
            (&[("D", 66)][..], 64, 65, "Hydrogen takes"),
            (&[("Fe", 66)][..], 64, 65, "not an upper-case symbol"),
            (&[("", 66)][..], 64, 65, "not an upper-case symbol"),
        ] {
            let error = table(elements, other, hydrogen).expect_err(message);
            assert!(error.contains(message), "{error}");
        }
    }
}

//! Caller-owned atom slots and allowed residue states.

use serde::Deserialize;

#[derive(Clone, Deserialize)]
pub(super) struct State {
    pub present: Vec<bool>,
    pub n_terminal: bool,
    pub c_terminal: bool,
    pub linked: bool,
}

impl State {
    pub fn topology(&self) -> [i64; 3] {
        [
            i64::from(self.n_terminal),
            i64::from(self.c_terminal),
            i64::from(self.linked),
        ]
    }
}

#[derive(Clone, Deserialize)]
pub(super) struct Residue {
    pub residue: String,
    pub names: Vec<String>,
    pub parents: Vec<i64>,
    pub states: Vec<State>,
}

#[derive(Clone, Deserialize)]
pub(crate) struct Schema {
    pub(super) slots: usize,
    pub(super) hydrogen_start: usize,
    pub(super) residues: Vec<Residue>,
}

impl Schema {
    pub fn read(json: &str) -> Result<Self, String> {
        let value: Self = serde_json::from_str(json).map_err(|error| error.to_string())?;
        if value.slots != 32 || value.hydrogen_start != 15 || value.residues.len() != 20 {
            return Err(
                "Observation packing requires 20 residues with 32 slots and hydrogen_start 15"
                    .into(),
            );
        }
        for residue in &value.residues {
            if residue.names.len() != value.slots
                || residue.parents.len() != value.slots
                || residue.states.is_empty()
                || residue
                    .parents
                    .iter()
                    .any(|&parent| !(-1..15).contains(&parent))
                || residue
                    .states
                    .iter()
                    .any(|state| state.present.len() != value.slots)
            {
                return Err(format!(
                    "Invalid observation schema for {}",
                    residue.residue
                ));
            }
        }
        Ok(value)
    }

    pub fn identity(&self, name: &str) -> usize {
        self.residues
            .iter()
            .position(|residue| residue.residue == name)
            .unwrap_or(20)
    }
}

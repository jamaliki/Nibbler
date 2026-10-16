//! Explicit linked-product preparation. Ordinary model construction stays lossless.

mod input;
pub(crate) mod library;
mod links;
mod product;
mod rearrange;
mod valence;

use super::{ComponentRegistry, PdbxModel, SemanticError};
use serde::Serialize;

/// An atom in the product graph, including coordinate-missing template sites.
#[derive(Clone, Debug, Serialize)]
pub struct ProductAtom {
    /// Deposited model number.
    pub model_number: i64,
    /// Chain instance in the supplied document.
    pub asym_id: String,
    /// Deposited sequence index, when supplied.
    pub sequence_id: Option<i64>,
    /// Author residue number.
    pub auth_sequence_id: String,
    /// Deposited insertion code, when supplied.
    pub insertion_code: Option<String>,
    /// Deposited component identifier.
    pub component_id: String,
    /// Component atom identifier.
    pub atom_id: String,
    /// Atom identifier in the component template before name mapping.
    pub template_atom_id: Option<String>,
    /// Selected alternate identifier.
    pub alternate_id: Option<String>,
    /// Zero-based source atom_site row; None for an unobserved template site.
    pub source_row: Option<usize>,
    /// Supplied coordinates; None for a missing template site.
    pub coordinates: Option<[f64; 3]>,
    /// Element from the component definition or source.
    pub element: String,
    /// Product formal charge, or None without a known definition.
    pub formal_charge: Option<i32>,
    /// Formal charge before product preparation.
    pub template_charge: Option<i32>,
    /// Charge explicitly supplied in atom_site, separate from template chemistry.
    pub deposited_charge: Option<i32>,
    /// Expected attached hydrogen count, or None without a known definition.
    pub hydrogen_count: Option<usize>,
    /// False for displaced atoms, which keep their source identity.
    pub active: bool,
    #[serde(skip)]
    residue: usize,
    #[serde(skip)]
    leaving: bool,
    #[serde(skip)]
    impacted: bool,
}

impl ProductAtom {
    fn is_hydrogen(&self) -> bool {
        matches!(self.element.as_str(), "H" | "D" | "T")
    }
}

/// A bond in the prepared graph. Order zero denotes metal coordination.
#[derive(Clone, Debug, Serialize)]
pub struct ProductBond {
    /// Endpoints index the product atom list, including inactive atoms.
    pub atom_indices: [usize; 2],
    /// Covalent order 1, 2 or 3; zero for coordination.
    pub order: u8,
    /// Component bond order before preparation, when applicable.
    pub template_order: Option<u8>,
    /// Component, deposited struct_conn identifier, or polymer sequence evidence.
    pub source: BondSource,
    /// Deposited connection identifiers supporting this bond.
    pub connection_ids: Vec<String>,
    /// Bond order explicitly supplied by a connection, when known.
    pub declared_order: Option<u8>,
    /// False for removed bonds.
    pub active: bool,
}

/// The operation that first inserted a product bond.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum BondSource {
    /// The chemical component definition supplied the bond.
    Component,
    /// A deposited connection supplied the bond.
    StructConn,
    /// Consecutive sequence sites supplied the bond.
    PolymerSequence,
    /// Product formation transferred this hydrogen.
    HydrogenTransfer,
}

/// A chemical edit, separate from the deposited source document.
#[derive(Clone, Debug, Serialize)]
pub struct ChemicalEdit {
    /// Affected product atom indices.
    pub atom_indices: Vec<usize>,
    /// Stable repair rule name.
    pub reason: &'static str,
}

/// An immutable product graph. Preparation never changes its source model.
#[derive(Clone, Debug, Serialize)]
pub struct PreparedChemistry {
    schema_version: &'static str,
    atoms: Vec<ProductAtom>,
    bonds: Vec<ProductBond>,
    edits: Vec<ChemicalEdit>,
    diagnostics: Vec<String>,
    #[serde(skip)]
    adjacency: Vec<Vec<usize>>,
    #[serde(skip)]
    residues: Vec<input::Residue>,
}

impl PreparedChemistry {
    /// Version of the graph and chemical repair rules.
    pub fn schema_version(&self) -> &'static str {
        self.schema_version
    }
    /// Atoms retain stable indices after leaving-group removal.
    pub fn atoms(&self) -> &[ProductAtom] {
        &self.atoms
    }
    /// Bonds include provenance and pre-repair orders.
    pub fn bonds(&self) -> &[ProductBond] {
        &self.bonds
    }
    /// Explicit chemical changes made during preparation.
    pub fn edits(&self) -> &[ChemicalEdit] {
        &self.edits
    }
    /// Records excluded links and long authored glycosylation links.
    pub fn diagnostics(&self) -> &[String] {
        &self.diagnostics
    }

    fn error(&self, atom: usize, message: &str) -> SemanticError {
        let a = &self.atoms[atom];
        SemanticError::new(
            "PDBX_LINK_PRODUCT",
            format!(
                "{message} at {}/{}/{}/{}",
                a.asym_id, a.auth_sequence_id, a.component_id, a.atom_id
            ),
            vec!["atom_site".into(), a.atom_id.clone()],
        )
    }

    fn neighbors(&self, atom: usize) -> impl Iterator<Item = (usize, u8)> + '_ {
        self.adjacency[atom].iter().filter_map(move |&k| {
            let b = &self.bonds[k];
            (b.active && b.order > 0).then(|| {
                (
                    if b.atom_indices[0] == atom {
                        b.atom_indices[1]
                    } else {
                        b.atom_indices[0]
                    },
                    b.order,
                )
            })
        })
    }

    fn degree(&self, atom: usize) -> usize {
        self.neighbors(atom)
            .map(|(_, order)| usize::from(order))
            .sum()
    }

    fn hydrogens(&self, atom: usize) -> Vec<usize> {
        let mut h: Vec<_> = self
            .neighbors(atom)
            .filter_map(|(other, _)| self.atoms[other].is_hydrogen().then_some(other))
            .collect();
        h.sort_by_key(|&a| self.atoms[a].source_row.is_some());
        h
    }

    fn bond_index(&self, a: usize, b: usize) -> Option<usize> {
        self.adjacency[a]
            .iter()
            .copied()
            .find(|&k| self.bonds[k].active && self.bonds[k].atom_indices.contains(&b))
    }

    fn add_bond(
        &mut self,
        a: usize,
        b: usize,
        order: u8,
        source: BondSource,
    ) -> Result<(), SemanticError> {
        if a == b {
            return Err(self.error(a, "Self bond"));
        }
        if let Some(k) = self.bond_index(a, b) {
            return if self.bonds[k].order == order {
                Ok(())
            } else {
                Err(self.error(a, "Conflicting bond orders"))
            };
        }
        let k = self.bonds.len();
        self.bonds.push(ProductBond {
            atom_indices: [a, b],
            order,
            template_order: (source == BondSource::Component).then_some(order),
            source,
            connection_ids: Vec::new(),
            declared_order: None,
            active: true,
        });
        self.adjacency[a].push(k);
        self.adjacency[b].push(k);
        Ok(())
    }

    fn change_order(&mut self, a: usize, b: usize, order: u8, reason: &'static str) {
        if let Some(k) = self.bond_index(a, b) {
            self.bonds[k].order = order;
            self.atoms[a].impacted = true;
            self.atoms[b].impacted = true;
            self.edits.push(ChemicalEdit {
                atom_indices: vec![a, b],
                reason,
            });
        }
    }

    fn remove(&mut self, atoms: &[usize], reason: &'static str) {
        for &a in atoms {
            self.atoms[a].active = false;
            for &k in &self.adjacency[a] {
                self.bonds[k].active = false;
            }
        }
        self.edits.push(ChemicalEdit {
            atom_indices: atoms.to_vec(),
            reason,
        });
    }

    fn set_charge(&mut self, a: usize, charge: i32, reason: &'static str) {
        if self.atoms[a].formal_charge != Some(charge) {
            self.atoms[a].formal_charge = Some(charge);
            self.atoms[a].impacted = true;
            self.edits.push(ChemicalEdit {
                atom_indices: vec![a],
                reason,
            });
        }
    }
}

/// Resolve linked chemistry for one deposited model and one conformer per residue.
///
/// `alternate_id=None` selects the highest summed occupancy within each residue.
/// A supplied alternate is used where present; other residues keep their highest
/// occupancy conformer. This is not a clash-aware ensemble selection.
///
/// # Errors
///
/// Rejects ambiguous identities, unsupported orders, absent linked templates,
/// nonidentity symmetry links and products without a supported radical-free valence.
pub fn prepare_chemistry(
    model: &PdbxModel,
    registry: Option<&ComponentRegistry>,
    model_number: i64,
    alternate_id: Option<&str>,
) -> Result<PreparedChemistry, SemanticError> {
    prepare_product(model, registry, model_number, alternate_id, true)
}

pub(crate) fn prepare_product(model: &PdbxModel, registry: Option<&ComponentRegistry>,
    model_number: i64, alternate_id: Option<&str>, use_observed_hydrogens: bool,
) -> Result<PreparedChemistry, SemanticError> {
    let mut graph = input::build(model, registry, model_number, alternate_id, use_observed_hydrogens)?;
    links::add(model, &mut graph)?;
    product::resolve(&mut graph)?;
    Ok(graph)
}

//! Explicit linked-product preparation. Ordinary model construction stays lossless.

mod input;
pub(crate) mod library;
mod links;
mod product;
mod rearrange;
mod valence;

use super::{ComponentRegistry, PdbxModel, SemanticError, ShortText};
use serde::Serialize;

/// An atom in the product graph, including coordinate-missing template sites.
#[derive(Clone, Debug, Serialize)]
pub struct ProductAtom {
    /// Deposited model number.
    pub model_number: i64,
    /// Chain instance in the supplied document.
    pub asym_id: ShortText,
    /// Deposited sequence index, when supplied.
    pub sequence_id: Option<i64>,
    /// Author residue number.
    pub auth_sequence_id: ShortText,
    /// Deposited insertion code, when supplied.
    pub insertion_code: Option<ShortText>,
    /// Deposited component identifier.
    pub component_id: ShortText,
    /// Component atom identifier.
    pub atom_id: ShortText,
    /// Atom identifier in the component template before name mapping.
    pub template_atom_id: Option<ShortText>,
    /// Selected alternate identifier.
    pub alternate_id: Option<ShortText>,
    /// Zero-based source atom_site row; None for an unobserved template site.
    pub source_row: Option<usize>,
    /// Supplied coordinates; None for a missing template site.
    pub coordinates: Option<[f64; 3]>,
    /// Element from the component definition or source.
    pub element: ShortText,
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

/// Bond indices of one atom. Up to four are stored inline, which covers almost every atom
/// and saves one allocation per atom.
#[derive(Clone, Debug, Default)]
struct BondList {
    inline: [usize; 4],
    count: usize,
    spilled: Vec<usize>,
}

impl BondList {
    fn as_slice(&self) -> &[usize] {
        if self.spilled.is_empty() {
            &self.inline[..self.count]
        } else {
            &self.spilled
        }
    }

    fn iter(&self) -> std::slice::Iter<'_, usize> {
        self.as_slice().iter()
    }

    fn push(&mut self, bond: usize) {
        if !self.spilled.is_empty() {
            self.spilled.push(bond);
        } else if self.count < self.inline.len() {
            self.inline[self.count] = bond;
            self.count += 1;
        } else {
            self.spilled = self.inline.to_vec();
            self.spilled.push(bond);
        }
    }
}

impl<'a> IntoIterator for &'a BondList {
    type Item = &'a usize;
    type IntoIter = std::slice::Iter<'a, usize>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
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
    adjacency: Vec<BondList>,
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
            vec!["atom_site".into(), a.atom_id.as_str().to_owned()],
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

/// Map the product that [`prepare_product`] built from `model` without observed hydrogens
/// onto `other`, when `other` gives the same product.
///
/// The product reads the components, the entities, the `struct_conn` entries and the heavy
/// atom sites of `model_number`. Of their coordinates it reads only the distances between
/// the `CA`, `C`, `N`, `O3'` and `P` atoms of sequence links. So `other` gives the same
/// product when those inputs are equal, the residues of its heavy sites come in the same
/// order with the same sites, and the sites of those five atoms keep their coordinates.
/// Sites of one residue may come in another order when that order cannot change the
/// product (see [`order_free`]). Returns the row of `other` that holds each product atom,
/// or `None` when `other` may give another product.
pub(crate) fn product_rows_in(product: &PreparedChemistry, model: &PdbxModel, other: &PdbxModel,
    model_number: i64) -> Option<Vec<Option<usize>>> {
    if model.components != other.components || model.entities != other.entities
        || connection_entries(model) != connection_entries(other) {
        return None;
    }
    let sites = &other.atom_sites;
    let theirs: Vec<usize> = sites.iter().enumerate()
        .filter(|(_, site)| site.model_number == model_number && !hydrogen(&site.element))
        .map(|(row, _)| row)
        .collect();
    product_rows_with(product, model, model_number, &theirs, |row| &sites[row],
        |row| [sites[row].x, sites[row].y, sites[row].z])
}

/// Map the product as [`product_rows_in`] does, onto a document whose model gives the same
/// components, entities and `struct_conn` entries as `model`. `theirs` lists its heavy rows
/// of `model_number` in order; `identity` gives a site that holds the values of a row in
/// every field the product reads apart from coordinates, and `position` gives the row's
/// coordinates.
pub(crate) fn product_rows_with<'a>(product: &PreparedChemistry, model: &PdbxModel, model_number: i64,
    theirs: &[usize], identity: impl Fn(usize) -> &'a super::AtomSite, position: impl Fn(usize) -> [f64; 3],
) -> Option<Vec<Option<usize>>> {
    let own: Vec<usize> = model.atom_sites.iter().enumerate()
        .filter(|(_, site)| site.model_number == model_number && !hydrogen(&site.element))
        .map(|(row, _)| row)
        .collect();
    if own.len() != theirs.len() {
        return None;
    }
    let mut rows = vec![usize::MAX; model.atom_sites.len()];
    // Rows whose identities are the sites of `model` themselves need no residue check.
    if own.iter().zip(theirs).all(|(&row, &their_row)| std::ptr::eq(&model.atom_sites[row], identity(their_row))) {
        for (&row, &their_row) in own.iter().zip(theirs) {
            rows[row] = their_row;
        }
        return mapped_rows(product, model, &rows, position);
    }
    let mut templated = vec![false; model.atom_sites.len()];
    for atom in &product.atoms {
        if let (Some(row), Some(_)) = (atom.source_row, &atom.template_atom_id) {
            *templated.get_mut(row)? = true;
        }
    }
    let mut start = 0;
    while let Some(&first) = own.get(start) {
        let key = residue_key(&model.atom_sites[first]);
        let end = own[start..].iter().position(|&row| residue_key(&model.atom_sites[row]) != key)
            .map_or(own.len(), |length| start + length);
        let (mine, yours) = (&own[start..end], &theirs[start..end]);
        if mine.iter().zip(yours).all(|(&row, &their_row)| same_identity(&model.atom_sites[row], identity(their_row))) {
            for (&row, &their_row) in mine.iter().zip(yours) {
                rows[row] = their_row;
            }
        } else {
            if !order_free(&model.atom_sites, mine, &templated) {
                return None;
            }
            fn sorted<'s>(run: &[usize], site: impl Fn(usize) -> &'s super::AtomSite) -> Vec<usize> {
                let mut run = run.to_vec();
                run.sort_by(|&left, &right| site(left).label_atom_id.cmp(&site(right).label_atom_id));
                run
            }
            let sites = &model.atom_sites;
            for (row, their_row) in sorted(mine, |row| &sites[row]).into_iter().zip(sorted(yours, &identity)) {
                if !same_identity(&model.atom_sites[row], identity(their_row)) {
                    return None;
                }
                rows[row] = their_row;
            }
        }
        start = end;
    }
    mapped_rows(product, model, &rows, position)
}

/// The row of `rows` that holds each product atom, when every source row of the product has
/// one and the sites of the `CA`, `C`, `N`, `O3'` and `P` atoms keep their coordinates.
fn mapped_rows(product: &PreparedChemistry, model: &PdbxModel, rows: &[usize],
    position: impl Fn(usize) -> [f64; 3]) -> Option<Vec<Option<usize>>> {
    let mut mapped = Vec::with_capacity(product.atoms.len());
    for atom in &product.atoms {
        let row = match atom.source_row {
            Some(row) => Some(*rows.get(row).filter(|&&their_row| their_row != usize::MAX)?),
            None => None,
        };
        if let (Some(own), Some(their_row)) = (atom.source_row, row)
            && ["CA", "C", "N", "O3'", "P"].iter().any(|name| atom.atom_id == *name)
        {
            let site = &model.atom_sites[own];
            if [site.x, site.y, site.z].map(f64::to_bits) != position(their_row).map(f64::to_bits) {
                return None;
            }
        }
        mapped.push(row);
    }
    Some(mapped)
}

/// Whether `element` names hydrogen, as [`input::is_hydrogen_element`] decides, without
/// checking the text as UTF-8.
pub(crate) fn hydrogen(element: &ShortText) -> bool {
    ["H", "D", "T", "h", "d", "t"].iter().any(|name| element == name)
}

/// The fields that group atom sites into residues, as [`input::build`] reads them.
fn residue_key(site: &super::AtomSite) -> (&ShortText, Option<i64>, &ShortText, Option<&ShortText>) {
    (&site.label_asym_id, site.label_seq_id, &site.auth_seq_id, site.insertion_code.as_ref())
}

/// Whether the product reads the sites `run` of one residue the same in any order.
///
/// A residue built from a component template pairs its sites with template atoms by name,
/// and its first site gives the residue fields. So the order cannot matter when every site
/// was paired with a template atom, none names an alternate, whose occupancies are summed in
/// order, and all share the residue fields.
fn order_free(sites: &[super::AtomSite], run: &[usize], templated: &[bool]) -> bool {
    let Some(first) = run.first().map(|&row| &sites[row]) else {
        return true;
    };
    run.iter().all(|&row| {
        let site = &sites[row];
        templated[row]
            && site.label_alt_id.is_none()
            && site.label_comp_id == first.label_comp_id
            && site.label_entity_id == first.label_entity_id
            && site.auth_asym_id == first.auth_asym_id
            && site.auth_comp_id == first.auth_comp_id
    })
}

/// The block entries that hold `struct_conn` items, and any save frame.
fn connection_entries(model: &PdbxModel) -> Vec<&crate::cif::CifEntry> {
    use crate::cif::{CifEntry, in_category};
    model.document.blocks().first().map_or_else(Vec::new, |block| {
        block.entries().iter().filter(|entry| match entry {
            CifEntry::Item(item) => in_category(item.tag(), "struct_conn"),
            CifEntry::Loop(values) => values.tags().iter().any(|tag| in_category(tag, "struct_conn")),
            CifEntry::Frame(_) => true,
        }).collect()
    })
}

/// Whether two atom sites agree in every field the product reads except coordinates.
fn same_identity(site: &super::AtomSite, other: &super::AtomSite) -> bool {
    std::ptr::eq(site, other) || site.model_number == other.model_number
        && site.element == other.element
        && site.label_atom_id == other.label_atom_id
        && site.label_alt_id == other.label_alt_id
        && site.label_comp_id == other.label_comp_id
        && site.label_asym_id == other.label_asym_id
        && site.label_entity_id == other.label_entity_id
        && site.label_seq_id == other.label_seq_id
        && site.insertion_code == other.insertion_code
        && site.occupancy.to_bits() == other.occupancy.to_bits()
        && site.formal_charge == other.formal_charge
        && site.auth_seq_id == other.auth_seq_id
        && site.auth_comp_id == other.auth_comp_id
        && site.auth_asym_id == other.auth_asym_id
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use super::{PreparedChemistry, prepare_product, product_rows_in};
    use crate::cif::parse;
    use crate::pdbx::{PdbxModel, build_model};

    /// A Gly-Asn-Gly peptide whose Asn binds a zinc ion through `struct_conn`.
    const PEPTIDE: &str = "data_shared
_entry.id SHARED
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'GNG peptide'
2 non-polymer 'ZINC ION'
loop_
_entity_poly.entity_id
_entity_poly.type
1 polypeptide(L)
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1 GLY n
1 2 ASN n
1 3 GLY n
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
loop_
_chem_comp.id
_chem_comp.type
GLY 'L-peptide linking'
ASN 'L-peptide linking'
ZN non-polymer
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_leaving_atom_flag
GLY N N 0 N
GLY CA C 0 N
GLY C C 0 N
GLY O O 0 N
GLY OXT O 0 Y
GLY H H 0 N
GLY H2 H 0 Y
GLY HA2 H 0 N
GLY HA3 H 0 N
GLY HXT H 0 Y
ASN N N 0 N
ASN CA C 0 N
ASN C C 0 N
ASN O O 0 N
ASN CB C 0 N
ASN CG C 0 N
ASN OD1 O 0 N
ASN ND2 N 0 N
ASN OXT O 0 Y
ASN H H 0 N
ASN H2 H 0 Y
ASN HA H 0 N
ASN HB2 H 0 N
ASN HB3 H 0 N
ASN HD21 H 0 N
ASN HD22 H 0 N
ASN HXT H 0 Y
ZN ZN Zn 2 N
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
GLY N CA SING
GLY N H SING
GLY N H2 SING
GLY CA C SING
GLY CA HA2 SING
GLY CA HA3 SING
GLY C O DOUB
GLY C OXT SING
GLY OXT HXT SING
ASN N CA SING
ASN N H SING
ASN N H2 SING
ASN CA C SING
ASN CA CB SING
ASN CA HA SING
ASN C O DOUB
ASN C OXT SING
ASN CB CG SING
ASN CB HB2 SING
ASN CB HB3 SING
ASN CG OD1 DOUB
ASN CG ND2 SING
ASN ND2 HD21 SING
ASN ND2 HD22 SING
ASN OXT HXT SING
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
ATOM 1 N N . GLY A 1 1 ? -1.200 0.800 0.000 1.00 10.00 ? 1 GLY A N 1
ATOM 2 C CA . GLY A 1 1 ? 0.000 0.000 0.000 1.00 10.00 ? 1 GLY A CA 1
ATOM 3 C C . GLY A 1 1 ? 1.400 0.500 0.000 1.00 10.00 ? 1 GLY A C 1
ATOM 4 O O . GLY A 1 1 ? 1.600 1.700 0.000 1.00 10.00 ? 1 GLY A O 1
ATOM 5 N N . ASN A 1 2 ? 2.400 -0.300 0.000 1.00 10.00 ? 2 ASN A N 1
ATOM 6 C CA . ASN A 1 2 ? 3.800 0.000 0.000 1.00 10.00 ? 2 ASN A CA 1
ATOM 7 C C . ASN A 1 2 ? 4.900 0.900 0.000 1.00 10.00 ? 2 ASN A C 1
ATOM 8 O O . ASN A 1 2 ? 4.700 2.100 0.000 1.00 10.00 ? 2 ASN A O 1
ATOM 9 C CB . ASN A 1 2 ? 4.200 -1.400 0.400 1.00 10.00 ? 2 ASN A CB 1
ATOM 10 C CG . ASN A 1 2 ? 5.600 -1.800 0.100 1.00 10.00 ? 2 ASN A CG 1
ATOM 11 O OD1 . ASN A 1 2 ? 6.400 -1.000 -0.300 1.00 10.00 ? 2 ASN A OD1 1
ATOM 12 N ND2 . ASN A 1 2 ? 5.900 -3.100 0.200 1.00 10.00 ? 2 ASN A ND2 1
ATOM 13 N N . GLY A 1 3 ? 6.100 0.400 0.000 1.00 10.00 ? 3 GLY A N 1
ATOM 14 C CA . GLY A 1 3 ? 7.600 0.400 0.000 1.00 10.00 ? 3 GLY A CA 1
ATOM 15 C C . GLY A 1 3 ? 8.300 1.700 0.000 1.00 10.00 ? 3 GLY A C 1
ATOM 16 O O . GLY A 1 3 ? 7.700 2.800 0.000 1.00 10.00 ? 3 GLY A O 1
ATOM 17 O OXT . GLY A 1 3 ? 9.500 1.700 0.000 1.00 10.00 ? 3 GLY A OXT 1
HETATM 18 Zn ZN . ZN B 2 . ? 8.000 -1.500 -0.500 1.00 10.00 2 101 ZN A ZN 1
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_dist_value
metalc1 metalc A ASN 2 OD1 B ZN . ZN 2.10
";

    const OD1: &str = "OD1 . ASN A 1 2 ? 6.400 -1.000 -0.300";
    const ND2: &str = "ND2 . ASN A 1 2 ? 5.900 -3.100 0.200";

    fn model(text: &str) -> PdbxModel {
        build_model(&parse(text.as_bytes()).expect("fixture parses")).expect("fixture builds")
    }

    fn product(model: &PdbxModel) -> PreparedChemistry {
        prepare_product(model, None, 1, None, false).expect("fixture gives a product")
    }

    /// The rows of `edited` that hold the product of [`PEPTIDE`], after checking that a
    /// product built from `edited` differs only in its coordinates and source rows.
    fn shared_rows(edited: &str) -> Vec<Option<usize>> {
        let (original, edited) = (model(PEPTIDE), model(edited));
        let base = product(&original);
        let rows = product_rows_in(&base, &original, &edited, 1).expect("the product is shared");
        let mut fresh = product(&edited);
        assert_eq!(fresh.atoms.iter().map(|atom| atom.source_row).collect::<Vec<_>>(), rows);
        let mut base = base;
        for atom in base.atoms.iter_mut().chain(fresh.atoms.iter_mut()) {
            atom.coordinates = None;
            atom.source_row = None;
        }
        assert_eq!(format!("{base:?}"), format!("{fresh:?}"));
        rows
    }

    /// Check that `edited` gives another product and that [`product_rows_in`] sees it.
    fn not_shared(edited: &str) {
        let (original, edited) = (model(PEPTIDE), model(edited));
        let base = product(&original);
        let mut fresh = product(&edited);
        for atom in &mut fresh.atoms {
            atom.coordinates = None;
        }
        let mut stripped = base.clone();
        for atom in &mut stripped.atoms {
            atom.coordinates = None;
        }
        assert_ne!(format!("{stripped:?}"), format!("{fresh:?}"));
        assert!(product_rows_in(&base, &original, &edited, 1).is_none());
    }

    #[test]
    fn unchanged_document_shares_every_row() {
        let rows = shared_rows(PEPTIDE);
        let base = product(&model(PEPTIDE));
        assert_eq!(rows, base.atoms.iter().map(|atom| atom.source_row).collect::<Vec<_>>());
    }

    #[test]
    fn flipped_asparagine_shares_the_product() {
        let flipped = PEPTIDE
            .replace(OD1, "OD1 . ASN A 1 2 ? 5.900 -3.100 0.200")
            .replace(ND2, "ND2 . ASN A 1 2 ? 6.400 -1.000 -0.300");
        assert_ne!(flipped, PEPTIDE);
        shared_rows(&flipped);
    }

    #[test]
    fn reordered_residue_with_added_hydrogens_shares_the_product() {
        let site = |line: &str| PEPTIDE.lines().find(|row| row.contains(line)).expect("site").to_owned();
        let (od1, nd2) = (site(OD1), site(ND2));
        let hydrogen = "ATOM 19 H HD21 . ASN A 1 2 ? 6.700 -3.600 0.000 1.00 10.00 ? 2 ASN A HD21 1";
        let reordered = PEPTIDE
            .replace(&format!("{od1}\n{nd2}\n"), &format!("{nd2}\n{hydrogen}\n{od1}\n"));
        assert_ne!(reordered, PEPTIDE);
        let rows = shared_rows(&reordered);
        let base = product(&model(PEPTIDE));
        let row_of = |name: &str| {
            let atom = base.atoms.iter().position(|atom| atom.component_id == "ASN" && atom.atom_id == name);
            rows[atom.expect("atom")]
        };
        assert_eq!((row_of("ND2"), row_of("OD1")), (Some(10), Some(12)));
    }

    #[test]
    fn moved_alpha_carbon_does_not_share_the_product() {
        // The second and third alpha carbons end up 5.1 A apart, which breaks their link.
        not_shared(&PEPTIDE.replace("CA . GLY A 1 3 ? 7.600", "CA . GLY A 1 3 ? 8.900"));
    }

    #[test]
    fn changed_connection_does_not_share_the_product() {
        not_shared(&PEPTIDE.replace("ASN 2 OD1 B ZN", "ASN 2 ND2 B ZN"));
    }
}

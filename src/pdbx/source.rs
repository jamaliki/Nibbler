//! Construction of one immutable semantic model from generic CIF category rows.

use std::collections::{BTreeMap, BTreeSet};

use crate::cif::{BlockKind, CategoryView, CifDocument, FastMap};

use super::component_source::{component_definitions, definitions_conflict};
use super::error::SemanticError;
use super::fields::{
    Table, case_key, category_rows, category_tables, duplicate, optional_bool, optional_integer,
    optional_text, required_integer, required_single, required_text, row_context,
};
use super::model::{
    AsymUnit, AtomSite, BranchSchemeRow, BranchedLink, BranchedNode, ComponentResolution,
    Connection, ConnectionEndpoint, Entity, EntityKind, NonPolySchemeRow, PdbxModel, PolySchemeRow,
    PolymerMonomer,
};
use super::registry::{ComponentRegistry, resolve_builtin, unresolved};

/// Build a source-backed semantic PDBx model without performing profile validation.
///
/// # Errors
///
/// Returns [`SemanticError`] when required structural identifiers or typed values are invalid.
pub fn build_model(document: &CifDocument) -> Result<PdbxModel, SemanticError> {
    build_model_with_registry(document, None)
}

/// Build a source-backed PDBx model using an optional local CCD registry.
///
/// # Errors
///
/// Returns [`SemanticError`] for structural failures or incompatible component definitions.
pub fn build_model_with_registry(
    document: &CifDocument,
    registry: Option<&ComponentRegistry>,
) -> Result<PdbxModel, SemanticError> {
    build(document, registry, None)
}

/// Build the registry-free model of `document` with atom sites the caller already decoded.
///
/// `sites` must equal what decoding the document's `_atom_site` rows would give. Every
/// other category is decoded from `document`.
#[cfg(feature = "python")]
pub(crate) fn build_model_with_sites(
    document: &CifDocument,
    sites: Vec<AtomSite>,
) -> Result<PdbxModel, SemanticError> {
    build(document, None, Some(sites))
}

fn build(
    document: &CifDocument,
    registry: Option<&ComponentRegistry>,
    sites: Option<Vec<AtomSite>>,
) -> Result<PdbxModel, SemanticError> {
    let [block] = document.blocks() else {
        return Err(SemanticError::new(
            "PDBX_MODEL_BLOCK_COUNT",
            "a PDBx model requires exactly one data block",
            Vec::new(),
        ));
    };
    if block.kind() != BlockKind::Data {
        return Err(SemanticError::new(
            "PDBX_MODEL_DATA_BLOCK",
            "a PDBx model cannot be constructed from a global block",
            Vec::new(),
        ));
    }
    let categories = CategoryView::new(block.entries());
    let entry_id = required_single(categories, "entry", "id")?;
    let audit_dictionary = audit_conform(categories)?;
    let mut components = component_definitions(categories, ComponentResolution::Embedded)?;
    let entities = entities(categories)?;
    let asym_units = asym_units(categories)?;
    let atom_sites = match sites {
        Some(sites) => sites,
        None => atom_sites(categories)?,
    };

    for component_id in referenced_components(&entities, &atom_sites) {
        let key = case_key(component_id);
        if let (Some(embedded), Some(cached)) = (
            components.get(&key),
            registry.and_then(|value| value.get(component_id)),
        ) && definitions_conflict(embedded, cached)
        {
            return Err(SemanticError::new(
                "PDBX_COMPONENT_CONFLICT",
                format!(
                    "embedded and local CCD definitions for component {component_id:?} conflict"
                ),
                vec![format!("component={component_id}")],
            ));
        }
        components.entry(key).or_insert_with(|| {
            registry
                .and_then(|value| value.get(component_id))
                .cloned()
                .or_else(|| resolve_builtin(component_id))
                .unwrap_or_else(|| unresolved(component_id.to_owned()))
        });
    }

    Ok(PdbxModel {
        document: document.clone(),
        entry_id,
        audit_dictionary,
        entities,
        asym_units,
        components: components.into_values().collect(),
        atom_sites,
        connections: connections(categories)?,
        poly_scheme: poly_scheme(categories)?,
        nonpoly_scheme: nonpoly_scheme(categories)?,
        branch_scheme: branch_scheme(categories)?,
    })
}

/// Load component definitions from an already-parsed local CCD artifact.
///
/// # Errors
///
/// Returns [`SemanticError`] when component identifiers or typed chemical fields are invalid.
pub fn build_component_registry(
    document: &CifDocument,
) -> Result<ComponentRegistry, SemanticError> {
    let mut definitions = BTreeMap::new();
    for block in document
        .blocks()
        .iter()
        .filter(|block| block.kind() == BlockKind::Data)
    {
        for (key, definition) in component_definitions(
            CategoryView::new(block.entries()),
            ComponentResolution::LocalCcd,
        )? {
            if definitions.insert(key, definition).is_some() {
                return Err(SemanticError::new(
                    "PDBX_COMPONENT_DUPLICATE",
                    "the local CCD cache defines a component more than once",
                    Vec::new(),
                ));
            }
        }
    }
    Ok(ComponentRegistry::new(definitions))
}

fn audit_conform(categories: CategoryView<'_>) -> Result<Option<(String, String)>, SemanticError> {
    let Some(row) =
        category_rows(categories, "audit_conform", ["dict_name", "dict_version"]).next()
    else {
        return Ok(None);
    };
    Ok(Some((
        required_text(&row, "audit_conform", "dict_name", 0)?,
        required_text(&row, "audit_conform", "dict_version", 1)?,
    )))
}

fn entities(categories: CategoryView<'_>) -> Result<Vec<Entity>, SemanticError> {
    let polymer_rows = keyed_optional_text(categories, "entity_poly", "entity_id", "type")?;
    let nonpoly_rows =
        keyed_optional_text(categories, "pdbx_entity_nonpoly", "entity_id", "comp_id")?;
    let mut sequences: BTreeMap<String, Vec<PolymerMonomer>> = BTreeMap::new();
    for row in category_rows(
        categories,
        "entity_poly_seq",
        ["entity_id", "num", "mon_id", "hetero"],
    ) {
        let entity_id = required_text(&row, "entity_poly_seq", "entity_id", 0)?;
        sequences
            .entry(case_key(&entity_id))
            .or_default()
            .push(PolymerMonomer {
                number: required_integer(&row, "entity_poly_seq", "num", 1)?,
                component_id: required_text(&row, "entity_poly_seq", "mon_id", 2)?,
                heterogeneous: optional_bool(&row, 3).unwrap_or(false),
            });
    }
    for sequence in sequences.values_mut() {
        sequence.sort_by_key(|monomer| monomer.number);
    }

    let mut branch_nodes: BTreeMap<String, Vec<BranchedNode>> = BTreeMap::new();
    for row in category_rows(
        categories,
        "pdbx_entity_branch_list",
        ["entity_id", "num", "comp_id"],
    ) {
        let entity_id = required_text(&row, "pdbx_entity_branch_list", "entity_id", 0)?;
        branch_nodes
            .entry(case_key(&entity_id))
            .or_default()
            .push(BranchedNode {
                number: required_integer(&row, "pdbx_entity_branch_list", "num", 1)?,
                component_id: required_text(&row, "pdbx_entity_branch_list", "comp_id", 2)?,
            });
    }
    let branch_links = branch_links(categories)?;

    category_rows(categories, "entity", ["id", "type", "pdbx_description"])
        .map(|row| {
            let id = required_text(&row, "entity", "id", 0)?;
            let kind_text = required_text(&row, "entity", "type", 1)?;
            let kind = match case_key(&kind_text).as_str() {
                "polymer" => EntityKind::Polymer,
                "non-polymer" => EntityKind::NonPolymer,
                "branched" => EntityKind::Branched,
                "water" => EntityKind::Water,
                _ => {
                    return Err(SemanticError::new(
                        "PDBX_ENTITY_TYPE",
                        format!("entity {id:?} has unsupported type {kind_text:?}"),
                        row_context("entity", "type", row.row_index()),
                    ));
                }
            };
            let key = case_key(&id);
            Ok(Entity {
                id,
                kind,
                description: optional_text(&row, 2),
                polymer_type: polymer_rows.get(&key).cloned().flatten(),
                sequence: sequences.remove(&key).unwrap_or_default(),
                component_id: nonpoly_rows.get(&key).cloned().flatten(),
                branch_nodes: branch_nodes.remove(&key).unwrap_or_default(),
                branch_links: branch_links.get(&key).cloned().unwrap_or_default(),
            })
        })
        .collect()
}

fn keyed_optional_text(
    categories: CategoryView<'_>,
    category: &str,
    key_item: &str,
    value_item: &str,
) -> Result<BTreeMap<String, Option<String>>, SemanticError> {
    let mut output = BTreeMap::new();
    for row in category_rows(categories, category, [key_item, value_item]) {
        let key = required_text(&row, category, key_item, 0)?;
        if output
            .insert(case_key(&key), optional_text(&row, 1))
            .is_some()
        {
            return Err(duplicate(category, key_item, &key, row.row_index()));
        }
    }
    Ok(output)
}

fn branch_links(
    categories: CategoryView<'_>,
) -> Result<BTreeMap<String, Vec<BranchedLink>>, SemanticError> {
    let category = "pdbx_entity_branch_link";
    let mut output: BTreeMap<String, Vec<BranchedLink>> = BTreeMap::new();
    for row in category_rows(
        categories,
        category,
        [
            "entity_id",
            "entity_branch_list_num_1",
            "entity_branch_list_num_2",
            "atom_id_1",
            "atom_id_2",
            "value_order",
        ],
    ) {
        let entity_id = required_text(&row, category, "entity_id", 0)?;
        output
            .entry(case_key(&entity_id))
            .or_default()
            .push(BranchedLink {
                first_node: required_integer(&row, category, "entity_branch_list_num_1", 1)?,
                second_node: required_integer(&row, category, "entity_branch_list_num_2", 2)?,
                first_atom_id: required_text(&row, category, "atom_id_1", 3)?,
                second_atom_id: required_text(&row, category, "atom_id_2", 4)?,
                order: optional_text(&row, 5),
            });
    }
    Ok(output)
}

fn asym_units(categories: CategoryView<'_>) -> Result<Vec<AsymUnit>, SemanticError> {
    category_rows(categories, "struct_asym", ["id", "entity_id", "details"])
        .map(|row| {
            Ok(AsymUnit {
                id: required_text(&row, "struct_asym", "id", 0)?,
                entity_id: required_text(&row, "struct_asym", "entity_id", 1)?,
                details: optional_text(&row, 2),
            })
        })
        .collect()
}

/// The `_atom_site` items a model decodes, in [`AtomSite`] field order.
const ATOM_SITE_ITEMS: [&str; 21] = [
    "group_pdb",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_pdb_ins_code",
    "cartn_x",
    "cartn_y",
    "cartn_z",
    "occupancy",
    "b_iso_or_equiv",
    "pdbx_formal_charge",
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "pdbx_pdb_model_num",
];

/// Decoded sites whose values rows of the `_atom_site` table repeat.
trait RepeatedSites {
    /// Return the decoded site whose values logical row `row` repeats, and a mask whose
    /// bit `f` is set when the row holds that site's value of item `f` of
    /// [`ATOM_SITE_ITEMS`]. Without a site the row decodes every value.
    fn site(&self, row: usize) -> (Option<&AtomSite>, u32);
}

/// Rows that repeat no decoded site.
struct NoRepeats;

impl RepeatedSites for NoRepeats {
    #[inline(always)]
    fn site(&self, _row: usize) -> (Option<&AtomSite>, u32) {
        (None, 0)
    }
}

pub(super) fn atom_sites(categories: CategoryView<'_>) -> Result<Vec<AtomSite>, SemanticError> {
    decode_atom_sites(categories, &NoRepeats)
}

/// Decode the `_atom_site` rows, copying the values that `repeats` already holds.
///
/// A copied value equals the decode of the row's text, and that decode cannot fail, so
/// the result and the first error match a decode of every value.
fn decode_atom_sites(
    categories: CategoryView<'_>,
    repeats: &impl RepeatedSites,
) -> Result<Vec<AtomSite>, SemanticError> {
    let category = "atom_site";
    let tables = category_tables(categories, category, ATOM_SITE_ITEMS);
    let mut sites = Vec::with_capacity(tables.iter().map(Table::rows).sum());
    for table in &tables {
        let mut group_pdb = table.column(0);
        let mut id = table.column(1);
        let mut element = table.column(2);
        let mut label_atom_id = table.column(3);
        let mut label_alt_id = table.column(4);
        let mut label_comp_id = table.column(5);
        let mut label_asym_id = table.column(6);
        let mut label_entity_id = table.column(7);
        let mut label_seq_id = table.column(8);
        let mut insertion_code = table.column(9);
        let mut x = table.column(10);
        let mut y = table.column(11);
        let mut z = table.column(12);
        let mut occupancy = table.column(13);
        let mut b_iso = table.column(14);
        let mut formal_charge = table.column(15);
        let mut auth_seq_id = table.column(16);
        let mut auth_comp_id = table.column(17);
        let mut auth_asym_id = table.column(18);
        let mut auth_atom_id = table.column(19);
        let mut model_number = table.column(20);
        for row in 0..table.rows() {
            // Errors number rows across every occurrence of the category.
            let n = sites.len();
            let (repeated, held) = repeats.site(n);
            macro_rules! value {
                ($field:literal, $name:ident, $decode:expr) => {
                    match repeated {
                        Some(site) if held >> $field & 1 != 0 => site.$name.clone(),
                        _ => $decode,
                    }
                };
            }
            sites.push(AtomSite {
                group_pdb: value!(
                    0,
                    group_pdb,
                    group_pdb.required_short(row, n, category, "group_pdb")?
                ),
                id: value!(1, id, id.required_short(row, n, category, "id")?),
                element: value!(
                    2,
                    element,
                    element.required_short(row, n, category, "type_symbol")?
                ),
                label_atom_id: value!(
                    3,
                    label_atom_id,
                    label_atom_id.required_short(row, n, category, "label_atom_id")?
                ),
                label_alt_id: value!(4, label_alt_id, label_alt_id.optional_short(row)),
                label_comp_id: value!(
                    5,
                    label_comp_id,
                    label_comp_id.required_short(row, n, category, "label_comp_id")?
                ),
                label_asym_id: value!(
                    6,
                    label_asym_id,
                    label_asym_id.required_short(row, n, category, "label_asym_id")?
                ),
                label_entity_id: value!(
                    7,
                    label_entity_id,
                    label_entity_id.required_short(row, n, category, "label_entity_id")?
                ),
                label_seq_id: value!(
                    8,
                    label_seq_id,
                    label_seq_id.optional_integer(row, n, category, "label_seq_id")?
                ),
                insertion_code: value!(9, insertion_code, insertion_code.optional_short(row)),
                x: value!(10, x, x.required_float(row, n, category, "cartn_x")?),
                y: value!(11, y, y.required_float(row, n, category, "cartn_y")?),
                z: value!(12, z, z.required_float(row, n, category, "cartn_z")?),
                occupancy: value!(
                    13,
                    occupancy,
                    occupancy.required_float(row, n, category, "occupancy")?
                ),
                b_iso: value!(
                    14,
                    b_iso,
                    b_iso.optional_float(row, n, category, "b_iso_or_equiv")?
                ),
                formal_charge: value!(
                    15,
                    formal_charge,
                    formal_charge.optional_charge(row, n, category, "pdbx_formal_charge")?
                ),
                auth_seq_id: value!(
                    16,
                    auth_seq_id,
                    auth_seq_id.required_short(row, n, category, "auth_seq_id")?
                ),
                auth_comp_id: value!(
                    17,
                    auth_comp_id,
                    auth_comp_id.required_short(row, n, category, "auth_comp_id")?
                ),
                auth_asym_id: value!(
                    18,
                    auth_asym_id,
                    auth_asym_id.required_short(row, n, category, "auth_asym_id")?
                ),
                auth_atom_id: value!(
                    19,
                    auth_atom_id,
                    auth_atom_id.required_short(row, n, category, "auth_atom_id")?
                ),
                model_number: value!(
                    20,
                    model_number,
                    model_number.required_integer(row, n, category, "pdbx_pdb_model_num")?
                ),
            });
        }
    }
    Ok(sites)
}

/// Rows of a Reduce3 output that copy values of the decoded sites of its source.
#[cfg(feature = "python")]
struct CopiedSites<'a> {
    sites: &'a [AtomSite],
    rows: &'a [u32],
    made: &'a [u64],
    // The bit of `made` for the column of each item, or `None` when rows decode the item.
    bits: [Option<u32>; ATOM_SITE_ITEMS.len()],
    // The last mask of `made` and the items it leaves copied. Most rows share a mask.
    last: std::cell::Cell<(u64, u32)>,
}

#[cfg(feature = "python")]
impl CopiedSites<'_> {
    fn held(&self, made: u64) -> u32 {
        let (last_made, last_held) = self.last.get();
        if last_made == made {
            return last_held;
        }
        let held = self.bits.iter().rev().fold(0, |held, bit| {
            held << 1 | u32::from(bit.is_some_and(|bit| made >> bit & 1 == 0))
        });
        self.last.set((made, held));
        held
    }
}

#[cfg(feature = "python")]
impl RepeatedSites for CopiedSites<'_> {
    #[inline(always)]
    fn site(&self, row: usize) -> (Option<&AtomSite>, u32) {
        let site = self
            .rows
            .get(row)
            .and_then(|source| self.sites.get(usize::try_from(*source).ok()?));
        match (site, self.made.get(row)) {
            (Some(site), Some(made)) => (Some(site), self.held(*made)),
            _ => (None, 0),
        }
    }
}

/// Build the registry-free model of `document`, which Reduce3 wrote from the document of
/// `source`, without decoding again the `_atom_site` values it copied.
///
/// Row `i` of the output `_atom_site` loop holds the text of source row `rows[i]` in every
/// column `c` below 64 whose bit `made[i] >> c & 1` is clear, or copies nothing when
/// `rows[i]` names no source row. Both loops have the same tags. When `source` or
/// `document` does not keep one `_atom_site` loop whose sites hold its text, this decodes
/// every row.
///
/// # Errors
///
/// Returns the error of [`build_model`].
#[cfg(feature = "python")]
pub(crate) fn build_model_with_copies(
    document: &CifDocument,
    source: &PdbxModel,
    rows: &[u32],
    made: &[u64],
) -> Result<PdbxModel, SemanticError> {
    let (Some((_, sites)), [block]) = (source.loop_sites(), document.blocks()) else {
        return build(document, None, None);
    };
    let Some(output) = super::model::only_atom_site_loop(block)
        .filter(|output| output.row_count() == rows.len() && rows.len() == made.len())
    else {
        return build(document, None, None);
    };
    let bits = ATOM_SITE_ITEMS.map(|item| {
        output
            .tags()
            .iter()
            .position(|tag| {
                crate::cif::split_tag(tag).is_some_and(|(category, name)| {
                    category.eq_ignore_ascii_case("atom_site") && name.eq_ignore_ascii_case(item)
                })
            })
            .and_then(|column| u32::try_from(column).ok())
            .filter(|bit| *bit < u64::BITS)
    });
    let copies = CopiedSites {
        sites,
        rows,
        made,
        bits,
        last: std::cell::Cell::new((u64::MAX, 0)),
    };
    let sites = decode_atom_sites(CategoryView::new(block.entries()), &copies)?;
    build(document, None, Some(sites))
}

fn poly_scheme(categories: CategoryView<'_>) -> Result<Vec<PolySchemeRow>, SemanticError> {
    let category = "pdbx_poly_seq_scheme";
    category_rows(
        categories,
        category,
        ["asym_id", "entity_id", "seq_id", "mon_id"],
    )
    .map(|row| {
        Ok(PolySchemeRow {
            asym_id: required_text(&row, category, "asym_id", 0)?,
            entity_id: required_text(&row, category, "entity_id", 1)?,
            seq_id: required_integer(&row, category, "seq_id", 2)?,
            component_id: required_text(&row, category, "mon_id", 3)?,
        })
    })
    .collect()
}

fn nonpoly_scheme(categories: CategoryView<'_>) -> Result<Vec<NonPolySchemeRow>, SemanticError> {
    let category = "pdbx_nonpoly_scheme";
    category_rows(
        categories,
        category,
        ["asym_id", "entity_id", "mon_id", "auth_seq_num"],
    )
    .map(|row| {
        Ok(NonPolySchemeRow {
            asym_id: required_text(&row, category, "asym_id", 0)?,
            entity_id: required_text(&row, category, "entity_id", 1)?,
            component_id: required_text(&row, category, "mon_id", 2)?,
            auth_seq_id: required_text(&row, category, "auth_seq_num", 3)?,
        })
    })
    .collect()
}

fn branch_scheme(categories: CategoryView<'_>) -> Result<Vec<BranchSchemeRow>, SemanticError> {
    let category = "pdbx_branch_scheme";
    category_rows(
        categories,
        category,
        [
            "asym_id",
            "entity_id",
            "num",
            "mon_id",
            "auth_seq_num",
            "pdb_seq_num",
        ],
    )
    .map(|row| {
        Ok(BranchSchemeRow {
            asym_id: required_text(&row, category, "asym_id", 0)?,
            entity_id: required_text(&row, category, "entity_id", 1)?,
            number: required_integer(&row, category, "num", 2)?,
            component_id: required_text(&row, category, "mon_id", 3)?,
            auth_seq_id: match optional_text(&row, 4) {
                Some(value) => value,
                None => required_text(&row, category, "pdb_seq_num", 5)?,
            },
        })
    })
    .collect()
}

pub(super) fn connections(categories: CategoryView<'_>) -> Result<Vec<Connection>, SemanticError> {
    let category = "struct_conn";
    category_rows(
        categories,
        category,
        [
            "id",
            "conn_type_id",
            "ptnr1_label_asym_id",
            "ptnr1_label_comp_id",
            "ptnr1_label_seq_id",
            "ptnr1_label_atom_id",
            "ptnr2_label_asym_id",
            "ptnr2_label_comp_id",
            "ptnr2_label_seq_id",
            "ptnr2_label_atom_id",
        ],
    )
    .map(|row| {
        Ok(Connection {
            id: required_text(&row, category, "id", 0)?,
            kind: required_text(&row, category, "conn_type_id", 1)?,
            first: ConnectionEndpoint {
                asym_id: required_text(&row, category, "ptnr1_label_asym_id", 2)?,
                component_id: optional_text(&row, 3),
                sequence_id: optional_integer(&row, category, "ptnr1_label_seq_id", 4)?,
                atom_id: required_text(&row, category, "ptnr1_label_atom_id", 5)?,
            },
            second: ConnectionEndpoint {
                asym_id: required_text(&row, category, "ptnr2_label_asym_id", 6)?,
                component_id: optional_text(&row, 7),
                sequence_id: optional_integer(&row, category, "ptnr2_label_seq_id", 8)?,
                atom_id: required_text(&row, category, "ptnr2_label_atom_id", 9)?,
            },
        })
    })
    .collect()
}

fn referenced_components<'a>(entities: &'a [Entity], atoms: &'a [AtomSite]) -> BTreeSet<&'a str> {
    // A few identifiers repeat thousands of times. A hash set drops the repeats, so the
    // ordered set sorts each identifier once.
    let mut output: FastMap<&str, ()> = FastMap::default();
    for entity in entities {
        output.extend(
            entity
                .sequence
                .iter()
                .map(|monomer| (monomer.component_id.as_str(), ())),
        );
        output.extend(entity.component_id.as_deref().map(|id| (id, ())));
        output.extend(
            entity
                .branch_nodes
                .iter()
                .map(|node| (node.component_id.as_str(), ())),
        );
    }
    // Atoms repeat a few component identifiers, so skip each repeat of the previous one
    // before the set lookup.
    let mut previous = None;
    for atom in atoms {
        let component = atom.label_comp_id.as_str();
        if previous != Some(component) {
            output.insert(component, ());
            previous = Some(component);
        }
    }
    output.into_keys().collect()
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::cif::parse;

    use super::*;

    const CHEMISTRY: &[u8] = include_bytes!("../../tests/fixtures/chemistry/ligand_ion_water.cif");
    const GLYCAN: &[u8] = include_bytes!("../../tests/fixtures/chemistry/branched_glycan.cif");

    #[test]
    fn retains_coordinate_and_chemistry_records_after_decoding() {
        let source = String::from_utf8(CHEMISTRY.to_vec())
            .expect("fixture is UTF-8")
            .replace("PG A ATP B 2 . ?", "PG A ATP B 2 . A");
        let document = parse(source.as_bytes()).expect("fixture parses");
        let model = build_model(&document).expect("fixture builds");

        assert_eq!(
            model.audit_dictionary,
            Some(("mmcif_pdbx.dic".to_owned(), "5.416".to_owned()))
        );
        assert_eq!(
            model.entities[0],
            Entity {
                id: "1".to_owned(),
                kind: EntityKind::Polymer,
                description: Some("ALA-MSE test peptide".to_owned()),
                polymer_type: Some("polypeptide(L)".to_owned()),
                sequence: vec![
                    PolymerMonomer {
                        number: 1,
                        component_id: "ALA".to_owned(),
                        heterogeneous: false,
                    },
                    PolymerMonomer {
                        number: 2,
                        component_id: "MSE".to_owned(),
                        heterogeneous: false,
                    },
                ],
                component_id: None,
                branch_nodes: Vec::new(),
                branch_links: Vec::new(),
            }
        );
        assert_eq!(model.entities[3].component_id.as_deref(), Some("HOH"));
        assert_eq!(
            model.asym_units[1],
            AsymUnit {
                id: "B".to_owned(),
                entity_id: "2".to_owned(),
                details: Some(
                    "ligand with alternate locations and negative author residue number".to_owned()
                ),
            }
        );
        assert_eq!(
            model.atom_sites[2],
            AtomSite {
                group_pdb: "HETATM".into(),
                id: "3".into(),
                element: "P".into(),
                label_atom_id: "PG".into(),
                label_alt_id: Some("A".into()),
                label_comp_id: "ATP".into(),
                label_asym_id: "B".into(),
                label_entity_id: "2".into(),
                label_seq_id: None,
                insertion_code: Some("A".into()),
                x: 7.0,
                y: 1.0,
                z: 0.0,
                occupancy: 0.6,
                b_iso: Some(20.0),
                formal_charge: None,
                auth_seq_id: "-1".into(),
                auth_comp_id: "ATP".into(),
                auth_asym_id: "L".into(),
                auth_atom_id: "PG".into(),
                model_number: 1,
            }
        );
        assert_eq!(model.atom_sites[4].formal_charge, Some(2));
        assert_eq!(
            model.connections[0],
            Connection {
                id: "ZN1".to_owned(),
                kind: "metalc".to_owned(),
                first: ConnectionEndpoint {
                    asym_id: "B".to_owned(),
                    component_id: Some("ATP".to_owned()),
                    sequence_id: None,
                    atom_id: "PG".to_owned(),
                },
                second: ConnectionEndpoint {
                    asym_id: "C".to_owned(),
                    component_id: Some("ZN".to_owned()),
                    sequence_id: None,
                    atom_id: "ZN".to_owned(),
                },
            }
        );
        assert_eq!(
            model.poly_scheme[0],
            PolySchemeRow {
                asym_id: "A".to_owned(),
                entity_id: "1".to_owned(),
                seq_id: 1,
                component_id: "ALA".to_owned(),
            }
        );
        assert_eq!(
            model.nonpoly_scheme[0],
            NonPolySchemeRow {
                asym_id: "B".to_owned(),
                entity_id: "2".to_owned(),
                component_id: "ATP".to_owned(),
                auth_seq_id: "-1".to_owned(),
            }
        );

        let ccd = parse(
            br#"data_TST
_chem_comp.id TST
_chem_comp.name 'TEST COMPONENT'
_chem_comp.type non-polymer
_chem_comp.formula 'C N'
_chem_comp.formula_weight 26.02
_chem_comp.pdbx_formal_charge 1
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
TST C1 C 1
TST N1 N 0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
TST C1 N1 doub y
"#,
        )
        .expect("component fixture parses");
        let registry = build_component_registry(&ccd).expect("component fixture builds");
        let component = registry.get("TST").expect("component exists");
        assert_eq!(component.id, "TST");
        assert_eq!(component.name.as_deref(), Some("TEST COMPONENT"));
        assert_eq!(component.component_type.as_deref(), Some("non-polymer"));
        assert_eq!(component.formula.as_deref(), Some("C N"));
        assert_eq!(component.formula_weight, Some(26.02));
        assert_eq!(component.formal_charge, Some(1));
        assert_eq!(component.resolution, ComponentResolution::LocalCcd);
        assert_eq!(
            component.atoms,
            [
                super::super::model::ComponentAtom {
                    atom_id: "C1".to_owned(),
                    alternate_id: None,
                    element: "C".to_owned(),
                    formal_charge: Some(1),
                    leaving: None,
                },
                super::super::model::ComponentAtom {
                    atom_id: "N1".to_owned(),
                    alternate_id: None,
                    element: "N".to_owned(),
                    formal_charge: Some(0),
                    leaving: None,
                },
            ]
        );
        assert_eq!(
            component.bonds,
            [super::super::model::ComponentBond {
                first_atom_id: "C1".to_owned(),
                second_atom_id: "N1".to_owned(),
                order: "doub".to_owned(),
                aromatic: Some(true),
            }]
        );
    }

    #[test]
    fn retains_branched_entity_and_scheme_records_after_decoding() {
        let document = parse(GLYCAN).expect("fixture parses");
        let model = build_model(&document).expect("fixture builds");
        let entity = &model.entities[0];

        assert_eq!(entity.kind, EntityKind::Branched);
        assert_eq!(
            entity.branch_nodes,
            [
                BranchedNode {
                    number: 1,
                    component_id: "NAG".to_owned(),
                },
                BranchedNode {
                    number: 2,
                    component_id: "MAN".to_owned(),
                },
            ]
        );
        assert_eq!(
            entity.branch_links,
            [BranchedLink {
                first_node: 1,
                second_node: 2,
                first_atom_id: "O4".to_owned(),
                second_atom_id: "C1".to_owned(),
                order: Some("sing".to_owned()),
            }]
        );
        assert_eq!(
            model.branch_scheme[1],
            BranchSchemeRow {
                asym_id: "G".to_owned(),
                entity_id: "1".to_owned(),
                number: 2,
                component_id: "MAN".to_owned(),
                auth_seq_id: "2".to_owned(),
            }
        );
    }
}

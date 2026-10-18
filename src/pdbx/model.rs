use crate::cif::CifDocument;
#[cfg(feature = "python")]
use crate::cif::{CifBlock, CifEntry, CifLoop, CifValueRef, in_category, split_tag};

use super::short_text::ShortText;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum EntityKind {
    Polymer,
    NonPolymer,
    Branched,
    Water,
}

impl EntityKind {
    pub(super) const fn as_str(self) -> &'static str {
        match self {
            Self::Polymer => "polymer",
            Self::NonPolymer => "non-polymer",
            Self::Branched => "branched",
            Self::Water => "water",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct PolymerMonomer {
    pub number: i64,
    pub component_id: String,
    pub heterogeneous: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct BranchedNode {
    pub number: i64,
    pub component_id: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct BranchedLink {
    pub first_node: i64,
    pub second_node: i64,
    pub first_atom_id: String,
    pub second_atom_id: String,
    pub order: Option<String>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Entity {
    pub id: String,
    pub kind: EntityKind,
    pub description: Option<String>,
    pub polymer_type: Option<String>,
    pub sequence: Vec<PolymerMonomer>,
    pub component_id: Option<String>,
    pub branch_nodes: Vec<BranchedNode>,
    pub branch_links: Vec<BranchedLink>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct AsymUnit {
    pub id: String,
    pub entity_id: String,
    pub details: Option<String>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum ComponentResolution {
    Embedded,
    MinimalRegistry,
    LocalCcd,
    Unresolved,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ComponentAtom {
    pub atom_id: String,
    pub alternate_id: Option<String>,
    pub element: String,
    pub formal_charge: Option<i32>,
    pub leaving: Option<bool>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ComponentBond {
    pub first_atom_id: String,
    pub second_atom_id: String,
    pub order: String,
    pub aromatic: Option<bool>,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct ComponentDefinition {
    pub id: String,
    pub name: Option<String>,
    pub component_type: Option<String>,
    pub formula: Option<String>,
    pub formula_weight: Option<f64>,
    pub formal_charge: Option<i32>,
    pub atoms: Vec<ComponentAtom>,
    pub bonds: Vec<ComponentBond>,
    pub resolution: ComponentResolution,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct AtomSite {
    pub group_pdb: ShortText,
    pub id: ShortText,
    pub element: ShortText,
    pub label_atom_id: ShortText,
    pub label_alt_id: Option<ShortText>,
    pub label_comp_id: ShortText,
    pub label_asym_id: ShortText,
    pub label_entity_id: ShortText,
    pub label_seq_id: Option<i64>,
    pub insertion_code: Option<ShortText>,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f64,
    pub b_iso: Option<f64>,
    pub formal_charge: Option<i32>,
    pub auth_seq_id: ShortText,
    pub auth_comp_id: ShortText,
    pub auth_asym_id: ShortText,
    pub auth_atom_id: ShortText,
    pub model_number: i64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ConnectionEndpoint {
    pub asym_id: String,
    pub component_id: Option<String>,
    pub sequence_id: Option<i64>,
    pub atom_id: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Connection {
    pub id: String,
    pub kind: String,
    pub first: ConnectionEndpoint,
    pub second: ConnectionEndpoint,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct PolySchemeRow {
    pub asym_id: String,
    pub entity_id: String,
    pub seq_id: i64,
    pub component_id: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct NonPolySchemeRow {
    pub asym_id: String,
    pub entity_id: String,
    pub component_id: String,
    pub auth_seq_id: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct BranchSchemeRow {
    pub asym_id: String,
    pub entity_id: String,
    pub number: i64,
    pub component_id: String,
    pub auth_seq_id: String,
}

/// An immutable, source-backed PDBx coordinate model.
#[derive(Clone, Debug, PartialEq)]
pub struct PdbxModel {
    pub(super) document: CifDocument,
    pub(super) entry_id: String,
    pub(super) audit_dictionary: Option<(String, String)>,
    pub(super) entities: Vec<Entity>,
    pub(super) asym_units: Vec<AsymUnit>,
    pub(super) components: Vec<ComponentDefinition>,
    pub(super) atom_sites: Vec<AtomSite>,
    pub(super) connections: Vec<Connection>,
    pub(super) poly_scheme: Vec<PolySchemeRow>,
    pub(super) nonpoly_scheme: Vec<NonPolySchemeRow>,
    pub(super) branch_scheme: Vec<BranchSchemeRow>,
}

/// Return the loop of `block` that holds every `_atom_site` item of the block, when it
/// holds only `_atom_site` items.
#[cfg(feature = "python")]
pub(crate) fn only_atom_site_loop(block: &CifBlock) -> Option<&CifLoop> {
    let mut holders = block.entries().iter().filter(|entry| match entry {
        CifEntry::Item(item) => in_category(item.tag(), "atom_site"),
        CifEntry::Loop(values) => values
            .tags()
            .iter()
            .any(|tag| in_category(tag, "atom_site")),
        CifEntry::Frame(_) => false,
    });
    let (Some(CifEntry::Loop(only)), None) = (holders.next(), holders.next()) else {
        return None;
    };
    only.tags()
        .iter()
        .all(|tag| in_category(tag, "atom_site"))
        .then_some(only)
}

/// One `_atom_site` item for every row, laid out for numpy.
#[cfg(feature = "python")]
pub(crate) enum SiteArray {
    /// `width` little-endian UCS-4 code points for each row, padded with zeros.
    Text { width: usize, codes: Vec<u8> },
    /// Little-endian `int64` values, 0 where missing, and the missing kind of each row:
    /// 0 when present, 1 for `?` and 2 for `.`.
    Integer { values: Vec<u8>, kinds: Vec<u8> },
}

/// Return `items` of the `_atom_site` rows of `document`, read from `model`, the model
/// built from it. Present text is the stored text, a missing value reads `?` or `.`, and
/// `label_seq_id` and `pdbx_PDB_model_num` are the parsed integers.
///
/// A model is built only when every required item of every row parses, so the result
/// equals a projection of the same items. Return `None` when the rows do not come from
/// the one `_atom_site` loop of the document, or an item is absent or none of the text
/// items and two integer items of [`AtomSite`].
#[cfg(feature = "python")]
pub(crate) fn atom_site_arrays(
    document: &CifDocument,
    model: &PdbxModel,
    items: &[String],
) -> Option<Vec<SiteArray>> {
    let [block] = document.blocks() else {
        return None;
    };
    let rows = only_atom_site_loop(block)?;
    let sites = model.atom_sites.as_slice();
    if sites.len() != rows.row_count() {
        return None;
    }
    let mut columns = items
        .iter()
        .map(|item| {
            let column = rows.tags().iter().position(|tag| {
                split_tag(tag).is_some_and(|(category, name)| {
                    category.eq_ignore_ascii_case("atom_site") && name.eq_ignore_ascii_case(item)
                })
            })?;
            Some(SiteColumn {
                field: SiteField::named(item)?,
                column,
                // `?` and `.` are one character wide.
                width: 1,
                kinds: vec![0; sites.len()],
            })
        })
        .collect::<Option<Vec<_>>>()?;
    // A site is several hundred bytes, so each pass reads all items of a site at once
    // instead of reading the sites once per item.
    for (row, site) in sites.iter().enumerate() {
        for column in &mut columns {
            let missing = match column.field.text(site) {
                // Every UTF-8 character has one byte that is not a continuation byte.
                Ok(Some(text)) => {
                    let count = text
                        .as_bytes()
                        .iter()
                        .filter(|&&byte| byte & 0xC0 != 0x80)
                        .count();
                    column.width = column.width.max(count);
                    false
                }
                Ok(None) => true,
                Err(()) => column.field.integer(site).is_none(),
            };
            // A site holds an optional value as `None`; the loop tells `?` from `.`.
            if missing {
                column.kinds[row] = match rows.value(row, column.column)? {
                    CifValueRef::Unknown => 1,
                    CifValueRef::NotApplicable => 2,
                    _ => return None,
                };
            }
        }
    }
    let mut arrays = columns
        .iter()
        .map(|column| match column.field {
            SiteField::LabelSeqId | SiteField::ModelNumber => SiteArray::Integer {
                values: Vec::with_capacity(sites.len() * 8),
                kinds: Vec::new(),
            },
            _ => SiteArray::Text {
                width: column.width,
                codes: vec![0; sites.len() * column.width * 4],
            },
        })
        .collect::<Vec<_>>();
    for (row, site) in sites.iter().enumerate() {
        for (column, array) in columns.iter().zip(&mut arrays) {
            match array {
                SiteArray::Text { width, codes } => {
                    let cells = codes.get_mut(row * *width * 4..(row + 1) * *width * 4)?;
                    match column.field.text(site) {
                        // The buffer is zeroed, so an ASCII character needs only its low byte.
                        Ok(Some(text)) if text.as_bytes().is_ascii() => {
                            for (cell, &byte) in cells.chunks_exact_mut(4).zip(text.as_bytes()) {
                                cell[0] = byte;
                            }
                        }
                        Ok(Some(text)) => {
                            for (cell, character) in
                                cells.chunks_exact_mut(4).zip(text.as_str().chars())
                            {
                                cell.copy_from_slice(&u32::from(character).to_le_bytes());
                            }
                        }
                        _ => cells[0] = if column.kinds[row] == 1 { b'?' } else { b'.' },
                    }
                }
                SiteArray::Integer { values, .. } => {
                    let value = column.field.integer(site).unwrap_or(0);
                    values.extend_from_slice(&value.to_le_bytes());
                }
            }
        }
    }
    for (column, array) in columns.into_iter().zip(&mut arrays) {
        if let SiteArray::Integer { kinds, .. } = array {
            *kinds = column.kinds;
        }
    }
    Some(arrays)
}

/// One requested `_atom_site` item: its field, loop column, widest text in characters,
/// and the missing kind of each row.
#[cfg(feature = "python")]
struct SiteColumn {
    field: SiteField,
    column: usize,
    width: usize,
    kinds: Vec<u8>,
}

/// The [`AtomSite`] fields that [`atom_site_arrays`] reads.
#[cfg(feature = "python")]
#[derive(Clone, Copy)]
enum SiteField {
    GroupPdb,
    Id,
    TypeSymbol,
    LabelAtomId,
    LabelAltId,
    LabelCompId,
    LabelAsymId,
    LabelEntityId,
    InsertionCode,
    AuthSeqId,
    AuthCompId,
    AuthAsymId,
    AuthAtomId,
    LabelSeqId,
    ModelNumber,
}

#[cfg(feature = "python")]
impl SiteField {
    fn named(item: &str) -> Option<Self> {
        Some(match item.to_ascii_lowercase().as_str() {
            "group_pdb" => Self::GroupPdb,
            "id" => Self::Id,
            "type_symbol" => Self::TypeSymbol,
            "label_atom_id" => Self::LabelAtomId,
            "label_alt_id" => Self::LabelAltId,
            "label_comp_id" => Self::LabelCompId,
            "label_asym_id" => Self::LabelAsymId,
            "label_entity_id" => Self::LabelEntityId,
            "pdbx_pdb_ins_code" => Self::InsertionCode,
            "auth_seq_id" => Self::AuthSeqId,
            "auth_comp_id" => Self::AuthCompId,
            "auth_asym_id" => Self::AuthAsymId,
            "auth_atom_id" => Self::AuthAtomId,
            "label_seq_id" => Self::LabelSeqId,
            "pdbx_pdb_model_num" => Self::ModelNumber,
            _ => return None,
        })
    }

    /// Return the text of a text field, or `Err` for an integer field.
    #[inline]
    fn text(self, site: &AtomSite) -> Result<Option<&ShortText>, ()> {
        Ok(Some(match self {
            Self::GroupPdb => &site.group_pdb,
            Self::Id => &site.id,
            Self::TypeSymbol => &site.element,
            Self::LabelAtomId => &site.label_atom_id,
            Self::LabelAltId => return Ok(site.label_alt_id.as_ref()),
            Self::LabelCompId => &site.label_comp_id,
            Self::LabelAsymId => &site.label_asym_id,
            Self::LabelEntityId => &site.label_entity_id,
            Self::InsertionCode => return Ok(site.insertion_code.as_ref()),
            Self::AuthSeqId => &site.auth_seq_id,
            Self::AuthCompId => &site.auth_comp_id,
            Self::AuthAsymId => &site.auth_asym_id,
            Self::AuthAtomId => &site.auth_atom_id,
            Self::LabelSeqId | Self::ModelNumber => return Err(()),
        }))
    }

    /// Return the value of an integer field, or `None` when it is missing or the field
    /// holds text.
    #[inline]
    const fn integer(self, site: &AtomSite) -> Option<i64> {
        match self {
            Self::LabelSeqId => site.label_seq_id,
            Self::ModelNumber => Some(site.model_number),
            _ => None,
        }
    }
}

impl PdbxModel {
    /// Return `_entry.id`.
    #[must_use]
    pub fn entry_id(&self) -> &str {
        &self.entry_id
    }

    /// Return the number of coordinate entities.
    #[must_use]
    pub fn entity_count(&self) -> usize {
        self.entities.len()
    }

    /// Return the number of asymmetric-unit instances.
    #[must_use]
    pub fn asym_unit_count(&self) -> usize {
        self.asym_units.len()
    }

    /// Return the number of component definitions and unresolved placeholders.
    #[must_use]
    pub fn component_count(&self) -> usize {
        self.components.len()
    }

    /// Return the number of coordinate atom sites.
    #[must_use]
    pub fn atom_site_count(&self) -> usize {
        self.atom_sites.len()
    }

    /// Return the number of explicit inter-site connections.
    #[must_use]
    pub fn connection_count(&self) -> usize {
        self.connections.len()
    }

    /// Iterate over canonical `_entity.type` values in source order.
    pub fn entity_kinds(&self) -> impl ExactSizeIterator<Item = &'static str> + '_ {
        self.entities.iter().map(|entity| entity.kind.as_str())
    }

    /// Iterate over resolved component identifiers in stable order.
    pub fn component_ids(&self) -> impl ExactSizeIterator<Item = &str> {
        self.components
            .iter()
            .map(|component| component.id.as_str())
    }

    pub(crate) const fn source_document(&self) -> &CifDocument {
        &self.document
    }

    /// The document's only `_atom_site` loop with the atom sites, one per row, when the
    /// sites hold the loop's text: every tag of the loop is an `_atom_site` item, the block
    /// has no other `_atom_site` items, and no value is a typed number with a source
    /// spelling, whose text the decode does not keep.
    #[cfg(feature = "python")]
    pub(crate) fn loop_sites(&self) -> Option<(&CifLoop, &[AtomSite])> {
        let [block] = self.document.blocks() else {
            return None;
        };
        let only = only_atom_site_loop(block)?;
        let fits = !only.holds_number_spellings() && self.atom_sites.len() == only.row_count();
        fits.then_some((only, self.atom_sites.as_slice()))
    }

    pub(crate) fn coordinate_entity_ids(&self) -> impl Iterator<Item = &str> {
        self.entities.iter().map(|entity| entity.id.as_str())
    }

    pub(crate) fn asym_entity_ids(&self) -> impl Iterator<Item = (&str, &str)> {
        self.asym_units
            .iter()
            .map(|asym| (asym.id.as_str(), asym.entity_id.as_str()))
    }

    pub(crate) fn coordinate_model_numbers(&self) -> impl Iterator<Item = i64> + '_ {
        self.atom_sites.iter().map(|atom| atom.model_number)
    }

    pub(crate) fn polymer_sites(&self) -> impl Iterator<Item = (i64, &str, i64, &str)> {
        self.atom_sites.iter().filter_map(|atom| {
            atom.label_seq_id.map(|sequence_id| {
                (
                    atom.model_number,
                    atom.label_asym_id.as_str(),
                    sequence_id,
                    atom.label_comp_id.as_str(),
                )
            })
        })
    }
}

#[cfg(all(test, feature = "python"))]
#[allow(clippy::expect_used)]
mod tests {
    use crate::cif::parse;
    use crate::pdbx::build_model;

    use super::*;

    const CHEMISTRY: &str = include_str!("../../tests/fixtures/chemistry/ligand_ion_water.cif");
    const ITEMS: [&str; 15] = [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "pdbx_PDB_ins_code",
        "auth_seq_id",
        "auth_comp_id",
        "auth_asym_id",
        "auth_atom_id",
        "label_seq_id",
        "pdbx_PDB_model_num",
    ];

    /// Return row `row` of a text array without its padding.
    fn text_at(width: usize, codes: &[u8], row: usize) -> String {
        codes[row * width * 4..(row + 1) * width * 4]
            .chunks_exact(4)
            .map(|code| u32::from_le_bytes([code[0], code[1], code[2], code[3]]))
            .take_while(|&code| code != 0)
            .map(|code| char::from_u32(code).expect("a code point"))
            .collect()
    }

    #[test]
    fn site_arrays_hold_the_loop_values() {
        let source = CHEMISTRY
            .replace("ATOM 1 C CA . ALA A 1 1 ?", "ATOM 1 C CA '?' ALA A 1 1 '.'")
            .replace(
                "ATOM 2 Se SE . MSE A 1 2 ?",
                "ATOM 2 Se SE ? MSE A 1 +2 \u{e9}",
            )
            .replace("? 900 HOH W O 1", "? '9 0' HOH W O 1");
        let document = parse(source.as_bytes()).expect("the edited fixture parses");
        let model = build_model(&document).expect("the edited fixture builds");
        let items = ITEMS.map(str::to_owned);
        let arrays =
            atom_site_arrays(&document, &model, &items).expect("the model holds every item");
        let rows = document
            .blocks()
            .first()
            .and_then(only_atom_site_loop)
            .expect("one atom-site loop");
        for (item, array) in items.iter().zip(&arrays) {
            let tag = format!("_atom_site.{item}");
            let column = rows
                .tags()
                .iter()
                .position(|name| *name == tag)
                .expect("a column");
            for row in 0..rows.row_count() {
                let value = rows.value(row, column).expect("a value");
                match array {
                    SiteArray::Text { width, codes } => {
                        let expected = match value {
                            CifValueRef::Text(text) => text.as_str().to_owned(),
                            CifValueRef::Unknown => "?".to_owned(),
                            other => {
                                assert!(matches!(other, CifValueRef::NotApplicable), "{other:?}");
                                ".".to_owned()
                            }
                        };
                        assert_eq!(text_at(*width, codes, row), expected, "{item} row {row}");
                    }
                    SiteArray::Integer { values, kinds } => {
                        let expected = match value {
                            CifValueRef::Text(text) => {
                                (text.as_str().parse().expect("an integer"), 0)
                            }
                            CifValueRef::Unknown => (0, 1),
                            other => {
                                assert!(matches!(other, CifValueRef::NotApplicable), "{other:?}");
                                (0, 2)
                            }
                        };
                        let bytes = values[row * 8..row * 8 + 8]
                            .try_into()
                            .expect("eight bytes");
                        assert_eq!(
                            (i64::from_le_bytes(bytes), kinds[row]),
                            expected,
                            "{item} row {row}"
                        );
                    }
                }
            }
        }
        let SiteArray::Text { width, codes } = &arrays[4] else {
            unreachable!("label_alt_id is text");
        };
        assert_eq!(
            [0, 1, 2].map(|row| text_at(*width, codes, row)),
            ["?", "?", "A"]
        );
        let SiteArray::Text { width, codes } = &arrays[8] else {
            unreachable!("pdbx_PDB_ins_code is text");
        };
        assert_eq!(
            [0, 1, 2].map(|row| text_at(*width, codes, row)),
            [".", "\u{e9}", "?"]
        );
        let SiteArray::Text { width, codes } = &arrays[9] else {
            unreachable!("auth_seq_id is text");
        };
        assert_eq!((*width, text_at(*width, codes, 5)), (3, "9 0".to_owned()));
        assert!(atom_site_arrays(&document, &model, &["Cartn_x".to_owned()]).is_none());
    }
}

//! Biological assemblies written out as explicit coordinates ([`assembly_document`]).

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};

use crate::cif::numeric::parse_float;
use crate::cif::{
    BlockCategory, BlockKind, CategoryView, CifBlock, CifDocument, CifEntry, CifLoop, CifValueRef,
    TextCell, TextColumnBuilder, in_category,
};

use super::PdbxModel;
use super::error::SemanticError;
use super::fields::{category_rows, required_float, required_text, row_context};

/// The most atom sites an assembly may have; a larger one fails with
/// `PDBX_ASSEMBLY_TOO_LARGE` rather than exhausting memory.
pub const MAX_ASSEMBLY_ATOM_SITES: usize = 50_000_000;

/// The most operator combinations one expression may produce.
const MAX_COMBINATIONS: usize = 1_000_000;

/// Copies of one atom closer than this (an atom on a symmetry axis of the assembly) are
/// written once, as gemmi's `transform_to_assembly` does by default.
const MERGE_DISTANCE: f64 = 0.2;

const UNCHANGED_BEFORE_SCHEMES: &[&str] = &[
    "entry",
    "audit_conform",
    "entity",
    "entity_poly",
    "entity_poly_seq",
    "pdbx_entity_nonpoly",
    "pdbx_entity_branch",
    "pdbx_entity_branch_list",
    "pdbx_entity_branch_link",
    "chem_comp",
    "chem_comp_atom",
    "chem_comp_bond",
];

/// A rigid transform: `x' = rotation · x + translation`.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Transform {
    rotation: [[f64; 3]; 3],
    translation: [f64; 3],
}

impl Transform {
    const IDENTITY: Self = Self {
        rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        translation: [0.0; 3],
    };

    /// The transform that applies `inner` first, then `self`.
    fn after(&self, inner: &Self) -> Self {
        let rotation = std::array::from_fn(|i| {
            std::array::from_fn(|k| {
                (0..3)
                    .map(|j| self.rotation[i][j] * inner.rotation[j][k])
                    .sum()
            })
        });
        Self {
            rotation,
            translation: self.apply(inner.translation),
        }
    }

    fn rotate(&self, point: [f64; 3]) -> [f64; 3] {
        std::array::from_fn(|i| (0..3).map(|j| self.rotation[i][j] * point[j]).sum())
    }

    fn apply(&self, point: [f64; 3]) -> [f64; 3] {
        let rotated = self.rotate(point);
        std::array::from_fn(|i| rotated[i] + self.translation[i])
    }

    fn rotates(&self) -> bool {
        self.rotation != Self::IDENTITY.rotation
    }

    /// `R · U · Rᵀ` for a symmetric tensor stored as 11, 22, 33, 12, 13, 23.
    fn rotate_tensor(&self, u: [f64; 6]) -> [f64; 6] {
        let [u11, u22, u33, u12, u13, u23] = u;
        let full = [[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]];
        let r = &self.rotation;
        let ru: [[f64; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|k| (0..3).map(|j| r[i][j] * full[j][k]).sum())
        });
        let out = |i: usize, k: usize| (0..3).map(|j| ru[i][j] * r[k][j]).sum::<f64>();
        [
            out(0, 0),
            out(1, 1),
            out(2, 2),
            out(0, 1),
            out(0, 2),
            out(1, 2),
        ]
    }
}

/// One copy of part of the asymmetric unit.
#[derive(Debug)]
struct AssemblyCopy {
    number: usize,
    transform: Transform,
    asym_ids: Vec<String>,
}

impl AssemblyCopy {
    fn has(&self, asym_id: &str) -> bool {
        self.asym_ids.iter().any(|id| id == asym_id)
    }

    /// `source` renamed for this copy: unchanged for copy 1, `source-n` otherwise.
    fn rename<'a>(&self, source: &'a str) -> Cow<'a, str> {
        if self.number == 1 {
            Cow::Borrowed(source)
        } else {
            Cow::Owned(format!("{source}-{}", self.number))
        }
    }
}

fn assembly_error(code: &'static str, message: impl Into<String>) -> SemanticError {
    SemanticError::new(code, message, Vec::new())
}

/// Expand `assembly_id` (by default, the first `_pdbx_struct_assembly`) of `model` into a
/// new one-block document in which every copy of every asymmetric unit the assembly uses
/// is written out. A program that reads coordinates (Reduce3, for one) then sees the
/// contacts between copies without applying any symmetry itself.
///
/// Copies are the distinct operator combinations of the assembly's
/// `_pdbx_struct_assembly_gen` rows, numbered from 1 in order of first appearance. Copy 1
/// keeps the source `label_asym_id` and `auth_asym_id` values; copy *n* > 1 appends `-n`
/// to both (`A` becomes `A-2`). An operator expression is a list (`1,2`, `1-4`), a
/// parenthesized list, or a product of parenthesized lists (`(1-60)(61)`); a combination
/// applies its rightmost operator first. A copy of an atom that lands within 0.2 A of an
/// earlier copy of the same atom (an atom on a symmetry axis of the assembly) is written
/// once, as gemmi's `transform_to_assembly` does; an asymmetric-unit copy all of whose
/// atoms are written by one earlier copy (an ion on an axis, say) is left out, and the
/// connections to it name that copy instead.
///
/// The result keeps only what describes the assembly's coordinates:
///
/// - unchanged: `_entry`, `_audit_conform`, the entity categories (except that
///   `_entity_poly.pdbx_strand_id` lists the copies' chains), the `_chem_comp` categories,
///   `_atom_type`, and `_struct_conn_type`;
/// - once per copy: `_struct_asym`, the sequence schemes, `_atom_site` (coordinates
///   transformed, atom ids renumbered, model by model), `_atom_site_anisotrop` (tensors
///   rotated; their uncertainties become unknown when the copy is rotated), and the
///   `_struct_conn` rows whose partners are both in the copy and in the same crystal copy
///   (symmetry `1_555` or absent);
/// - left out: everything else, including the crystal cell and symmetry, which no longer
///   apply, and the assembly definitions themselves.
///
/// # Errors
///
/// Returns `PDBX_ASSEMBLY_ABSENT` when the model defines no assembly,
/// `PDBX_ASSEMBLY_UNKNOWN` for an assembly id it does not define,
/// `PDBX_ASSEMBLY_EXPRESSION` for an operator expression that cannot be read,
/// `PDBX_ASSEMBLY_OPERATOR` for an operator missing from `_pdbx_struct_oper_list`,
/// `PDBX_ASSEMBLY_ASYM` for an asymmetric unit the model does not have,
/// `PDBX_ASSEMBLY_ID_COLLISION` when a renamed chain would take an existing chain's name,
/// `PDBX_ASSEMBLY_TOO_LARGE` beyond [`MAX_ASSEMBLY_ATOM_SITES`] atom sites, and the
/// `PDBX_ITEM_*` codes for missing or untyped operator and coordinate values.
pub fn assembly_document(
    model: &PdbxModel,
    assembly_id: Option<&str>,
) -> Result<CifDocument, SemanticError> {
    let block = model
        .source_document()
        .blocks()
        .iter()
        .find(|block| block.kind() == BlockKind::Data)
        .ok_or_else(|| assembly_error("PDBX_ASSEMBLY_ABSENT", "the model has no data block"))?;
    let view = CategoryView::new(block.entries());
    let copies = copies(model, view, assembly_id)?;
    let entries = Expansion::new(block, &copies)?.entries()?;
    Ok(CifDocument::new(vec![CifBlock::data(
        block.code().unwrap_or_default().to_owned(),
        entries,
    )]))
}

fn copies(
    model: &PdbxModel,
    view: CategoryView<'_>,
    assembly_id: Option<&str>,
) -> Result<Vec<AssemblyCopy>, SemanticError> {
    let selected = match assembly_id {
        Some(id) => id.to_owned(),
        None => category_rows(view, "pdbx_struct_assembly", ["id"])
            .next()
            .and_then(|row| row.text(0))
            .or_else(|| {
                category_rows(view, "pdbx_struct_assembly_gen", ["assembly_id"])
                    .next()
                    .and_then(|row| row.text(0))
            })
            .ok_or_else(|| {
                assembly_error("PDBX_ASSEMBLY_ABSENT", "the model defines no assembly")
            })?,
    };
    let operators = operators(view)?;
    let known_asyms: HashSet<&str> = model
        .asym_units
        .iter()
        .map(|asym| asym.id.as_str())
        .collect();

    let mut copies: Vec<AssemblyCopy> = Vec::new();
    let mut by_combination: HashMap<Vec<String>, usize> = HashMap::new();
    let generators = category_rows(
        view,
        "pdbx_struct_assembly_gen",
        ["assembly_id", "oper_expression", "asym_id_list"],
    );
    for row in generators {
        if row.text(0).as_deref() != Some(selected.as_str()) {
            continue;
        }
        let expression = required_text(&row, "pdbx_struct_assembly_gen", "oper_expression", 1)?;
        let asym_list = required_text(&row, "pdbx_struct_assembly_gen", "asym_id_list", 2)?;
        let asym_ids: Vec<&str> = asym_list
            .split(',')
            .map(str::trim)
            .filter(|id| !id.is_empty())
            .collect();
        if let Some(unknown) = asym_ids.iter().find(|id| !known_asyms.contains(**id)) {
            return Err(SemanticError::new(
                "PDBX_ASSEMBLY_ASYM",
                format!(
                    "assembly {selected:?} uses asymmetric unit {unknown:?}, which the model does not have"
                ),
                row_context("pdbx_struct_assembly_gen", "asym_id_list", row.row_index()),
            ));
        }
        let combinations = operator_combinations(&expression).map_err(|message| {
            SemanticError::new(
                "PDBX_ASSEMBLY_EXPRESSION",
                format!("operator expression {expression:?}: {message}"),
                row_context(
                    "pdbx_struct_assembly_gen",
                    "oper_expression",
                    row.row_index(),
                ),
            )
        })?;
        for combination in combinations {
            let index = match by_combination.get(&combination) {
                Some(&index) => index,
                None => {
                    let transform =
                        combination
                            .iter()
                            .rev()
                            .try_fold(Transform::IDENTITY, |applied, id| {
                                operators
                                    .get(id)
                                    .map(|operator| operator.after(&applied))
                                    .ok_or_else(|| {
                                        SemanticError::new(
                                            "PDBX_ASSEMBLY_OPERATOR",
                                            format!(
                                                "operator {id:?} is not in _pdbx_struct_oper_list"
                                            ),
                                            row_context(
                                                "pdbx_struct_assembly_gen",
                                                "oper_expression",
                                                row.row_index(),
                                            ),
                                        )
                                    })
                            })?;
                    copies.push(AssemblyCopy {
                        number: copies.len() + 1,
                        transform,
                        asym_ids: Vec::new(),
                    });
                    by_combination.insert(combination, copies.len() - 1);
                    copies.len() - 1
                }
            };
            let Some(copy) = copies.get_mut(index) else {
                continue;
            };
            for id in &asym_ids {
                if !copy.has(id) {
                    copy.asym_ids.push((*id).to_owned());
                }
            }
        }
    }
    if copies.is_empty() {
        let defined = category_rows(view, "pdbx_struct_assembly", ["id"])
            .chain(category_rows(
                view,
                "pdbx_struct_assembly_gen",
                ["assembly_id"],
            ))
            .next()
            .is_some();
        return Err(if defined {
            assembly_error(
                "PDBX_ASSEMBLY_UNKNOWN",
                format!("the model defines no assembly {selected:?}"),
            )
        } else {
            assembly_error("PDBX_ASSEMBLY_ABSENT", "the model defines no assembly")
        });
    }
    Ok(copies)
}

fn operators(view: CategoryView<'_>) -> Result<HashMap<String, Transform>, SemanticError> {
    const CATEGORY: &str = "pdbx_struct_oper_list";
    const ITEMS: [&str; 13] = [
        "id",
        "matrix[1][1]",
        "matrix[1][2]",
        "matrix[1][3]",
        "vector[1]",
        "matrix[2][1]",
        "matrix[2][2]",
        "matrix[2][3]",
        "vector[2]",
        "matrix[3][1]",
        "matrix[3][2]",
        "matrix[3][3]",
        "vector[3]",
    ];
    let mut operators = HashMap::new();
    for row in category_rows(view, CATEGORY, ITEMS) {
        let id = required_text(&row, CATEGORY, "id", 0)?;
        let mut values = [0.0; 12];
        for (field, value) in values.iter_mut().enumerate() {
            let item = ITEMS.get(field + 1).copied().unwrap_or_default();
            *value = required_float(&row, CATEGORY, item, field + 1)?;
        }
        let [r11, r12, r13, t1, r21, r22, r23, t2, r31, r32, r33, t3] = values;
        operators.insert(
            id,
            Transform {
                rotation: [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]],
                translation: [t1, t2, t3],
            },
        );
    }
    Ok(operators)
}

/// The operator-id combinations of `expression`, leftmost group outermost.
fn operator_combinations(expression: &str) -> Result<Vec<Vec<String>>, String> {
    let expression = expression.trim();
    let groups: Vec<&str> = if expression.contains('(') {
        let mut groups = Vec::new();
        let mut rest = expression;
        while !rest.is_empty() {
            let inner = rest
                .strip_prefix('(')
                .and_then(|after| after.split_once(')'))
                .ok_or("expected a parenthesized list")?;
            groups.push(inner.0);
            rest = inner.1.trim_start();
        }
        groups
    } else {
        vec![expression]
    };
    let mut combinations: Vec<Vec<String>> = vec![Vec::new()];
    for group in groups {
        let ids = operator_ids(group)?;
        if combinations.len().saturating_mul(ids.len()) > MAX_COMBINATIONS {
            return Err(format!(
                "more than {MAX_COMBINATIONS} operator combinations"
            ));
        }
        combinations = combinations
            .iter()
            .flat_map(|prefix| {
                ids.iter().map(move |id| {
                    let mut combination = prefix.clone();
                    combination.push(id.clone());
                    combination
                })
            })
            .collect();
    }
    Ok(combinations)
}

/// The ids of one list: comma-separated ids and numeric ranges (`1-4`).
fn operator_ids(list: &str) -> Result<Vec<String>, String> {
    let mut ids = Vec::new();
    for item in list.split(',').map(str::trim) {
        if item.is_empty() || item.contains(char::is_whitespace) {
            return Err(format!("empty or malformed operator id {item:?}"));
        }
        let range = item.split_once('-').and_then(|(first, last)| {
            Some((first.parse::<u64>().ok()?, last.parse::<u64>().ok()?))
        });
        match range {
            Some((first, last)) if first <= last => {
                if last - first >= MAX_COMBINATIONS as u64 {
                    return Err(format!("range {item:?} is too long"));
                }
                ids.extend((first..=last).map(|number| number.to_string()));
            }
            Some(_) => return Err(format!("range {item:?} runs backwards")),
            None => ids.push(item.to_owned()),
        }
    }
    Ok(ids)
}

/// The source block and the copies being written out.
struct Expansion<'a> {
    block: &'a CifBlock,
    copies: &'a [AssemblyCopy],
    atom_sites: BlockCategory<'a>,
    /// `label_asym_id` to the `auth_asym_id` values its atoms use, in source order.
    auth_ids: HashMap<&'a str, Vec<Cow<'a, str>>>,
    /// Source coordinates of each `_atom_site` row.
    positions: Vec<[f64; 3]>,
    /// (copy index, `_atom_site` row) of atoms an earlier copy already places.
    merged: HashSet<(usize, usize)>,
    /// (copy index, `label_asym_id`) of asymmetric-unit copies whose every atom an
    /// earlier copy places, to that copy's index.
    absorbed: HashMap<(usize, String), usize>,
}

impl<'a> Expansion<'a> {
    fn new(block: &'a CifBlock, copies: &'a [AssemblyCopy]) -> Result<Self, SemanticError> {
        let atom_sites = BlockCategory::find(block.entries(), "atom_site").ok_or_else(|| {
            assembly_error(
                "PDBX_CATEGORY_REQUIRED",
                "required category \"atom_site\" is absent",
            )
        })?;
        let mut auth_ids: HashMap<&'a str, Vec<Cow<'a, str>>> = HashMap::new();
        if let (Some(label), Some(auth)) = (
            atom_sites.column("label_asym_id"),
            atom_sites.column("auth_asym_id"),
        ) {
            for row in 0..atom_sites.row_count() {
                let (Some(Cow::Borrowed(label)), Some(auth)) =
                    (atom_sites.text(row, label), atom_sites.text(row, auth))
                else {
                    continue;
                };
                let known = auth_ids.entry(label).or_default();
                if !known.contains(&auth) {
                    known.push(auth);
                }
            }
        }
        let positions = source_positions(&atom_sites)?;
        let (merged, absorbed) = coincident_copies(&atom_sites, &positions, copies)?;
        let expansion = Self {
            block,
            copies,
            atom_sites,
            auth_ids,
            positions,
            merged,
            absorbed,
        };
        expansion.check_names()?;
        Ok(expansion)
    }

    fn is_absorbed(&self, copy_index: usize, asym_id: &str) -> bool {
        self.absorbed
            .contains_key(&(copy_index, asym_id.to_owned()))
    }

    /// Every renamed chain identifier must be new.
    fn check_names(&self) -> Result<(), SemanticError> {
        for (kind, names) in [
            ("label_asym_id", self.label_names()),
            ("auth_asym_id", self.auth_names()),
        ] {
            let mut seen: HashMap<String, (usize, &str)> = HashMap::new();
            for (copy, source, name) in names {
                if let Some(&(other_copy, other_source)) = seen.get(&name)
                    && other_source != source
                {
                    return Err(assembly_error(
                        "PDBX_ASSEMBLY_ID_COLLISION",
                        format!(
                            "{kind} {name:?} of copy {copy} is already {other_source:?} of copy {other_copy}"
                        ),
                    ));
                }
                seen.entry(name).or_insert((copy, source));
            }
        }
        Ok(())
    }

    fn label_names(&self) -> Vec<(usize, &str, String)> {
        self.copies
            .iter()
            .flat_map(|copy| {
                copy.asym_ids
                    .iter()
                    .map(move |id| (copy.number, id.as_str(), copy.rename(id).into_owned()))
            })
            .collect()
    }

    fn auth_names(&self) -> Vec<(usize, &str, String)> {
        self.copies
            .iter()
            .flat_map(|copy| {
                copy.asym_ids.iter().flat_map(move |id| {
                    self.auth_ids
                        .get(id.as_str())
                        .into_iter()
                        .flatten()
                        .map(move |auth| {
                            (copy.number, auth.as_ref(), copy.rename(auth).into_owned())
                        })
                })
            })
            .collect()
    }

    fn entries(&self) -> Result<Vec<CifEntry>, SemanticError> {
        let mut entries = Vec::new();
        for category in UNCHANGED_BEFORE_SCHEMES {
            if *category == "entity_poly" {
                entries.extend(self.entity_poly()?);
            } else {
                entries.extend(self.unchanged(category));
            }
        }
        entries.extend(self.per_copy("struct_asym", &[("id", &[])], false)?);
        let strand: &[&str] = &["pdb_strand_id"];
        entries.extend(self.per_copy("pdbx_poly_seq_scheme", &[("asym_id", strand)], false)?);
        entries.extend(self.per_copy("pdbx_nonpoly_scheme", &[("asym_id", strand)], false)?);
        let branch: &[&str] = &["pdb_asym_id", "auth_asym_id"];
        entries.extend(self.per_copy("pdbx_branch_scheme", &[("asym_id", branch)], false)?);
        entries.extend(self.unchanged("atom_type"));
        let anisotrop = BlockCategory::find(self.block.entries(), "atom_site_anisotrop");
        let (atom_site, new_ids) = self.atom_site(anisotrop.is_some())?;
        entries.push(atom_site);
        if let Some(anisotrop) = anisotrop {
            entries.extend(self.anisotrop(&anisotrop, &new_ids)?);
        }
        entries.extend(self.unchanged("struct_conn_type"));
        let partners: &[(&str, &[&str])] = &[
            ("ptnr1_label_asym_id", &["ptnr1_auth_asym_id"]),
            ("ptnr2_label_asym_id", &["ptnr2_auth_asym_id"]),
            ("pdbx_ptnr3_label_asym_id", &["pdbx_ptnr3_auth_asym_id"]),
        ];
        entries.extend(self.per_copy("struct_conn", partners, true)?);
        Ok(entries)
    }

    fn unchanged(&self, category: &str) -> Vec<CifEntry> {
        self.block
            .entries()
            .iter()
            .filter(|entry| match entry {
                CifEntry::Item(item) => in_category(item.tag(), category),
                CifEntry::Loop(cif_loop) => cif_loop
                    .tags()
                    .first()
                    .is_some_and(|tag| in_category(tag, category)),
                CifEntry::Frame(_) => false,
            })
            .cloned()
            .collect()
    }

    /// `_entity_poly` with `pdbx_strand_id` listing each polymer's chains in every copy.
    fn entity_poly(&self) -> Result<Vec<CifEntry>, SemanticError> {
        let Some(table) = BlockCategory::find(self.block.entries(), "entity_poly") else {
            return Ok(Vec::new());
        };
        let (Some(entity), Some(strands)) =
            (table.column("entity_id"), table.column("pdbx_strand_id"))
        else {
            return Ok(self.unchanged("entity_poly"));
        };
        let asym_entity = self.asym_entities();
        let tags = table.tags();
        let mut columns = new_columns(tags.len(), table.row_count());
        for row in 0..table.row_count() {
            let entity_id = table.text(row, entity);
            let mut chains: Vec<String> = Vec::new();
            for (copy_index, copy) in self.copies.iter().enumerate() {
                for asym in &copy.asym_ids {
                    if asym_entity.get(asym.as_str()).map(|id| Cow::Borrowed(*id)) != entity_id
                        || self.is_absorbed(copy_index, asym)
                    {
                        continue;
                    }
                    for auth in self.auth_ids.get(asym.as_str()).into_iter().flatten() {
                        let name = copy.rename(auth).into_owned();
                        if !chains.contains(&name) {
                            chains.push(name);
                        }
                    }
                }
            }
            let joined = chains.join(",");
            for (column_index, column) in columns.iter_mut().enumerate() {
                let pushed = if column_index == strands {
                    column.push(if chains.is_empty() {
                        TextCell::Unknown
                    } else {
                        TextCell::Text(&joined)
                    })
                } else {
                    push_source(column, &table, row, column_index)
                };
                pushed.ok_or_else(too_large)?;
            }
        }
        Ok(vec![finish_loop(&tags, columns, table.row_count())])
    }

    fn asym_entities(&self) -> HashMap<&'a str, &'a str> {
        let mut entities = HashMap::new();
        if let Some(table) = BlockCategory::find(self.block.entries(), "struct_asym")
            && let (Some(id), Some(entity)) = (table.column("id"), table.column("entity_id"))
        {
            for row in 0..table.row_count() {
                if let (Some(Cow::Borrowed(id)), Some(Cow::Borrowed(entity))) =
                    (table.text(row, id), table.text(row, entity))
                {
                    entities.insert(id, entity);
                }
            }
        }
        entities
    }

    /// One row per copy of every row of `category` whose `endpoints` (a `label_asym_id`
    /// column and the `auth_asym_id` columns that go with it) name asymmetric units of
    /// the copy, renamed for the copy. Rows naming an asymmetric-unit copy an earlier copy
    /// absorbed are left out, except in `_struct_conn` (`connection`): there such an
    /// endpoint names the absorbing copy instead, and only a row whose every endpoint is
    /// absorbed is left out. In `_struct_conn` the `id` gains the copy suffix too, and only
    /// rows whose partners are in the same crystal copy are kept.
    fn per_copy(
        &self,
        category: &str,
        endpoints: &[(&str, &[&str])],
        connection: bool,
    ) -> Result<Option<CifEntry>, SemanticError> {
        let Some(table) = BlockCategory::find(self.block.entries(), category) else {
            return Ok(None);
        };
        let endpoints: Vec<(usize, Vec<usize>)> = endpoints
            .iter()
            .filter_map(|(label, auths)| {
                Some((
                    table.column(label)?,
                    auths.iter().filter_map(|item| table.column(item)).collect(),
                ))
            })
            .collect();
        let symmetry: Vec<usize> = if connection {
            ["ptnr1_symmetry", "ptnr2_symmetry"]
                .iter()
                .filter_map(|item| table.column(item))
                .collect()
        } else {
            Vec::new()
        };
        let id_column = if connection { table.column("id") } else { None };
        let tags = table.tags();
        let mut columns = new_columns(tags.len(), table.row_count());
        let mut rows = 0;
        for (copy_index, copy) in self.copies.iter().enumerate() {
            for row in 0..table.row_count() {
                // the copy each present endpoint is written as, by label column
                let mut written_as: Vec<(usize, &AssemblyCopy, bool)> = Vec::new();
                for (label, _) in &endpoints {
                    let Some(asym) = table.text(row, *label) else {
                        continue;
                    };
                    if !copy.has(&asym) {
                        written_as.clear();
                        break;
                    }
                    let target = self.absorbed.get(&(copy_index, asym.into_owned())).copied();
                    let shown = target
                        .and_then(|index| self.copies.get(index))
                        .unwrap_or(copy);
                    written_as.push((*label, shown, target.is_some()));
                }
                let absorbed = written_as
                    .iter()
                    .filter(|(_, _, absorbed)| *absorbed)
                    .count();
                if written_as.is_empty()
                    || (connection && absorbed == written_as.len())
                    || (!connection && absorbed > 0)
                {
                    continue;
                }
                if symmetry
                    .iter()
                    .any(|&column| table.text(row, column).is_some_and(|code| code != "1_555"))
                {
                    continue;
                }
                for (column_index, column) in columns.iter_mut().enumerate() {
                    let shown = endpoints
                        .iter()
                        .find(|(label, auths)| {
                            *label == column_index || auths.contains(&column_index)
                        })
                        .and_then(|(label, _)| written_as.iter().find(|(l, _, _)| l == label))
                        .map(|(_, shown, _)| *shown)
                        .or((id_column == Some(column_index)).then_some(copy));
                    let pushed = match (shown, table.text(row, column_index)) {
                        (Some(shown), Some(source)) => {
                            column.push(TextCell::Text(&shown.rename(&source)))
                        }
                        _ => push_source(column, &table, row, column_index),
                    };
                    pushed.ok_or_else(too_large)?;
                }
                rows += 1;
            }
        }
        Ok((rows > 0).then(|| finish_loop(&tags, columns, rows)))
    }

    /// `_atom_site` for every copy, model by model, and (with `track_ids`) each copy's new
    /// atom id per source row, 0 for rows outside the copy.
    fn atom_site(&self, track_ids: bool) -> Result<(CifEntry, Vec<Vec<u32>>), SemanticError> {
        let table = &self.atom_sites;
        let required = |item: &str| {
            table.column(item).ok_or_else(|| {
                assembly_error(
                    "PDBX_ITEM_REQUIRED",
                    format!("required item _atom_site.{item} is absent"),
                )
            })
        };
        let label = required("label_asym_id")?;
        let coordinates = [
            required("Cartn_x")?,
            required("Cartn_y")?,
            required("Cartn_z")?,
        ];
        let merged = &self.merged;
        let id = table.column("id");
        let auth = table.column("auth_asym_id");
        let model = table.column("pdbx_PDB_model_num");
        let tensors = [
            Tensor::find(table, "aniso_U"),
            Tensor::find(table, "aniso_B"),
        ];
        let row_count = table.row_count();

        let total = self
            .copies
            .iter()
            .enumerate()
            .map(|(copy_index, copy)| {
                (0..row_count)
                    .filter(|&row| {
                        !merged.contains(&(copy_index, row))
                            && table.text(row, label).is_some_and(|asym| copy.has(&asym))
                    })
                    .count()
            })
            .try_fold(0usize, usize::checked_add)
            .filter(|total| *total <= MAX_ASSEMBLY_ATOM_SITES)
            .ok_or_else(too_large)?;

        let mut models: Vec<Option<Cow<'_, str>>> = Vec::new();
        for row in 0..row_count {
            let number = model.and_then(|column| table.text(row, column));
            if !models.contains(&number) {
                models.push(number);
            }
        }
        let tags = table.tags();
        let mut columns = new_columns(tags.len(), total);
        let mut new_ids = if track_ids {
            vec![vec![0u32; row_count]; self.copies.len()]
        } else {
            Vec::new()
        };
        let mut next_id = 0usize;
        for number in &models {
            for (copy_index, copy) in self.copies.iter().enumerate() {
                for row in 0..row_count {
                    if model.and_then(|column| table.text(row, column)) != *number {
                        continue;
                    }
                    let Some(asym) = table.text(row, label) else {
                        continue;
                    };
                    if !copy.has(&asym) || merged.contains(&(copy_index, row)) {
                        continue;
                    }
                    let point = self.positions.get(row).copied().unwrap_or_default();
                    let moved = copy.transform.apply(point);
                    next_id += 1;
                    if let Some(slot) = new_ids.get_mut(copy_index).and_then(|ids| ids.get_mut(row))
                    {
                        *slot = u32::try_from(next_id).map_err(|_| too_large())?;
                    }
                    let tensor_cells: Vec<(usize, Option<String>)> = tensors
                        .iter()
                        .flatten()
                        .flat_map(|tensor| tensor.cells(table, row, &copy.transform))
                        .collect();
                    for (column_index, column) in columns.iter_mut().enumerate() {
                        let pushed = if let Some(axis) =
                            coordinates.iter().position(|&c| c == column_index)
                        {
                            column.push(TextCell::Text(&coordinate_text(
                                moved.get(axis).copied().unwrap_or_default(),
                            )))
                        } else if Some(column_index) == id {
                            column.push(TextCell::Text(&next_id.to_string()))
                        } else if column_index == label || Some(column_index) == auth {
                            match table.text(row, column_index) {
                                Some(source) => column.push(TextCell::Text(&copy.rename(&source))),
                                None => push_source(column, table, row, column_index),
                            }
                        } else if let Some(cell) = tensor_cell(&tensor_cells, column_index) {
                            column.push(cell)
                        } else {
                            push_source(column, table, row, column_index)
                        };
                        pushed.ok_or_else(too_large)?;
                    }
                }
            }
        }
        Ok((finish_loop(&tags, columns, total), new_ids))
    }

    /// `table` (`_atom_site_anisotrop`) for every copy's atoms, with the copies' atom ids.
    fn anisotrop(
        &self,
        table: &BlockCategory<'_>,
        new_ids: &[Vec<u32>],
    ) -> Result<Option<CifEntry>, SemanticError> {
        let (Some(atom_id), Some(source_id)) = (table.column("id"), self.atom_sites.column("id"))
        else {
            return Ok(None);
        };
        let mut source_rows: HashMap<Cow<'_, str>, usize> = HashMap::new();
        for row in 0..self.atom_sites.row_count() {
            if let Some(id) = self.atom_sites.text(row, source_id) {
                source_rows.entry(id).or_insert(row);
            }
        }
        let renamed: Vec<usize> = ["pdbx_label_asym_id", "pdbx_auth_asym_id"]
            .iter()
            .filter_map(|item| table.column(item))
            .collect();
        let tensors = [Tensor::find(table, "U"), Tensor::find(table, "B")];
        let tags = table.tags();
        let mut columns = new_columns(tags.len(), table.row_count());
        let mut rows = 0;
        for (copy_index, copy) in self.copies.iter().enumerate() {
            for row in 0..table.row_count() {
                let Some(new_id) = table
                    .text(row, atom_id)
                    .and_then(|id| source_rows.get(&id).copied())
                    .and_then(|source_row| new_ids.get(copy_index)?.get(source_row).copied())
                    .filter(|&id| id != 0)
                else {
                    continue;
                };
                let tensor_cells: Vec<(usize, Option<String>)> = tensors
                    .iter()
                    .flatten()
                    .flat_map(|tensor| tensor.cells(table, row, &copy.transform))
                    .collect();
                for (column_index, column) in columns.iter_mut().enumerate() {
                    let pushed = if column_index == atom_id {
                        column.push(TextCell::Text(&new_id.to_string()))
                    } else if let (true, Some(source)) = (
                        renamed.contains(&column_index),
                        table.text(row, column_index),
                    ) {
                        column.push(TextCell::Text(&copy.rename(&source)))
                    } else if let Some(cell) = tensor_cell(&tensor_cells, column_index) {
                        column.push(cell)
                    } else {
                        push_source(column, table, row, column_index)
                    };
                    pushed.ok_or_else(too_large)?;
                }
                rows += 1;
            }
        }
        Ok((rows > 0).then(|| finish_loop(&tags, columns, rows)))
    }
}

/// The coordinates of every `_atom_site` row.
fn source_positions(atom_sites: &BlockCategory<'_>) -> Result<Vec<[f64; 3]>, SemanticError> {
    let mut columns = [0; 3];
    for (column, item) in columns.iter_mut().zip(["Cartn_x", "Cartn_y", "Cartn_z"]) {
        *column = atom_sites.column(item).ok_or_else(|| {
            assembly_error(
                "PDBX_ITEM_REQUIRED",
                format!("required item _atom_site.{item} is absent"),
            )
        })?;
    }
    (0..atom_sites.row_count())
        .map(|row| {
            let mut point = [0.0; 3];
            for (value, &column) in point.iter_mut().zip(&columns) {
                *value = number_value(atom_sites, row, column).ok_or_else(|| {
                    SemanticError::new(
                        "PDBX_ITEM_TYPE",
                        "an _atom_site coordinate is missing or not a finite float",
                        row_context("atom_site", "Cartn_x", row),
                    )
                })?;
            }
            Ok(point)
        })
        .collect()
}

type MergedAtoms = HashSet<(usize, usize)>;
type AbsorbedCopies = HashMap<(usize, String), usize>;

/// The atoms whose copy lands within [`MERGE_DISTANCE`] of an earlier copy of the same
/// atom, and the asymmetric-unit copies all of whose atoms do so onto one earlier copy.
/// Each atom's copies are hashed into cells of the merge distance, so the work grows with
/// the number of atoms written, not with its square.
fn coincident_copies(
    atom_sites: &BlockCategory<'_>,
    positions: &[[f64; 3]],
    copies: &[AssemblyCopy],
) -> Result<(MergedAtoms, AbsorbedCopies), SemanticError> {
    let label = atom_sites.column("label_asym_id").ok_or_else(|| {
        assembly_error(
            "PDBX_ITEM_REQUIRED",
            "required item _atom_site.label_asym_id is absent",
        )
    })?;
    let mut merged = MergedAtoms::new();
    // (copy, asym) -> the one earlier copy holding all its atoms so far, or None
    let mut fate: HashMap<(usize, String), Option<usize>> = HashMap::new();
    let mut cells: HashMap<[i64; 3], Vec<(usize, [f64; 3])>> = HashMap::new();
    for (row, &point) in positions.iter().enumerate() {
        let Some(asym) = atom_sites.text(row, label) else {
            continue;
        };
        cells.clear();
        for (copy_index, copy) in copies.iter().enumerate() {
            if !copy.has(&asym) {
                continue;
            }
            let moved = copy.transform.apply(point);
            let cell = moved.map(|value| (value / MERGE_DISTANCE).floor() as i64);
            let earlier = neighbor_cells(cell)
                .filter_map(|key| cells.get(&key))
                .flatten()
                .find(|(_, other)| {
                    distance_squared(moved, *other) <= MERGE_DISTANCE * MERGE_DISTANCE
                })
                .map(|(index, _)| *index);
            let key = (copy_index, asym.clone().into_owned());
            let state = fate.entry(key).or_insert(earlier);
            if *state != earlier {
                *state = None;
            }
            match earlier {
                Some(_) => {
                    merged.insert((copy_index, row));
                }
                None => cells.entry(cell).or_default().push((copy_index, moved)),
            }
        }
    }
    let absorbed = fate
        .into_iter()
        .filter_map(|(key, state)| state.map(|index| (key, index)))
        .collect();
    Ok((merged, absorbed))
}

fn neighbor_cells(cell: [i64; 3]) -> impl Iterator<Item = [i64; 3]> {
    let [x, y, z] = cell;
    (-1..=1).flat_map(move |dx| {
        (-1..=1).flat_map(move |dy| (-1..=1).map(move |dz| [x + dx, y + dy, z + dz]))
    })
}

fn distance_squared(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(&b).map(|(p, q)| (p - q) * (p - q)).sum()
}

/// The columns of one anisotropic tensor (`U[1][1]` ... `U[2][3]`) and their
/// uncertainties.
struct Tensor {
    components: [usize; 6],
    uncertainties: Vec<usize>,
}

/// The text to write for `column` among one row's tensor `cells`: a rotated value, or
/// unknown for an uncertainty; `None` when `column` is not a rotated tensor column.
fn tensor_cell(cells: &[(usize, Option<String>)], column: usize) -> Option<TextCell<'_>> {
    cells
        .iter()
        .find(|(index, _)| *index == column)
        .map(|(_, value)| value.as_deref().map_or(TextCell::Unknown, TextCell::Text))
}

impl Tensor {
    const INDICES: [&'static str; 6] = ["[1][1]", "[2][2]", "[3][3]", "[1][2]", "[1][3]", "[2][3]"];

    fn find(table: &BlockCategory<'_>, prefix: &str) -> Option<Self> {
        let mut components = [0; 6];
        for (slot, index) in components.iter_mut().zip(Self::INDICES) {
            *slot = table.column(&format!("{prefix}{index}"))?;
        }
        let uncertainties = Self::INDICES
            .iter()
            .filter_map(|index| table.column(&format!("{prefix}{index}_esd")))
            .collect();
        Some(Self {
            components,
            uncertainties,
        })
    }

    /// The row's tensor rotated by `transform`, or left unchanged when it does not rotate
    /// or a component is missing; uncertainties of a rotated tensor become unknown.
    fn cells(
        &self,
        table: &BlockCategory<'_>,
        row: usize,
        transform: &Transform,
    ) -> Vec<(usize, Option<String>)> {
        if !transform.rotates() {
            return Vec::new();
        }
        let mut values = [0.0; 6];
        for (value, &column) in values.iter_mut().zip(&self.components) {
            match number_value(table, row, column) {
                Some(number) => *value = number,
                None => return Vec::new(),
            }
        }
        let rotated = transform.rotate_tensor(values);
        let mut cells: Vec<(usize, Option<String>)> = self
            .components
            .iter()
            .zip(rotated)
            .map(|(&column, value)| (column, Some(format!("{value:.4}"))))
            .collect();
        cells.extend(self.uncertainties.iter().map(|&column| (column, None)));
        cells
    }
}

fn number_value(table: &BlockCategory<'_>, row: usize, column: usize) -> Option<f64> {
    let number = match table.value(row, column)? {
        CifValueRef::Text(text) => parse_float(text.as_str())?,
        CifValueRef::Integer(number, _) => number as f64,
        CifValueRef::Float(number, _, _) => number,
        CifValueRef::Unknown | CifValueRef::NotApplicable => return None,
    };
    number.is_finite().then_some(number)
}

fn coordinate_text(value: f64) -> String {
    let text = format!("{value:.3}");
    if text == "-0.000" {
        "0.000".to_owned()
    } else {
        text
    }
}

fn new_columns(count: usize, rows: usize) -> Vec<TextColumnBuilder> {
    (0..count)
        .map(|_| TextColumnBuilder::with_capacity(rows))
        .collect()
}

fn push_source(
    column: &mut TextColumnBuilder,
    table: &BlockCategory<'_>,
    row: usize,
    index: usize,
) -> Option<()> {
    match table.value(row, index) {
        Some(CifValueRef::NotApplicable) => column.push(TextCell::NotApplicable),
        Some(CifValueRef::Unknown) | None => column.push(TextCell::Unknown),
        Some(_) => match table.text(row, index) {
            Some(text) => column.push(TextCell::Text(&text)),
            None => column.push(TextCell::Unknown),
        },
    }
}

fn finish_loop(tags: &[&str], columns: Vec<TextColumnBuilder>, rows: usize) -> CifEntry {
    CifEntry::Loop(CifLoop::from_columns(
        tags.iter().map(|tag| (*tag).to_owned()).collect(),
        columns.into_iter().map(TextColumnBuilder::finish).collect(),
        rows,
    ))
}

fn too_large() -> SemanticError {
    assembly_error(
        "PDBX_ASSEMBLY_TOO_LARGE",
        format!(
            "the assembly exceeds {MAX_ASSEMBLY_ATOM_SITES} atom sites or 4 GiB of text in one column"
        ),
    )
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::cif::{BlockCategory, CifDocument, parse};
    use crate::pdbx::{build_model, validate_model};

    use super::*;

    const CHEMISTRY: &str = include_str!("../../tests/fixtures/chemistry/ligand_ion_water.cif");

    /// Two assemblies: 1 is the whole asymmetric unit under the identity and a two-fold
    /// about z shifted 10 A along x; 2 is chain A alone.
    const ASSEMBLIES: &str = "
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
1 author_defined_assembly
2 'chain A alone'

loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 (1,2) A,B,C,D
2 1 A

loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 'identity operation' 1 0 0 0 0 1 0 0 0 0 1 0
2 'crystal symmetry operation' -1 0 0 10 0 -1 0 0 0 0 1 0

_cell.length_a 40.0
_symmetry.space_group_name_H-M 'P 2'

loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.type_symbol
_atom_site_anisotrop.U[1][1]
_atom_site_anisotrop.U[2][2]
_atom_site_anisotrop.U[3][3]
_atom_site_anisotrop.U[1][2]
_atom_site_anisotrop.U[1][3]
_atom_site_anisotrop.U[2][3]
_atom_site_anisotrop.U[1][1]_esd
1 C 0.1 0.2 0.3 0.01 0.02 0.03 0.001
";

    const CONNECTIONS: &str = "loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_symmetry
ZN1 metalc ATP B . PG 1_555 ZN C . ZN 1_555
ZN2 metalc ATP B . PG 1_555 ZN C . ZN 2_555
";

    fn source(with_assemblies: bool) -> String {
        let start = CHEMISTRY
            .find("loop_\n_struct_conn.id")
            .expect("fixture has struct_conn");
        let mut text = format!("{}{CONNECTIONS}", &CHEMISTRY[..start]);
        text = text.replace(
            "_entity_poly.pdbx_seq_one_letter_code_can\n1 polypeptide(L) no yes 'A(MSE)' AM",
            "_entity_poly.pdbx_seq_one_letter_code_can\n_entity_poly.pdbx_strand_id\n1 polypeptide(L) no yes 'A(MSE)' AM X",
        );
        if with_assemblies {
            text.push_str(ASSEMBLIES);
        }
        text
    }

    fn expand(text: &str, assembly_id: Option<&str>) -> Result<CifDocument, SemanticError> {
        let document = parse(text.as_bytes()).expect("fixture parses");
        let model = build_model(&document).expect("fixture builds");
        assembly_document(&model, assembly_id)
    }

    fn column(document: &CifDocument, category: &str, item: &str) -> Vec<Option<String>> {
        let block = document.blocks().first().expect("one block");
        let table = BlockCategory::find(block.entries(), category).expect("category present");
        let index = table.column(item).expect("item present");
        (0..table.row_count())
            .map(|row| table.text(row, index).map(Cow::into_owned))
            .collect()
    }

    fn texts(values: &[&str]) -> Vec<Option<String>> {
        values
            .iter()
            .map(|value| Some((*value).to_owned()))
            .collect()
    }

    #[test]
    fn writes_every_copy_with_transformed_coordinates_and_renamed_chains() {
        let document = expand(&source(true), None).expect("assembly 1 expands");

        assert_eq!(
            column(&document, "atom_site", "label_asym_id"),
            texts(&[
                "A", "A", "B", "B", "C", "D", "A-2", "A-2", "B-2", "B-2", "C-2", "D-2"
            ])
        );
        assert_eq!(
            column(&document, "atom_site", "auth_asym_id"),
            texts(&[
                "X", "X", "L", "L", "Z", "W", "X-2", "X-2", "L-2", "L-2", "Z-2", "W-2"
            ])
        );
        assert_eq!(
            column(&document, "atom_site", "id"),
            texts(&[
                "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"
            ])
        );
        // copy 2: (x, y, z) -> (10 - x, -y, z)
        assert_eq!(
            column(&document, "atom_site", "Cartn_x")[7],
            Some("6.200".to_owned())
        );
        assert_eq!(
            column(&document, "atom_site", "Cartn_y")[8],
            Some("-1.000".to_owned())
        );
        assert_eq!(
            column(&document, "atom_site", "Cartn_y")[6],
            Some("0.000".to_owned())
        );
        assert_eq!(
            column(&document, "atom_site", "label_alt_id")[8],
            Some("A".to_owned())
        );
        assert_eq!(
            column(&document, "struct_asym", "id"),
            texts(&["A", "B", "C", "D", "A-2", "B-2", "C-2", "D-2"])
        );
        assert_eq!(
            column(&document, "pdbx_nonpoly_scheme", "pdb_strand_id"),
            texts(&["L", "Z", "W", "L-2", "Z-2", "W-2"])
        );
        assert_eq!(
            column(&document, "entity_poly", "pdbx_strand_id"),
            texts(&["X,X-2"])
        );
    }

    #[test]
    fn the_result_is_a_valid_model_without_crystal_or_assembly_categories() {
        let document = expand(&source(true), None).expect("assembly 1 expands");
        let block = document.blocks().first().expect("one block");
        for category in [
            "cell",
            "symmetry",
            "pdbx_struct_assembly",
            "pdbx_struct_assembly_gen",
            "pdbx_struct_oper_list",
        ] {
            assert!(
                BlockCategory::find(block.entries(), category).is_none(),
                "{category} kept"
            );
        }
        let model = build_model(&document).expect("the assembly builds a model");
        assert_eq!(model.atom_site_count(), 12);
        assert_eq!(model.asym_unit_count(), 8);
        let report = validate_model(&model);
        assert!(
            report.diagnostics().is_empty(),
            "{:?}",
            report.diagnostics()
        );
    }

    #[test]
    fn keeps_connections_within_a_crystal_copy_once_per_copy() {
        let document = expand(&source(true), None).expect("assembly 1 expands");
        assert_eq!(
            column(&document, "struct_conn", "id"),
            texts(&["ZN1", "ZN1-2"])
        );
        assert_eq!(
            column(&document, "struct_conn", "ptnr2_label_asym_id"),
            texts(&["C", "C-2"])
        );
    }

    #[test]
    fn rotates_anisotropic_tensors_and_drops_their_uncertainties() {
        let document = expand(&source(true), None).expect("assembly 1 expands");
        assert_eq!(
            column(&document, "atom_site_anisotrop", "id"),
            texts(&["1", "7"])
        );
        // R = diag(-1, -1, 1): U13 and U23 change sign
        assert_eq!(
            column(&document, "atom_site_anisotrop", "U[1][3]")[1],
            Some("-0.0200".to_owned())
        );
        assert_eq!(
            column(&document, "atom_site_anisotrop", "U[2][3]")[1],
            Some("-0.0300".to_owned())
        );
        assert_eq!(
            column(&document, "atom_site_anisotrop", "U[1][2]")[1],
            Some("0.0100".to_owned())
        );
        assert_eq!(
            column(&document, "atom_site_anisotrop", "U[1][1]_esd"),
            vec![Some("0.001".to_owned()), None]
        );
    }

    #[test]
    fn selects_an_assembly_by_id() {
        let document = expand(&source(true), Some("2")).expect("assembly 2 expands");
        assert_eq!(
            column(&document, "atom_site", "label_asym_id"),
            texts(&["A", "A"])
        );
        assert_eq!(
            column(&document, "atom_site", "Cartn_x")[1],
            Some("3.800".to_owned())
        );
        let block = document.blocks().first().expect("one block");
        assert!(BlockCategory::find(block.entries(), "struct_conn").is_none());
    }

    #[test]
    fn writes_copies_model_by_model() {
        let text = source(true).replace(
            "HETATM 6 O O . HOH D 4 . ? 12.000 0.000 0.000 1.00 25.00 ? 900 HOH W O 1",
            "HETATM 6 O O . HOH D 4 . ? 12.000 0.000 0.000 1.00 25.00 ? 900 HOH W O 1\n\
             ATOM 7 C CA . ALA A 1 1 ? 0.500 0.000 0.000 1.00 10.00 ? 10 ALA X CA 2",
        );
        let document = parse(text.as_bytes()).expect("fixture parses");
        let model = build_model(&document).expect("fixture builds");
        let document = assembly_document(&model, Some("2")).expect("assembly 2 expands");
        assert_eq!(
            column(&document, "atom_site", "pdbx_PDB_model_num"),
            texts(&["1", "1", "2"])
        );
    }

    #[test]
    fn reports_assemblies_it_cannot_build() {
        let code = |text: &str, id: Option<&str>| {
            expand(text, id).map(|_| ()).map_err(|error| error.code())
        };
        assert_eq!(code(&source(false), None), Err("PDBX_ASSEMBLY_ABSENT"));
        assert_eq!(code(&source(true), Some("9")), Err("PDBX_ASSEMBLY_UNKNOWN"));
        let broken = |from: &str, to: &str| source(true).replace(from, to);
        assert_eq!(
            code(&broken("1 (1,2) A,B,C,D", "1 (1,2 A,B,C,D"), None),
            Err("PDBX_ASSEMBLY_EXPRESSION")
        );
        assert_eq!(
            code(&broken("1 (1,2) A,B,C,D", "1 (1,3) A,B,C,D"), None),
            Err("PDBX_ASSEMBLY_OPERATOR")
        );
        assert_eq!(
            code(&broken("1 (1,2) A,B,C,D", "1 (1,2) A,E"), None),
            Err("PDBX_ASSEMBLY_ASYM")
        );
        assert_eq!(
            code(
                &broken(
                    "2 'crystal symmetry operation' -1 0 0 10",
                    "2 'crystal symmetry operation' -1 0 0 ?"
                ),
                None
            ),
            Err("PDBX_ITEM_REQUIRED")
        );
    }

    #[test]
    fn rejects_a_renamed_chain_that_already_exists() {
        let text = source(true)
            .replace("A,B,C,D", "A,B,C,D-2")
            .replace("\nD 4 .", "\nD-2 4 .")
            .replace("HOH D 4", "HOH D-2 4")
            .replace("D 4 HOH 1", "D-2 4 HOH 1");
        let document = parse(text.as_bytes()).expect("fixture parses");
        let model = build_model(&document).expect("fixture builds");
        let renamed = assembly_document(&model, None)
            .map(|_| ())
            .map_err(|error| error.code());
        assert_eq!(renamed, Ok(()), "D-2 alone does not collide");
        let text = text.replace("A,B,C,D-2", "A,B,C,D-2,D").replace(
            "loop_\n_struct_asym.id\n_struct_asym.entity_id\n_struct_asym.details\n",
            "loop_\n_struct_asym.id\n_struct_asym.entity_id\n_struct_asym.details\nD 4 .\n",
        );
        let document = parse(text.as_bytes()).expect("fixture parses");
        let model = build_model(&document).expect("fixture builds");
        let collided = assembly_document(&model, None)
            .map(|_| ())
            .map_err(|error| error.code());
        assert_eq!(collided, Err("PDBX_ASSEMBLY_ID_COLLISION"));
    }

    #[test]
    fn writes_an_atom_on_a_symmetry_axis_once() {
        // the two-fold maps (5, 0, z) onto itself: the zinc's copies coincide
        let text = source(true).replace(
            "HETATM 5 Zn ZN . ZN C 3 . ? 9.000 0.000 0.000",
            "HETATM 5 Zn ZN . ZN C 3 . ? 5.000 0.000 0.000",
        );
        let document = expand(&text, None).expect("assembly 1 expands");

        assert_eq!(
            column(&document, "atom_site", "label_asym_id"),
            texts(&[
                "A", "A", "B", "B", "C", "D", "A-2", "A-2", "B-2", "B-2", "D-2"
            ])
        );
        assert_eq!(
            column(&document, "struct_asym", "id"),
            texts(&["A", "B", "C", "D", "A-2", "B-2", "D-2"])
        );
        assert_eq!(
            column(&document, "pdbx_nonpoly_scheme", "pdb_strand_id"),
            texts(&["L", "Z", "W", "L-2", "W-2"])
        );
        // the second copy's ATP coordinates the one zinc
        assert_eq!(
            column(&document, "struct_conn", "id"),
            texts(&["ZN1", "ZN1-2"])
        );
        assert_eq!(
            column(&document, "struct_conn", "ptnr1_label_asym_id"),
            texts(&["B", "B-2"])
        );
        assert_eq!(
            column(&document, "struct_conn", "ptnr2_label_asym_id"),
            texts(&["C", "C"])
        );
        let model = build_model(&document).expect("the assembly builds a model");
        let report = validate_model(&model);
        assert!(
            report.diagnostics().is_empty(),
            "{:?}",
            report.diagnostics()
        );
    }

    #[test]
    fn expands_operator_expressions_leftmost_group_outermost() {
        let ids = |expression: &str| {
            operator_combinations(expression).map(|combinations| {
                combinations
                    .into_iter()
                    .map(|combination| combination.join("x"))
                    .collect::<Vec<_>>()
            })
        };
        assert_eq!(ids("1"), Ok(vec!["1".to_owned()]));
        assert_eq!(
            ids("1,2,X0"),
            Ok(vec!["1".to_owned(), "2".to_owned(), "X0".to_owned()])
        );
        assert_eq!(
            ids("(1-3)"),
            Ok(vec!["1".to_owned(), "2".to_owned(), "3".to_owned()])
        );
        assert_eq!(
            ids("(1,2) (3-4)"),
            Ok(vec![
                "1x3".to_owned(),
                "1x4".to_owned(),
                "2x3".to_owned(),
                "2x4".to_owned()
            ])
        );
        assert!(ids("(1,2").is_err());
        assert!(ids("(3-1)").is_err());
        assert!(ids("1, ,2").is_err());
    }

    #[test]
    fn a_combination_applies_its_rightmost_operator_first() {
        let shift = Transform {
            rotation: Transform::IDENTITY.rotation,
            translation: [1.0, 0.0, 0.0],
        };
        let turn = Transform {
            rotation: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0; 3],
        };
        // turn after shift: (0,0,0) -> (1,0,0) -> (0,1,0)
        assert_eq!(turn.after(&shift).apply([0.0; 3]), [0.0, 1.0, 0.0]);
    }
}

//! Atom84 molecular targets from one source table and one linked chemical product.
//!
//! The view reads the source and the product once. Python perceives the bounded
//! ligand domains with RDKit; `finish` then builds the legacy Atom32 tokens, the
//! source snapshot and the typed Atom84 rows without another document read.

mod arrays;
mod json;
mod order;
mod pack;
mod product;
mod snapshot;

use super::product::ProductTarget;
use super::uncertainty::TargetUncertainty;
use super::{HydrogenParents, LigandFeatures, PackOptions, PreparedSource, Schema, TokenTable};
use crate::cif::{BlockCategory, CifDocument};
use crate::pdbx::{ConformerPolicy, PdbxModel, PreparedChemistry, SemanticError, ShortText};
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

pub(crate) use arrays::Array;
pub(crate) use product::Columns;

/// Coordinate slots in one typed molecular token.
pub(crate) const SLOTS: usize = 84;

/// A failure carries the Python exception class of the step it replaces.
pub(crate) enum Failure {
    Value(String),
    Chemistry(SemanticError),
    Index(String),
}

impl From<String> for Failure {
    fn from(message: String) -> Self {
        Self::Value(message)
    }
}

impl From<&str> for Failure {
    fn from(message: &str) -> Self {
        Self::Value(message.to_owned())
    }
}

impl From<SemanticError> for Failure {
    fn from(error: SemanticError) -> Self {
        Self::Chemistry(error)
    }
}

/// Ligand-domain perception rules and the cached templates of canonical residues.
#[derive(Deserialize)]
pub(crate) struct Perception {
    pub protein: HashSet<String>,
    pub nucleotides: HashSet<String>,
    pub unknown: HashSet<String>,
    pub water: HashSet<String>,
    pub numbers: HashMap<String, i16>,
    pub templates: HashMap<String, HashMap<String, [i64; 4]>>,
    pub columns: Columns,
}

/// Fixed names and tables of the versioned Atom84 contract.
#[derive(Deserialize)]
pub(crate) struct MolecularSchema {
    pub atom_names: Vec<String>,
    pub hydrogen_start: usize,
    pub index: Vec<[usize; 32]>,
    pub defined: Vec<[bool; 32]>,
    pub required_heavy: Vec<[bool; 14]>,
    pub nucleotides: Vec<String>,
    pub unknown_components: HashSet<String>,
    pub backbone: Vec<String>,
    pub bases: HashMap<String, Vec<String>>,
    pub sugar_hydrogens: Vec<String>,
    pub atomic_numbers: HashMap<String, i16>,
    pub exclusions: Vec<String>,
    pub schema_json: String,
    pub schema_sha256: String,
    pub perception: Perception,
}

impl MolecularSchema {
    pub(crate) fn read(text: &str) -> Result<Self, String> {
        let schema: Self = serde_json::from_str(text).map_err(|error| error.to_string())?;
        if schema.atom_names.len() != SLOTS
            || schema.index.len() != 22
            || schema.defined.len() != 22
            || schema.required_heavy.len() != 20
            || schema.backbone.len() > SLOTS
            || schema.hydrogen_start + 3 > SLOTS
            || schema.exclusions.len() != 9
            || schema.index.iter().flatten().any(|&slot| slot >= SLOTS)
        {
            return Err("The molecular schema does not describe 84 typed slots".into());
        }
        Ok(schema)
    }
}

/// The packer tables used by the legacy Atom32 pass.
pub(crate) struct Packer<'a> {
    pub schema: &'a Schema,
    pub features: &'a LigandFeatures,
    pub hydrogen_parents: &'a HydrogenParents,
    pub table: &'a TokenTable,
    pub excluded: &'a [String],
    pub molecular: &'a MolecularSchema,
}

/// Options that differ between the heavy index and the protonated preparation.
pub(crate) struct Options<'a> {
    pub complete: bool,
    pub sasa: bool,
    pub uncertainties: Option<&'a str>,
    pub expected: Option<&'a str>,
    /// Return both key lists as compact JSON text.
    pub keys_json: bool,
    /// Return only the token layout the heavy index stores and hashes.
    pub index: bool,
    /// Also return the component, element and name labels as `S8` arrays, which Python
    /// reads without decoding the label JSON.
    pub reference_labels: bool,
}

/// The atoms of one perception graph and its bonds as (first, second, order).
pub(crate) type Graph<'a> = (Vec<GraphAtom<'a>>, &'a [(usize, usize, u8)]);

/// One source table, one linked product and the ligand domains to perceive.
pub(crate) struct View {
    pub(crate) source: PreparedSource,
    pub(crate) document: Arc<CifDocument>,
    chem_data: PathBuf,
    /// A view of Reduce3's output shares the product of the view of its input.
    product: Arc<Product>,
    /// The row of this document that holds each product atom.
    rows: Vec<Option<usize>>,
    observed: Vec<ShortText>,
}

/// One linked product, read without observed hydrogens, and the ligand domains drawn from it.
struct Product {
    /// The model the product was read from, when its atoms name rows of the whole model.
    /// Their coordinates and source rows belong to that model.
    model: Option<Arc<PdbxModel>>,
    /// The hash of the molecular schema whose perception rules drew the domains.
    schema: String,
    chains: HashMap<String, String>,
    chemistry: PreparedChemistry,
    residues: product::Residues,
    features: Vec<[i16; 16]>,
    domains: Vec<product::Domain>,
    graphs: Vec<product::LocalGraph>,
}

/// One perception atom as Python's RDKit graph builder reads it.
pub(crate) struct GraphAtom<'a> {
    pub atom: &'a crate::pdbx::ProductAtom,
    pub observed: &'a str,
    pub source_row: Option<usize>,
    pub cap: bool,
}

/// A list of token keys, as compact JSON text or as owned lists.
pub(crate) enum Keys {
    Text(String),
    Lists(Vec<Vec<String>>),
}

/// Typed outputs in the key order of the Python packer.
pub(crate) struct Finished {
    pub arrays: Vec<(&'static str, Array)>,
    pub keys: Keys,
    pub required: Keys,
    pub stats: pack::Stats,
}

/// Legacy keys as compact JSON, written as `json.dumps` writes them.
fn legacy_key_text<'a>(
    keys: impl Iterator<Item = &'a [crate::pdbx::ShortText; 5]>,
    count: usize,
) -> String {
    let mut out = String::with_capacity(count * 28);
    out.push('[');
    // Tokens of one residue share the first four fields, so keep their text while they repeat.
    let mut previous: Option<&[crate::pdbx::ShortText; 5]> = None;
    let mut prefix = String::new();
    for (index, key) in keys.enumerate() {
        if index > 0 {
            out.push(',');
        }
        if previous.is_none_or(|last| !super::record::same_residue(last, key)) {
            prefix.clear();
            prefix.push('[');
            for field in &key[..4] {
                json::string(&mut prefix, field);
                prefix.push(',');
            }
        }
        previous = Some(key);
        out.push_str(&prefix);
        json::string(&mut out, &key[4]);
        out.push(']');
    }
    out.push(']');
    out
}

/// The `label_atom_id` text of the row that holds each product atom, or the product name.
fn observed(
    atoms: &[crate::pdbx::ProductAtom],
    rows: &[Option<usize>],
    table: Option<&BlockCategory<'_>>,
    column: Option<usize>,
) -> Vec<ShortText> {
    atoms
        .iter()
        .zip(rows)
        .map(|(atom, row)| observed_name(atom, *row, table, column))
        .collect()
}

/// The `label_atom_id` text of `row`, or the name of `atom` without a row.
fn observed_name(
    atom: &crate::pdbx::ProductAtom,
    row: Option<usize>,
    table: Option<&BlockCategory<'_>>,
    column: Option<usize>,
) -> ShortText {
    match (row, table, column) {
        (Some(row), Some(table), Some(column)) => table
            .text(row, column)
            .map_or_else(|| ShortText::new(""), |text| ShortText::new(&text)),
        _ => atom.atom_id.clone(),
    }
}

impl View {
    /// Read a view; `built` is the registry-free PDBx model of `document` when the caller
    /// already has it. A model of other blocks is ignored. When `base` is a view whose product
    /// this document gives again, as Reduce3 gives the product of its input, the view shares it.
    pub(crate) fn read_with_model(
        document: Arc<CifDocument>,
        built: Option<Arc<PdbxModel>>,
        model: i64,
        chem_data: &Path,
        schema: &MolecularSchema,
        base: Option<&Self>,
    ) -> Result<Self, Failure> {
        let built =
            built.filter(|full| std::ptr::eq(full.source_document().blocks(), document.blocks()));
        let full = match built {
            Some(full) => full,
            None => match crate::pdbx::build_model(&document) {
                Ok(full) => Arc::new(full),
                Err(error) => {
                    // The source read comes first, so its error wins.
                    PreparedSource::read(&document, model)?;
                    return Err(error.into());
                }
            },
        };
        // The model has decoded the atom rows, so the source reads them from its sites.
        let source = PreparedSource {
            source: super::source::read_with_sites(&document, model, full.loop_sites())?,
            model,
        };
        // The full model already holds the decoded sites and connections.
        let selection = crate::pdbx::select_conformer_with_policy(
            &full,
            Some(model),
            0,
            ConformerPolicy::Occupancy,
        )?;
        let table = document
            .blocks()
            .first()
            .and_then(|block| BlockCategory::find(block.entries(), "atom_site"));
        let count = table.as_ref().map_or(0, BlockCategory::row_count);
        let column = table
            .as_ref()
            .and_then(|table| table.column("label_atom_id"));
        let whole = selection.atom_rows.iter().copied().eq(0..count);
        let shared = base.filter(|base| {
            whole
                && base.source.model == model
                && base.chem_data == chem_data
                && base.product.schema == schema.schema_sha256
        });
        // A failure leaves the product to the full read, which reports errors in its own order.
        if let Some(base) = shared
            && let Some(own) = &base.product.model
            && product::polymer_chains(&document).is_ok_and(|chains| chains == base.product.chains)
            && let Some(rows) =
                crate::pdbx::product_rows_in(&base.product.chemistry, own, &full, model)
        {
            let observed = observed(
                base.product.chemistry.atoms(),
                &rows,
                table.as_ref(),
                column,
            );
            return Ok(Self {
                source,
                document,
                chem_data: chem_data.to_owned(),
                product: Arc::clone(&base.product),
                rows,
                observed,
            });
        }
        let selected = if whole {
            full
        } else {
            Arc::new(crate::pdbx::build_model(
                &crate::pdbx::selected_atom_document(&document, &selection.atom_rows)?,
            )?)
        };
        let registry = crate::pdbx::chemical_library_registry(&selected, chem_data)
            .map_err(|message| SemanticError::new("PDBX_LINK_PRODUCT", message, Vec::new()))?;
        let chemistry =
            crate::pdbx::prepare_product(&selected, Some(&registry), model, None, false)?;
        let atoms = chemistry.atoms();
        let rows: Vec<_> = atoms
            .iter()
            .map(|atom| {
                atom.source_row
                    .and_then(|row| selection.atom_rows.get(row).copied())
            })
            .collect();
        let observed = observed(atoms, &rows, table.as_ref(), column);
        let chains = product::polymer_chains(&document)?;
        let perception = &schema.perception;
        let residues = product::residues(atoms, &chains, perception);
        let features = product::native_features(atoms, &residues, perception, perception.columns);
        let domains = product::domains(atoms, chemistry.bonds(), &residues, &chains, perception);
        let incident = product::incident(atoms, chemistry.bonds());
        let mut local = vec![u32::MAX; atoms.len()];
        let graphs = domains
            .iter()
            .map(|domain| {
                product::local_graph(
                    atoms,
                    chemistry.bonds(),
                    &domain.atoms,
                    &incident,
                    &mut local,
                )
            })
            .collect();
        let product = Product {
            model: whole.then_some(selected),
            schema: schema.schema_sha256.clone(),
            chains,
            chemistry,
            residues,
            features,
            domains,
            graphs,
        };
        Ok(Self {
            source,
            document,
            chem_data: chem_data.to_owned(),
            product: Arc::new(product),
            rows,
            observed,
        })
    }

    /// Read the view of `document`, which Reduce3 wrote from the document of `source`, from
    /// the view `base` of that document and the values Reduce3 wrote (see
    /// [`super::source::Copied`]), without building the model of `document`. Returns `None`
    /// when [`Self::read_with_model`] might give another view or an error.
    #[cfg(feature = "reduce3")]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn read_copied(
        document: Arc<CifDocument>,
        source: &Arc<PdbxModel>,
        rows: &[u32],
        made: &[u64],
        model: i64,
        chem_data: &Path,
        schema: &MolecularSchema,
        base: &Self,
    ) -> Option<Self> {
        let own = base
            .product
            .model
            .as_ref()
            .filter(|own| Arc::ptr_eq(own, source))?;
        if base.source.model != model
            || base.chem_data != chem_data
            || base.product.schema != schema.schema_sha256
        {
            return None;
        }
        let copied =
            super::source::Copied::new(&document, own, rows, made, model, &base.source.source)?;
        if product::polymer_chains(&document).ok()? != base.product.chains {
            return None;
        }
        let atoms = base.product.chemistry.atoms();
        let rows = crate::pdbx::product_rows_with(
            &base.product.chemistry,
            own,
            model,
            copied.heavy(),
            |row| copied.identity(row),
            |row| copied.position(row),
        )?;
        let table = document
            .blocks()
            .first()
            .and_then(|block| BlockCategory::find(block.entries(), "atom_site"));
        let column = table
            .as_ref()
            .and_then(|table| table.column("label_atom_id"));
        // A row that copies the name of the row of `base` holds the same text.
        let observed = atoms
            .iter()
            .zip(&rows)
            .zip(base.rows.iter().zip(&base.observed))
            .map(|((atom, row), (base_row, name))| match (*row, *base_row) {
                (Some(row), Some(base_row)) if copied.copies_name(row, base_row) => name.clone(),
                _ => observed_name(atom, *row, table.as_ref(), column),
            })
            .collect();
        let source = PreparedSource {
            source: copied.source,
            model,
        };
        Some(Self {
            source,
            document,
            chem_data: chem_data.to_owned(),
            product: Arc::clone(&base.product),
            rows,
            observed,
        })
    }

    /// Atoms and bonds of each perception graph; caps follow the product atoms.
    pub(crate) fn perception_graphs(&self) -> Vec<Graph<'_>> {
        let atoms = self.product.chemistry.atoms();
        self.product
            .graphs
            .iter()
            .map(|graph| {
                let mut result: Vec<_> = graph
                    .ordered
                    .iter()
                    .map(|&old| GraphAtom {
                        atom: &atoms[old],
                        observed: &self.observed[old],
                        source_row: self.rows[old],
                        cap: false,
                    })
                    .collect();
                result.extend(graph.caps.iter().map(|&old| GraphAtom {
                    atom: &atoms[old],
                    observed: &self.observed[old],
                    source_row: None,
                    cap: true,
                }));
                (result, graph.bonds.as_slice())
            })
            .collect()
    }

    pub(crate) fn ligand_components(&self) -> Vec<String> {
        self.source.components()
    }

    pub(crate) fn missing_hydrogen_dictionaries(
        &self,
        chem_data: &Path,
    ) -> Result<Vec<String>, String> {
        self.source.missing_hydrogen_dictionaries(chem_data)
    }

    /// Merge RDKit domain features; caps never change native H counts or charges.
    fn merge(
        &self,
        perceived: &[Vec<u8>],
        columns: Columns,
    ) -> Result<(Vec<[i16; 16]>, String), Failure> {
        if perceived.len() != self.product.graphs.len() || columns.count != 16 {
            return Err("Perceived features must describe every ligand domain".into());
        }
        let product = &self.product;
        let atoms = product.chemistry.atoms();
        let mut features = product.features.clone();
        let mut audit = String::new();
        for (index, (graph, values)) in product.graphs.iter().zip(perceived).enumerate() {
            let total = graph.ordered.len() + graph.caps.len();
            if values.len() != total * 32 {
                return Err("Perceived features must describe every ligand domain atom".into());
            }
            for (new, &old) in graph.ordered.iter().enumerate() {
                let row = &mut features[old];
                for (column, value) in row.iter_mut().enumerate() {
                    let at = (new * 16 + column) * 2;
                    *value = i16::from_ne_bytes([values[at], values[at + 1]]);
                }
                if graph.boundary.binary_search(&old).is_ok() && row[columns.atomic_number] != 0 {
                    row[columns.total_hydrogen] = atoms[old].hydrogen_count.unwrap_or(0) as i16;
                    row[columns.formal_charge] = atoms[old].formal_charge.unwrap_or(0) as i16;
                }
                if graph.uncertain.binary_search(&old).is_ok() {
                    row[columns.hybridization] = 0;
                }
            }
            if index > 0 {
                audit.push(',');
            }
            audit.push_str("{\"boundary_atoms\":");
            audit.push_str(&graph.boundary.len().to_string());
            audit.push_str(",\"boundary_hybridization_unknown_atoms\":");
            audit.push_str(&graph.uncertain.len().to_string());
            audit.push_str(",\"host_residues\":[");
            for (position, &residue) in product.domains[index].hosts.iter().enumerate() {
                let atom = &atoms[product.residues.first[residue as usize]];
                if position > 0 {
                    audit.push(',');
                }
                audit.push('[');
                json::string(&mut audit, &atom.asym_id);
                audit.push(',');
                json::string(&mut audit, &atom.auth_sequence_id);
                audit.push(',');
                match &atom.insertion_code {
                    Some(code) => json::string(&mut audit, code),
                    None => audit.push_str("null"),
                }
                audit.push(',');
                json::string(&mut audit, &atom.component_id);
                audit.push(']');
            }
            audit.push_str("],\"product_atoms\":");
            audit.push_str(&graph.ordered.len().to_string());
            audit.push_str(",\"rdkit_atoms\":");
            audit.push_str(&total.to_string());
            audit.push('}');
        }
        let provenance = format!(
            "{{\"boundary\":\"perception-only-H-saturation;native-H-charge-retained;\
            boundary-hybridization-unknown\",\"domains\":[{audit}],\"policy\":\"ligand-and-direct-host-residue.v1\",\
            \"polymer_targets\":\"native-linked-product-H-charge-and-unchanged-CCD-hybridization\"}}"
        );
        Ok((features, provenance))
    }

    /// Product targets for heavy atoms, in product order.
    ///
    /// Atoms of one product residue share the first four key fields, so the targets of a
    /// residue borrow them from its first atom. `product::apply` then finds a repeated
    /// residue by address without reading the text.
    fn targets(
        &self,
        features: &[[i16; 16]],
        complete: &[bool],
        modified: &[bool],
        columns: Columns,
    ) -> Vec<ProductTarget<&str>> {
        let atoms = self.product.chemistry.atoms();
        let mut previous: Option<(u32, [&str; 4])> = None;
        // Reserve for every atom: a filtered iterator gives no length, and growing a list this
        // large copies it each time.
        let mut targets = Vec::with_capacity(atoms.len());
        for (index, atom) in atoms.iter().enumerate() {
            if !atom.active
                || atom.source_row.is_none()
                || atom.element == "H"
                || atom.element == "D"
                || atom.element == "T"
            {
                continue;
            }
            let id = self.product.residues.of_atom[index];
            let prefix = match previous {
                Some((last, prefix)) if last == id => prefix,
                _ => {
                    let prefix = [
                        atom.asym_id.as_str(),
                        atom.auth_sequence_id.as_str(),
                        atom.insertion_code
                            .as_deref()
                            .filter(|code| !code.is_empty())
                            .unwrap_or("?"),
                        atom.component_id.as_str(),
                    ];
                    previous = Some((id, prefix));
                    prefix
                }
            };
            let residue = id as usize;
            targets.push(ProductTarget::new(
                [
                    prefix[0],
                    prefix[1],
                    prefix[2],
                    prefix[3],
                    self.observed[index].as_str(),
                ],
                atom.sequence_id,
                atom.formal_charge,
                i64::from(features[index][columns.hybridization]),
                atom.hydrogen_count,
                complete[residue],
                modified[residue],
            ));
        }
        targets
    }
}

/// Build legacy tokens, the source snapshot and typed Atom84 rows.
pub(crate) fn finish(
    view: &View,
    perceived: &[Vec<u8>],
    packer: &Packer<'_>,
    options: &Options<'_>,
) -> Result<Finished, Failure> {
    let columns = packer.molecular.perception.columns;
    let (features, provenance) = view.merge(perceived, columns)?;
    let atoms = view.product.chemistry.atoms();
    let bonds = view.product.chemistry.bonds();
    let (complete, modified) =
        product::completeness(atoms, bonds, &view.product.residues, &features, columns);
    let targets = view.targets(&features, &complete, &modified, columns);
    let mut legacy = view.source.pack(
        packer.schema,
        packer.features,
        packer.hydrogen_parents,
        packer.table,
        packer.excluded,
        &PackOptions {
            model: view.source.model,
            complete: options.complete,
            topology_policy: "observed_connectivity",
            protein_only: false,
            chem_data: Some(&view.chem_data),
        },
    )?;
    super::product::apply(&mut legacy, packer.schema, packer.table, &targets)?;
    let uncertainties: Vec<TargetUncertainty> = match options.uncertainties {
        Some(text) => serde_json::from_str(text).map_err(|error| error.to_string())?,
        None => Vec::new(),
    };
    super::uncertainty::apply(&mut legacy, packer.schema, &uncertainties)?;
    legacy.protein_gaps();
    legacy.validate(packer.schema, packer.table)?;
    let mut source = snapshot::Snapshot::read(
        &view.source.source,
        atoms,
        bonds,
        &view.rows,
        &features,
        &view.product.residues,
        &complete,
        columns,
        provenance,
        &packer.molecular.unknown_components,
    );
    // Disabled SASA changes only its provenance text, which the index output omits.
    if options.sasa || !options.index {
        source.sasa(
            &view.source.source,
            &view.document,
            view.source.model,
            options.sasa,
            &packer.molecular.perception.protein,
        );
    }
    source.mask_uncertain_hydrogens(&view.source.source, &uncertainties);
    source.mask_unknown_chemistry();
    let mut packed = pack::pack(
        &legacy,
        &view.source.source,
        &source,
        packer.molecular,
        packer.table,
        options.index,
    )?;
    let ordered = match options.expected {
        Some(text) => Some(order::indexed(&mut packed, text)?),
        None => None,
    };
    let retained = || {
        legacy
            .keys
            .iter()
            .zip(&packed.legacy_token_index)
            .filter(|(_, token)| **token >= 0)
            .map(|(key, _)| key)
    };
    let (keys, required) = if options.keys_json {
        (
            Keys::Text(ordered.unwrap_or_else(|| order::key_text(&packed.keys))),
            Keys::Text(legacy_key_text(retained(), legacy.keys.len())),
        )
    } else {
        (
            Keys::Lists(
                packed
                    .keys
                    .iter()
                    .map(|key| key.map(str::to_owned).to_vec())
                    .collect(),
            ),
            Keys::Lists(
                retained()
                    .map(|key| key.iter().map(String::from).collect())
                    .collect(),
            ),
        )
    };
    let arrays = if options.index {
        packed.index_arrays()
    } else {
        let mut arrays = packed.row_arrays();
        arrays.extend(source.arrays(
            &view.source.source,
            &packer.molecular.exclusions,
            options.reference_labels,
        ));
        arrays.extend(packed.target_arrays(packer.molecular));
        arrays
    };
    Ok(Finished {
        arrays,
        keys,
        required,
        stats: packed.stats,
    })
}

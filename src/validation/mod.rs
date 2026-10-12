//! Native chemical interpretation, restraint checks and Probe contacts.
//!
//! Reports retain unsupported coverage. Callers must run the reference backend
//! for those checks; findings from a partial graph do not establish completion.

mod atoms;
mod connections;
mod dictionaries;
mod report;
mod source;

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

use reduce3::cifsource::CifSource;
use reduce3::interp::{FlatAtoms, validation};
use reduce3::monlib::MonLib;
use serde_json::Value;

use crate::cif::CifDocument;
use crate::geometry::RestraintGraph;
use crate::reduce::BlockSource;

/// Options shared by chemical targets and contact checks.
#[derive(Clone, Debug)]
pub struct Options {
    /// Use nuclear hydrogen targets and radii.
    pub nuclear: bool,
    /// Apply the conformation-dependent backbone library.
    pub use_cdl: bool,
    /// Strict standardized-deviation cutoff.
    pub sigma_cutoff: f64,
    /// Include the independently built graph for parity tests.
    pub include_graph: bool,
    /// Link IDs and two source atom rows that select their residue sides.
    pub restraint_links: Vec<(String, u32, u32)>,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            nuclear: false,
            use_cdl: true,
            sigma_cutoff: 5.0,
            include_graph: false,
            restraint_links: Vec::new(),
        }
    }
}

type Libraries = HashMap<PathBuf, Arc<MonLib>>;
static LIBRARIES: OnceLock<Mutex<Libraries>> = OnceLock::new();

/// Load and cache validation dictionaries by canonical root.
///
/// The caller must keep dictionary files immutable during a process.
pub fn load_library(root: &Path) -> Result<Arc<MonLib>, String> {
    let root = root.canonicalize().map_err(|error| error.to_string())?;
    let cache = LIBRARIES.get_or_init(|| Mutex::new(HashMap::new()));
    let mut libraries = cache
        .lock()
        .map_err(|_| "validation library cache is poisoned")?;
    if let Some(library) = libraries.get(&root) {
        return Ok(Arc::clone(library));
    }
    let library = Arc::new(MonLib::load_validation(&root)?);
    libraries.insert(root, Arc::clone(&library));
    Ok(library)
}

/// Cache private user dictionaries without changing the shared chemistry tree.
pub fn load_with_dictionaries(
    root: &Path,
    documents: &[Arc<CifDocument>],
) -> Result<Arc<MonLib>, String> {
    dictionaries::overlay(load_library(root)?, documents)
}

/// Validate one parsed model without text serialization or external processes.
///
/// This checks existing sites. Hydrogen addition remains a separate Reduce3
/// operation so validation cannot repair a defect before reporting it.
pub fn run(document: &CifDocument, library: &MonLib, options: &Options) -> Result<Value, String> {
    if !options.sigma_cutoff.is_finite() || options.sigma_cutoff <= 0.0 {
        return Err("sigma_cutoff must be finite and positive".into());
    }
    let mut blocks = document
        .blocks()
        .iter()
        .filter(|block| BlockSource::new(block).table("atom_site").is_some());
    let block = blocks.next().ok_or("no block has an atom_site table")?;
    if blocks.next().is_some() {
        return Err("Multiple coordinate blocks need separate validation views".into());
    }
    source::validate_charges(block)?;
    let inline_library;
    let library = if source::has_restraints(document) {
        let text = dictionaries::dictionary_text(document)?;
        inline_library =
            library.with_validation_cif(&[(&text, Path::new("<inline-dictionary>"))])?;
        &inline_library
    } else {
        library
    };
    let original = reduce3::mmcif::structure_from_cif(&BlockSource::new(block))?;
    let (structure, normalization) = validation::normalized_copy(&original);
    let flat = FlatAtoms::from_structure(&structure);
    if flat.pos.is_empty() {
        return Err("the model has no atoms".into());
    }
    let dictionary_links = options
        .restraint_links
        .iter()
        .map(|(id, first, second)| {
            let find = |row| {
                flat.path
                    .iter()
                    .position(|&p| structure.atom(p).src == row)
                    .map(|i| i as u32)
                    .ok_or("A dictionary link atom row is outside the source model")
            };
            Ok(validation::DictionaryLink {
                id: id.clone(),
                atom_indices: [find(*first)?, find(*second)?],
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let mut graph = validation::build(
        &structure,
        &flat,
        library,
        &validation::ValidationOptions {
            nuclear: options.nuclear,
            use_cdl: options.use_cdl,
            dictionary_links,
            ..Default::default()
        },
    )?;
    for reason in source::unsupported(block)
        .into_iter()
        .chain(connections::uncovered(block, &structure, &flat, &graph))
    {
        for coverage in &mut graph.coverage {
            coverage.unsupported.push(reason.clone());
        }
    }
    let atoms = atoms::build(&structure, &flat, options.include_graph);
    let restraints = report::restraints(&graph);
    let evaluator = RestraintGraph::new(
        atoms.ids.clone(),
        graph.hydrogen_flags.clone(),
        restraints.bonds,
        restraints.angles,
        restraints.chiralities,
        restraints.planes,
    )?;
    let sites: Vec<_> = graph.sites.iter().map(|p| [p.x, p.y, p.z]).collect();
    let mut evaluation = evaluator.evaluate(&sites, options.sigma_cutoff)?;
    evaluation.findings.retain(|finding| match finding.check {
        "angles" => graph.angles[finding.proxy_index].origin == 0,
        "chiralities" => graph.chiralities[finding.proxy_index].origin == 0,
        _ => true,
    });
    if !graph.symmetry_bonds.is_empty() {
        let (cell, _) = reduce3::cell::crystal_operators(&structure)
            .ok_or("Symmetry bonds require a valid crystal cell and space group")?;
        for bond in &graph.symmetry_bonds {
            let first = graph.sites[bond.indices[0] as usize];
            let second = cell.orthogonalize(bond.image.apply(cell.fractionalize(graph.sites[bond.indices[1] as usize])));
            let delta = bond.ideal - first.dist(second);
            if !delta.is_finite() || !(delta * delta * bond.weight).is_finite() {
                return Err("Symmetry bond evaluation overflowed".into());
            }
        }
        if let Some((_, count)) = evaluation.counts.iter_mut().find(|(name, _)| *name == "bonds") {
            *count += graph.symmetry_bonds.len();
        }
    }
    let contacts = report::contacts(&structure, &flat, library, &graph, options);
    let mut result = report::finish(&graph, atoms, &evaluator, evaluation, contacts, options);
    result["normalization"] = serde_json::json!(normalization.iter().map(|change| {
        serde_json::json!({"kind": "symmetric_label_sites", "source_rows": change.source_rows})
    }).collect::<Vec<_>>());
    result["requested_links"] = serde_json::json!(options.restraint_links);
    Ok(result)
}

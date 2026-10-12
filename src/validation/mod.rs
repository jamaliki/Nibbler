//! Native chemical interpretation, restraint checks and Probe contacts.
//!
//! Reports retain unsupported coverage. Callers must run the reference backend
//! for those checks; findings from a partial graph do not establish completion.

mod atoms;
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
}

impl Default for Options {
    fn default() -> Self {
        Self {
            nuclear: false,
            use_cdl: true,
            sigma_cutoff: 5.0,
            include_graph: false,
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
    let structure = reduce3::mmcif::structure_from_cif(&BlockSource::new(block))?;
    let flat = FlatAtoms::from_structure(&structure);
    if flat.pos.is_empty() {
        return Err("the model has no atoms".into());
    }
    let mut graph = validation::build(
        &structure,
        &flat,
        library,
        &validation::ValidationOptions {
            nuclear: options.nuclear,
            use_cdl: options.use_cdl,
            ..Default::default()
        },
    )?;
    for reason in source::unsupported(document, block) {
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
    let evaluation = evaluator.evaluate(&sites, options.sigma_cutoff)?;
    let contacts = report::contacts(&structure, &flat, library, &graph, options);
    Ok(report::finish(
        &graph, atoms, &evaluator, evaluation, contacts, options,
    ))
}

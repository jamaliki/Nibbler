//! Private Python binding for running Reduce3 on documents.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyModule;

use crate::python::PyCifDocument;
use crate::reduce::{
    Approach, NTermCharge, OptParams, Params, ProbeParams, ReduceError, load_monomer_library, run, source_model_number,
};

/// Options already validated by `nibbler.reduce.run`.
#[derive(FromPyObject)]
#[pyo3(from_item_all)]
struct Options {
    chem_data: Option<PathBuf>,
    approach: String,
    compat: bool,
    add_flip_movers: bool,
    n_terminal_charge: String,
    keep_existing_h: bool,
    exclude_water: bool,
    use_neutron_distances: bool,
    preference_magnitude: f64,
    non_flip_preference: f64,
    planar_hydroxyl_preference: f64,
    acid_syn_preference: f64,
    skip_bond_fix_up: bool,
    set_flip_states: Option<String>,
    model_id: Option<usize>,
    alt_id: Option<String>,
    bonded_neighbor_depth: usize,
    stop_on_any_missing_hydrogen: bool,
    ignore_missing_restraints: bool,
    allow_missing_restraints_for: Vec<String>,
    required_hydrogen_target_keys: Option<String>,
    allow_undefined_ligand_targets: bool,
    allow_unresolved_bond_targets: bool,
    verbosity: i32,
    probe: HashMap<String, f64>,
    probe_flags: HashMap<String, bool>,
}

fn invalid(message: String) -> PyErr {
    PyValueError::new_err(message)
}

fn set_probe(probe: &mut ProbeParams, name: &str, value: f64) -> PyResult<()> {
    let field = match name {
        "probe_radius" => &mut probe.probe_radius,
        "density" => &mut probe.density,
        "worse_clash_cutoff" => &mut probe.worse_clash_cutoff,
        "clash_cutoff" => &mut probe.clash_cutoff,
        "contact_cutoff" => &mut probe.contact_cutoff,
        "uncharged_hydrogen_cutoff" => &mut probe.uncharged_hydrogen_cutoff,
        "charged_hydrogen_cutoff" => &mut probe.charged_hydrogen_cutoff,
        "bump_weight" => &mut probe.bump_weight,
        "hydrogen_bond_weight" => &mut probe.hydrogen_bond_weight,
        "gap_weight" => &mut probe.gap_weight,
        _ => return Err(invalid(format!("unknown probe parameter {name:?}"))),
    };
    *field = value;
    Ok(())
}

fn set_probe_flag(probe: &mut ProbeParams, name: &str, value: bool) -> PyResult<()> {
    let field = match name {
        "allow_weak_hydrogen_bonds" => &mut probe.allow_weak_hydrogen_bonds,
        "ignore_ion_interactions" => &mut probe.ignore_ion_interactions,
        "set_polar_hydrogen_radius" => &mut probe.set_polar_hydrogen_radius,
        _ => return Err(invalid(format!("unknown probe flag {name:?}"))),
    };
    *field = value;
    Ok(())
}

fn params(options: &Options) -> PyResult<Params> {
    let approach = match options.approach.as_str() {
        "add" => Approach::Add,
        "remove" => Approach::Remove,
        "optimize" => Approach::Optimize,
        other => return Err(invalid(format!("unknown approach {other:?}"))),
    };
    let n_terminal_charge = match options.n_terminal_charge.as_str() {
        "residue_one" => NTermCharge::ResidueOne,
        "first_in_chain" => NTermCharge::FirstInChain,
        "no_charge" => NTermCharge::NoCharge,
        other => return Err(invalid(format!("unknown n_terminal_charge {other:?}"))),
    };
    let mut probe = ProbeParams::reduce2_defaults();
    for (name, value) in &options.probe {
        set_probe(&mut probe, name, *value)?;
    }
    for (name, value) in &options.probe_flags {
        set_probe_flag(&mut probe, name, *value)?;
    }
    Ok(Params {
        approach,
        compat: options.compat,
        keep_existing_h: options.keep_existing_h,
        n_terminal_charge,
        exclude_water: options.exclude_water,
        model_id: options.model_id,
        ignore_missing_restraints: options.ignore_missing_restraints,
        stop_on_any_missing_hydrogen: options.stop_on_any_missing_hydrogen,
        allow_missing_restraints_for: options.allow_missing_restraints_for.clone(),
        required_hydrogen_parent_sources: None,
        allowed_unplaced_hydrogen_parent_sources: Vec::new(),
        allowed_unplaced_bridge_hydrogens: Vec::new(),
        allow_unresolved_bond_targets: options.allow_unresolved_bond_targets,
        opt: OptParams {
            probe,
            add_flip_movers: options.add_flip_movers,
            alt_id: options.alt_id.clone(),
            bonded_neighbor_depth: options.bonded_neighbor_depth,
            use_neutron_distances: options.use_neutron_distances,
            preference_magnitude: options.preference_magnitude,
            non_flip_preference: options.non_flip_preference,
            planar_hydroxyl_preference: options.planar_hydroxyl_preference,
            acid_syn_preference: options.acid_syn_preference,
            skip_bond_fixup: options.skip_bond_fix_up,
            flip_states: options.set_flip_states.clone().unwrap_or_default(),
            verbosity: options.verbosity,
            compat: options.compat,
            ..OptParams::default()
        },
    })
}

fn reduce_error_to_python(error: ReduceError) -> PyErr {
    PyValueError::new_err((
        error.code().as_str().to_owned(),
        error.message().to_owned(),
        Vec::<String>::new(),
    ))
}

/// Run Reduce3 on a document and return the output document and report.
#[pyfunction]
fn reduce_document(
    py: Python<'_>,
    document: &PyCifDocument,
    options: Options,
) -> PyResult<(PyCifDocument, String, String)> {
    let mut params = params(&options)?;
    let source = Arc::clone(document.document());
    let chem_data = options.chem_data;
    let target_keys = options.required_hydrogen_target_keys;
    let allow_undefined = options.allow_undefined_ligand_targets;
    if (allow_undefined || params.allow_unresolved_bond_targets) && target_keys.is_none() {
        return Err(invalid("Target masking requires retained hydrogen target keys".into()));
    }
    if target_keys.is_some() && (params.compat || params.approach != Approach::Add
        || !params.stop_on_any_missing_hydrogen || params.model_id.is_none()) {
        return Err(PyValueError::new_err(("REDUCE_INVALID_MODEL",
            "Required hydrogen target keys need an explicit model_id and strict fixed-mode addition", Vec::<String>::new())));
    }
    let reduction = py
        .detach(move || {
            if let Some(keys) = target_keys {
                let model = params.model_id.map(source_model_number).transpose()
                    .map_err(reduce_error_to_python)?
                    .ok_or_else(|| invalid("model_id is missing".into()))?;
                let retained = crate::observations::RetainedParents::read(&source, &keys, model, chem_data.as_deref())
                    .map_err(|message| PyValueError::new_err(("REDUCE_INVALID_MODEL", message, Vec::<String>::new())))?;
                if allow_undefined {
                    params.allowed_unplaced_hydrogen_parent_sources = retained.ligand_sources();
                    params.allowed_unplaced_bridge_hydrogens = retained.bridge_sources(chem_data.as_deref())
                        .map_err(invalid)?;
                }
                params.required_hydrogen_parent_sources = Some(retained.required);
            }
            let monomers = load_monomer_library(chem_data.as_deref()).map_err(reduce_error_to_python)?;
            run(&source, &monomers, &params).map_err(reduce_error_to_python)
        })?;
    let statuses = reduction.target_uncertainties_json.clone();
    let (output, report) = reduction.into_parts();
    Ok((
        PyCifDocument::new(Arc::new(output), document.schema_name().map(str::to_owned)),
        report,
        statuses,
    ))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(reduce_document, module)?)?;
    Ok(())
}

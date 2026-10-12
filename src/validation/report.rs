use reduce3::atominfo::{bonded_lists, extra_atom_info_validation};
use reduce3::clashcheck;
use reduce3::interp::{
    FlatAtoms,
    validation::{Coverage, ValidationGraph},
};
use reduce3::model::Structure;
use reduce3::monlib::MonLib;
use reduce3::world::World;
use serde_json::{Value, json};

use super::Options;
use super::atoms::AtomMap;
use crate::geometry::{Angle, Bond, Chirality, Evaluation, Plane, RestraintGraph};

pub(super) struct Restraints {
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
    pub chiralities: Vec<Chirality>,
    pub planes: Vec<Plane>,
}

pub(super) fn restraints(graph: &ValidationGraph) -> Restraints {
    Restraints {
        bonds: graph
            .bonds
            .iter()
            .map(|b| {
                (
                    b.indices[0] as usize,
                    b.indices[1] as usize,
                    b.ideal,
                    b.weight,
                    0.0,
                    b.origin == 0,
                )
            })
            .collect(),
        angles: graph
            .angles
            .iter()
            .map(|a| {
                (
                    a.indices[0] as usize,
                    a.indices[1] as usize,
                    a.indices[2] as usize,
                    a.ideal,
                    a.weight,
                )
            })
            .collect(),
        chiralities: graph
            .chiralities
            .iter()
            .map(|c| {
                (
                    c.indices[0] as usize,
                    c.indices[1] as usize,
                    c.indices[2] as usize,
                    c.indices[3] as usize,
                    c.ideal,
                    c.weight,
                    c.both_signs,
                    c.pseudochiral,
                )
            })
            .collect(),
        planes: graph
            .planes
            .iter()
            .map(|p| {
                (
                    p.indices.iter().map(|&i| i as usize).collect(),
                    p.weights.clone(),
                )
            })
            .collect(),
    }
}

fn geometry_results(coverage: &[Coverage], evaluation: Evaluation) -> (Value, Vec<Value>, bool) {
    let mut checks = serde_json::Map::new();
    for (name, count) in &evaluation.counts {
        let coverage = coverage.iter().find(|c| c.check == *name);
        let unsupported = coverage
            .map(|c| c.unsupported.clone())
            .unwrap_or_else(|| vec!["Chemical coverage is missing".into()]);
        let unresolved = coverage.map_or(0, |c| c.unresolved);
        checks.insert(
            (*name).into(),
            json!({"status": if unsupported.is_empty() && unresolved == 0 {
            "checked"
        } else { "incomplete" }, "count": count, "unsupported": unsupported,
            "unresolved": unresolved}),
        );
    }
    let mut failures: Vec<_> = evaluation
        .findings
        .into_iter()
        .map(|f| {
            json!({
                "check": f.check, "proxy_index": f.proxy_index, "indices": f.indices,
                "model": f.model, "target": f.target, "delta": f.delta, "sigma": f.sigma,
                "score": f.score, "residual": f.residual
            })
        })
        .collect();
    failures.extend(
        evaluation
            .hydrogen_failures
            .into_iter()
            .map(|(i, parents)| {
                json!({
                    "check": "hydrogens", "indices": [i], "parents": parents,
                    "reason": "hydrogen_parent_count"
                })
            }),
    );
    (Value::Object(checks), failures, evaluation.near_cutoff)
}

pub(super) fn contacts(
    structure: &Structure,
    flat: &FlatAtoms,
    library: &MonLib,
    graph: &ValidationGraph,
    options: &Options,
) -> Result<(clashcheck::Report, Option<Value>), String> {
    for name in ["atom_types", "bonds", "hydrogens"] {
        let coverage = graph
            .coverage
            .iter()
            .find(|c| c.check == name)
            .ok_or_else(|| format!("Chemical coverage is missing for {name}"))?;
        if !coverage.complete() || coverage.unresolved > 0 {
            return Err(format!(
                "{} coverage: {} ({} unresolved)",
                coverage.check,
                coverage.unsupported.join("; "),
                coverage.unresolved
            ));
        }
    }
    let bonded = bonded_lists(
        flat.pos.len(),
        graph.bonds.iter().map(|b| (b.indices[0], b.indices[1]))
            .chain(graph.symmetry_bonds.iter().map(|b| (b.indices[0], b.indices[1]))),
    );
    let mut info = extra_atom_info_validation(
        structure,
        flat,
        &graph.topology.etype,
        library,
        &bonded,
        options.nuclear,
    );
    if graph.typed_charges.len() != info.info.len() {
        return Err("Chemical charges must cover every Probe atom".into());
    }
    for (atom, &charge) in info.info.iter_mut().zip(&graph.typed_charges) {
        atom.charge = i8::try_from(charge).map_err(|_| "Chemical charge exceeds Probe range")?;
    }
    if info
        .warnings
        .lines()
        .any(|line| line.starts_with("Warning:"))
    {
        return Err(info.warnings);
    }
    let mut world = World::from_structure(structure);
    let initial_info = options.include_graph.then(|| {
        json!(
            info.info
                .iter()
                .enumerate()
                .map(|(i, a)| json!({"index":i,"vdw_radius":a.vdw_radius,
            "is_acceptor":a.is_acceptor,"is_donor":a.is_donor,"is_ion":a.is_ion,
            "is_dummy_hydrogen":a.is_dummy_hydrogen,"charge":a.charge,"alt":a.alt}))
                .collect::<Vec<_>>()
        )
    });
    world.info = info.info;
    world.bonded = bonded;
    let probe_options = clashcheck::Options::for_clashscore(&world, options.nuclear);
    let report = clashcheck::run(world, &probe_options)?;
    Ok((report, initial_info))
}

fn contact_results(
    contacts: Result<(clashcheck::Report, Option<Value>), String>,
    include_graph: bool,
    checks: &mut Value,
    failures: &mut Vec<Value>,
) -> Value {
    match contacts {
        Ok((report, initial_info)) => {
            failures.extend(report.clashes.iter().filter(|c| c.overlap.abs() > 0.4).map(|c| json!({
                "check": "clashes", "indices": [c.first, c.second], "overlap": c.overlap,
                "location": [c.location.x, c.location.y, c.location.z], "max_b_factor": c.max_b_factor
            })));
            checks["clashes"] = json!({"status": "checked", "count": report.clashes.len(),
                "clashscore": report.clashscore, "phantom_count": report.phantom_count});
            let mut output = json!({"clashes": report.clashes.iter().map(|c| json!({
                "indices": [c.first,c.second], "overlap": c.overlap,
                "location": [c.location.x,c.location.y,c.location.z]
            })).collect::<Vec<_>>(), "clashscore": report.clashscore,
                "initial_atom_info":initial_info});
            if include_graph {
                output["directed"] = json!(report.contacts.iter().map(|c|json!({
                    "source":c.source,"target":c.target,"class":format!("{:?}",c.interaction_type),
                    "dot_count":c.dot_count,"atom_gap":c.atom_gap,"min_dot_gap":c.min_dot_gap,
                    "score_over_density":c.score_over_density,
                    "location":[c.location.x,c.location.y,c.location.z]
                })).collect::<Vec<_>>());
                output["contact_atoms"] =
                    json!(report.atoms.iter().enumerate().map(|(i,a)|json!({
                    "index":i,"model":a.model_id.as_ref(),"chain":a.chain.as_ref(),"resseq":a.resseq_raw.as_ref(),
                    "icode":a.icode.as_ref(),"component":a.resname.as_ref(),"altloc":a.altloc.as_ref(),"name":a.name,
                    "element":a.element,"phantom":i>=report.n_real
                })).collect::<Vec<_>>());
            }
            output
        }
        Err(error) => {
            checks["clashes"] = json!({"status": "incomplete", "reason": error});
            Value::Null
        }
    }
}

pub(super) fn finish(
    graph: &ValidationGraph,
    atoms: AtomMap,
    restraints: &RestraintGraph,
    evaluation: Evaluation,
    contacts: Result<(clashcheck::Report, Option<Value>), String>,
    options: &Options,
) -> Value {
    let (mut checks, mut failures, near_cutoff) = geometry_results(&graph.coverage, evaluation);
    let contact_report =
        contact_results(contacts, options.include_graph, &mut checks, &mut failures);
    atoms.label_findings(&mut failures);
    let complete = graph.complete()
        && !near_cutoff
        && [
            "bonds",
            "angles",
            "chiralities",
            "planarities",
            "hydrogens",
            "clashes",
        ]
        .iter()
        .all(|name| checks[name]["status"] == "checked");
    let mut result = json!({"schema": "nibbler.native-validation.v1", "complete": complete,
        "checks": checks, "failures": failures, "near_cutoff": near_cutoff,
        "nuclear": options.nuclear, "use_cdl": options.use_cdl, "contacts": contact_report});
    if options.include_graph {
        result["atoms"] = json!(atoms.details);
        result["graph"] = json!({"hydrogen_flags": graph.hydrogen_flags,
            "bonds": restraints.bonds, "angles": restraints.angles,
            "chiralities": restraints.chiralities, "planes": restraints.planes,
            "angle_origins": graph.angles.iter().map(|a| a.origin).collect::<Vec<_>>(),
            "chirality_origins": graph.chiralities.iter().map(|c| c.origin).collect::<Vec<_>>(),
            "bonds_asu": graph.symmetry_bonds.iter().map(|b| json!({
                "i_seq": b.indices[0], "j_seq": b.indices[1],
                "rotation": b.image.r, "translation": b.image.t,
                "distance_ideal": b.ideal, "weight": b.weight, "slack":0.0,
                "origin_id": b.origin,
            })).collect::<Vec<_>>()});
        result["coverage"] = json!(
            graph
                .coverage
                .iter()
                .map(|c| json!({
                    "check": c.check, "unsupported": c.unsupported, "unresolved": c.unresolved
                }))
                .collect::<Vec<_>>()
        );
        result["dictionaries"] = json!(
            graph
                .dictionaries
                .iter()
                .map(|d| json!({
                    "component": d.residue, "source": d.source, "modifications": d.modifications
                }))
                .collect::<Vec<_>>()
        );
    }
    result
}

#[cfg(test)]
mod tests {
    use super::geometry_results;
    use crate::geometry::Evaluation;

    #[test]
    fn a_missing_coverage_record_is_incomplete() {
        let evaluation = Evaluation {
            near_cutoff: false,
            counts: vec![("bonds", 0)],
            findings: Vec::new(),
            hydrogen_failures: Vec::new(),
        };
        let (checks, _, _) = geometry_results(&[], evaluation);
        assert_eq!(checks["bonds"]["status"], "incomplete");
    }
}

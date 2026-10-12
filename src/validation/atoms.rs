use reduce3::interp::FlatAtoms;
use reduce3::model::Structure;
use serde_json::{Value, json};

pub(super) struct AtomMap {
    pub ids: Vec<String>,
    sources: Vec<u32>,
    pub details: Option<Vec<Value>>,
}

pub(super) fn build(structure: &Structure, flat: &FlatAtoms, detailed: bool) -> AtomMap {
    let mut map = AtomMap {
        ids: Vec::with_capacity(flat.path.len()),
        sources: Vec::with_capacity(flat.path.len()),
        details: detailed.then(Vec::new),
    };
    for (index, &p) in flat.path.iter().enumerate() {
        let atom = structure.atom(p);
        let group = structure.atom_group(p);
        let residue = structure.residue_group(p);
        let chain = structure.chain(p);
        let model = &structure.models[p.model as usize];
        let id = json!([
            model.id,
            chain.id,
            residue.resseq,
            residue.icode,
            group.resname,
            group.altloc,
            atom.name
        ])
        .to_string();
        if let Some(details) = &mut map.details {
            details.push(json!({
            "index": index, "id": id, "source_index": atom.src,
            "model": model.id, "chain": chain.id, "resseq": residue.resseq,
            "icode": residue.icode.trim(), "component": group.resname.trim(),
            "altloc": group.altloc.trim(), "name": atom.name.trim(), "element": atom.elem(),
            "xyz": [atom.xyz.x, atom.xyz.y, atom.xyz.z], "occupancy": atom.occ,
            "b_iso": atom.b, "hetero": atom.hetero, "charge": atom.charge.trim()}));
        }
        map.ids.push(id);
        map.sources.push(atom.src);
    }
    map
}

impl AtomMap {
    pub fn label_findings(&self, failures: &mut [Value]) {
        for failure in failures {
            if let Some(indices) = failure["indices"].as_array() {
                let indices: Vec<_> = indices
                    .iter()
                    .filter_map(Value::as_u64)
                    .map(|i| i as usize)
                    .collect();
                failure["atom_ids"] =
                    json!(indices.iter().map(|&i| &self.ids[i]).collect::<Vec<_>>());
                failure["source_rows"] =
                    json!(indices.iter().map(|&i| self.sources[i]).collect::<Vec<_>>());
            }
            if let Some(parents) = failure["parents"].as_array() {
                let parents: Vec<_> = parents
                    .iter()
                    .filter_map(Value::as_u64)
                    .map(|i| i as usize)
                    .collect();
                failure["parent_ids"] =
                    json!(parents.iter().map(|&i| &self.ids[i]).collect::<Vec<_>>());
            }
        }
    }
}

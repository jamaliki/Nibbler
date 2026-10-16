//! Read immutable CCD definitions from the selected Reduce3 chemical library.

use super::super::{ComponentRegistry, PdbxModel, build_component_registry};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock, PoisonError};

static DEFINITIONS: OnceLock<Mutex<BTreeMap<PathBuf, ComponentRegistry>>> = OnceLock::new();

pub(crate) fn registry(model: &PdbxModel, root: &Path) -> Result<ComponentRegistry, String> {
    let root = root.canonicalize().map_err(|error| error.to_string())?;
    let cache = DEFINITIONS.get_or_init(Mutex::default);
    let mut cache = cache.lock().unwrap_or_else(PoisonError::into_inner);
    let mut result = BTreeMap::new();
    for component in &model.components {
        let id = component.id.to_ascii_uppercase();
        if !id.bytes().all(|c| c.is_ascii_alphanumeric()) || id.is_empty() {
            continue;
        }
        let file = root
            .join("chemical_components")
            .join(id[..1].to_ascii_lowercase())
            .join(format!("data_{id}.cif"));
        if !cache.contains_key(&file) {
            let definitions = match std::fs::read(&file) {
                Ok(bytes) => {
                    let document = crate::cif::parse(&bytes).map_err(|error| error.to_string())?;
                    build_component_registry(&document).map_err(|error| error.to_string())?
                }
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                    ComponentRegistry::default()
                }
                Err(error) => return Err(format!("Cannot read CCD {}: {error}", file.display())),
            };
            cache.insert(file.clone(), definitions);
        }
        if let Some(definition) = cache.get(&file).and_then(|registry| registry.get(&id)) {
            result.insert(id.to_ascii_lowercase(), definition.clone());
        }
    }
    Ok(ComponentRegistry::new(result))
}

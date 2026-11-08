//! PyO3 bindings for the private `nibbler._core` module.
//!
//! Each submodule serves one Python facade and registers its own classes and functions;
//! the Python-visible names and `module=` paths do not depend on this layout.

pub(crate) mod document;
mod geometry;
mod mmcif;
pub(crate) mod observations;
#[cfg(feature = "reduce3")]
mod reduce;
mod scan;
pub(crate) mod table;
#[cfg(feature = "reduce3")]
mod validation;

use pyo3::prelude::*;
use pyo3::types::PyModule;

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    document::register(module)?;
    table::register(module)?;
    mmcif::register(module)?;
    scan::register(module)?;
    geometry::register(module)?;
    observations::register(module)?;
    #[cfg(feature = "reduce3")]
    reduce::register(module)?;
    #[cfg(feature = "reduce3")]
    validation::register(module)?;
    Ok(())
}

//! Native core for the Nibbler CIF toolkit.
//!
//! The [`cif`] module contains the shared CIF 1.1 syntax implementation. Higher-level
//! projection, schema, and Python APIs build on that one production path. With the
//! `reduce3` feature, `reduce` adds hydrogens to parsed models in memory.

pub mod cif;
pub mod modelcif;
pub mod pdbx;
#[cfg(feature = "reduce3")]
pub mod reduce;

#[cfg(feature = "python")]
mod python;
#[cfg(all(feature = "python", feature = "reduce3"))]
mod python_reduce;
#[cfg(feature = "python")]
mod python_scan;

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Version of the public Python/Rust contract.
pub const CONTRACT_VERSION: u16 = 1;

/// Return the native crate version.
#[must_use]
#[cfg_attr(feature = "python", pyfunction)]
pub const fn native_version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// Return the public contract version.
#[must_use]
#[cfg_attr(feature = "python", pyfunction)]
pub const fn contract_version() -> u16 {
    CONTRACT_VERSION
}

/// Return the standard Cargo profile class inferred from debug assertions.
#[must_use]
#[cfg_attr(feature = "python", pyfunction)]
pub const fn build_profile() -> &'static str {
    if cfg!(debug_assertions) {
        "debug"
    } else {
        "release"
    }
}

/// Initialize the private native Python module.
#[cfg(feature = "python")]
#[pymodule]
fn _core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(native_version, module)?)?;
    module.add_function(wrap_pyfunction!(contract_version, module)?)?;
    module.add_function(wrap_pyfunction!(build_profile, module)?)?;
    python::register(module)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{CONTRACT_VERSION, build_profile, contract_version, native_version};

    #[test]
    fn exposes_contract_metadata() {
        assert_eq!(contract_version(), CONTRACT_VERSION);
        assert_eq!(native_version(), env!("CARGO_PKG_VERSION"));
        assert_eq!(
            build_profile(),
            if cfg!(debug_assertions) {
                "debug"
            } else {
                "release"
            }
        );
    }
}

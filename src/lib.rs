//! Native core for the Nibbler CIF toolkit.
//!
//! Phase 0 exposes build and contract metadata only. CIF parsing will be added through
//! the single production path specified in `DESIGN.md`.

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Version of the public Python/Rust contract established by Phase 0.
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

/// Initialize the private native Python module.
#[cfg(feature = "python")]
#[pymodule]
fn _core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(native_version, module)?)?;
    module.add_function(wrap_pyfunction!(contract_version, module)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{CONTRACT_VERSION, contract_version, native_version};

    #[test]
    fn exposes_contract_metadata() {
        assert_eq!(contract_version(), CONTRACT_VERSION);
        assert_eq!(native_version(), env!("CARGO_PKG_VERSION"));
    }
}

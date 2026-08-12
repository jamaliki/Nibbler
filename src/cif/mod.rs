//! Strict CIF 1.1 parsing and canonical serialization.
//!
//! [`parse`] accepts immutable bytes and returns an order-preserving logical document.
//! Parsed text values retain spans into a shared source buffer, so materializing a
//! document does not allocate a separate string for every value.
//!
//! # Example
//!
//! ```
//! use _core::cif::{parse, write_canonical};
//!
//! let document = parse(b"data_example\n_entry.id example\n")?;
//! let canonical = write_canonical(&document)?;
//! assert_eq!(parse(canonical.as_bytes())?, document);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod arrow;
mod binary;
mod dictionary;
mod dictionary_compiler;
mod dictionary_source;
mod document;
mod error;
mod input;
mod lexer;
pub(crate) mod numeric;
mod parallel;
mod parallel_index;
mod parallel_probe;
mod parser;
mod projection;
mod projection_sink;
mod schema;
mod schema_artifact;
mod sink;
mod source;
mod table;
mod token;
mod validation;
mod validation_engine;
mod writer;

pub use dictionary::{
    CategoryDefinition, Dictionary, DictionaryError, DictionaryMetadata, ItemDefinition, ItemRange,
    TypeDefinition,
};
pub use dictionary_compiler::compile_dictionary;
pub use document::{
    BlockKind, CifBlock, CifDocument, CifEntry, CifFrame, CifItem, CifLoop, CifRow, CifValue,
    CifValueRef, CifValues, OriginalLexeme, QuoteStyle, StandardUncertainty, TextValue,
    TextValueRef,
};
pub use error::{ParseError, ParseErrorCode, SourceSpan, WriteError, WriteErrorCode};
pub use input::{InputError, InputErrorCode, InputLimits, decode_source, read_source_file};
#[cfg(feature = "python")]
pub(crate) use input::{LoadedBytes, decode_input, read_input_file};
pub use parser::{
    Limits, ParseOptions, parse, parse_source, parse_source_with_options, parse_with_options,
};
pub use projection::{
    Predicate, ProjectionError, ProjectionErrorCode, ProjectionPlan, project, project_source,
    project_source_with_options, project_with_options,
};
pub use schema::{SchemaError, SchemaName};
pub use schema_artifact::{SchemaArtifactError, decode_dictionary, encode_dictionary};
pub use source::SourceBuffer;
pub use table::{CifCell, CifCellRef, CifColumn, CifTable, ColumnType, MissingKind, RowProvenance};
pub use validation::{Diagnostic, Severity, ValidationReport, validate_document};
pub use writer::{format_text, write_canonical, write_preserving};

pub use arrow::{ArrowExportError, ArrowMissingPolicy, export_arrow_stream};
pub use binary::{
    BinaryCifError, BinaryCifErrorCode, decode_binary, encode_binary, project_binary,
};
#[cfg(feature = "python")]
pub(crate) use binary::{is_binary_container, project_binary_named};
#[cfg(feature = "fuzzing")]
#[doc(hidden)]
pub use lexer::fuzz_lexer;

//! Hydrogen addition and optimization of parsed models with Reduce3.
//!
//! [`run`] takes the first data block of a [`CifDocument`] that has an `_atom_site`
//! loop, runs [Reduce3](https://github.com/jamaliki/reduce3) on it in memory, and
//! returns the result as a new document. No text is written or parsed again: Reduce3
//! reads the parsed values through [`BlockSource`] and streams its output model into a
//! [`DocumentSink`].
//!
//! The result block keeps the source block code and, by default, every category of the
//! source block in its order: `_atom_site` is rebuilt with the source's items and label
//! identifiers (new hydrogens take their residue's), atom ids are renumbered,
//! `_atom_site_anisotrop` follows the new ids, and `_atom_type` gains the elements it
//! lacks; every other entry is shared unchanged, so `_struct_conn`, the entities and the
//! sequence schemes stay valid ([`block_document`]). With [`Params::compat`] the block
//! instead has the layout Reduce2 writes (cell, space group, `_struct_asym`,
//! `_chem_comp`, `_atom_site`, and `_atom_site_anisotrop`, with regenerated label
//! identifiers; [`structure_document`]). Either way the entries equal a parse of the
//! file the `reduce3` program writes in that mode.
//!
//! Reduce3 needs the cctbx `chem_data` monomer library; [`load_monomer_library`] finds
//! and caches it.
//!
//! # Example
//!
//! ```no_run
//! use _core::cif::parse;
//! use _core::reduce::{Params, load_monomer_library, run};
//!
//! let document = parse(&std::fs::read("model.cif")?)?;
//! let monomers = load_monomer_library(None)?;
//! let reduction = run(&document, &monomers, &Params::default())?;
//! println!("{}", reduction.report());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::borrow::Cow;
use std::error::Error;
use std::fmt::{self, Display, Formatter};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock, PoisonError};

use reduce3::cifsource::{CifCell, CifSink, CifSource, CifTable, ColumnValue};
use reduce3::model::Structure;

pub use reduce3::hplace::NTermCharge;
pub use reduce3::optimizer::OptParams;
pub use reduce3::probe::ProbeParams;
pub use reduce3::{Approach, MonLib, Params};

use crate::cif::{
    BlockCategory, CifBlock, CifDocument, CifEntry, CifItem, CifLoop, CifValue, CifValueRef,
    LoopColumn, SourceColumnBuilder, TextCell, TextColumn, TextColumnBuilder, in_category,
    split_tag,
};

/// The class of a [`ReduceError`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReduceErrorCode {
    /// No block has an `_atom_site` loop, or the loop cannot be read as a model.
    InvalidModel,
    /// The `chem_data` monomer library was not found.
    MonomerLibraryNotFound,
    /// The monomer library could not be read.
    MonomerLibraryInvalid,
    /// Reduce3 rejected the model, for example for missing restraints.
    Failed,
    /// The result does not fit the in-memory document representation.
    ResultTooLarge,
}

impl ReduceErrorCode {
    /// Return the stable machine-readable code.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::InvalidModel => "REDUCE_INVALID_MODEL",
            Self::MonomerLibraryNotFound => "REDUCE_MONOMER_LIBRARY_NOT_FOUND",
            Self::MonomerLibraryInvalid => "REDUCE_MONOMER_LIBRARY_INVALID",
            Self::Failed => "REDUCE_FAILED",
            Self::ResultTooLarge => "REDUCE_RESULT_TOO_LARGE",
        }
    }
}

/// A failure while running Reduce3 on a document.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReduceError {
    code: ReduceErrorCode,
    message: String,
}

impl ReduceError {
    fn new(code: ReduceErrorCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }

    /// Return the error class.
    #[must_use]
    pub const fn code(&self) -> ReduceErrorCode {
        self.code
    }

    /// Return the human-readable failure detail.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }
}

impl Display for ReduceError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {}", self.code.as_str(), self.message)
    }
}

impl Error for ReduceError {}

/// The outcome of [`run`]: the new document and Reduce3's report.
#[derive(Clone, Debug)]
pub struct Reduction {
    document: CifDocument,
    report: String,
    #[cfg(feature = "python")]
    pub(crate) target_uncertainties_json: String,
}

impl Reduction {
    /// Return the document holding the output model.
    #[must_use]
    pub const fn document(&self) -> &CifDocument {
        &self.document
    }

    /// Return the report Reduce3 writes to its description file: scores, the
    /// configuration of every movable group, and deleted hydrogens.
    #[must_use]
    pub fn report(&self) -> &str {
        &self.report
    }

    /// Return the document and the report.
    #[must_use]
    pub fn into_parts(self) -> (CifDocument, String) {
        (self.document, self.report)
    }
}

/// Run Reduce3 on the first data block of `document` that has an `_atom_site` loop.
/// In fixed mode, `model_id` names the source `_atom_site.pdbx_PDB_model_num`.
/// Compatibility mode keeps Reduce2's ordinal selection.
///
/// # Errors
///
/// Returns [`ReduceErrorCode::InvalidModel`] when no block holds a readable model,
/// [`ReduceErrorCode::Failed`] when Reduce3 rejects it, and
/// [`ReduceErrorCode::ResultTooLarge`] when the output cannot be stored.
pub fn run(
    document: &CifDocument,
    monomers: &MonLib,
    params: &Params,
) -> Result<Reduction, ReduceError> {
    run_block(model_block(document)?, monomers, params)
}

/// [`run`], which with `record` also returns the source rows that the output
/// `_atom_site` rows copy.
///
/// # Errors
///
/// Returns the same errors as [`run`].
#[cfg(feature = "python")]
pub(crate) fn run_recording(
    document: &CifDocument,
    monomers: &MonLib,
    params: &Params,
    record: bool,
) -> Result<(Reduction, Option<SourceCopies>), ReduceError> {
    run_block_with(model_block(document)?, monomers, params, record)
}

/// Return the first data block of `document` that has an `_atom_site` loop.
fn model_block(document: &CifDocument) -> Result<&CifBlock, ReduceError> {
    document
        .blocks()
        .iter()
        .find(|block| BlockSource::new(block).table("atom_site").is_some())
        .ok_or_else(|| {
            ReduceError::new(
                ReduceErrorCode::InvalidModel,
                "no block has an _atom_site loop",
            )
        })
}

/// Run Reduce3 on one data block.
///
/// # Errors
///
/// Returns the same errors as [`run`].
pub fn run_block(
    block: &CifBlock,
    monomers: &MonLib,
    params: &Params,
) -> Result<Reduction, ReduceError> {
    run_block_with(block, monomers, params, false).map(|(reduction, _)| reduction)
}

/// [`run_block`], which with `record` also returns the source rows that the output
/// `_atom_site` rows copy.
fn run_block_with(
    block: &CifBlock,
    monomers: &MonLib,
    params: &Params,
    record: bool,
) -> Result<(Reduction, Option<SourceCopies>), ReduceError> {
    let source = BlockSource::new(block);
    // The writer compares output coordinates with these numbers instead of parsing
    // the source again.
    let (mut structure, numbers) = reduce3::mmcif::structure_and_numbers_from_cif(&source)
        .map_err(|message| ReduceError::new(ReduceErrorCode::InvalidModel, message))?;
    if structure.atoms_size() == 0 {
        return Err(ReduceError::new(
            ReduceErrorCode::InvalidModel,
            "the _atom_site loop has no atoms",
        ));
    }
    if !params.compat {
        restore_single_model_id(&mut structure, &source);
    }
    let mut pipeline_params = params.clone();
    if !params.compat
        && let Some(model_id) = params.model_id
    {
        pipeline_params.model_id = Some(source_model_index(&structure, model_id)? + 1);
    }
    let output = reduce3::pipeline::run(structure, monomers, &pipeline_params)
        .map_err(|message| ReduceError::new(ReduceErrorCode::Failed, message))?;
    #[cfg(feature = "python")]
    let target_uncertainties_json = if !params.allowed_unplaced_hydrogen_parent_sources.is_empty()
        || !params.allowed_unplaced_bridge_hydrogens.is_empty()
        || params.allow_unresolved_bond_targets
    {
        let model = params
            .model_id
            .map(source_model_number)
            .transpose()?
            .ok_or_else(|| {
                ReduceError::new(
                    ReduceErrorCode::InvalidModel,
                    "Target masking requires an explicit source model ID",
                )
            })?;
        let statuses = crate::observations::uncertainty::from_output(
            block,
            model,
            &output,
            params.required_hydrogen_parent_sources.as_deref(),
            params.allow_unresolved_bond_targets,
        )
        .map_err(|message| ReduceError::new(ReduceErrorCode::InvalidModel, message))?;
        serde_json::to_string(&statuses)
            .map_err(|error| ReduceError::new(ReduceErrorCode::ResultTooLarge, error.to_string()))?
    } else {
        "[]".to_owned()
    };
    let (document, copies) = if params.compat {
        (
            structure_document(&output.structure, block.code().unwrap_or("default"))?,
            None,
        )
    } else {
        block_document_with(&output.structure, block, Some(&numbers), record)?
    };
    let reduction = Reduction {
        document,
        report: output.description,
        #[cfg(feature = "python")]
        target_uncertainties_json,
    };
    Ok((reduction, copies))
}

pub(crate) fn source_model_number(model_id: usize) -> Result<i64, ReduceError> {
    i64::try_from(model_id)
        .ok()
        .filter(|number| *number > 0)
        .ok_or_else(|| {
            ReduceError::new(
                ReduceErrorCode::InvalidModel,
                "Source model ID must be a positive signed 64-bit integer",
            )
        })
}

fn restore_single_model_id(structure: &mut Structure, source: &BlockSource<'_>) {
    if structure.models.len() == 1
        && structure.models[0].id.is_empty()
        && let Some(table) = source.table("atom_site")
    {
        structure.models[0].id = table.column("pdbx_PDB_model_num").map_or_else(
            || "1".to_owned(),
            |column| table.cell(0, column).into_owned(),
        );
    }
}

fn source_model_index(structure: &Structure, model_id: usize) -> Result<usize, ReduceError> {
    let number = source_model_number(model_id)?;
    let mut matches = structure
        .models
        .iter()
        .enumerate()
        .filter(|(_, model)| model.id.trim().parse::<i64>() == Ok(number));
    let index = matches.next().map(|(index, _)| index).ok_or_else(|| {
        ReduceError::new(
            ReduceErrorCode::InvalidModel,
            format!("Source model ID {number} is not present"),
        )
    })?;
    if matches.next().is_some() {
        return Err(ReduceError::new(
            ReduceErrorCode::InvalidModel,
            format!("Source model ID {number} is ambiguous"),
        ));
    }
    Ok(index)
}

#[cfg(test)]
mod source_model_tests {
    use super::{ReduceErrorCode, source_model_index, source_model_number};
    use reduce3::model::{Model, Structure};

    fn structure(ids: &[&str]) -> Structure {
        Structure {
            models: ids
                .iter()
                .map(|id| Model {
                    id: (*id).into(),
                    chains: Vec::new(),
                })
                .collect(),
            ..Structure::default()
        }
    }

    #[test]
    fn source_model_ids_map_without_renumbering() {
        let input = structure(&["35", "18"]);
        assert_eq!(source_model_index(&input, 18).ok(), Some(1));
        assert_eq!(source_model_index(&input, 35).ok(), Some(0));
        assert_eq!(input.models[0].id, "35");
        assert_eq!(input.models[1].id, "18");
    }

    #[test]
    fn missing_source_ids_do_not_fall_back_to_ordinals() {
        assert_eq!(
            source_model_index(&structure(&["18"]), 1)
                .err()
                .map(|error| error.code()),
            Some(ReduceErrorCode::InvalidModel)
        );
    }

    #[test]
    fn repeated_and_numeric_alias_ids_are_ambiguous() {
        for ids in [["18", "35", "18"], ["18", "018", "35"]] {
            assert!(
                source_model_index(&structure(&ids), 18)
                    .err()
                    .is_some_and(|error| error.message().contains("ambiguous"))
            );
        }
    }

    #[test]
    fn source_model_numbers_reject_zero_and_overflow() {
        assert!(source_model_number(0).is_err());
        assert!(source_model_number(usize::MAX).is_err());
    }

    #[test]
    fn parsed_single_model_ids_match_source_values() {
        use super::{BlockSource, restore_single_model_id};
        for text in [
            include_bytes!("../tests/reduce3/linked_8zje.cif").as_slice(),
            include_bytes!("../tests/reduce3/linked_9hj3.cif").as_slice(),
        ] {
            let document = crate::cif::parse(text).unwrap();
            let source = BlockSource::new(&document.blocks()[0]);
            let mut structure = reduce3::mmcif::structure_from_cif(&source).unwrap();
            restore_single_model_id(&mut structure, &source);
            assert_eq!(
                structure
                    .models
                    .iter()
                    .map(|model| &model.id)
                    .collect::<Vec<_>>(),
                vec!["1"]
            );
            assert_eq!(source_model_index(&structure, 1).ok(), Some(0));
        }
    }
}

#[cfg(all(test, feature = "python"))]
mod copy_tests {
    use super::{Params, SourceCopies, load_monomer_library, run_recording};

    #[test]
    fn copied_sites_build_the_model_that_a_full_decode_builds() {
        let monomers = load_monomer_library(None).expect("Reduce tests require REDUCE3_CHEM_DATA");
        for text in [
            include_bytes!("../tests/reduce3/1crn.cif").as_slice(),
            include_bytes!("../tests/fixtures/chemistry/ligand_ion_water.cif").as_slice(),
        ] {
            let source = crate::cif::parse(text).expect("The fixture must parse");
            let model = crate::pdbx::build_model(&source).expect("The fixture must form a model");
            assert!(model.loop_sites().is_some());
            let (reduction, copies) = run_recording(&source, &monomers, &Params::default(), true)
                .expect("Reduce3 must protonate the fixture");
            let copies = copies.expect("The sink must rewrite the source _atom_site loop");
            assert!(copies.rows.iter().any(|row| *row != SourceCopies::NONE));
            assert!(
                copies.rows.iter().any(|row| *row == SourceCopies::NONE)
                    || copies.made.iter().any(|made| *made != 0)
            );
            let full = crate::pdbx::build_model(reduction.document())
                .expect("The output must form a model");
            let copied = crate::pdbx::build_model_with_copies(
                reduction.document(),
                &model,
                &copies.rows,
                &copies.made,
            )
            .expect("The copied output must form a model");
            // Debug text also tells apart the signs of zero.
            assert_eq!(format!("{copied:?}"), format!("{full:?}"));
        }
    }
}

/// Build a one-block document holding `structure` in the block it was read from:
/// every category of `block` in order, with the atom tables rebuilt (see the module
/// documentation). Entries Reduce3 does not change are shared with `block`.
///
/// # Errors
///
/// Returns [`ReduceErrorCode::InvalidModel`] when `block` has no `_atom_site` loop and
/// [`ReduceErrorCode::ResultTooLarge`] when a column exceeds the in-memory representation.
pub fn block_document(structure: &Structure, block: &CifBlock) -> Result<CifDocument, ReduceError> {
    block_document_with(structure, block, None, false).map(|(document, _)| document)
}

/// [`block_document`], with the numbers read from `block` with its model. With `record`,
/// it also returns the source rows that the `_atom_site` rows copy, when the sink
/// rewrote the source loop.
fn block_document_with(
    structure: &Structure,
    block: &CifBlock,
    numbers: Option<&reduce3::mmcif::SourceNumbers>,
    record: bool,
) -> Result<(CifDocument, Option<SourceCopies>), ReduceError> {
    let code = block.code().unwrap_or("default");
    let mut sink = DocumentSink::with_source(block);
    sink.record = record;
    reduce3::mmcif::write_cif_preserving_with(
        structure,
        &BlockSource::new(block),
        numbers,
        code,
        &mut sink,
    )
    .map_err(|message| ReduceError::new(ReduceErrorCode::InvalidModel, message))?;
    let copies = sink.copies.take();
    sink.finish(code).map(|document| (document, copies))
}

/// Build a one-block document holding `structure` as Reduce3 writes it.
///
/// # Errors
///
/// Returns [`ReduceErrorCode::ResultTooLarge`] when a column exceeds the in-memory
/// representation.
pub fn structure_document(structure: &Structure, code: &str) -> Result<CifDocument, ReduceError> {
    let mut sink = DocumentSink::new();
    reduce3::mmcif::write_cif(structure, &mut sink);
    sink.finish(code)
}

/// Monomer libraries loaded so far, by canonical directory.
type LoadedLibraries = Mutex<Vec<(PathBuf, Arc<MonLib>)>>;

static MONOMER_LIBRARIES: OnceLock<LoadedLibraries> = OnceLock::new();

/// Load the cctbx `chem_data` monomer library, once per directory and process.
///
/// With `chem_data`, that directory is used. Otherwise the library is looked up the way
/// the `reduce3` program does: `$REDUCE3_CHEM_DATA`, `$CHEM_DATA`, the parent of
/// `$MMTBX_CCP4_MONOMER_LIB` or `$CLIBD_MON`, and the active conda environment.
/// Reduce3 keeps reading residue definitions from the directory as it needs them.
///
/// # Errors
///
/// Returns [`ReduceErrorCode::MonomerLibraryNotFound`] when no directory with `geostd`
/// and `mon_lib` is found, and [`ReduceErrorCode::MonomerLibraryInvalid`] when its
/// index files cannot be read.
pub fn load_monomer_library(chem_data: Option<&Path>) -> Result<Arc<MonLib>, ReduceError> {
    let root = match chem_data {
        Some(path) if path.join("geostd").is_dir() && path.join("mon_lib").is_dir() => {
            path.to_path_buf()
        }
        Some(path) => {
            return Err(ReduceError::new(
                ReduceErrorCode::MonomerLibraryNotFound,
                format!(
                    "{} is not a chem_data directory (it needs geostd and mon_lib)",
                    path.display()
                ),
            ));
        }
        None => MonLib::locate(None).ok_or_else(|| {
            ReduceError::new(
                ReduceErrorCode::MonomerLibraryNotFound,
                "could not find chem_data (the monomer library); pass its directory or set \
                 REDUCE3_CHEM_DATA",
            )
        })?,
    };
    let root = std::fs::canonicalize(&root).unwrap_or(root);
    let mut loaded = MONOMER_LIBRARIES
        .get_or_init(Mutex::default)
        .lock()
        .unwrap_or_else(PoisonError::into_inner);
    if let Some((_, library)) = loaded.iter().find(|(path, _)| *path == root) {
        return Ok(Arc::clone(library));
    }
    let library =
        Arc::new(MonLib::load(&root).map_err(|message| {
            ReduceError::new(ReduceErrorCode::MonomerLibraryInvalid, message)
        })?);
    loaded.push((root, Arc::clone(&library)));
    Ok(library)
}

/// A parsed data block as Reduce3 input.
#[derive(Clone, Copy, Debug)]
pub struct BlockSource<'a> {
    block: &'a CifBlock,
}

impl<'a> BlockSource<'a> {
    /// Wrap one block.
    #[must_use]
    pub const fn new(block: &'a CifBlock) -> Self {
        Self { block }
    }
}

impl CifSource for BlockSource<'_> {
    type Table<'t>
        = BlockTable<'t>
    where
        Self: 't;

    fn table(&self, category: &str) -> Option<BlockTable<'_>> {
        BlockCategory::find(self.block.entries(), category).map(BlockTable::new)
    }

    fn categories(&self) -> Vec<String> {
        let mut names: Vec<String> = Vec::new();
        for entry in self.block.entries() {
            let tag = match entry {
                CifEntry::Loop(cif_loop) => cif_loop.tags().first().map(String::as_str),
                CifEntry::Item(item) => Some(item.tag()),
                CifEntry::Frame(_) => None,
            };
            if let Some((name, _)) = tag.and_then(split_tag)
                && !names.iter().any(|n| n.eq_ignore_ascii_case(name))
            {
                names.push(name.to_owned());
            }
        }
        names
    }
}

/// One category of a [`BlockSource`]: a loop, or the category's scalar items as one row.
#[derive(Clone, Debug)]
pub struct BlockTable<'a> {
    category: BlockCategory<'a>,
    /// Direct access to the columns a columnar loop stores as text.
    texts: Vec<Option<TextColumn<'a>>>,
}

impl<'a> BlockTable<'a> {
    fn new(category: BlockCategory<'a>) -> Self {
        let texts = match &category {
            BlockCategory::Loop(cif_loop) => (0..cif_loop.column_count())
                .map(|column| cif_loop.text_column(column))
                .collect(),
            BlockCategory::Items(_) => Vec::new(),
        };
        Self { category, texts }
    }

    fn value(&self, row: usize, column: usize) -> Option<CifValueRef<'a>> {
        match self.texts.get(column) {
            Some(Some(text)) => text.value(row),
            _ => self.category.value(row, column),
        }
    }
}

impl CifTable for BlockTable<'_> {
    fn row_count(&self) -> usize {
        self.category.row_count()
    }

    fn column(&self, tag: &str) -> Option<usize> {
        self.category.column(tag)
    }

    fn cell(&self, row: usize, column: usize) -> Cow<'_, str> {
        match self.value(row, column) {
            Some(CifValueRef::Text(text)) => Cow::Borrowed(text.as_str()),
            Some(CifValueRef::Integer(number, original)) => original.map_or_else(
                || Cow::Owned(number.to_string()),
                |text| Cow::Borrowed(text.as_str()),
            ),
            Some(CifValueRef::Float(number, _, original)) => original.map_or_else(
                || Cow::Owned(number.to_string()),
                |text| Cow::Borrowed(text.as_str()),
            ),
            Some(CifValueRef::NotApplicable) => Cow::Borrowed("."),
            Some(CifValueRef::Unknown) | None => Cow::Borrowed("?"),
        }
    }

    fn number(&self, row: usize, column: usize) -> Option<f64> {
        match self.value(row, column)? {
            CifValueRef::Text(text) => reduce3::cif::parse_f64(text.as_str()),
            CifValueRef::Integer(number, _) => Some(number as f64),
            CifValueRef::Float(number, _, _) => Some(number),
            CifValueRef::Unknown | CifValueRef::NotApplicable => None,
        }
    }

    fn tags(&self) -> Vec<Cow<'_, str>> {
        self.category
            .tags()
            .into_iter()
            .map(Cow::Borrowed)
            .collect()
    }

    fn is_loop(&self) -> bool {
        self.category.is_loop()
    }

    fn missing(&self, row: usize, column: usize) -> Option<CifCell<'static>> {
        match self.value(row, column) {
            Some(CifValueRef::Unknown) | None => Some(CifCell::Unknown),
            Some(CifValueRef::NotApplicable) => Some(CifCell::NotApplicable),
            Some(CifValueRef::Text(_) | CifValueRef::Integer(..) | CifValueRef::Float(..)) => None,
        }
    }

    fn cell_with_missing(
        &self,
        row: usize,
        column: usize,
    ) -> (Option<CifCell<'static>>, Cow<'_, str>) {
        match self.value(row, column) {
            Some(CifValueRef::Text(text)) => (None, Cow::Borrowed(text.as_str())),
            Some(CifValueRef::Integer(value, original)) => (
                None,
                original.map_or_else(
                    || Cow::Owned(value.to_string()),
                    |text| Cow::Borrowed(text.as_str()),
                ),
            ),
            Some(CifValueRef::Float(value, _, original)) => (
                None,
                original.map_or_else(
                    || Cow::Owned(value.to_string()),
                    |text| Cow::Borrowed(text.as_str()),
                ),
            ),
            Some(CifValueRef::NotApplicable) => (Some(CifCell::NotApplicable), Cow::Borrowed(".")),
            Some(CifValueRef::Unknown) | None => (Some(CifCell::Unknown), Cow::Borrowed("?")),
        }
    }
}

/// Builds a native one-block document from Reduce3's output.
///
/// Loops Reduce3 writes are stored column-wise with one string column per item, the
/// representation BinaryCIF decoding uses; present values are unquoted text. With a
/// source block, categories Reduce3 leaves unchanged are copied from it as they are.
#[derive(Debug, Default)]
pub struct DocumentSink<'a> {
    source: Option<&'a CifBlock>,
    entries: Vec<CifEntry>,
    open: Option<OpenLoop>,
    open_source: Option<OpenSourceLoop<'a>>,
    failure: Option<&'static str>,
    // Whether to record the source rows of a rewritten `_atom_site` loop in `copies`.
    record: bool,
    copies: Option<SourceCopies>,
}

/// The source row that each row of a rewritten `_atom_site` loop copies, and the columns
/// in which it does not hold that row's value.
///
/// A row's source row is the row named by its first value that Reduce3 copies from the
/// source loop. Bit `c` of the row's mask is set when column `c` holds any other value.
/// Columns from 64 on have no bit, and readers treat them as not copied.
#[derive(Clone, Debug, Default)]
pub(crate) struct SourceCopies {
    /// The source row of each output row, or [`SourceCopies::NONE`].
    pub(crate) rows: Vec<u32>,
    /// The columns of each output row that do not copy its source row.
    pub(crate) made: Vec<u64>,
}

impl SourceCopies {
    /// The source row of a row that copies no value from the source loop.
    pub(crate) const NONE: u32 = u32::MAX;

    fn new(rows: usize) -> Self {
        Self {
            rows: vec![Self::NONE; rows],
            made: vec![0; rows],
        }
    }

    /// Record the values of column `column`, one per row.
    fn note(&mut self, column: usize, values: &[ColumnValue]) {
        let Some(bit) = u32::try_from(column).ok().filter(|bit| *bit < u64::BITS) else {
            return;
        };
        let mask = 1_u64 << bit;
        for ((row, made), value) in self.rows.iter_mut().zip(&mut self.made).zip(values) {
            match *value {
                ColumnValue::Source(source) if *row == Self::NONE => {
                    match u32::try_from(source) {
                        // No earlier column copied a value, so all of them are made.
                        Ok(source) if source != Self::NONE => {
                            *row = source;
                            *made |= mask - 1;
                        }
                        _ => *made |= mask,
                    }
                }
                ColumnValue::Source(source) if u32::try_from(source) == Ok(*row) => {}
                _ => *made |= mask,
            }
        }
    }
}

/// Return the loop of `block` that a sink rewrites when Reduce3 sends `category` column
/// by column.
pub(crate) fn source_loop<'a>(block: &'a CifBlock, category: &str) -> Option<&'a CifLoop> {
    match BlockCategory::find(block.entries(), category) {
        Some(BlockCategory::Loop(source)) => Some(source),
        _ => None,
    }
}

#[derive(Debug)]
struct OpenLoop {
    tags: Vec<String>,
    columns: Vec<TextColumnBuilder>,
    rows: usize,
}

/// A loop that rewrites a source loop, received column by column.
#[derive(Debug)]
struct OpenSourceLoop<'a> {
    source: &'a CifLoop,
    tags: Vec<String>,
    columns: Vec<LoopColumn>,
    rows: usize,
    copies: Option<SourceCopies>,
}

const TOO_LARGE: &str = "a column of the output model exceeds 4 GiB of text";

/// Build one output column from the values of `column` of `source` that `values` name.
fn source_column(
    source: &CifLoop,
    column: usize,
    values: &[ColumnValue],
    text: &str,
) -> Result<LoopColumn, &'static str> {
    let made = |start: usize, end: usize| {
        text.get(start..end)
            .map(TextCell::Text)
            .ok_or("Reduce3 named text outside the text of a loop column")
    };
    // A text column shares its entries; any other column is read value by value.
    let made_count = values
        .iter()
        .filter(|value| matches!(value, ColumnValue::Made(..)))
        .count();
    if let Some(mut builder) = source.text_column(column).and_then(|text_column| {
        SourceColumnBuilder::new(text_column, values.len(), made_count, text.len())
    }) {
        for &value in values {
            match value {
                ColumnValue::Source(row) => builder.push_source(row),
                ColumnValue::Made(start, end) => {
                    builder.push(made(start, end)?).ok_or(TOO_LARGE)?
                }
                ColumnValue::Unknown => builder.push(TextCell::Unknown).ok_or(TOO_LARGE)?,
                ColumnValue::NotApplicable => {
                    builder.push(TextCell::NotApplicable).ok_or(TOO_LARGE)?
                }
            }
        }
        return Ok(builder.finish());
    }
    use std::fmt::Write as _;
    let mut builder = TextColumnBuilder::with_capacity(values.len());
    let mut number = String::new();
    for &value in values {
        let cell = match value {
            ColumnValue::Source(row) => match source.value(row, column) {
                Some(CifValueRef::Text(text)) => TextCell::Text(text.as_str()),
                Some(CifValueRef::Integer(value, original)) => match original {
                    Some(text) => TextCell::Text(text.as_str()),
                    None => {
                        number.clear();
                        let _ = write!(number, "{value}");
                        TextCell::Text(&number)
                    }
                },
                Some(CifValueRef::Float(value, _, original)) => match original {
                    Some(text) => TextCell::Text(text.as_str()),
                    None => {
                        number.clear();
                        let _ = write!(number, "{value}");
                        TextCell::Text(&number)
                    }
                },
                Some(CifValueRef::NotApplicable) => TextCell::NotApplicable,
                Some(CifValueRef::Unknown) | None => TextCell::Unknown,
            },
            ColumnValue::Made(start, end) => made(start, end)?,
            ColumnValue::Unknown => TextCell::Unknown,
            ColumnValue::NotApplicable => TextCell::NotApplicable,
        };
        builder.push(cell).ok_or(TOO_LARGE)?;
    }
    Ok(builder.finish())
}

const fn text_cell(value: CifCell<'_>) -> TextCell<'_> {
    match value {
        CifCell::Text(text) => TextCell::Text(text),
        CifCell::NotApplicable => TextCell::NotApplicable,
        CifCell::Unknown => TextCell::Unknown,
    }
}

fn owned_value(value: CifCell<'_>) -> CifValue {
    match value {
        CifCell::Text(text) => CifValue::text(text),
        CifCell::Unknown => CifValue::Unknown,
        CifCell::NotApplicable => CifValue::NotApplicable,
    }
}

impl<'a> DocumentSink<'a> {
    /// An empty sink.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// An empty sink that copies unchanged categories from `source`.
    #[must_use]
    pub fn with_source(source: &'a CifBlock) -> Self {
        Self {
            source: Some(source),
            ..Self::default()
        }
    }

    fn close_loop(&mut self) {
        if let Some(open) = self.open_source.take()
            && open.rows > 0
        {
            if open.columns.len() == open.tags.len() {
                if open.copies.is_some() {
                    self.copies = open.copies;
                }
                self.entries.push(CifEntry::Loop(CifLoop::from_columns(
                    open.tags,
                    open.columns,
                    open.rows,
                )));
            } else {
                self.failure = Some("Reduce3 wrote a loop without all of its columns");
            }
        }
        if let Some(open) = self.open.take()
            && open.rows > 0
        {
            let columns = open
                .columns
                .into_iter()
                .map(TextColumnBuilder::finish)
                .collect();
            self.entries.push(CifEntry::Loop(CifLoop::from_columns(
                open.tags, columns, open.rows,
            )));
        }
    }

    /// Return the document, with `code` as its data block code. Loops without rows are
    /// left out.
    ///
    /// # Errors
    ///
    /// Returns [`ReduceErrorCode::ResultTooLarge`] when a column outgrew 32-bit offsets
    /// or the rows did not match the loops.
    pub fn finish(mut self, code: &str) -> Result<CifDocument, ReduceError> {
        self.close_loop();
        if let Some(failure) = self.failure {
            return Err(ReduceError::new(ReduceErrorCode::ResultTooLarge, failure));
        }
        Ok(CifDocument::new(vec![CifBlock::data(
            code.to_owned(),
            self.entries,
        )]))
    }
}

impl CifSink for DocumentSink<'_> {
    fn begin_block(&mut self, _code: &str) {}

    fn copy_category(&mut self, category: &str) -> bool {
        let Some(source) = self.source else {
            return false;
        };
        self.close_loop();
        let before = self.entries.len();
        for entry in source.entries() {
            let copy = match entry {
                CifEntry::Loop(cif_loop) => cif_loop
                    .tags()
                    .first()
                    .is_some_and(|tag| in_category(tag, category)),
                CifEntry::Item(item) => in_category(item.tag(), category),
                CifEntry::Frame(_) => false,
            };
            if copy {
                self.entries.push(entry.clone());
            }
        }
        self.entries.len() > before
    }

    fn item(&mut self, tag: &str, value: CifCell<'_>) {
        self.close_loop();
        self.entries.push(CifEntry::Item(CifItem::new(
            tag.to_owned(),
            owned_value(value),
        )));
    }

    fn begin_loop(&mut self, tags: &[&str], rows: usize) {
        self.close_loop();
        self.open = Some(OpenLoop {
            tags: tags.iter().map(|tag| (*tag).to_owned()).collect(),
            columns: tags
                .iter()
                .map(|_| TextColumnBuilder::with_capacity(rows))
                .collect(),
            rows: 0,
        });
    }

    fn row(&mut self, values: &[CifCell<'_>]) {
        let Some(open) = self.open.as_mut() else {
            self.failure = Some("Reduce3 wrote a loop row outside a loop");
            return;
        };
        if values.len() != open.columns.len() {
            self.failure = Some("Reduce3 wrote a loop row with the wrong number of values");
            return;
        }
        for (column, &value) in open.columns.iter_mut().zip(values) {
            if column.push(text_cell(value)).is_none() {
                self.failure = Some(TOO_LARGE);
            }
        }
        open.rows += 1;
    }

    fn end_loop(&mut self) {
        self.close_loop();
    }

    fn begin_source_loop(&mut self, category: &str, tags: &[&str], rows: usize) -> bool {
        let Some(source) = self.source.and_then(|block| source_loop(block, category)) else {
            return false;
        };
        if source.tags().len() != tags.len()
            || source
                .tags()
                .iter()
                .zip(tags)
                .any(|(source_tag, tag)| source_tag != tag)
        {
            return false;
        }
        self.close_loop();
        self.open_source = Some(OpenSourceLoop {
            source,
            tags: tags.iter().map(|tag| (*tag).to_owned()).collect(),
            columns: Vec::with_capacity(tags.len()),
            rows,
            copies: (self.record && category == "atom_site").then(|| SourceCopies::new(rows)),
        });
        true
    }

    fn source_column(&mut self, column: usize, values: &[ColumnValue], text: &str) {
        let Some(open) = self.open_source.as_mut() else {
            self.failure = Some("Reduce3 wrote a loop column outside a loop");
            return;
        };
        if column != open.columns.len() || column >= open.tags.len() || values.len() != open.rows {
            self.failure = Some("Reduce3 wrote a loop column out of place");
            return;
        }
        match source_column(open.source, column, values, text) {
            Ok(built) => open.columns.push(built),
            Err(failure) => self.failure = Some(failure),
        }
        if let Some(copies) = open.copies.as_mut() {
            copies.note(column, values);
        }
    }
}

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

use reduce3::cifsource::{CifCell, CifSink, CifSource, CifTable};
use reduce3::model::Structure;

pub use reduce3::hplace::NTermCharge;
pub use reduce3::optimizer::OptParams;
pub use reduce3::probe::ProbeParams;
pub use reduce3::{Approach, MonLib, Params};

use crate::cif::{
    BlockCategory, CifBlock, CifDocument, CifEntry, CifItem, CifLoop, CifValue, CifValueRef,
    TextCell, TextColumnBuilder, in_category, split_tag,
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
    let block = document
        .blocks()
        .iter()
        .find(|block| BlockSource::new(block).table("atom_site").is_some())
        .ok_or_else(|| {
            ReduceError::new(
                ReduceErrorCode::InvalidModel,
                "no block has an _atom_site loop",
            )
        })?;
    run_block(block, monomers, params)
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
    let structure = reduce3::mmcif::structure_from_cif(&BlockSource::new(block))
        .map_err(|message| ReduceError::new(ReduceErrorCode::InvalidModel, message))?;
    if structure.atoms_size() == 0 {
        return Err(ReduceError::new(
            ReduceErrorCode::InvalidModel,
            "the _atom_site loop has no atoms",
        ));
    }
    let output = reduce3::pipeline::run(structure, monomers, params)
        .map_err(|message| ReduceError::new(ReduceErrorCode::Failed, message))?;
    let document = if params.compat {
        structure_document(&output.structure, block.code().unwrap_or("default"))?
    } else {
        block_document(&output.structure, block)?
    };
    Ok(Reduction {
        document,
        report: output.description,
    })
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
    let code = block.code().unwrap_or("default");
    let mut sink = DocumentSink::with_source(block);
    reduce3::mmcif::write_cif_preserving(structure, &BlockSource::new(block), code, &mut sink)
        .map_err(|message| ReduceError::new(ReduceErrorCode::InvalidModel, message))?;
    sink.finish(code)
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
        BlockCategory::find(self.block.entries(), category).map(BlockTable)
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
pub struct BlockTable<'a>(BlockCategory<'a>);

impl CifTable for BlockTable<'_> {
    fn row_count(&self) -> usize {
        self.0.row_count()
    }

    fn column(&self, tag: &str) -> Option<usize> {
        self.0.column(tag)
    }

    fn cell(&self, row: usize, column: usize) -> Cow<'_, str> {
        self.0
            .text(row, column)
            .unwrap_or(Cow::Borrowed(match self.0.value(row, column) {
                Some(CifValueRef::NotApplicable) => ".",
                _ => "?",
            }))
    }

    fn number(&self, row: usize, column: usize) -> Option<f64> {
        match self.0.value(row, column)? {
            CifValueRef::Text(text) => reduce3::cif::parse_f64(text.as_str()),
            CifValueRef::Integer(number, _) => Some(number as f64),
            CifValueRef::Float(number, _, _) => Some(number),
            CifValueRef::Unknown | CifValueRef::NotApplicable => None,
        }
    }

    fn tags(&self) -> Vec<Cow<'_, str>> {
        self.0.tags().into_iter().map(Cow::Borrowed).collect()
    }

    fn is_loop(&self) -> bool {
        self.0.is_loop()
    }

    fn missing(&self, row: usize, column: usize) -> Option<CifCell<'static>> {
        match self.0.value(row, column) {
            Some(CifValueRef::Unknown) | None => Some(CifCell::Unknown),
            Some(CifValueRef::NotApplicable) => Some(CifCell::NotApplicable),
            Some(CifValueRef::Text(_) | CifValueRef::Integer(..) | CifValueRef::Float(..)) => None,
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
    failure: Option<&'static str>,
}

#[derive(Debug)]
struct OpenLoop {
    tags: Vec<String>,
    columns: Vec<TextColumnBuilder>,
    rows: usize,
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
                self.failure = Some("a column of the output model exceeds 4 GiB of text");
            }
        }
        open.rows += 1;
    }

    fn end_loop(&mut self) {
        self.close_loop();
    }
}

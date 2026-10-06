//! Hydrogen addition and optimization of parsed models with Reduce3.
//!
//! [`run`] takes the first data block of a [`CifDocument`] that has an `_atom_site`
//! loop, runs [Reduce3](https://github.com/jamaliki/reduce3) on it in memory, and
//! returns the result as a new document. No text is written or parsed again: Reduce3
//! reads the parsed values through [`BlockSource`] and streams its output model into a
//! [`DocumentSink`].
//!
//! The result block keeps the source block code. Its items and loops are the ones the
//! `reduce3` program writes to an mmCIF file (cell, space group, `_struct_asym`,
//! `_chem_comp`, `_atom_site`, and `_atom_site_anisotrop`), held as text values, so the
//! entries equal a parse of that file. Other source categories are not carried over,
//! because Reduce3 renumbers atoms and reassigns label asym identifiers.
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
    CifBlock, CifDocument, CifEntry, CifItem, CifLoop, CifValue, CifValueRef, ColumnValues,
    LoopColumn, StringColumn, split_tag,
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
    Ok(Reduction {
        document: structure_document(&output.structure, block.code().unwrap_or("default"))?,
        report: output.description,
    })
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

fn in_category(tag: &str, category: &str) -> bool {
    split_tag(tag).is_some_and(|(name, _)| name.eq_ignore_ascii_case(category))
}

fn item_is(tag: &str, item: &str) -> bool {
    split_tag(tag).is_some_and(|(_, name)| name.eq_ignore_ascii_case(item))
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
        let mut items = Vec::new();
        for entry in self.block.entries() {
            match entry {
                CifEntry::Loop(cif_loop)
                    if cif_loop
                        .tags()
                        .first()
                        .is_some_and(|tag| in_category(tag, category)) =>
                {
                    return Some(BlockTable(TableKind::Loop(cif_loop)));
                }
                CifEntry::Item(item) if in_category(item.tag(), category) => items.push(item),
                CifEntry::Loop(_) | CifEntry::Item(_) | CifEntry::Frame(_) => {}
            }
        }
        (!items.is_empty()).then_some(BlockTable(TableKind::Items(items)))
    }
}

/// One category of a [`BlockSource`]: a loop, or the category's scalar items as one row.
#[derive(Clone, Debug)]
pub struct BlockTable<'a>(TableKind<'a>);

#[derive(Clone, Debug)]
enum TableKind<'a> {
    Loop(&'a CifLoop),
    Items(Vec<&'a CifItem>),
}

impl<'a> BlockTable<'a> {
    fn value(&self, row: usize, column: usize) -> Option<CifValueRef<'a>> {
        match &self.0 {
            TableKind::Loop(cif_loop) => cif_loop.value(row, column),
            TableKind::Items(items) => (row == 0)
                .then(|| items.get(column).map(|item| item.value().as_ref()))
                .flatten(),
        }
    }
}

impl CifTable for BlockTable<'_> {
    fn row_count(&self) -> usize {
        match &self.0 {
            TableKind::Loop(cif_loop) => cif_loop.row_count(),
            TableKind::Items(_) => 1,
        }
    }

    fn column(&self, tag: &str) -> Option<usize> {
        match &self.0 {
            TableKind::Loop(cif_loop) => cif_loop.tags().iter().position(|t| item_is(t, tag)),
            TableKind::Items(items) => items.iter().position(|item| item_is(item.tag(), tag)),
        }
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
}

/// Builds a native one-block document from Reduce3's output.
///
/// Loops are stored column-wise with one string column per item, the representation
/// BinaryCIF decoding uses; present values are unquoted text.
#[derive(Debug, Default)]
pub struct DocumentSink {
    entries: Vec<CifEntry>,
    open: Option<OpenLoop>,
    failure: Option<&'static str>,
}

#[derive(Debug)]
struct OpenLoop {
    tags: Vec<String>,
    columns: Vec<ColumnBuilder>,
    rows: usize,
}

#[derive(Debug)]
struct ColumnBuilder {
    data: String,
    offsets: Vec<u32>,
    indices: Vec<u32>,
    mask: Option<Vec<u8>>,
}

impl ColumnBuilder {
    fn with_capacity(rows: usize) -> Self {
        let mut offsets = Vec::with_capacity(rows.min(1 << 20) + 1);
        offsets.push(0);
        Self {
            data: String::new(),
            offsets,
            indices: Vec::with_capacity(rows),
            mask: None,
        }
    }

    fn last_entry(&self) -> Option<&str> {
        let [.., start, end] = self.offsets[..] else {
            return None;
        };
        self.data.get(start as usize..end as usize)
    }

    /// Append one value; `None` when the column outgrows 32-bit offsets.
    fn push(&mut self, value: CifCell<'_>) -> Option<()> {
        let row = self.indices.len();
        let (index, kind) = match value {
            CifCell::Text("") => (0, 0),
            CifCell::Text(text) => {
                if self.last_entry() != Some(text) {
                    self.data.push_str(text);
                    self.offsets.push(u32::try_from(self.data.len()).ok()?);
                }
                (u32::try_from(self.offsets.len() - 1).ok()?, 0)
            }
            CifCell::NotApplicable => (0, 1),
            CifCell::Unknown => (0, 2),
        };
        if kind != 0 || self.mask.is_some() {
            self.mask.get_or_insert_with(|| vec![0; row]).push(kind);
        }
        self.indices.push(index);
        Some(())
    }

    fn finish(self) -> LoopColumn {
        LoopColumn {
            values: ColumnValues::Strings(StringColumn::new(self.data, self.offsets, self.indices)),
            mask: self.mask,
        }
    }
}

fn owned_value(value: CifCell<'_>) -> CifValue {
    match value {
        CifCell::Text(text) => CifValue::text(text),
        CifCell::Unknown => CifValue::Unknown,
        CifCell::NotApplicable => CifValue::NotApplicable,
    }
}

impl DocumentSink {
    /// An empty sink.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    fn close_loop(&mut self) {
        if let Some(open) = self.open.take()
            && open.rows > 0
        {
            let columns = open
                .columns
                .into_iter()
                .map(ColumnBuilder::finish)
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

impl CifSink for DocumentSink {
    fn begin_block(&mut self, _code: &str) {}

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
                .map(|_| ColumnBuilder::with_capacity(rows))
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
            if column.push(value).is_none() {
                self.failure = Some("a column of the output model exceeds 4 GiB of text");
            }
        }
        open.rows += 1;
    }

    fn end_loop(&mut self) {
        self.close_loop();
    }
}

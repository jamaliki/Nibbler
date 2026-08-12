//! Explicit confidence-to-B-factor compatibility views.

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{self, Display, Formatter};

use crate::cif::{
    BlockKind, CifBlock, CifDocument, CifEntry, CifLoop, CifRow, CifValue, CifValueRef,
};

use super::aggregate::ModelCifModel;
use super::qa::QaValue;
use super::writer::canonical_document;

/// Explicit policy for the viewer-compatibility B-factor view.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MirrorPolicy {
    /// Preserve `_atom_site.B_iso_or_equiv` exactly as supplied.
    Disabled,
    /// Copy one explicitly selected local QA metric into B-factor output values.
    LocalMetric(i64),
}

/// A failure to construct the requested ModelCIF compatibility view.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct WriteError {
    code: &'static str,
    message: String,
}

impl WriteError {
    fn new(code: &'static str, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }

    /// Return the stable machine-readable error code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.code
    }

    /// Return the human-readable diagnostic message.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }
}

impl Display for WriteError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {}", self.code, self.message)
    }
}

impl Error for WriteError {}

/// Construct canonical ModelCIF output with an explicit compatibility policy.
///
/// The selected metric remains present in `_ma_qa_metric*`; mirroring never changes its
/// definition or claims that confidence is an experimental displacement parameter.
///
/// # Errors
///
/// Returns an error if the selected metric is absent, is not local, or has duplicate
/// values for one model/residue identity.
pub fn canonical_document_with_mirror(
    model: &ModelCifModel,
    policy: MirrorPolicy,
) -> Result<CifDocument, WriteError> {
    let document = canonical_document(model);
    match policy {
        MirrorPolicy::Disabled => Ok(document),
        MirrorPolicy::LocalMetric(metric_id) => mirror_local_metric(document, model, metric_id),
    }
}

fn mirror_local_metric(
    document: CifDocument,
    model: &ModelCifModel,
    metric_id: i64,
) -> Result<CifDocument, WriteError> {
    let Some(metric) = model
        .qa_metrics()
        .iter()
        .find(|metric| metric.id() == metric_id)
    else {
        return Err(WriteError::new(
            "MODELCIF_MIRROR_METRIC",
            format!("QA metric {metric_id} does not exist"),
        ));
    };
    if !metric.mode().eq_ignore_ascii_case("local") {
        return Err(WriteError::new(
            "MODELCIF_MIRROR_MODE",
            format!(
                "QA metric {metric_id} has mode {:?}, not local",
                metric.mode()
            ),
        ));
    }
    let mut values = BTreeMap::new();
    for qa in model.qa_values() {
        let QaValue::Local {
            model_id,
            site,
            metric_id: value_metric_id,
            value,
            ..
        } = qa
        else {
            continue;
        };
        if *value_metric_id != metric_id {
            continue;
        }
        let key = (
            *model_id,
            fold(site.asym_id()),
            site.sequence_id(),
            fold(site.component_id()),
        );
        if values.insert(key, *value).is_some() {
            return Err(WriteError::new(
                "MODELCIF_MIRROR_DUPLICATE",
                format!("metric {metric_id} has duplicate values for one model/residue"),
            ));
        }
    }
    if values.is_empty() {
        return Err(WriteError::new(
            "MODELCIF_MIRROR_VALUES",
            format!("QA metric {metric_id} has no local values"),
        ));
    }
    let blocks = document
        .blocks()
        .iter()
        .map(|block| mirror_block(block, &values))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(CifDocument::new(blocks))
}

type ResidueValues = BTreeMap<(i64, String, i64, String), f64>;

fn mirror_block(block: &CifBlock, values: &ResidueValues) -> Result<CifBlock, WriteError> {
    let entries = block
        .entries()
        .iter()
        .map(|entry| match entry {
            CifEntry::Loop(cif_loop)
                if loop_category(cif_loop)
                    .is_some_and(|value| value.eq_ignore_ascii_case("atom_site")) =>
            {
                mirror_atom_loop(cif_loop, values).map(CifEntry::Loop)
            }
            _ => Ok(entry.clone()),
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(match block.kind() {
        BlockKind::Data => CifBlock::data(block.code().unwrap_or_default().to_owned(), entries),
        BlockKind::Global => CifBlock::global(entries),
    })
}

fn mirror_atom_loop(cif_loop: &CifLoop, values: &ResidueValues) -> Result<CifLoop, WriteError> {
    let asym = column(cif_loop, "label_asym_id")?;
    let sequence = column(cif_loop, "label_seq_id")?;
    let component = column(cif_loop, "label_comp_id")?;
    let model = column(cif_loop, "pdbx_pdb_model_num")?;
    let b_iso = optional_column(cif_loop, "b_iso_or_equiv");
    let mut tags = cif_loop.tags().to_vec();
    let output_b_iso = b_iso.unwrap_or_else(|| {
        tags.push("_atom_site.B_iso_or_equiv".to_owned());
        tags.len() - 1
    });
    let mut output = Vec::with_capacity(cif_loop.row_count() * tags.len());
    for row_index in 0..cif_loop.row_count() {
        let row = cif_loop.row(row_index).ok_or_else(|| {
            WriteError::new(
                "MODELCIF_MIRROR_ATOM_ROW",
                "atom-site loop contains an incomplete row",
            )
        })?;
        let mirrored = match optional_integer(row_value(row, sequence)?, "_atom_site.label_seq_id")?
        {
            Some(sequence_id) => {
                let key = (
                    integer(row_value(row, model)?, "_atom_site.pdbx_PDB_model_num")?,
                    fold(text_value(
                        row_value(row, asym)?,
                        "_atom_site.label_asym_id",
                    )?),
                    sequence_id,
                    fold(text_value(
                        row_value(row, component)?,
                        "_atom_site.label_comp_id",
                    )?),
                );
                values.get(&key).copied()
            }
            None => None,
        };
        for output_column in 0..tags.len() {
            if output_column == output_b_iso {
                let value = if let Some(value) = mirrored {
                    CifValue::Float(value, None, None)
                } else if let Some(source) = b_iso {
                    row_value(row, source)?.to_owned()
                } else {
                    CifValue::Unknown
                };
                output.push(value);
            } else {
                output.push(row_value(row, output_column)?.to_owned());
            }
        }
    }
    Ok(CifLoop::new(tags, output))
}

fn column(cif_loop: &CifLoop, item: &str) -> Result<usize, WriteError> {
    optional_column(cif_loop, item).ok_or_else(|| {
        WriteError::new(
            "MODELCIF_MIRROR_ATOM_ITEM",
            format!("atom-site loop lacks _atom_site.{item}"),
        )
    })
}

fn optional_column(cif_loop: &CifLoop, item: &str) -> Option<usize> {
    cif_loop
        .tags()
        .iter()
        .position(|tag| item_name(tag).is_some_and(|value| value.eq_ignore_ascii_case(item)))
}

fn row_value(row: CifRow<'_>, column: usize) -> Result<CifValueRef<'_>, WriteError> {
    row.get(column).ok_or_else(|| {
        WriteError::new(
            "MODELCIF_MIRROR_ATOM_ROW",
            "atom-site loop contains an incomplete row",
        )
    })
}

fn integer(value: CifValueRef<'_>, tag: &str) -> Result<i64, WriteError> {
    match value {
        CifValueRef::Integer(value, _) => Ok(value),
        CifValueRef::Text(value) => value.as_str().parse().map_err(|_| {
            WriteError::new(
                "MODELCIF_MIRROR_ATOM_TYPE",
                format!("{tag} is not an integer"),
            )
        }),
        CifValueRef::Float(..) | CifValueRef::Unknown | CifValueRef::NotApplicable => {
            Err(WriteError::new(
                "MODELCIF_MIRROR_ATOM_TYPE",
                format!("{tag} is not a present integer"),
            ))
        }
    }
}

fn optional_integer(value: CifValueRef<'_>, tag: &str) -> Result<Option<i64>, WriteError> {
    match value {
        CifValueRef::Unknown | CifValueRef::NotApplicable => Ok(None),
        _ => integer(value, tag).map(Some),
    }
}

fn text_value<'a>(value: CifValueRef<'a>, tag: &str) -> Result<&'a str, WriteError> {
    value.as_text().ok_or_else(|| {
        WriteError::new(
            "MODELCIF_MIRROR_ATOM_TYPE",
            format!("{tag} is not present text"),
        )
    })
}

fn loop_category(cif_loop: &CifLoop) -> Option<&str> {
    category(cif_loop.tags().first()?.as_str())
}

fn category(tag: &str) -> Option<&str> {
    tag.strip_prefix('_')?.split_once('.').map(|(name, _)| name)
}

fn item_name(tag: &str) -> Option<&str> {
    tag.strip_prefix('_')?.split_once('.').map(|(_, item)| item)
}

fn fold(value: &str) -> String {
    value.to_ascii_lowercase()
}

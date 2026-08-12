use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{self, Display, Formatter};
use std::sync::Arc;

use super::SchemaName;
use super::error::ParseError;
use super::parser::{ParseOptions, parse_source_into};
use super::projection_sink::TableSink;
use super::schema::loaded_schema;
use super::source::SourceBuffer;
use super::table::{CifCell, CifTable, ColumnType, MissingKind};

/// One native predicate evaluated while parsing projected rows.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Predicate {
    /// Retain rows whose present text equals the operand.
    Equal {
        /// Item name to inspect.
        column: String,
        /// Required present text.
        value: String,
    },
    /// Retain rows whose present text differs from the operand.
    NotEqual {
        /// Item name to inspect.
        column: String,
        /// Present text to reject.
        value: String,
    },
    /// Retain rows whose present text belongs to the operand set.
    In {
        /// Item name to inspect.
        column: String,
        /// Accepted present-text values.
        values: Vec<String>,
    },
    /// Retain rows with either or one precise missing state.
    Missing {
        /// Item name to inspect.
        column: String,
        /// Precise missing state, or `None` to accept either.
        kind: Option<MissingKind>,
    },
}

impl Predicate {
    pub(super) fn column(&self) -> &str {
        match self {
            Self::Equal { column, .. }
            | Self::NotEqual { column, .. }
            | Self::In { column, .. }
            | Self::Missing { column, .. } => column,
        }
    }
}

/// An immutable, validated category projection.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ProjectionPlan {
    pub(super) category: String,
    pub(super) columns: Option<Vec<ColumnSpec>>,
    pub(super) predicates: Vec<PlannedPredicate>,
    pub(super) schema_types: Option<Arc<BTreeMap<String, ColumnType>>>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ColumnSpec {
    pub(super) name: String,
    pub(super) key: String,
    pub(super) column_type: ColumnType,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct PlannedPredicate {
    pub(super) predicate: Predicate,
    pub(super) column: ColumnSpec,
}

pub(super) struct ProjectionColumns {
    pub(super) columns: Vec<ColumnSpec>,
    pub(super) predicate_columns: Vec<usize>,
}

impl ColumnSpec {
    pub(super) fn new(name: &str, column_type: ColumnType) -> Self {
        Self {
            name: name.to_owned(),
            key: name.to_ascii_lowercase(),
            column_type,
        }
    }
}

impl ProjectionPlan {
    /// Select one CIF category. A leading underscore is accepted and removed.
    ///
    /// # Errors
    ///
    /// Returns [`ProjectionErrorCode::InvalidCategory`] for an empty category or a
    /// value containing an item separator (`.`).
    pub fn new(category: impl Into<String>) -> Result<Self, ProjectionError> {
        let category = normalize_category(&category.into())?;
        Ok(Self {
            category,
            columns: None,
            predicates: Vec::new(),
            schema_types: None,
        })
    }

    /// Restrict output to item names in exactly this order.
    ///
    /// # Errors
    ///
    /// Rejects empty, duplicate, or cross-category selectors.
    pub fn with_columns<I, S>(mut self, columns: I) -> Result<Self, ProjectionError>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let mut specs = Vec::new();
        for column in columns {
            let mut spec = normalize_column(&self.category, &column.into())?;
            if let Some(types) = &self.schema_types {
                spec.column_type = schema_column_type(types, &self.category, &spec, "item")?;
            }
            if specs
                .iter()
                .any(|existing: &ColumnSpec| existing.key == spec.key)
            {
                return Err(ProjectionError::new(
                    ProjectionErrorCode::DuplicateColumn,
                    format!("duplicate projected column {:?}", spec.name),
                ));
            }
            specs.push(spec);
        }
        self.columns = Some(specs);
        Ok(self)
    }

    /// Add a predicate evaluated natively before a row is retained.
    ///
    /// # Errors
    ///
    /// Rejects a predicate column outside this projection's category.
    pub fn with_predicate(mut self, predicate: Predicate) -> Result<Self, ProjectionError> {
        let column = normalize_column(&self.category, predicate.column())?;
        if let Some(types) = &self.schema_types {
            schema_column_type(types, &self.category, &column, "predicate item")?;
        }
        self.predicates.push(PlannedPredicate { predicate, column });
        Ok(self)
    }

    /// Apply one lock-pinned dictionary to selected output columns.
    ///
    /// # Errors
    ///
    /// Rejects an unknown category or explicitly selected item, and propagates embedded
    /// schema artifact failures.
    pub fn with_schema(mut self, schema: SchemaName) -> Result<Self, ProjectionError> {
        let loaded = loaded_schema(schema).map_err(|error| {
            ProjectionError::new(ProjectionErrorCode::Schema, error.to_string())
        })?;
        let dictionary = loaded.dictionary();
        let Some(category) = dictionary.category(&self.category) else {
            return Err(ProjectionError::new(
                ProjectionErrorCode::SchemaCategory,
                format!(
                    "category {:?} is not defined by {}",
                    self.category,
                    dictionary.metadata().dictionary_name()
                ),
            ));
        };
        let mut types = BTreeMap::new();
        for item_name in category.canonical_items() {
            let Some(item) = dictionary.item(item_name) else {
                continue;
            };
            let Some((_, column)) = item
                .name()
                .strip_prefix('_')
                .and_then(|name| name.split_once('.'))
            else {
                continue;
            };
            types.insert(
                column.to_ascii_lowercase(),
                column_type(dictionary, item.type_code()),
            );
        }
        if let Some(columns) = &mut self.columns {
            for column in columns {
                column.column_type = schema_column_type(&types, &self.category, column, "item")?;
            }
        }
        for predicate in &self.predicates {
            schema_column_type(&types, &self.category, &predicate.column, "predicate item")?;
        }
        self.schema_types = Some(Arc::new(types));
        Ok(self)
    }

    pub(super) fn column_type(&self, item: &str) -> Option<ColumnType> {
        self.schema_types
            .as_ref()
            .and_then(|types| types.get(&item.to_ascii_lowercase()).copied())
            .or_else(|| self.schema_types.is_none().then_some(ColumnType::Text))
    }

    pub(super) fn prepare_columns(&self, outputs: Vec<ColumnSpec>) -> ProjectionColumns {
        let mut columns = outputs;
        let mut predicate_columns = Vec::with_capacity(self.predicates.len());
        for predicate in &self.predicates {
            let column_index = match columns
                .iter()
                .position(|candidate| candidate.key == predicate.column.key)
            {
                Some(index) => index,
                None => {
                    columns.push(predicate.column.clone());
                    columns.len() - 1
                }
            };
            predicate_columns.push(column_index);
        }
        ProjectionColumns {
            columns,
            predicate_columns,
        }
    }
}

pub(super) fn predicate_matches(predicate: &Predicate, cell: &CifCell) -> bool {
    match predicate {
        Predicate::Equal { value, .. } => cell.as_text() == Some(value),
        Predicate::NotEqual { value, .. } => cell.as_text().is_some_and(|text| text != value),
        Predicate::In { values, .. } => cell
            .as_text()
            .is_some_and(|text| values.iter().any(|value| value == text)),
        Predicate::Missing { kind, .. } => cell
            .missing_kind()
            .is_some_and(|actual| kind.is_none_or(|expected| actual == expected)),
    }
}

pub(super) fn normalize_category(category: &str) -> Result<String, ProjectionError> {
    let category = category.strip_prefix('_').unwrap_or(category);
    if category.is_empty()
        || category.contains('.')
        || category.bytes().any(|byte| byte.is_ascii_whitespace())
    {
        return Err(ProjectionError::new(
            ProjectionErrorCode::InvalidCategory,
            format!("invalid CIF category {category:?}"),
        ));
    }
    Ok(category.to_ascii_lowercase())
}

pub(super) fn normalize_column(
    category: &str,
    selector: &str,
) -> Result<ColumnSpec, ProjectionError> {
    let selector = selector.strip_prefix('_').unwrap_or(selector);
    let item = if let Some((selected_category, item)) = selector.split_once('.') {
        if !selected_category.eq_ignore_ascii_case(category) {
            return Err(ProjectionError::new(
                ProjectionErrorCode::InvalidColumn,
                format!("column {selector:?} does not belong to category {category:?}"),
            ));
        }
        item
    } else {
        selector
    };
    if item.is_empty() || item.bytes().any(|byte| byte.is_ascii_whitespace()) {
        return Err(ProjectionError::new(
            ProjectionErrorCode::InvalidColumn,
            format!("invalid CIF column {selector:?}"),
        ));
    }
    Ok(ColumnSpec::new(item, ColumnType::Text))
}

fn schema_column_type(
    types: &BTreeMap<String, ColumnType>,
    category: &str,
    column: &ColumnSpec,
    role: &str,
) -> Result<ColumnType, ProjectionError> {
    types.get(&column.key).copied().ok_or_else(|| {
        ProjectionError::new(
            ProjectionErrorCode::SchemaItem,
            format!(
                "{role} _{category}.{} is not defined by the selected schema",
                column.name
            ),
        )
    })
}

fn column_type(dictionary: &super::Dictionary, type_code: &str) -> ColumnType {
    let Some(item_type) = dictionary.value_type(type_code) else {
        return ColumnType::Text;
    };
    match item_type.code().to_ascii_lowercase().as_str() {
        "int" | "positive_int" => ColumnType::Integer,
        "float" => ColumnType::Float,
        _ => ColumnType::Text,
    }
}

/// Stable categories for projection failures.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProjectionErrorCode {
    /// The selected category is empty or malformed.
    InvalidCategory,
    /// A selected item name is empty or belongs to a different category.
    InvalidColumn,
    /// A selected output column occurs more than once.
    DuplicateColumn,
    /// A native predicate operator or operand list is malformed.
    InvalidPredicate,
    /// A matching category occurrence lacks a required output or predicate column.
    MissingColumn,
    /// Matching occurrences expose incompatible implicit column layouts.
    IncompatibleColumns,
    /// The underlying CIF source is malformed.
    Parse,
    /// The embedded schema could not be loaded.
    Schema,
    /// The projected category is absent from the selected schema.
    SchemaCategory,
    /// An explicitly selected item is absent from the selected schema.
    SchemaItem,
    /// A projected value does not match its dictionary-selected physical type.
    SchemaType,
}

impl ProjectionErrorCode {
    /// Return the stable machine-readable diagnostic code.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::InvalidCategory => "CIF_PROJECTION_INVALID_CATEGORY",
            Self::InvalidColumn => "CIF_PROJECTION_INVALID_COLUMN",
            Self::DuplicateColumn => "CIF_PROJECTION_DUPLICATE_COLUMN",
            Self::InvalidPredicate => "CIF_PROJECTION_INVALID_PREDICATE",
            Self::MissingColumn => "CIF_PROJECTION_MISSING_COLUMN",
            Self::IncompatibleColumns => "CIF_PROJECTION_INCOMPATIBLE_COLUMNS",
            Self::Parse => "CIF_PROJECTION_PARSE",
            Self::Schema => "CIF_SCHEMA_ARTIFACT",
            Self::SchemaCategory => "CIF_SCHEMA_CATEGORY_UNKNOWN",
            Self::SchemaItem => "CIF_SCHEMA_ITEM_UNKNOWN",
            Self::SchemaType => "CIF_SCHEMA_TYPE",
        }
    }
}

/// A structured projection-plan, shape, or syntax failure.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ProjectionError {
    code: ProjectionErrorCode,
    message: String,
    parse_error: Option<ParseError>,
}

impl ProjectionError {
    pub(super) fn new(code: ProjectionErrorCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
            parse_error: None,
        }
    }

    #[cfg(feature = "python")]
    pub(crate) fn invalid_predicate(message: impl Into<String>) -> Self {
        Self::new(ProjectionErrorCode::InvalidPredicate, message)
    }

    #[cfg(feature = "python")]
    pub(crate) fn schema(code: &str, message: &str) -> Self {
        let code = match code {
            "CIF_SCHEMA_CATEGORY_UNKNOWN" => ProjectionErrorCode::SchemaCategory,
            "CIF_SCHEMA_ITEM_UNKNOWN" => ProjectionErrorCode::SchemaItem,
            _ => ProjectionErrorCode::Schema,
        };
        Self::new(code, message)
    }

    /// Return the stable projection error category.
    #[must_use]
    pub const fn code(&self) -> ProjectionErrorCode {
        self.code
    }

    /// Return the human-readable error detail.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }

    /// Return a located parser failure when [`ProjectionErrorCode::Parse`] applies.
    #[must_use]
    pub const fn parse_error(&self) -> Option<&ParseError> {
        self.parse_error.as_ref()
    }
}

impl From<ParseError> for ProjectionError {
    fn from(error: ParseError) -> Self {
        Self {
            code: ProjectionErrorCode::Parse,
            message: error.to_string(),
            parse_error: Some(error),
        }
    }
}

impl Display for ProjectionError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.message)
    }
}

impl Error for ProjectionError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        self.parse_error
            .as_ref()
            .map(|error| error as &(dyn Error + 'static))
    }
}

/// Parse and project one named source with strict defaults.
///
/// Only selected and predicate columns allocate per-value strings.
///
/// # Errors
///
/// Returns a syntax, resource-limit, plan, or projection-shape failure.
pub fn project_source(
    source: SourceBuffer,
    plan: ProjectionPlan,
) -> Result<CifTable, ProjectionError> {
    project_source_with_options(source, plan, ParseOptions::default())
}

/// Parse and project in-memory UTF-8 CIF bytes using strict defaults.
///
/// # Errors
///
/// Returns a syntax, resource-limit, plan, or projection-shape failure.
pub fn project(bytes: &[u8], plan: ProjectionPlan) -> Result<CifTable, ProjectionError> {
    let source = SourceBuffer::from_bytes("<memory>", bytes)?;
    project_source(source, plan)
}

/// Parse and project in-memory UTF-8 CIF bytes with explicit resource limits.
///
/// # Errors
///
/// Returns the same failures as [`project`].
pub fn project_with_options(
    bytes: &[u8],
    plan: ProjectionPlan,
    options: ParseOptions,
) -> Result<CifTable, ProjectionError> {
    if bytes.len() > options.limits.source_bytes {
        return Err(ParseError::new(
            super::error::ParseErrorCode::ResourceLimit,
            format!(
                "source exceeds the configured limit of {} bytes",
                options.limits.source_bytes
            ),
            "<memory>",
            super::error::SourceSpan::new(0, 0, 1, 1),
        )
        .into());
    }
    let source = SourceBuffer::from_bytes("<memory>", bytes)?;
    project_source_with_options(source, plan, options)
}

/// Parse and project one named source with explicit resource limits.
///
/// # Errors
///
/// Returns the same failures as [`project_source`].
pub fn project_source_with_options(
    source: SourceBuffer,
    plan: ProjectionPlan,
    options: ParseOptions,
) -> Result<CifTable, ProjectionError> {
    let mut sink = TableSink::new(source.name(), plan);
    parse_source_into(source, options, &mut sink)?;
    sink.finish()
}

//! Constraint evaluation over one logical CIF document and one loaded schema.

use std::borrow::Cow;
use std::collections::{BTreeMap, HashMap, HashSet};

use super::dictionary::{CategoryDefinition, Dictionary, ItemDefinition};
use super::dictionary_source::{CategoryColumn, CategoryView, split_tag};
use super::document::{BlockKind, CifDocument, CifEntry, CifValueRef};
use super::numeric::parse_float;
use super::schema::LoadedSchema;
use super::validation::{Diagnostic, MAX_DIAGNOSTICS, Severity};

pub(super) fn validate(document: &CifDocument, schema: &LoadedSchema) -> (Vec<Diagnostic>, bool) {
    let mut validator = Validator {
        schema,
        diagnostics: Vec::new(),
        truncated: false,
    };
    for block in document.blocks() {
        let block_context = block.code().unwrap_or("<global>");
        let enforce_categories = block.kind() == BlockKind::Data;
        validator.validate_entries(block.entries(), block_context, None, enforce_categories);
        for entry in block.entries() {
            if let CifEntry::Frame(frame) = entry {
                validator.validate_entries(
                    frame.entries(),
                    block_context,
                    Some(frame.code()),
                    false,
                );
            }
        }
    }
    (validator.diagnostics, validator.truncated)
}

struct Validator<'a> {
    schema: &'a LoadedSchema,
    diagnostics: Vec<Diagnostic>,
    truncated: bool,
}

struct Occurrence<'a> {
    category: &'a CategoryDefinition,
    columns: BTreeMap<String, Column<'a>>,
    row_count: usize,
}

struct Column<'a> {
    item: &'a ItemDefinition,
    source: CategoryColumn<'a>,
}

impl<'a> Validator<'a> {
    fn validate_entries(
        &mut self,
        entries: &'a [CifEntry],
        block: &str,
        frame: Option<&str>,
        enforce_categories: bool,
    ) {
        let occurrences = self.collect_occurrences(entries, block, frame);
        if enforce_categories {
            self.validate_mandatory_categories(&occurrences, block);
        }
        let parent_values = value_index(&occurrences, self.schema.dictionary());
        for occurrence in &occurrences {
            self.validate_occurrence(occurrence, &parent_values, block, frame);
            if self.truncated {
                return;
            }
        }
    }

    fn collect_occurrences(
        &mut self,
        entries: &'a [CifEntry],
        block: &str,
        frame: Option<&str>,
    ) -> Vec<Occurrence<'a>> {
        let dictionary = self.schema.dictionary();
        let mut scalar_columns: BTreeMap<String, BTreeMap<String, Column<'a>>> = BTreeMap::new();
        let mut occurrences = Vec::new();
        for source in CategoryView::new(entries).occurrences() {
            let mut by_category: BTreeMap<String, BTreeMap<String, Column<'a>>> = BTreeMap::new();
            for column in source.columns() {
                let Some(definition) = dictionary.item(column.tag()) else {
                    self.unknown_item(column.tag(), block, frame);
                    continue;
                };
                let columns = if source.is_scalar() {
                    scalar_columns
                        .entry(definition.category().to_ascii_lowercase())
                        .or_default()
                } else {
                    by_category
                        .entry(definition.category().to_ascii_lowercase())
                        .or_default()
                };
                self.insert_column(columns, definition, column, block, frame);
            }
            if source.is_scalar() {
                continue;
            }
            if by_category.len() > 1 {
                self.push(
                    "CIF_SCHEMA_LOOP_CATEGORY",
                    Severity::Error,
                    "a dictionary loop must contain items from one category".to_owned(),
                    context(block, frame, None, None, None, dictionary),
                );
            }
            for (category_key, columns) in by_category {
                if let Some(category) = dictionary.category(&category_key) {
                    occurrences.push(Occurrence {
                        category,
                        columns,
                        row_count: source.row_count(),
                    });
                }
            }
        }
        for (category_key, columns) in scalar_columns {
            if let Some(category) = dictionary.category(&category_key) {
                occurrences.push(Occurrence {
                    category,
                    columns,
                    row_count: 1,
                });
            }
        }
        occurrences
    }

    #[allow(clippy::too_many_arguments)]
    fn insert_column(
        &mut self,
        columns: &mut BTreeMap<String, Column<'a>>,
        item: &'a ItemDefinition,
        source: CategoryColumn<'a>,
        block: &str,
        frame: Option<&str>,
    ) {
        let key = item.name().to_ascii_lowercase();
        if columns.insert(key, Column { item, source }).is_some() {
            self.push(
                "CIF_SCHEMA_ALIAS_DUPLICATE",
                Severity::Error,
                format!(
                    "item {:?} is present through more than one name",
                    item.name()
                ),
                context(
                    block,
                    frame,
                    Some(item.category()),
                    Some(item.name()),
                    None,
                    self.schema.dictionary(),
                ),
            );
        }
    }

    fn unknown_item(&mut self, tag: &str, block: &str, frame: Option<&str>) {
        self.push(
            "CIF_SCHEMA_ITEM_UNKNOWN",
            Severity::Error,
            format!("item {tag:?} is not defined by the selected dictionary"),
            context(
                block,
                frame,
                tag_category(tag),
                Some(tag),
                None,
                self.schema.dictionary(),
            ),
        );
    }

    fn validate_mandatory_categories(&mut self, occurrences: &[Occurrence<'a>], block: &str) {
        for category in self
            .schema
            .dictionary()
            .categories()
            .filter(|category| category.mandatory())
        {
            if !occurrences
                .iter()
                .any(|occurrence| occurrence.category.id().eq_ignore_ascii_case(category.id()))
            {
                self.push(
                    "CIF_SCHEMA_CATEGORY_REQUIRED",
                    Severity::Error,
                    format!("mandatory category {:?} is absent", category.id()),
                    context(
                        block,
                        None,
                        Some(category.id()),
                        None,
                        None,
                        self.schema.dictionary(),
                    ),
                );
            }
        }
    }

    fn validate_occurrence(
        &mut self,
        occurrence: &Occurrence<'a>,
        parent_values: &HashMap<String, HashSet<String>>,
        block: &str,
        frame: Option<&str>,
    ) {
        for item_name in occurrence.category.canonical_items() {
            let Some(item) = self.schema.dictionary().item(item_name) else {
                continue;
            };
            if item.mandatory() && !occurrence.columns.contains_key(item_name) {
                self.push(
                    "CIF_SCHEMA_ITEM_REQUIRED",
                    Severity::Error,
                    format!("mandatory item {:?} is absent", item.name()),
                    context(
                        block,
                        frame,
                        Some(occurrence.category.id()),
                        Some(item.name()),
                        None,
                        self.schema.dictionary(),
                    ),
                );
            }
        }
        for column in occurrence.columns.values() {
            self.validate_column(column, block, frame);
            self.validate_parents(column, parent_values, block, frame);
        }
        self.validate_keys(occurrence, block, frame);
    }

    fn validate_column(&mut self, column: &Column<'a>, block: &str, frame: Option<&str>) {
        let dictionary = self.schema.dictionary();
        let Some(item_type) = dictionary.value_type(column.item.type_code()) else {
            return;
        };
        let pattern = self.schema.pattern(item_type.code());
        for (row_index, value) in column.source.values().enumerate() {
            let Some(text) = present_text(value) else {
                continue;
            };
            if pattern.is_some_and(|pattern| !pattern.is_match(&text)) {
                self.push(
                    "CIF_SCHEMA_TYPE",
                    Severity::Error,
                    format!(
                        "value {:?} does not match DDL2 type {:?}",
                        text,
                        item_type.code()
                    ),
                    context(
                        block,
                        frame,
                        Some(column.item.category()),
                        Some(column.source.tag()),
                        Some(row_index),
                        dictionary,
                    ),
                );
                continue;
            }
            if !column.item.enumerations().is_empty()
                && !column.item.enumerations().iter().any(|allowed| {
                    if item_type.primitive().eq_ignore_ascii_case("uchar") {
                        allowed.eq_ignore_ascii_case(&text)
                    } else {
                        allowed == text.as_ref()
                    }
                })
            {
                self.push(
                    "CIF_SCHEMA_ENUMERATION",
                    Severity::Error,
                    format!("value {:?} is not an allowed enumeration", text),
                    context(
                        block,
                        frame,
                        Some(column.item.category()),
                        Some(column.source.tag()),
                        Some(row_index),
                        dictionary,
                    ),
                );
            }
            self.validate_range(column, &text, row_index, block, frame);
        }
    }

    fn validate_range(
        &mut self,
        column: &Column<'a>,
        text: &str,
        row_index: usize,
        block: &str,
        frame: Option<&str>,
    ) {
        if column.item.ranges().is_empty() {
            return;
        }
        let Some(number) = parse_float(text) else {
            return;
        };
        let allowed = column.item.ranges().iter().any(|range| {
            let minimum = range.minimum().and_then(parse_float);
            let maximum = range.maximum().and_then(parse_float);
            minimum.is_none_or(|minimum| number >= minimum)
                && maximum.is_none_or(|maximum| number <= maximum)
        });
        if !allowed {
            self.push(
                "CIF_SCHEMA_RANGE",
                Severity::Error,
                format!("value {text:?} is outside every allowed range"),
                context(
                    block,
                    frame,
                    Some(column.item.category()),
                    Some(column.source.tag()),
                    Some(row_index),
                    self.schema.dictionary(),
                ),
            );
        }
    }

    fn validate_keys(&mut self, occurrence: &Occurrence<'a>, block: &str, frame: Option<&str>) {
        if occurrence.category.keys().is_empty() {
            return;
        }
        let mut seen = HashSet::new();
        for row_index in 0..occurrence.row_count {
            let key = occurrence
                .category
                .keys()
                .iter()
                .map(|name| {
                    let column = occurrence.columns.get(name)?;
                    let value = column.source.value(row_index)?;
                    comparison_text(value, column.item, self.schema.dictionary())
                })
                .collect::<Option<Vec<_>>>();
            let Some(key) = key else {
                self.push(
                    "CIF_SCHEMA_KEY_MISSING",
                    Severity::Error,
                    "category key contains an absent or missing value".to_owned(),
                    context(
                        block,
                        frame,
                        Some(occurrence.category.id()),
                        None,
                        Some(row_index),
                        self.schema.dictionary(),
                    ),
                );
                continue;
            };
            if !seen.insert(key) {
                self.push(
                    "CIF_SCHEMA_KEY_DUPLICATE",
                    Severity::Error,
                    "category key is not unique".to_owned(),
                    context(
                        block,
                        frame,
                        Some(occurrence.category.id()),
                        None,
                        Some(row_index),
                        self.schema.dictionary(),
                    ),
                );
            }
        }
    }

    fn validate_parents(
        &mut self,
        column: &Column<'a>,
        parent_values: &HashMap<String, HashSet<String>>,
        block: &str,
        frame: Option<&str>,
    ) {
        if column.item.parents().is_empty() {
            return;
        }
        let available_parents = column
            .item
            .parents()
            .iter()
            .filter_map(|parent| parent_values.get(parent))
            .collect::<Vec<_>>();
        for (row_index, value) in column.source.values().enumerate() {
            let Some(display_text) = present_text(value) else {
                continue;
            };
            let Some(comparison_text) =
                comparison_text(value, column.item, self.schema.dictionary())
            else {
                continue;
            };
            let found = available_parents
                .iter()
                .any(|values| values.contains(&comparison_text));
            if !found {
                self.push(
                    "CIF_SCHEMA_PARENT_MISSING",
                    Severity::Error,
                    format!("value {display_text:?} has no matching parent item value"),
                    context(
                        block,
                        frame,
                        Some(column.item.category()),
                        Some(column.source.tag()),
                        Some(row_index),
                        self.schema.dictionary(),
                    ),
                );
            }
        }
    }

    fn push(
        &mut self,
        code: &'static str,
        severity: Severity,
        message: String,
        context: Vec<String>,
    ) {
        if self.diagnostics.len() >= MAX_DIAGNOSTICS {
            self.truncated = true;
            return;
        }
        self.diagnostics.push(Diagnostic {
            code,
            severity,
            message,
            context,
        });
    }
}

fn value_index(
    occurrences: &[Occurrence<'_>],
    dictionary: &Dictionary,
) -> HashMap<String, HashSet<String>> {
    let mut index: HashMap<String, HashSet<String>> = HashMap::new();
    for column in occurrences
        .iter()
        .flat_map(|occurrence| occurrence.columns.values())
    {
        let values = index
            .entry(column.item.name().to_ascii_lowercase())
            .or_default();
        values.extend(
            column
                .source
                .values()
                .filter_map(|value| comparison_text(value, column.item, dictionary)),
        );
    }
    index
}

fn comparison_text(
    value: CifValueRef<'_>,
    item: &ItemDefinition,
    dictionary: &Dictionary,
) -> Option<String> {
    let text = present_text(value)?;
    let case_insensitive = dictionary
        .value_type(item.type_code())
        .is_some_and(|item_type| item_type.primitive().eq_ignore_ascii_case("uchar"));
    Some(if case_insensitive {
        text.to_ascii_lowercase()
    } else {
        text.into_owned()
    })
}

fn present_text(value: CifValueRef<'_>) -> Option<Cow<'_, str>> {
    match value {
        CifValueRef::Text(text) => Some(Cow::Borrowed(text.as_str())),
        CifValueRef::Integer(number, _) => Some(Cow::Owned(number.to_string())),
        CifValueRef::Float(number, _, _) => Some(Cow::Owned(number.to_string())),
        CifValueRef::Unknown | CifValueRef::NotApplicable => None,
    }
}

fn context(
    block: &str,
    frame: Option<&str>,
    category: Option<&str>,
    item: Option<&str>,
    row: Option<usize>,
    dictionary: &Dictionary,
) -> Vec<String> {
    let mut output = vec![format!("block={block}")];
    if let Some(frame) = frame {
        output.push(format!("frame={frame}"));
    }
    if let Some(category) = category {
        output.push(format!("category={category}"));
    }
    if let Some(item) = item {
        output.push(format!("item={item}"));
    }
    if let Some(row) = row {
        output.push(format!("row={}", row + 1));
    }
    output.push(format!(
        "dictionary={}@{}",
        dictionary.metadata().dictionary_name(),
        dictionary.metadata().version()
    ));
    output
}

fn tag_category(tag: &str) -> Option<&str> {
    split_tag(tag).map(|(category, _)| category)
}

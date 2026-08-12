use std::thread;

use super::document::BlockKind;
use super::error::ParseErrorCode;
use super::lexer::Lexer;
use super::parallel::{LoopKernel, LoopRange, ParallelLoop, ProjectedChunk, scan_projected_loop};
use super::projection::{
    ColumnSpec, Predicate, ProjectionError, ProjectionErrorCode, ProjectionPlan,
};
use super::sink::{ParseSink, ParsedValue};
use super::source::SourceBuffer;
use super::table::{CifCell, CifTable, ColumnBuilder, ColumnType, TableBuilder};
use super::token::TokenKind;

pub(super) struct TableSink {
    plan: ProjectionPlan,
    table: TableBuilder,
    segments: Vec<CifTable>,
    source_name: String,
    block_code: Option<String>,
    frame_code: Option<String>,
    scalar_values: Vec<(ColumnSpec, CifCell)>,
    current_loop: Option<LoopProjection>,
    pending_parallel: Option<Vec<CifTable>>,
    error: Option<ProjectionError>,
}

struct LoopProjection {
    target_by_tag: Vec<Option<usize>>,
    needed: Vec<ColumnSpec>,
    values: Vec<Option<CifCell>>,
    predicate_free: bool,
}

impl TableSink {
    pub(super) fn new(source_name: &str, plan: ProjectionPlan) -> Self {
        let columns = plan.columns.as_ref().map_or_else(Vec::new, |specs| {
            specs
                .iter()
                .map(|spec| ColumnBuilder::new(spec.name.clone(), spec.column_type))
                .collect()
        });
        Self {
            table: TableBuilder {
                category: plan.category.clone(),
                columns,
                provenance: Vec::new(),
            },
            segments: Vec::new(),
            plan,
            source_name: source_name.to_owned(),
            block_code: None,
            frame_code: None,
            scalar_values: Vec::new(),
            current_loop: None,
            pending_parallel: None,
            error: None,
        }
    }

    pub(super) fn finish(mut self) -> Result<CifTable, ProjectionError> {
        self.flush_scalars();
        if let Some(error) = self.error {
            return Err(error);
        }
        let tail = self.table.finish();
        if self.segments.is_empty() {
            return Ok(tail);
        }
        let mut segments = self.segments.into_iter();
        let Some(mut table) = segments.next() else {
            return Ok(tail);
        };
        for segment in segments {
            table.append(segment).map_err(|message| {
                ProjectionError::new(ProjectionErrorCode::IncompatibleColumns, message)
            })?;
        }
        if tail.row_count() != 0 {
            table.append(tail).map_err(|message| {
                ProjectionError::new(ProjectionErrorCode::IncompatibleColumns, message)
            })?;
        }
        Ok(table)
    }

    fn take_table(&mut self) -> CifTable {
        let replacement = self.table.empty_like();
        std::mem::replace(&mut self.table, replacement).finish()
    }

    fn flush_scalars(&mut self) {
        if self.scalar_values.is_empty() || self.error.is_some() {
            self.scalar_values.clear();
            return;
        }
        let available: Vec<ColumnSpec> = self
            .scalar_values
            .iter()
            .map(|(spec, _)| spec.clone())
            .collect();
        let Some(needed) = self.prepare_columns(&available) else {
            self.scalar_values.clear();
            return;
        };
        let mut values = needed
            .iter()
            .map(|spec| {
                self.scalar_values
                    .iter()
                    .find(|(candidate, _)| candidate.key == spec.key)
                    .map(|(_, value)| value.clone())
            })
            .collect::<Vec<_>>();
        self.scalar_values.clear();
        self.retain_row(&needed, &mut values);
    }

    fn prepare_columns(&mut self, available: &[ColumnSpec]) -> Option<Vec<ColumnSpec>> {
        let outputs = self
            .plan
            .columns
            .clone()
            .unwrap_or_else(|| available.to_vec());
        if self.plan.columns.is_none() && self.table.columns.is_empty() {
            self.table.columns = outputs
                .iter()
                .map(|spec| ColumnBuilder::new(spec.name.clone(), spec.column_type))
                .collect();
        } else if !same_column_keys(&outputs, &self.table.columns) {
            self.set_error(
                ProjectionErrorCode::IncompatibleColumns,
                "matching category occurrences have incompatible column layouts",
            );
            return None;
        }

        let mut needed = outputs;
        for predicate in &self.plan.predicates {
            let Ok(spec) = normalize_column(&self.plan.category, predicate.column()) else {
                debug_assert!(false, "predicate was validated when the plan was built");
                continue;
            };
            if !needed.iter().any(|candidate| candidate.key == spec.key) {
                needed.push(spec);
            }
        }
        if let Some(missing) = needed.iter().find(|required| {
            !available
                .iter()
                .any(|candidate| candidate.key == required.key)
        }) {
            self.set_error(
                ProjectionErrorCode::MissingColumn,
                format!(
                    "category {:?} does not contain required column {:?}",
                    self.plan.category, missing.name
                ),
            );
            return None;
        }
        Some(needed)
    }

    fn retain_row(&mut self, needed: &[ColumnSpec], values: &mut [Option<CifCell>]) {
        if self.error.is_some() || values.iter().any(Option::is_none) {
            return;
        }
        if !self.predicates_match(needed, values) {
            return;
        }
        let mut append_error = None;
        for (column_index, column) in self.table.columns.iter_mut().enumerate() {
            let Some(value) = values.get_mut(column_index).and_then(Option::take) else {
                debug_assert!(false, "every output value was populated");
                return;
            };
            if let Err(message) = column.append(value) {
                append_error = Some(message);
                break;
            }
        }
        if let Some(message) = append_error {
            self.set_error(ProjectionErrorCode::SchemaType, message);
            return;
        }
        self.append_provenance();
    }

    fn predicates_match(&self, needed: &[ColumnSpec], values: &[Option<CifCell>]) -> bool {
        self.plan.predicates.iter().all(|predicate| {
            let Ok(spec) = normalize_column(&self.plan.category, predicate.column()) else {
                return false;
            };
            let Some(index) = needed
                .iter()
                .position(|candidate| candidate.key == spec.key)
            else {
                return false;
            };
            let Some(cell) = values.get(index).and_then(Option::as_ref) else {
                return false;
            };
            predicate_matches(predicate, cell)
        })
    }

    fn set_error(&mut self, code: ProjectionErrorCode, message: impl Into<String>) {
        if self.error.is_none() {
            self.error = Some(ProjectionError::new(code, message));
        }
    }

    fn append_provenance(&mut self) {
        self.table.append_provenance(
            &self.source_name,
            self.block_code.as_deref(),
            self.frame_code.as_deref(),
        );
    }

    fn column_spec(&mut self, item: &str) -> Option<ColumnSpec> {
        let Some(column_type) = self.plan.column_type(item) else {
            self.set_error(
                ProjectionErrorCode::SchemaItem,
                format!(
                    "item _{}.{} is not defined by the selected schema",
                    self.plan.category, item
                ),
            );
            return None;
        };
        Some(ColumnSpec::from_item(item, column_type))
    }
}

impl ParseSink for TableSink {
    fn start_block(&mut self, kind: BlockKind, code: Option<&str>) {
        self.block_code = match kind {
            BlockKind::Data => code.map(str::to_owned),
            BlockKind::Global => None,
        };
    }

    fn end_block(&mut self) {
        self.flush_scalars();
        self.block_code = None;
    }

    fn start_frame(&mut self, code: &str) {
        self.flush_scalars();
        self.frame_code = Some(code.to_owned());
    }

    fn end_frame(&mut self) {
        self.flush_scalars();
        self.frame_code = None;
    }

    fn item(&mut self, tag: &str, value: ParsedValue<'_>) {
        let Some(spec) = split_tag(tag)
            .filter(|(category, _)| category.eq_ignore_ascii_case(&self.plan.category))
        else {
            return;
        };
        let Some(spec) = self.column_spec(spec.1) else {
            return;
        };
        self.scalar_values.push((spec, cell_from_value(&value)));
    }

    fn start_loop(&mut self, tags: &[String]) {
        debug_assert!(self.pending_parallel.is_none());
        let mut available = Vec::new();
        for tag in tags {
            let Some((_, item)) = split_tag(tag)
                .filter(|(category, _)| category.eq_ignore_ascii_case(&self.plan.category))
            else {
                continue;
            };
            let Some(spec) = self.column_spec(item) else {
                continue;
            };
            available.push(spec);
        }
        if available.is_empty() || self.error.is_some() {
            self.current_loop = None;
            return;
        }
        let Some(needed) = self.prepare_columns(&available) else {
            return;
        };
        let target_by_tag = tags
            .iter()
            .map(|tag| {
                split_tag(tag).and_then(|(category, item)| {
                    if !category.eq_ignore_ascii_case(&self.plan.category) {
                        return None;
                    }
                    let key = item.to_ascii_lowercase();
                    needed.iter().position(|spec| spec.key == key)
                })
            })
            .collect();
        let values = vec![None; needed.len()];
        self.current_loop = Some(LoopProjection {
            target_by_tag,
            needed,
            values,
            predicate_free: self.plan.predicates.is_empty(),
        });
    }

    fn wants_loop_value(&self, column_index: usize) -> bool {
        self.current_loop
            .as_ref()
            .is_some_and(|cif_loop| cif_loop.target_by_tag[column_index].is_some())
    }

    fn loop_value(&mut self, column_index: usize, value: ParsedValue<'_>) {
        if self
            .current_loop
            .as_ref()
            .is_some_and(|cif_loop| cif_loop.predicate_free)
        {
            let target_index = self
                .current_loop
                .as_ref()
                .and_then(|cif_loop| cif_loop.target_by_tag[column_index]);
            if let Some(target_index) = target_index {
                let result = self.table.columns[target_index].append_parsed(&value);
                if let Err(message) = result {
                    self.set_error(ProjectionErrorCode::SchemaType, message);
                    return;
                }
            }
            return;
        }
        let Some(cif_loop) = &mut self.current_loop else {
            return;
        };
        if column_index == 0 {
            cif_loop.values.fill(None);
        }
        if let Some(target_index) = cif_loop.target_by_tag[column_index] {
            cif_loop.values[target_index] = Some(cell_from_value(&value));
        }
    }

    fn loop_kernel(&self) -> Option<LoopKernel> {
        (self.current_loop.is_some() && self.scalar_values.is_empty() && self.error.is_none())
            .then_some(LoopKernel::Project)
    }

    fn try_parallel_loop(
        &mut self,
        source: &SourceBuffer,
        tags: &[String],
        value_start: usize,
        max_token_bytes: usize,
    ) -> Result<Option<ParallelLoop>, super::error::ParseError> {
        if self.current_loop.is_none() {
            return Ok(None);
        }
        let plan = &self.plan;
        let source_name = &self.source_name;
        let block_code = self.block_code.as_deref();
        let frame_code = self.frame_code.as_deref();
        let Some(projected) = scan_projected_loop(
            source,
            value_start,
            max_token_bytes,
            tags.len(),
            |source, range| {
                project_chunk_speculative(
                    source,
                    plan,
                    source_name,
                    block_code,
                    frame_code,
                    tags,
                    range,
                    max_token_bytes,
                )
            },
        )?
        else {
            return Ok(None);
        };

        let mut tables = Vec::with_capacity(projected.outputs.len());
        for output in projected.outputs {
            match output {
                Ok(table) => tables.push(table),
                Err(error) if self.error.is_none() => self.error = Some(error),
                Err(_) => {}
            }
        }
        self.pending_parallel = Some(tables);
        Ok(Some(projected.cif_loop))
    }

    fn parallel_loop(
        &mut self,
        source: &SourceBuffer,
        tags: &[String],
        cif_loop: &mut ParallelLoop,
    ) -> Result<bool, super::error::ParseError> {
        if self.current_loop.is_none() {
            return Ok(false);
        }
        if let Some(tables) = self.pending_parallel.take() {
            let prefix = self.take_table();
            if prefix.row_count() != 0 {
                self.segments.push(prefix);
            }
            self.segments
                .extend(tables.into_iter().filter(|table| table.row_count() != 0));
            return Ok(true);
        }
        let chunks = cif_loop.row_aligned_chunks(source, tags.len())?;
        if chunks.len() < 2 {
            return Ok(false);
        }

        let prefix = self.take_table();
        if prefix.row_count() != 0 {
            self.segments.push(prefix);
        }
        let plan = &self.plan;
        let source_name = &self.source_name;
        let block_code = self.block_code.as_deref();
        let frame_code = self.frame_code.as_deref();
        let max_token_bytes = cif_loop.max_token_bytes;
        let result = thread::scope(|scope| {
            let handles = chunks
                .iter()
                .copied()
                .map(|chunk| {
                    scope.spawn(move || {
                        project_chunk(
                            source,
                            plan,
                            source_name,
                            block_code,
                            frame_code,
                            tags,
                            chunk.start,
                            chunk.end,
                            max_token_bytes,
                        )
                    })
                })
                .collect::<Vec<_>>();
            handles
                .into_iter()
                .zip(chunks)
                .map(|(handle, chunk)| {
                    handle.join().map_err(|_| {
                        ProjectionError::from(source.error(
                            ParseErrorCode::UnexpectedControl,
                            "parallel projection worker terminated unexpectedly",
                            chunk.start,
                            chunk.start,
                        ))
                    })?
                })
                .collect::<Result<Vec<_>, ProjectionError>>()
        });
        match result {
            Ok(tables) => self
                .segments
                .extend(tables.into_iter().filter(|table| table.row_count() != 0)),
            Err(error) => self.error = Some(error),
        }
        Ok(true)
    }

    fn end_loop_row(&mut self) {
        if self
            .current_loop
            .as_ref()
            .is_some_and(|cif_loop| cif_loop.predicate_free)
        {
            self.append_provenance();
            return;
        }
        let Some(mut cif_loop) = self.current_loop.take() else {
            return;
        };
        self.retain_row(&cif_loop.needed, &mut cif_loop.values);
        cif_loop.values.fill(None);
        self.current_loop = Some(cif_loop);
    }

    fn end_loop(&mut self) {
        self.current_loop = None;
    }
}

#[allow(clippy::too_many_arguments)]
fn project_chunk_speculative(
    source: SourceBuffer,
    plan: &ProjectionPlan,
    source_name: &str,
    block_code: Option<&str>,
    frame_code: Option<&str>,
    tags: &[String],
    range: LoopRange,
    max_token_bytes: usize,
) -> ProjectedChunk<Result<CifTable, ProjectionError>> {
    let mut sink = TableSink::new(source_name, plan.clone());
    sink.block_code = block_code.map(str::to_owned);
    sink.frame_code = frame_code.map(str::to_owned);
    sink.start_loop(tags);
    let mut lexer = Lexer::bounded(source.clone(), range.start, range.end, max_token_bytes);
    let mut values = 0;
    let mut last_value = None;
    let mut control = None;
    let mut error = None;
    loop {
        match lexer.next_token() {
            Ok(Some(token)) if token.kind == TokenKind::Value => {
                let column_index = values % tags.len();
                if sink.wants_loop_value(column_index) {
                    sink.loop_value(column_index, ParsedValue::new(&source, &token));
                }
                values += 1;
                if values % tags.len() == 0 {
                    sink.end_loop_row();
                }
                last_value = Some(token);
            }
            Ok(Some(token)) => {
                control = Some(token);
                break;
            }
            Ok(None) => break,
            Err(parse_error) => {
                error = Some(parse_error);
                break;
            }
        }
    }
    sink.end_loop();
    ProjectedChunk {
        values,
        last_value,
        control,
        error,
        output: sink.finish(),
    }
}

#[allow(clippy::too_many_arguments)]
fn project_chunk(
    source: &SourceBuffer,
    plan: &ProjectionPlan,
    source_name: &str,
    block_code: Option<&str>,
    frame_code: Option<&str>,
    tags: &[String],
    start: usize,
    end: usize,
    max_token_bytes: usize,
) -> Result<CifTable, ProjectionError> {
    let mut sink = TableSink::new(source_name, plan.clone());
    sink.block_code = block_code.map(str::to_owned);
    sink.frame_code = frame_code.map(str::to_owned);
    sink.start_loop(tags);
    let mut lexer = Lexer::bounded(source.clone(), start, end, max_token_bytes);
    let mut column_index = 0;
    while let Some(token) = lexer.next_token()? {
        if token.kind != TokenKind::Value {
            return Err(source
                .error(
                    ParseErrorCode::UnexpectedControl,
                    "parallel projection chunk crossed the loop boundary",
                    token.span.start,
                    token.span.end,
                )
                .into());
        }
        if sink.wants_loop_value(column_index) {
            sink.loop_value(column_index, ParsedValue::new(source, &token));
        }
        column_index += 1;
        if column_index == tags.len() {
            sink.end_loop_row();
            column_index = 0;
        }
    }
    if column_index != 0 {
        return Err(source
            .error(
                ParseErrorCode::LoopValueCount,
                "parallel projection chunk ended inside a loop row",
                end,
                end,
            )
            .into());
    }
    sink.end_loop();
    sink.finish()
}

impl ColumnSpec {
    pub(super) fn from_item(item: &str, column_type: ColumnType) -> Self {
        Self {
            name: item.to_owned(),
            key: item.to_ascii_lowercase(),
            column_type,
        }
    }
}

fn cell_from_value(value: &ParsedValue<'_>) -> CifCell {
    match (value.is_unquoted(), value.as_str()) {
        (true, "?") => CifCell::Unknown,
        (true, ".") => CifCell::NotApplicable,
        _ => CifCell::Text(value.as_str().to_owned()),
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
    Ok(ColumnSpec::from_item(item, ColumnType::Text))
}

fn split_tag(tag: &str) -> Option<(&str, &str)> {
    tag.strip_prefix('_')?.split_once('.')
}

fn same_column_keys(specs: &[ColumnSpec], columns: &[ColumnBuilder]) -> bool {
    specs.len() == columns.len()
        && specs
            .iter()
            .zip(columns)
            .all(|(spec, column)| spec.key.eq_ignore_ascii_case(&column.name))
}

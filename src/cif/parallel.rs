use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use super::document::SourceCell;
use super::error::{ParseError, ParseErrorCode};
use super::lexer::Lexer;
use super::parallel_index::{RawChunk, safe_chunks};
use super::parallel_probe::loop_reaches;
use super::source::SourceBuffer;
use super::token::{Token, TokenKind};

#[cfg(not(feature = "robustness"))]
pub(crate) const PARALLEL_LOOP_MIN_BYTES: usize = 8 * 1024 * 1024;
#[cfg(feature = "robustness")]
pub(crate) const PARALLEL_LOOP_MIN_BYTES: usize = 1024;
#[cfg(not(feature = "robustness"))]
const BYTES_PER_WORKER: usize = 2 * 1024 * 1024;
#[cfg(feature = "robustness")]
const BYTES_PER_WORKER: usize = 256;
const CELLS_PER_SEGMENT: usize = 64 * 1024;

static ACTIVE_WORKERS: AtomicUsize = AtomicUsize::new(0);

#[derive(Clone, Copy, Eq, PartialEq)]
pub(crate) enum LoopKernel {
    Project,
    Retain,
}

#[derive(Clone, Copy)]
pub(crate) struct LoopChunk {
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) value_prefix: usize,
}

pub(crate) struct ProjectedChunk<T> {
    pub(crate) values: usize,
    pub(crate) last_value: Option<Token>,
    pub(crate) control: Option<Token>,
    pub(crate) error: Option<ParseError>,
    pub(crate) output: T,
}

pub(crate) struct ProjectedLoop<T> {
    pub(crate) cif_loop: ParallelLoop,
    pub(crate) outputs: Vec<T>,
}

pub(crate) struct ParallelLoop {
    pub(crate) chunks: Vec<LoopChunk>,
    pub(crate) loop_end: usize,
    pub(crate) value_count: usize,
    pub(crate) last_value: Option<Token>,
    pub(crate) control: Option<Token>,
    pub(crate) max_token_bytes: usize,
    retained: Vec<Vec<SourceCell>>,
    _lease: WorkerLease,
}

impl ParallelLoop {
    pub(crate) fn take_source_cells(&mut self) -> Vec<Vec<SourceCell>> {
        std::mem::take(&mut self.retained)
    }

    pub(crate) fn row_aligned_chunks(
        &self,
        source: &SourceBuffer,
        column_count: usize,
    ) -> Result<Vec<LoopChunk>, ParseError> {
        let mut starts = Vec::with_capacity(self.chunks.len());
        for (chunk_index, chunk) in self.chunks.iter().enumerate() {
            if chunk_index == 0 || chunk.value_prefix % column_count == 0 {
                starts.push(chunk.start);
                continue;
            }
            let skip = column_count - chunk.value_prefix % column_count;
            let mut lexer =
                Lexer::bounded(source.clone(), chunk.start, chunk.end, self.max_token_bytes);
            let mut aligned_start = chunk.start;
            for _ in 0..skip {
                let Some(token) = lexer.next_token()? else {
                    aligned_start = chunk.end;
                    break;
                };
                if token.kind != TokenKind::Value {
                    aligned_start = token.span.start;
                    break;
                }
                aligned_start = token.span.end;
            }
            if aligned_start < chunk.end {
                let mut boundary = Lexer::bounded(
                    source.clone(),
                    aligned_start,
                    chunk.end,
                    self.max_token_bytes,
                );
                aligned_start = boundary
                    .next_token()?
                    .map_or(chunk.end, |token| token.span.start);
            }
            starts.push(aligned_start);
        }
        let mut chunks = Vec::with_capacity(starts.len());
        for index in 0..starts.len() {
            let start = starts[index];
            let end = starts.get(index + 1).copied().unwrap_or(self.loop_end);
            if start < end {
                chunks.push(LoopChunk {
                    start,
                    end,
                    value_prefix: 0,
                });
            }
        }
        Ok(chunks)
    }
}

struct WorkerLease {
    count: usize,
}

impl Drop for WorkerLease {
    fn drop(&mut self) {
        ACTIVE_WORKERS.fetch_sub(self.count, Ordering::AcqRel);
    }
}

pub(crate) fn scan_loop(
    source: &SourceBuffer,
    value_start: usize,
    max_token_bytes: usize,
    kernel: LoopKernel,
) -> Result<Option<ParallelLoop>, ParseError> {
    let Some((raw_chunks, lease)) = prepare_chunks(source, value_start, max_token_bytes, kernel)?
    else {
        return Ok(None);
    };
    let scans = run_workers(
        source,
        &raw_chunks,
        "parallel lexer worker terminated unexpectedly",
        |source, chunk| scan_chunk(source, chunk, max_token_bytes, kernel),
    )?;
    let (loop_end, control) = loop_boundary(source, &scans)?;
    let mut chunks = Vec::new();
    let mut value_count = 0;
    let mut last_value = None;
    let mut retained: Vec<Vec<SourceCell>> = Vec::new();
    for (mut scan, raw_chunk) in scans.into_iter().zip(raw_chunks) {
        if raw_chunk.start >= loop_end {
            break;
        }
        chunks.push(LoopChunk {
            start: raw_chunk.start,
            end: raw_chunk.end.min(loop_end),
            value_prefix: value_count,
        });
        value_count += scan.values;
        if scan.last_value.is_some() {
            last_value = scan.last_value.take();
        }
        retained.append(&mut scan.output);
        if scan.control.is_some() {
            break;
        }
    }
    Ok(Some(ParallelLoop {
        chunks,
        loop_end,
        value_count,
        last_value,
        control,
        max_token_bytes,
        retained,
        _lease: lease,
    }))
}

pub(crate) fn scan_projected_loop<T, F>(
    source: &SourceBuffer,
    value_start: usize,
    max_token_bytes: usize,
    column_count: usize,
    project: F,
) -> Result<Option<ProjectedLoop<T>>, ParseError>
where
    T: Send,
    F: Fn(SourceBuffer, Range<usize>) -> ProjectedChunk<T> + Sync,
{
    let Some((raw_chunks, lease)) =
        prepare_chunks(source, value_start, max_token_bytes, LoopKernel::Project)?
    else {
        return Ok(None);
    };
    let scans = run_workers(
        source,
        &raw_chunks,
        "parallel projection worker terminated unexpectedly",
        |source, chunk| project(source, chunk.start..chunk.end),
    )?;
    let (loop_end, control) = loop_boundary(source, &scans)?;
    let mut chunks = Vec::new();
    let mut outputs = Vec::new();
    let mut value_count = 0;
    let mut last_value = None;
    for (scan, raw_chunk) in scans.into_iter().zip(raw_chunks) {
        if raw_chunk.start >= loop_end {
            break;
        }
        if value_count % column_count != 0 {
            return Ok(None);
        }
        chunks.push(LoopChunk {
            start: raw_chunk.start,
            end: raw_chunk.end.min(loop_end),
            value_prefix: value_count,
        });
        value_count += scan.values;
        if scan.last_value.is_some() {
            last_value = scan.last_value;
        }
        outputs.push(scan.output);
        if scan.control.is_some() {
            break;
        }
    }
    Ok(Some(ProjectedLoop {
        cif_loop: ParallelLoop {
            chunks,
            loop_end,
            value_count,
            last_value,
            control,
            max_token_bytes,
            retained: Vec::new(),
            _lease: lease,
        },
        outputs,
    }))
}

fn prepare_chunks(
    source: &SourceBuffer,
    value_start: usize,
    max_token_bytes: usize,
    kernel: LoopKernel,
) -> Result<Option<(Vec<RawChunk>, WorkerLease)>, ParseError> {
    let remaining = source.len().saturating_sub(value_start);
    if remaining < PARALLEL_LOOP_MIN_BYTES {
        return Ok(None);
    }
    if kernel == LoopKernel::Retain
        && !loop_reaches(
            source,
            value_start,
            value_start + PARALLEL_LOOP_MIN_BYTES,
            max_token_bytes,
        )
    {
        return Ok(None);
    }
    let requested = remaining.div_ceil(BYTES_PER_WORKER);
    let Some(lease) = acquire_workers(requested) else {
        return Ok(None);
    };
    let chunks = safe_chunks(source, value_start, lease.count)?;
    if chunks.len() < 2 {
        return Ok(None);
    }
    Ok(Some((chunks, lease)))
}

fn loop_boundary<T>(
    source: &SourceBuffer,
    scans: &[ProjectedChunk<T>],
) -> Result<(usize, Option<Token>), ParseError> {
    let control = scans
        .iter()
        .filter_map(|scan| scan.control.as_ref())
        .min_by_key(|token| token.span.start)
        .cloned();
    let loop_end = control
        .as_ref()
        .map_or(source.len(), |token| token.span.start);
    if let Some(error) = scans
        .iter()
        .filter_map(|scan| scan.error.as_ref())
        .filter(|error| error.span().byte_start() < loop_end)
        .min_by_key(|error| error.span().byte_start())
    {
        return Err(error.clone());
    }
    Ok((loop_end, control))
}

fn acquire_workers(requested: usize) -> Option<WorkerLease> {
    let available = thread::available_parallelism().map_or(1, usize::from);
    if available < 2 {
        return None;
    }
    let requested = requested.clamp(2, available);
    let mut active = ACTIVE_WORKERS.load(Ordering::Acquire);
    loop {
        let granted = requested.min(available.saturating_sub(active));
        if granted < 2 {
            return None;
        }
        match ACTIVE_WORKERS.compare_exchange_weak(
            active,
            active + granted,
            Ordering::AcqRel,
            Ordering::Acquire,
        ) {
            Ok(_) => return Some(WorkerLease { count: granted }),
            Err(updated) => active = updated,
        }
    }
}

fn run_workers<T, F>(
    source: &SourceBuffer,
    chunks: &[RawChunk],
    panic_message: &'static str,
    work: F,
) -> Result<Vec<T>, ParseError>
where
    T: Send,
    F: Fn(SourceBuffer, RawChunk) -> T + Sync,
{
    thread::scope(|scope| {
        let handles = chunks
            .iter()
            .copied()
            .map(|chunk| {
                let source = source.clone();
                let work = &work;
                scope.spawn(move || work(source, chunk))
            })
            .collect::<Vec<_>>();
        handles
            .into_iter()
            .zip(chunks)
            .map(|(handle, chunk)| {
                handle.join().map_err(|_| {
                    source.error(
                        ParseErrorCode::UnexpectedControl,
                        panic_message,
                        chunk.start,
                        chunk.start,
                    )
                })
            })
            .collect()
    })
}

fn scan_chunk(
    source: SourceBuffer,
    chunk: RawChunk,
    max_token_bytes: usize,
    kernel: LoopKernel,
) -> ProjectedChunk<Vec<Vec<SourceCell>>> {
    let mut lexer = Lexer::bounded(source.clone(), chunk.start, chunk.end, max_token_bytes);
    let mut values = 0;
    let mut last_value = None;
    let mut retained: Vec<Vec<SourceCell>> = Vec::new();
    loop {
        match lexer.next_token() {
            Ok(Some(token)) if token.kind == TokenKind::Value => {
                if kernel == LoopKernel::Retain {
                    push_retained_cell(&mut retained, SourceCell::from_token(&token));
                }
                values += 1;
                last_value = Some(token);
            }
            Ok(Some(control)) => {
                return ProjectedChunk {
                    values,
                    last_value,
                    control: Some(control),
                    error: None,
                    output: retained,
                };
            }
            Ok(None) => {
                return ProjectedChunk {
                    values,
                    last_value,
                    control: None,
                    error: None,
                    output: retained,
                };
            }
            Err(error) => {
                return ProjectedChunk {
                    values,
                    last_value,
                    control: None,
                    error: Some(error),
                    output: retained,
                };
            }
        }
    }
}

fn push_retained_cell(retained: &mut Vec<Vec<SourceCell>>, cell: SourceCell) {
    match retained.last_mut() {
        Some(segment) if segment.len() < CELLS_PER_SEGMENT => segment.push(cell),
        _ => {
            let mut segment = Vec::with_capacity(CELLS_PER_SEGMENT);
            segment.push(cell);
            retained.push(segment);
        }
    }
}

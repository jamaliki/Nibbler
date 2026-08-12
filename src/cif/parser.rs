use std::collections::HashSet;

use super::document::{BlockKind, CifDocument};
use super::error::{ParseError, ParseErrorCode, SourceSpan};
use super::lexer::Lexer;
use super::parallel::scan_loop;
use super::sink::{DocumentSink, ParseSink, ParsedValue};
use super::source::SourceBuffer;
use super::token::{Token, TokenKind};

/// Hard limits applied by strict CIF parsing.
///
/// Counts are document-wide. A limit is checked before the corresponding value is
/// retained, preventing allocations based on untrusted declared sizes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Limits {
    /// Maximum UTF-8 source size.
    pub source_bytes: usize,
    /// Maximum bytes in one lexical token.
    pub token_bytes: usize,
    /// Maximum data and global blocks.
    pub blocks: usize,
    /// Maximum save frames.
    pub frames: usize,
    /// Maximum distinct loops.
    pub loops: usize,
    /// Maximum scalar and loop tags combined.
    pub tags: usize,
    /// Maximum columns in any one loop.
    pub loop_columns: usize,
    /// Maximum loop rows combined across the document.
    pub rows: usize,
    /// Maximum scalar and loop values combined.
    pub values: usize,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            source_bytes: 2 * 1024 * 1024 * 1024,
            token_bytes: 64 * 1024 * 1024,
            blocks: 100_000,
            frames: 1_000_000,
            loops: 10_000_000,
            tags: 10_000_000,
            loop_columns: 100_000,
            rows: 1_000_000_000,
            values: 1_000_000_000,
        }
    }
}

/// Options for strict CIF parsing.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ParseOptions {
    /// Resource limits enforced before document allocation.
    pub limits: Limits,
}

/// Parse UTF-8 CIF bytes using strict defaults.
///
/// Schema-less ordinary values remain text. Unknown (`?`) and not-applicable (`.`)
/// values remain distinct variants.
///
/// # Errors
///
/// Returns a located [`ParseError`] for malformed UTF-8, lexical errors, grammar
/// errors, duplicate tags, incomplete loops, or resource-limit failures.
pub fn parse(bytes: &[u8]) -> Result<CifDocument, ParseError> {
    parse_with_options(bytes, ParseOptions::default())
}

/// Parse UTF-8 CIF bytes with explicit strict options.
///
/// # Errors
///
/// Returns the same failures as [`parse`].
pub fn parse_with_options(bytes: &[u8], options: ParseOptions) -> Result<CifDocument, ParseError> {
    if bytes.len() > options.limits.source_bytes {
        return Err(ParseError::new(
            ParseErrorCode::ResourceLimit,
            format!(
                "source exceeds the configured limit of {} bytes",
                options.limits.source_bytes
            ),
            "<memory>",
            SourceSpan::new(0, 0, 1, 1),
        ));
    }
    let source = SourceBuffer::from_bytes("<memory>", bytes)?;
    parse_source_with_options(source, options)
}

/// Parse an immutable named source using strict defaults.
///
/// # Errors
///
/// Returns a located [`ParseError`] for grammar or resource-limit failures.
pub fn parse_source(source: SourceBuffer) -> Result<CifDocument, ParseError> {
    parse_source_with_options(source, ParseOptions::default())
}

/// Parse an immutable named source with explicit strict options.
///
/// # Errors
///
/// Returns a located [`ParseError`] for grammar or resource-limit failures.
pub fn parse_source_with_options(
    source: SourceBuffer,
    options: ParseOptions,
) -> Result<CifDocument, ParseError> {
    if source.len() > u32::MAX as usize {
        return Err(source.error(
            ParseErrorCode::ResourceLimit,
            "full CIF documents are limited to 4 GiB source buffers",
            0,
            0,
        ));
    }
    let mut sink = DocumentSink::new(source.clone());
    parse_source_into(source, options, &mut sink)?;
    Ok(sink.finish())
}

pub(crate) fn parse_source_into(
    source: SourceBuffer,
    options: ParseOptions,
    sink: &mut impl ParseSink,
) -> Result<(), ParseError> {
    if source.len() > options.limits.source_bytes {
        return Err(source.error(
            ParseErrorCode::ResourceLimit,
            format!(
                "source exceeds the configured limit of {} bytes",
                options.limits.source_bytes
            ),
            0,
            0,
        ));
    }
    Parser::new(source, options.limits).parse_document(sink)
}

#[derive(Default)]
struct Counts {
    blocks: usize,
    frames: usize,
    loops: usize,
    tags: usize,
    rows: usize,
    values: usize,
}

struct Parser {
    source: SourceBuffer,
    lexer: Lexer,
    lookahead: Option<Token>,
    reached_end: bool,
    limits: Limits,
    counts: Counts,
    #[cfg(feature = "robustness")]
    allow_parallel: bool,
}

impl Parser {
    fn new(source: SourceBuffer, limits: Limits) -> Self {
        Self {
            lexer: Lexer::new(source.clone(), limits.token_bytes),
            source,
            lookahead: None,
            reached_end: false,
            limits,
            counts: Counts::default(),
            #[cfg(feature = "robustness")]
            allow_parallel: true,
        }
    }

    #[cfg(feature = "robustness")]
    const fn serial(mut self) -> Self {
        self.allow_parallel = false;
        self
    }

    #[inline(always)]
    const fn parallel_enabled(&self) -> bool {
        #[cfg(feature = "robustness")]
        {
            self.allow_parallel
        }
        #[cfg(not(feature = "robustness"))]
        {
            true
        }
    }

    pub(crate) fn parse_document(mut self, sink: &mut impl ParseSink) -> Result<(), ParseError> {
        let mut block_codes = HashSet::new();
        let mut has_global_block = false;

        while let Some(header) = self.peek()?.cloned() {
            match header.kind {
                TokenKind::Data => {
                    self.charge(CountKind::Block, &header)?;
                    let header = self.next_required()?;
                    let code = self.control_suffix(&header, "data_").to_owned();
                    if code.is_empty() {
                        return Err(self.error_at(
                            ParseErrorCode::EmptyBlockCode,
                            "strict CIF requires a non-empty data block code",
                            &header,
                        ));
                    }
                    let folded_code = code.to_ascii_lowercase();
                    if !block_codes.insert(folded_code) {
                        return Err(self.error_at(
                            ParseErrorCode::DuplicateBlock,
                            format!("duplicate data block code {code:?}"),
                            &header,
                        ));
                    }
                    sink.start_block(BlockKind::Data, Some(&code));
                    self.parse_block_entries(sink)?;
                    sink.end_block();
                }
                TokenKind::Global => {
                    self.charge(CountKind::Block, &header)?;
                    let header = self.next_required()?;
                    if has_global_block {
                        return Err(self.error_at(
                            ParseErrorCode::DuplicateBlock,
                            "a document may contain at most one global block",
                            &header,
                        ));
                    }
                    has_global_block = true;
                    sink.start_block(BlockKind::Global, None);
                    self.parse_block_entries(sink)?;
                    sink.end_block();
                }
                _ => {
                    return Err(self.error_at(
                        ParseErrorCode::ExpectedBlock,
                        "expected data_ or global_ before document content",
                        &header,
                    ));
                }
            }
        }

        if self.counts.blocks == 0 {
            return Err(self.error_at_end(
                ParseErrorCode::ExpectedBlock,
                "CIF document contains no data_ or global_ block",
            ));
        }
        Ok(())
    }

    fn parse_block_entries(&mut self, sink: &mut impl ParseSink) -> Result<(), ParseError> {
        let mut tags = HashSet::new();

        while let Some(token) = self.peek()?.cloned() {
            match token.kind {
                TokenKind::Data | TokenKind::Global => break,
                TokenKind::Tag => {
                    self.parse_item(&mut tags, sink)?;
                }
                TokenKind::Loop => {
                    self.parse_loop(&mut tags, sink)?;
                }
                TokenKind::SaveStart => {
                    self.parse_frame(sink)?;
                }
                TokenKind::SaveEnd => {
                    return Err(self.error_at(
                        ParseErrorCode::UnexpectedSaveEnd,
                        "save_ closes no open save frame",
                        &token,
                    ));
                }
                TokenKind::Stop => {
                    return Err(self.error_at(
                        ParseErrorCode::UnexpectedStop,
                        "stop_ is only valid as an explicit loop terminator",
                        &token,
                    ));
                }
                TokenKind::Value => {
                    return Err(self.error_at(
                        ParseErrorCode::UnexpectedValue,
                        "bare value is not associated with a tag or loop",
                        &token,
                    ));
                }
            }
        }
        Ok(())
    }

    fn parse_frame(&mut self, sink: &mut impl ParseSink) -> Result<(), ParseError> {
        let opening = self.next_required()?;
        self.charge(CountKind::Frame, &opening)?;
        let code = self.control_suffix(&opening, "save_").to_owned();
        let mut tags = HashSet::new();
        sink.start_frame(&code);

        loop {
            let Some(token) = self.peek()?.cloned() else {
                return Err(self.error_at_end(
                    ParseErrorCode::UnterminatedSaveFrame,
                    format!("save frame {code:?} reaches end of input before save_"),
                ));
            };
            match token.kind {
                TokenKind::SaveEnd => {
                    self.next_required()?;
                    sink.end_frame();
                    return Ok(());
                }
                TokenKind::Data | TokenKind::Global => {
                    return Err(self.error_at(
                        ParseErrorCode::UnterminatedSaveFrame,
                        format!("save frame {code:?} reaches a block boundary before save_"),
                        &token,
                    ));
                }
                TokenKind::SaveStart => {
                    return Err(self.error_at(
                        ParseErrorCode::NestedSaveFrame,
                        "save frames cannot be nested",
                        &token,
                    ));
                }
                TokenKind::Tag => {
                    self.parse_item(&mut tags, sink)?;
                }
                TokenKind::Loop => {
                    self.parse_loop(&mut tags, sink)?;
                }
                TokenKind::Stop => {
                    return Err(self.error_at(
                        ParseErrorCode::UnexpectedStop,
                        "stop_ is only valid as an explicit loop terminator",
                        &token,
                    ));
                }
                TokenKind::Value => {
                    return Err(self.error_at(
                        ParseErrorCode::UnexpectedValue,
                        "bare value is not associated with a tag or loop",
                        &token,
                    ));
                }
            }
        }
    }

    fn parse_item(
        &mut self,
        seen_tags: &mut HashSet<String>,
        sink: &mut impl ParseSink,
    ) -> Result<(), ParseError> {
        let tag_token = self.next_required()?;
        let tag = self.register_tag(&tag_token, seen_tags)?;
        let Some(value_token) = self.peek()?.cloned() else {
            return Err(self.error_at(
                ParseErrorCode::MissingItemValue,
                format!("tag {tag:?} reaches end of input without a value"),
                &tag_token,
            ));
        };
        if value_token.kind != TokenKind::Value {
            return Err(self.error_at(
                ParseErrorCode::MissingItemValue,
                format!("tag {tag:?} is not followed by a value"),
                &tag_token,
            ));
        }
        let value_token = self.next_required()?;
        self.charge(CountKind::Value, &value_token)?;
        sink.item(&tag, ParsedValue::new(&self.source, &value_token));
        Ok(())
    }

    fn parse_loop(
        &mut self,
        seen_tags: &mut HashSet<String>,
        sink: &mut impl ParseSink,
    ) -> Result<(), ParseError> {
        let loop_token = self.next_required()?;
        self.charge(CountKind::Loop, &loop_token)?;
        let mut tags = Vec::new();
        while self
            .peek()?
            .is_some_and(|token| token.kind == TokenKind::Tag)
        {
            let tag_token = self.next_required()?;
            if tags.len() >= self.limits.loop_columns {
                return Err(self.error_at(
                    ParseErrorCode::ResourceLimit,
                    format!(
                        "loop exceeds the configured limit of {} columns",
                        self.limits.loop_columns
                    ),
                    &tag_token,
                ));
            }
            tags.push(self.register_tag(&tag_token, seen_tags)?);
        }
        if tags.is_empty() {
            return Err(self.error_at(
                ParseErrorCode::MissingLoopTag,
                "loop_ must be followed by at least one tag",
                &loop_token,
            ));
        }

        sink.start_loop(&tags);
        let value_start = self
            .peek()?
            .filter(|token| token.kind == TokenKind::Value)
            .map(|token| token.span.start);
        if let Some(value_start) = value_start
            && self.parallel_enabled()
            && let Some(kernel) = sink.loop_kernel()
            && self.parallel_limits_are_unreachable(value_start)
        {
            let direct =
                sink.try_parallel_loop(&self.source, &tags, value_start, self.limits.token_bytes)?;
            let parallel = match direct {
                Some(parallel) => Some(parallel),
                None => scan_loop(&self.source, value_start, self.limits.token_bytes, kernel)?,
            };
            let Some(mut parallel) = parallel else {
                return self.parse_loop_serial(loop_token, tags, sink);
            };
            if parallel.value_count == 0 {
                return Err(self.error_at(
                    ParseErrorCode::MissingLoopValue,
                    "strict CIF rejects a loop with no values",
                    &loop_token,
                ));
            }
            if parallel.value_count % tags.len() != 0 {
                let location = parallel.last_value.as_ref().unwrap_or(&loop_token);
                return Err(self.error_at(
                    ParseErrorCode::LoopValueCount,
                    format!(
                        "loop has {} tags but {} values; the final row is incomplete",
                        tags.len(),
                        parallel.value_count
                    ),
                    location,
                ));
            }
            if sink.parallel_loop(&self.source, &tags, &mut parallel)? {
                self.counts.values += parallel.value_count;
                self.counts.rows += parallel.value_count / tags.len();
                if let Some(control) = parallel.control {
                    self.lexer.set_offset(control.span.end);
                    if control.kind == TokenKind::Stop {
                        self.lookahead = None;
                    } else {
                        self.lookahead = Some(control);
                    }
                    self.reached_end = false;
                } else {
                    self.lexer.set_offset(parallel.loop_end);
                    self.lookahead = None;
                    self.reached_end = true;
                }
                sink.end_loop();
                return Ok(());
            }
        }
        self.parse_loop_serial(loop_token, tags, sink)
    }

    fn parse_loop_serial(
        &mut self,
        loop_token: Token,
        tags: Vec<String>,
        sink: &mut impl ParseSink,
    ) -> Result<(), ParseError> {
        let mut value_count = 0_usize;
        let mut column_index = 0_usize;
        let mut last_value = None;
        while self
            .peek()?
            .is_some_and(|token| token.kind == TokenKind::Value)
        {
            let value_token = self.next_required()?;
            if column_index == 0 {
                self.charge(CountKind::Row, &value_token)?;
            }
            self.charge(CountKind::Value, &value_token)?;
            if sink.wants_loop_value(column_index) {
                sink.loop_value(column_index, ParsedValue::new(&self.source, &value_token));
            }
            value_count += 1;
            column_index += 1;
            if column_index == tags.len() {
                column_index = 0;
                sink.end_loop_row();
            }
            last_value = Some(value_token);
        }
        if self
            .peek()?
            .is_some_and(|token| token.kind == TokenKind::Stop)
        {
            self.next_required()?;
        }
        if value_count == 0 {
            return Err(self.error_at(
                ParseErrorCode::MissingLoopValue,
                "strict CIF rejects a loop with no values",
                &loop_token,
            ));
        }
        if column_index != 0 {
            let location = last_value.as_ref().unwrap_or(&loop_token);
            return Err(self.error_at(
                ParseErrorCode::LoopValueCount,
                format!(
                    "loop has {} tags but {} values; the final row is incomplete",
                    tags.len(),
                    value_count
                ),
                location,
            ));
        }
        sink.end_loop();
        Ok(())
    }

    fn register_tag(
        &mut self,
        token: &Token,
        seen_tags: &mut HashSet<String>,
    ) -> Result<String, ParseError> {
        if self.token_text(token).len() == 1 {
            return Err(self.error_at(
                ParseErrorCode::InvalidTag,
                "a CIF tag must contain characters after its leading underscore",
                token,
            ));
        }
        self.charge(CountKind::Tag, token)?;
        let tag = self.token_text(token).to_owned();
        let folded_tag = tag.to_ascii_lowercase();
        if !seen_tags.insert(folded_tag) {
            return Err(self.error_at(
                ParseErrorCode::DuplicateTag,
                format!("duplicate tag {tag:?} in one block or save frame"),
                token,
            ));
        }
        Ok(tag)
    }

    fn parallel_limits_are_unreachable(&self, value_start: usize) -> bool {
        let remaining_bytes = self.source.len().saturating_sub(value_start);
        self.limits.values.saturating_sub(self.counts.values) >= remaining_bytes
            && self.limits.rows.saturating_sub(self.counts.rows) >= remaining_bytes
    }

    fn charge(&mut self, kind: CountKind, token: &Token) -> Result<(), ParseError> {
        let (count, limit, label) = match kind {
            CountKind::Block => (&mut self.counts.blocks, self.limits.blocks, "blocks"),
            CountKind::Frame => (&mut self.counts.frames, self.limits.frames, "save frames"),
            CountKind::Loop => (&mut self.counts.loops, self.limits.loops, "loops"),
            CountKind::Tag => (&mut self.counts.tags, self.limits.tags, "tags"),
            CountKind::Row => (&mut self.counts.rows, self.limits.rows, "loop rows"),
            CountKind::Value => (&mut self.counts.values, self.limits.values, "values"),
        };
        if *count >= limit {
            return Err(self.source.error(
                ParseErrorCode::ResourceLimit,
                format!("document exceeds the configured limit of {limit} {label}"),
                token.span.start,
                token.span.end,
            ));
        }
        *count += 1;
        Ok(())
    }

    fn peek(&mut self) -> Result<Option<&Token>, ParseError> {
        if self.lookahead.is_none() && !self.reached_end {
            self.lookahead = self.lexer.next_token()?;
            self.reached_end = self.lookahead.is_none();
        }
        Ok(self.lookahead.as_ref())
    }

    fn next_required(&mut self) -> Result<Token, ParseError> {
        if let Some(token) = self.lookahead.take() {
            return Ok(token);
        }
        if !self.reached_end {
            if let Some(token) = self.lexer.next_token()? {
                return Ok(token);
            }
            self.reached_end = true;
        }
        Err(self.error_at_end(
            ParseErrorCode::UnexpectedControl,
            "internal parser request reached end of input",
        ))
    }

    fn token_text(&self, token: &Token) -> &str {
        self.source.slice(token.content.clone())
    }

    fn control_suffix<'a>(&'a self, token: &Token, prefix: &str) -> &'a str {
        let text = self.token_text(token);
        debug_assert!(text.len() >= prefix.len());
        &text[prefix.len()..]
    }

    fn error_at(
        &self,
        code: ParseErrorCode,
        message: impl Into<String>,
        token: &Token,
    ) -> ParseError {
        self.source
            .error(code, message, token.span.start, token.span.end)
    }

    fn error_at_end(&self, code: ParseErrorCode, message: impl Into<String>) -> ParseError {
        self.source
            .error(code, message, self.source.len(), self.source.len())
    }
}

#[cfg(feature = "robustness")]
pub(crate) fn parse_source_into_serial_reference(
    source: SourceBuffer,
    options: ParseOptions,
    sink: &mut impl ParseSink,
) -> Result<(), ParseError> {
    if source.len() > options.limits.source_bytes {
        return Err(source.error(
            ParseErrorCode::ResourceLimit,
            format!(
                "source exceeds the configured limit of {} bytes",
                options.limits.source_bytes
            ),
            0,
            0,
        ));
    }
    Parser::new(source, options.limits)
        .serial()
        .parse_document(sink)
}

#[cfg(feature = "robustness")]
/// Parse a document with the parallel loop kernel disabled.
///
/// This entry point exists only for robustness tests against [`parse_with_options`].
#[doc(hidden)]
pub fn parse_serial_reference(
    bytes: &[u8],
    options: ParseOptions,
) -> Result<CifDocument, ParseError> {
    if bytes.len() > options.limits.source_bytes {
        return Err(ParseError::new(
            ParseErrorCode::ResourceLimit,
            format!(
                "source exceeds the configured limit of {} bytes",
                options.limits.source_bytes
            ),
            "<memory>",
            SourceSpan::new(0, 0, 1, 1),
        ));
    }
    let source = SourceBuffer::from_bytes("<memory>", bytes)?;
    let mut sink = DocumentSink::new(source.clone());
    parse_source_into_serial_reference(source, options, &mut sink)?;
    Ok(sink.finish())
}

#[derive(Clone, Copy)]
enum CountKind {
    Block,
    Frame,
    Loop,
    Tag,
    Row,
    Value,
}

use super::document::{
    BlockKind, CifBlock, CifDocument, CifEntry, CifFrame, CifItem, CifLoop, CifValue, SourceCell,
    TextValue,
};
use super::error::ParseError;
use super::parallel::{LoopKernel, ParallelLoop};
use super::source::SourceBuffer;
use super::token::Token;

/// A decoded parser value that still refers to the immutable source buffer.
pub(crate) struct ParsedValue<'a> {
    source: &'a SourceBuffer,
    token: &'a Token,
}

impl<'a> ParsedValue<'a> {
    pub(crate) const fn new(source: &'a SourceBuffer, token: &'a Token) -> Self {
        Self { source, token }
    }

    pub(crate) fn as_str(&self) -> &str {
        self.source.slice(self.token.content.clone())
    }

    pub(crate) fn is_unquoted(&self) -> bool {
        self.token.quote_style == super::document::QuoteStyle::Unquoted
    }

    pub(crate) fn to_cif_value(&self) -> CifValue {
        let text = self.as_str();
        if self.token.quote_style == super::document::QuoteStyle::Unquoted {
            if text == "?" {
                return CifValue::Unknown;
            }
            if text == "." {
                return CifValue::NotApplicable;
            }
        }
        CifValue::Text(TextValue::from_source(
            self.source.clone(),
            self.token.content.clone(),
            self.token.quote_style,
        ))
    }

    pub(crate) fn to_source_cell(&self) -> SourceCell {
        SourceCell::from_token(self.token)
    }
}

/// Events emitted by the one CIF grammar implementation.
pub(crate) trait ParseSink {
    fn start_block(&mut self, kind: BlockKind, code: Option<&str>);
    fn end_block(&mut self);
    fn start_frame(&mut self, code: &str);
    fn end_frame(&mut self);
    fn item(&mut self, tag: &str, value: ParsedValue<'_>);
    fn start_loop(&mut self, tags: &[String]);
    fn wants_loop_value(&self, _column_index: usize) -> bool {
        true
    }
    fn loop_value(&mut self, column_index: usize, value: ParsedValue<'_>);
    fn loop_kernel(&self) -> Option<LoopKernel> {
        None
    }
    fn try_parallel_loop(
        &mut self,
        _source: &SourceBuffer,
        _tags: &[String],
        _value_start: usize,
        _max_token_bytes: usize,
    ) -> Result<Option<ParallelLoop>, ParseError> {
        Ok(None)
    }
    fn parallel_loop(
        &mut self,
        _source: &SourceBuffer,
        _tags: &[String],
        _cif_loop: &mut ParallelLoop,
    ) -> Result<bool, ParseError> {
        Ok(false)
    }
    fn end_loop_row(&mut self) {}
    fn end_loop(&mut self);
}

pub(crate) struct DocumentSink {
    source: SourceBuffer,
    blocks: Vec<CifBlock>,
    block: Option<BlockBuilder>,
    frame: Option<FrameBuilder>,
    cif_loop: Option<LoopBuilder>,
}

struct BlockBuilder {
    kind: BlockKind,
    code: Option<String>,
    entries: Vec<CifEntry>,
}

struct FrameBuilder {
    code: String,
    entries: Vec<CifEntry>,
}

struct LoopBuilder {
    tags: Vec<String>,
    cell_chunks: Vec<Vec<SourceCell>>,
}

impl DocumentSink {
    pub(crate) fn new(source: SourceBuffer) -> Self {
        Self {
            source,
            blocks: Vec::new(),
            block: None,
            frame: None,
            cif_loop: None,
        }
    }

    pub(crate) fn finish(self) -> CifDocument {
        debug_assert!(self.block.is_none());
        debug_assert!(self.frame.is_none());
        debug_assert!(self.cif_loop.is_none());
        CifDocument::new(self.blocks)
    }

    fn push_entry(&mut self, entry: CifEntry) {
        if let Some(frame) = &mut self.frame {
            frame.entries.push(entry);
        } else if let Some(block) = &mut self.block {
            block.entries.push(entry);
        } else {
            debug_assert!(false, "parser emitted an entry outside a block");
        }
    }
}

impl ParseSink for DocumentSink {
    fn start_block(&mut self, kind: BlockKind, code: Option<&str>) {
        debug_assert!(self.block.is_none());
        self.block = Some(BlockBuilder {
            kind,
            code: code.map(str::to_owned),
            entries: Vec::new(),
        });
    }

    fn end_block(&mut self) {
        let Some(block) = self.block.take() else {
            debug_assert!(false, "parser ended a block that was not open");
            return;
        };
        let block = match block.kind {
            BlockKind::Data => CifBlock::data(block.code.unwrap_or_default(), block.entries),
            BlockKind::Global => CifBlock::global(block.entries),
        };
        self.blocks.push(block);
    }

    fn start_frame(&mut self, code: &str) {
        debug_assert!(self.frame.is_none());
        self.frame = Some(FrameBuilder {
            code: code.to_owned(),
            entries: Vec::new(),
        });
    }

    fn end_frame(&mut self) {
        let Some(frame) = self.frame.take() else {
            debug_assert!(false, "parser ended a frame that was not open");
            return;
        };
        self.push_entry(CifEntry::Frame(CifFrame::new(frame.code, frame.entries)));
    }

    fn item(&mut self, tag: &str, value: ParsedValue<'_>) {
        self.push_entry(CifEntry::Item(CifItem::new(
            tag.to_owned(),
            value.to_cif_value(),
        )));
    }

    fn start_loop(&mut self, tags: &[String]) {
        debug_assert!(self.cif_loop.is_none());
        self.cif_loop = Some(LoopBuilder {
            tags: tags.to_vec(),
            cell_chunks: vec![Vec::new()],
        });
    }

    fn loop_value(&mut self, _column_index: usize, value: ParsedValue<'_>) {
        if let Some(cif_loop) = &mut self.cif_loop {
            if let Some(cells) = cif_loop.cell_chunks.last_mut() {
                cells.push(value.to_source_cell());
            }
        } else {
            debug_assert!(false, "parser emitted a loop value outside a loop");
        }
    }

    fn loop_kernel(&self) -> Option<LoopKernel> {
        self.cif_loop.as_ref().map(|_| LoopKernel::Retain)
    }

    fn parallel_loop(
        &mut self,
        _source: &SourceBuffer,
        _tags: &[String],
        cif_loop: &mut ParallelLoop,
    ) -> Result<bool, ParseError> {
        let Some(builder) = &mut self.cif_loop else {
            return Ok(false);
        };
        builder.cell_chunks = cif_loop.take_source_cells();
        Ok(true)
    }

    fn end_loop(&mut self) {
        let Some(cif_loop) = self.cif_loop.take() else {
            debug_assert!(false, "parser ended a loop that was not open");
            return;
        };
        self.push_entry(CifEntry::Loop(CifLoop::from_source_chunks(
            cif_loop.tags,
            self.source.clone(),
            cif_loop.cell_chunks,
        )));
    }
}

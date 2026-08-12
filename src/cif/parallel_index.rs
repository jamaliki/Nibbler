use std::thread;

use memchr::{memchr, memchr2};

use super::error::{ParseError, ParseErrorCode};
use super::source::SourceBuffer;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct RawChunk {
    pub(super) start: usize,
    pub(super) end: usize,
}

fn is_line_start(bytes: &[u8], byte_offset: usize) -> bool {
    byte_offset == 0 || matches!(bytes[byte_offset - 1], b'\n' | b'\r')
}

const fn is_whitespace(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r')
}

#[derive(Clone, Copy)]
struct BoundarySummary {
    end_state: [bool; 2],
    first_safe: [Option<usize>; 2],
}

pub(super) fn safe_chunks(
    source: &SourceBuffer,
    start: usize,
    workers: usize,
) -> Result<Vec<RawChunk>, ParseError> {
    let bytes = source.as_str().as_bytes();
    if workers <= 1 {
        return Ok(vec![RawChunk {
            start,
            end: bytes.len(),
        }]);
    }
    let targets = chunk_targets(bytes.len(), start, workers);
    let (scan_starts, target_blocks) = boundary_scan_starts(bytes, start, &targets);
    let summaries = thread::scope(|scope| {
        let handles = scan_starts
            .iter()
            .copied()
            .enumerate()
            .map(|(index, block_start)| {
                let block_end = scan_starts.get(index + 1).copied().unwrap_or(bytes.len());
                scope.spawn(move || summarize_boundaries(bytes, block_start, block_end))
            })
            .collect::<Vec<_>>();
        handles
            .into_iter()
            .zip(&scan_starts)
            .map(|(handle, &block_start)| {
                handle.join().map_err(|_| {
                    source.error(
                        ParseErrorCode::UnexpectedControl,
                        "parallel boundary worker terminated unexpectedly",
                        block_start,
                        block_start,
                    )
                })
            })
            .collect::<Result<Vec<_>, _>>()
    })?;

    let mut incoming = Vec::with_capacity(summaries.len());
    let mut in_text_field = false;
    for summary in &summaries {
        incoming.push(in_text_field);
        in_text_field = summary.end_state[usize::from(in_text_field)];
    }

    let mut boundaries = Vec::with_capacity(workers + 1);
    boundaries.push(start);
    for &first_block in &target_blocks {
        let mut block_index = first_block;
        let boundary = loop {
            if block_index == summaries.len() {
                break bytes.len();
            }
            let state_index = usize::from(incoming[block_index]);
            if let Some(candidate) = summaries[block_index].first_safe[state_index] {
                break candidate;
            }
            block_index += 1;
        };
        if boundary != bytes.len() && boundaries.last().copied() != Some(boundary) {
            boundaries.push(boundary);
        }
    }
    boundaries.push(bytes.len());
    Ok(boundaries
        .windows(2)
        .map(|pair| RawChunk {
            start: pair[0],
            end: pair[1],
        })
        .collect())
}

fn summarize_boundaries(bytes: &[u8], start: usize, end: usize) -> BoundarySummary {
    let mut states = [false, true];
    let mut first_safe = [None, None];
    let mut eligible = [Some(start), None];
    let mut cursor = start;
    while let Some(relative) = memchr(b';', &bytes[cursor..end]) {
        let byte_offset = cursor + relative;
        cursor = byte_offset + 1;
        if !is_line_start(bytes, byte_offset) {
            continue;
        }
        for state_index in 0..2 {
            if first_safe[state_index].is_none()
                && eligible[state_index].is_some_and(|candidate| candidate < byte_offset)
            {
                first_safe[state_index] = eligible[state_index];
            }
        }
        let closes = byte_offset + 1 == bytes.len() || is_whitespace(bytes[byte_offset + 1]);
        let next_line = memchr2(b'\n', b'\r', &bytes[byte_offset..end])
            .map(|terminator| byte_offset + terminator + 1)
            .filter(|&candidate| candidate < end);
        for state_index in 0..2 {
            if states[state_index] {
                if closes {
                    states[state_index] = false;
                    eligible[state_index] = next_line;
                }
            } else {
                states[state_index] = true;
                eligible[state_index] = None;
            }
        }
    }
    for state_index in 0..2 {
        if first_safe[state_index].is_none() {
            first_safe[state_index] = eligible[state_index].filter(|&candidate| candidate < end);
        }
    }
    BoundarySummary {
        end_state: states,
        first_safe,
    }
}

fn chunk_targets(len: usize, start: usize, workers: usize) -> Vec<usize> {
    (1..workers)
        .map(|index| start + (len - start) * index / workers)
        .collect()
}

fn boundary_scan_starts(bytes: &[u8], start: usize, targets: &[usize]) -> (Vec<usize>, Vec<usize>) {
    let mut starts = Vec::with_capacity(targets.len() + 1);
    let mut target_blocks = Vec::with_capacity(targets.len());
    starts.push(start);
    for &target in targets {
        let aligned = if target == start || is_line_start(bytes, target) {
            target
        } else {
            memchr2(b'\n', b'\r', &bytes[target..])
                .map_or(bytes.len(), |relative| target + relative + 1)
        };
        if aligned == bytes.len() {
            target_blocks.push(starts.len());
        } else if starts.last().copied() == Some(aligned) {
            target_blocks.push(starts.len() - 1);
        } else {
            starts.push(aligned);
            target_blocks.push(starts.len() - 1);
        }
    }
    (starts, target_blocks)
}

#[cfg(test)]
fn serial_safe_chunks(bytes: &[u8], start: usize, workers: usize) -> Vec<RawChunk> {
    let targets = chunk_targets(bytes.len(), start, workers);
    let mut boundaries = Vec::with_capacity(workers + 1);
    boundaries.push(start);
    let mut target_index = 0;
    let mut byte_offset = start;
    let mut in_text_field = false;
    while byte_offset < bytes.len() && target_index < targets.len() {
        if bytes[byte_offset] == b';' && is_line_start(bytes, byte_offset) {
            if in_text_field {
                if byte_offset + 1 == bytes.len() || is_whitespace(bytes[byte_offset + 1]) {
                    in_text_field = false;
                }
            } else {
                in_text_field = true;
            }
        } else if !in_text_field {
            while target_index < targets.len() && byte_offset >= targets[target_index] {
                if boundaries.last().copied() != Some(byte_offset) {
                    boundaries.push(byte_offset);
                }
                target_index += 1;
            }
        }
        let Some(relative) = memchr2(b'\n', b'\r', &bytes[byte_offset..]) else {
            break;
        };
        byte_offset += relative + 1;
    }
    boundaries.push(bytes.len());
    boundaries
        .windows(2)
        .map(|pair| RawChunk {
            start: pair[0],
            end: pair[1],
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{SourceBuffer, safe_chunks, serial_safe_chunks};

    fn assert_boundary_equivalence(text: &str, start: usize) {
        let source = SourceBuffer::from_text("<boundary-test>", text);
        for workers in [1, 2, 4, 8, 16] {
            let serial = serial_safe_chunks(text.as_bytes(), start, workers);
            assert_eq!(
                safe_chunks(&source, start, workers),
                Ok(serial),
                "workers={workers}, start={start}"
            );
        }
    }

    #[test]
    fn parallel_boundaries_match_serial_text_field_state_machine() {
        let cases = [
            "value\nnext\nlast\n",
            "value\rnext\rlast\r",
            "value\r\nnext\r\nlast\r\n",
            "value\n;\n;\nnext\n",
            "value\n;open\nbody # ' \"\n;\nnext\n",
            "value\r;open\rbody\r;\rnext\r",
            "value\r\n;open\r\nbody\r\n;\r\nnext\r\n",
            "value\n;open\n;not-close\nbody\n;\nnext\n",
            "'quoted ;' # comment ;\n\"quoted ;\"\nnext\n",
            "value\n;open\nbody\nbody2\n",
            "value\n;",
            "value\n;open\nbody\n;",
        ];
        for text in cases {
            assert_boundary_equivalence(text, 0);
        }
        assert_boundary_equivalence("prefix value\nnext\nlast\n", 7);
    }
}

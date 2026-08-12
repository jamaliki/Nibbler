use std::error::Error;
use std::fmt::{self, Display, Formatter};
use std::fs;
use std::io::{self, Read};
use std::path::Path;

use flate2::read::GzDecoder;

use super::error::ParseError;
use super::source::SourceBuffer;

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const GZIP_INITIAL_CAPACITY_MULTIPLIER: usize = 5;
const GZIP_MAX_INITIAL_CAPACITY: usize = 256 * 1024 * 1024;

pub(crate) struct LoadedBytes {
    pub(crate) source_name: String,
    pub(crate) bytes: Vec<u8>,
}

/// Hard limits enforced before retaining decoded source bytes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct InputLimits {
    /// Maximum bytes read from one compressed or uncompressed input.
    pub input_bytes: usize,
    /// Maximum bytes retained after decompression.
    pub decompressed_bytes: usize,
    /// Maximum integer decompressed-to-input size ratio for gzip.
    pub decompression_ratio: usize,
}

impl Default for InputLimits {
    fn default() -> Self {
        Self {
            input_bytes: 2 * 1024 * 1024 * 1024,
            decompressed_bytes: 2 * 1024 * 1024 * 1024,
            decompression_ratio: 1_000,
        }
    }
}

/// Stable categories for file loading and gzip failures.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InputErrorCode {
    /// The source could not be opened or read.
    Io,
    /// Input bytes exceed the configured limit.
    InputLimit,
    /// Decompressed bytes exceed the configured limit.
    DecompressedLimit,
    /// The gzip stream expands beyond the configured ratio.
    DecompressionRatio,
    /// A gzip header, member, checksum, or trailer is malformed or truncated.
    InvalidGzip,
    /// The decoded source is not valid UTF-8.
    InvalidUtf8,
}

impl InputErrorCode {
    /// Return the stable machine-readable diagnostic code.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Io => "CIF_INPUT_IO",
            Self::InputLimit => "CIF_INPUT_LIMIT",
            Self::DecompressedLimit => "CIF_DECOMPRESSED_LIMIT",
            Self::DecompressionRatio => "CIF_DECOMPRESSION_RATIO",
            Self::InvalidGzip => "CIF_GZIP_INVALID",
            Self::InvalidUtf8 => "CIF_INVALID_UTF8",
        }
    }
}

/// A structured failure while loading or decompressing one CIF source.
#[derive(Debug)]
pub struct InputError {
    code: InputErrorCode,
    source_name: String,
    message: String,
    io_error: Option<Box<io::Error>>,
    parse_error: Option<Box<ParseError>>,
}

impl InputError {
    fn new(
        code: InputErrorCode,
        source_name: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self {
            code,
            source_name: source_name.into(),
            message: message.into(),
            io_error: None,
            parse_error: None,
        }
    }

    fn from_io(source_name: &str, code: InputErrorCode, error: io::Error) -> Self {
        Self {
            code,
            source_name: source_name.to_owned(),
            message: error.to_string(),
            io_error: Some(Box::new(error)),
            parse_error: None,
        }
    }

    fn from_parse(error: ParseError) -> Self {
        Self {
            code: InputErrorCode::InvalidUtf8,
            source_name: error.source_name().to_owned(),
            message: error.message().to_owned(),
            io_error: None,
            parse_error: Some(Box::new(error)),
        }
    }

    /// Return the stable input error category.
    #[must_use]
    pub const fn code(&self) -> InputErrorCode {
        self.code
    }

    /// Return the stable machine-readable diagnostic code.
    #[must_use]
    pub const fn code_str(&self) -> &'static str {
        self.code.as_str()
    }

    /// Return the source name or filesystem location supplied by the caller.
    #[must_use]
    pub fn source_name(&self) -> &str {
        &self.source_name
    }

    /// Return the human-readable failure detail.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }
}

impl Display for InputError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{}: {}: {}",
            self.source_name,
            self.code.as_str(),
            self.message
        )
    }
}

impl Error for InputError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        self.io_error
            .as_deref()
            .map(|error| error as &(dyn Error + 'static))
            .or_else(|| {
                self.parse_error
                    .as_deref()
                    .map(|error| error as &(dyn Error + 'static))
            })
    }
}

/// Read, optionally decompress, and UTF-8 validate a filesystem source.
///
/// Compression is detected by gzip magic bytes; the filename suffix is advisory.
///
/// # Errors
///
/// Returns a structured I/O, resource-limit, gzip, or UTF-8 failure.
pub fn read_source_file(
    file: impl AsRef<Path>,
    limits: InputLimits,
) -> Result<SourceBuffer, InputError> {
    let loaded = read_input_file(file, limits)?;
    SourceBuffer::from_owned_bytes(loaded.source_name, loaded.bytes).map_err(InputError::from_parse)
}

pub(crate) fn read_input_file(
    file: impl AsRef<Path>,
    limits: InputLimits,
) -> Result<LoadedBytes, InputError> {
    let file = file.as_ref();
    let source_name = file.to_string_lossy().into_owned();
    let metadata = fs::metadata(file)
        .map_err(|error| InputError::from_io(&source_name, InputErrorCode::Io, error))?;
    let input_bytes = usize::try_from(metadata.len()).unwrap_or(usize::MAX);
    if input_bytes > limits.input_bytes {
        return Err(limit_error(
            InputErrorCode::InputLimit,
            &source_name,
            limits.input_bytes,
            "input bytes",
        ));
    }
    let bytes = fs::read(file)
        .map_err(|error| InputError::from_io(&source_name, InputErrorCode::Io, error))?;
    decode_input(source_name, bytes, limits)
}

/// Optionally decompress and UTF-8 validate one named in-memory source.
///
/// # Errors
///
/// Returns a structured resource-limit, gzip, or UTF-8 failure.
pub fn decode_source(
    source_name: &str,
    bytes: &[u8],
    limits: InputLimits,
) -> Result<SourceBuffer, InputError> {
    let loaded = decode_input(source_name.to_owned(), bytes.to_vec(), limits)?;
    SourceBuffer::from_owned_bytes(loaded.source_name, loaded.bytes).map_err(InputError::from_parse)
}

pub(crate) fn decode_input(
    source_name: String,
    bytes: Vec<u8>,
    limits: InputLimits,
) -> Result<LoadedBytes, InputError> {
    if bytes.len() > limits.input_bytes {
        return Err(limit_error(
            InputErrorCode::InputLimit,
            &source_name,
            limits.input_bytes,
            "input bytes",
        ));
    }
    if !bytes.starts_with(&GZIP_MAGIC) {
        if bytes.len() > limits.decompressed_bytes {
            return Err(limit_error(
                InputErrorCode::DecompressedLimit,
                &source_name,
                limits.decompressed_bytes,
                "decompressed bytes",
            ));
        }
        return Ok(LoadedBytes { source_name, bytes });
    }

    let mut decoder = GzDecoder::new(bytes.as_slice());
    let read_limit = limits.decompressed_bytes.saturating_add(1) as u64;
    let initial_capacity = bytes
        .len()
        .saturating_mul(GZIP_INITIAL_CAPACITY_MULTIPLIER)
        .min(GZIP_MAX_INITIAL_CAPACITY)
        .min(limits.decompressed_bytes)
        .min(bytes.len().saturating_mul(limits.decompression_ratio));
    let mut decompressed = Vec::with_capacity(initial_capacity);
    decoder
        .by_ref()
        .take(read_limit)
        .read_to_end(&mut decompressed)
        .map_err(|error| InputError::from_io(&source_name, InputErrorCode::InvalidGzip, error))?;
    if decompressed.len() > limits.decompressed_bytes {
        return Err(limit_error(
            InputErrorCode::DecompressedLimit,
            &source_name,
            limits.decompressed_bytes,
            "decompressed bytes",
        ));
    }
    if limits.decompression_ratio == 0
        || decompressed.len() > bytes.len().saturating_mul(limits.decompression_ratio)
    {
        return Err(InputError::new(
            InputErrorCode::DecompressionRatio,
            &source_name,
            format!(
                "gzip expands to {} bytes from {} input bytes, exceeding ratio {}",
                decompressed.len(),
                bytes.len(),
                limits.decompression_ratio
            ),
        ));
    }
    Ok(LoadedBytes {
        source_name,
        bytes: decompressed,
    })
}

fn limit_error(code: InputErrorCode, source_name: &str, limit: usize, label: &str) -> InputError {
    InputError::new(
        code,
        source_name,
        format!("source exceeds the configured limit of {limit} {label}"),
    )
}

/*!Module: shared error type used by lexing, parsing, runtime execution, and
file I/O layers.

 */
use std::{fmt, io};

use crate::span::Span;

/// Convenience alias for std::result::Result<T, Error>.
pub type Result<T> = std::result::Result<T, Error>;

/// Error category: io, lex, parse, runtime.
#[derive(Debug)]
pub enum ErrorKind {
    Io(io::Error),
    Lex(String),
    Parse(String),
    Runtime(String),
}

/// Unified error type with file context, source span, and kind discriminator.
#[derive(Debug)]
pub struct Error {
    pub kind: ErrorKind,
    pub span: Option<Span>,
    pub file: Option<String>,
}

impl Error {
    /// wraps an I/O error and optional file name for reporting.
    pub fn io(err: io::Error, file: Option<String>) -> Self {
        Self { kind: ErrorKind::Io(err), span: None, file }
    }
    /// builds a lexical error with optional source span metadata.

    pub fn lex(msg: impl Into<String>, file: Option<String>, span: Option<Span>) -> Self {
        Self { kind: ErrorKind::Lex(msg.into()), span, file }
    }
    /// builds a parser error with optional source span metadata.

    pub fn parse(msg: impl Into<String>, file: Option<String>, span: Option<Span>) -> Self {
        Self { kind: ErrorKind::Parse(msg.into()), span, file }
    }
    /// builds an execution-time error with optional source span metadata.

    pub fn runtime(msg: impl Into<String>, file: Option<String>, span: Option<Span>) -> Self {
        Self { kind: ErrorKind::Runtime(msg.into()), span, file }
    }
}

impl fmt::Display for Error {
    /// renders errors as stable `kind file:start..end: message` text. -
    /// From<io::Error>::from: converts raw I/O failures into the shared error type
    /// when no file/span context is available.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (k, msg) = match &self.kind {
            ErrorKind::Io(e) => ("io", e.to_string()),
            ErrorKind::Lex(s) => ("lex", s.clone()),
            ErrorKind::Parse(s) => ("parse", s.clone()),
            ErrorKind::Runtime(s) => ("runtime", s.clone()),
        };

        match (&self.file, &self.span) {
            (Some(file), Some(span)) => {
                write!(f, "{k} {file}:{}..{}: {msg}", span.start, span.end)
            }
            (Some(file), None) => write!(f, "{k} {file}: {msg}"),
            (None, Some(span)) => write!(f, "{k} {}..{}: {msg}", span.start, span.end),
            (None, None) => write!(f, "{k} {msg}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<io::Error> for Error {
    fn from(value: io::Error) -> Self {
        Error { kind: ErrorKind::Io(value), span: None, file: None }
    }
}

/*!Module: byte-span utilities used to attach source locations to tokens,
parsed commands, and errors.

 */
/// Half-open byte range [start, end) for error location tracking.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}

impl Span {
    /// constructs an immutable half-open source range from start and end byte
    /// offsets.
    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }
}

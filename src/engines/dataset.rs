/*!Module: parser and compact state-set representation for `xread` character
matrices.

 */
use std::collections::BTreeSet;

/// Parsed phylogenetic character matrix: taxa, character states, dimensions.
#[derive(Debug, Clone)]
pub struct Dataset {
    pub title: Option<String>,
    pub nchar: usize,
    pub ntax: usize,
    pub taxa: Vec<String>,     /* original names */
    pub matrix: Vec<StateSet>, /*  ntax * nchar */
}

/// Bit-packed set of up to 36 character states using the low 36 bits of a u64.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StateSet {
    bits: u64, /*  36 bits used */
}

impl StateSet {
    /// Constant bitmask with all 36 state bits set (value is !(!0 << 36)).
    pub const ALL36: StateSet = StateSet { bits: (1u64 << 36) - 1 };
    /// creates a state set with no possible states.

    pub fn empty() -> Self {
        StateSet { bits: 0 }
    }
    /// creates a one-state bitset for a concrete state index.

    pub fn singleton(idx: u8) -> Self {
        StateSet { bits: 1u64 << idx }
    }
    /// adds one concrete state to an existing bitset.

    pub fn insert(&mut self, idx: u8) {
        self.bits |= 1u64 << idx;
    }
    /// reports whether no state bits are set.

    pub fn is_empty(&self) -> bool {
        self.bits == 0
    }
    /// exposes the raw low-36-bit representation used by scoring.

    pub fn bits(&self) -> u64 {
        self.bits
    }
}

/// Parser error specific to xread block format.
#[derive(Debug, Clone)]
pub struct XReadParseError {
    pub message: String,
}

impl XReadParseError {
    /// wraps a parser message in the xread-specific error type.
    fn new(msg: impl Into<String>) -> Self {
        Self { message: msg.into() }
    }
}

/// Parse content AFTER the `xread` keyword up to before the final ';' (raw collected by parser).
/// Supports:
/// - optional title in single quotes: 'title'
/// - dims: nchar ntax  (nchar first!)
/// - ntax taxon rows: `name` `states`
/// - states: 0-9 a-z A-Z (A-Z == a-z), '?' '-' unknown, [ ... ] polymorphic
/// - spaces optional in matrix
pub fn parse_xread_block(raw: &str) -> Result<Dataset, XReadParseError> {
    let mut p = Cursor::new(raw);

    p.skip_ws();
    let title =
        if p.peek_char() == Some('\'') { Some(p.parse_single_quoted_string()?) } else { None };

    p.skip_ws();
    let nchar = p.parse_usize().map_err(|e| XReadParseError::new(e))?;
    p.require_ws("expected whitespace between nchar and ntax")?;
    let ntax = p.parse_usize().map_err(|e| XReadParseError::new(e))?;

    if ntax < 1 || nchar < 1 {
        return Err(XReadParseError::new("nchar and ntax must be positive"));
    }

    let mut taxa = Vec::with_capacity(ntax);
    let mut matrix = Vec::with_capacity(ntax * nchar);
    let mut used_states: BTreeSet<u8> = BTreeSet::new();

    for row in 0..ntax {
        p.skip_ws();

        /* taxon name */
        let name = p.parse_taxon_name()?;
        taxa.push(name);

        /* at least some space before states (allow multiple spaces/newlines) */
        p.skip_ws();

        for col in 0..nchar {
            p.skip_ws();
            let ss = p.parse_state_set(&mut used_states)?;
            if ss.is_empty() {
                return Err(XReadParseError::new(format!(
                    "empty state at taxon {row}, character {col}"
                )));
            }
            matrix.push(ss);
        }
    }

    /* constraint: at most 32 distinct concrete states used (ignore ? and - because they mean ALL) */
    if used_states.len() > 32 {
        return Err(XReadParseError::new(format!(
            "matrix uses {} distinct states (>32). Reduce alphabet usage or allow >32.",
            used_states.len()
        )));
    }

    Ok(Dataset { title, nchar, ntax, taxa, matrix })
}

/*small cursor parser*/

struct Cursor<'a> {
    s: &'a str,
    i: usize,
}

impl<'a> Cursor<'a> {
    /// creates a byte cursor over raw xread text.
    fn new(s: &'a str) -> Self {
        Self { s, i: 0 }
    }
    /// returns the current UTF-8 character without consuming it.

    fn peek_char(&self) -> Option<char> {
        self.s[self.i..].chars().next()
    }
    /// consumes one character and advances the cursor.

    fn bump(&mut self) -> Option<char> {
        let c = self.peek_char()?;
        self.i += c.len_utf8();
        Some(c)
    }
    /// skips whitespace between xread fields.

    fn skip_ws(&mut self) {
        while let Some(c) = self.peek_char() {
            if c.is_whitespace() {
                self.bump();
            } else {
                break;
            }
        }
    }
    /// validates mandatory whitespace between dimensions.

    fn require_ws(&mut self, msg: &str) -> Result<(), XReadParseError> {
        match self.peek_char() {
            Some(c) if c.is_whitespace() => Ok(()),
            _ => Err(XReadParseError::new(msg)),
        }
    }
    /// parses the optional title literal.

    fn parse_single_quoted_string(&mut self) -> Result<String, XReadParseError> {
        let Some('\'') = self.bump() else {
            return Err(XReadParseError::new("expected opening quote"));
        };
        let mut out = String::new();
        while let Some(c) = self.peek_char() {
            if c == '\'' {
                self.bump();
                return Ok(out);
            }
            out.push(c);
            self.bump();
        }
        Err(XReadParseError::new("unterminated single-quoted title"))
    }
    /// parses decimal dimensions.

    fn parse_usize(&mut self) -> Result<usize, String> {
        self.skip_ws();
        let start = self.i;
        while let Some(c) = self.peek_char() {
            if c.is_ascii_digit() {
                self.bump();
            } else {
                break;
            }
        }
        if self.i == start {
            return Err("expected integer".to_string());
        }
        self.s[start..self.i].parse::<usize>().map_err(|_| "invalid integer".to_string())
    }
    /// parses a TNT/Hennig-style taxon identifier.

    fn parse_taxon_name(&mut self) -> Result<String, XReadParseError> {
        /* first char must be alpha; rest alpha/digit/underscore */
        let Some(c0) = self.peek_char() else {
            return Err(XReadParseError::new("unexpected eof reading taxon name"));
        };
        if !c0.is_ascii_alphabetic() {
            return Err(XReadParseError::new("taxon name must start with an alphabetic character"));
        }

        let mut out = String::new();
        while let Some(c) = self.peek_char() {
            if c.is_ascii_alphanumeric() || c == '_' {
                out.push(c);
                self.bump();
            } else {
                break;
            }
        }
        Ok(out)
    }
    /// parses singleton, missing, unknown, or polymorphic state sets and records
    /// observed concrete states.

    fn parse_state_set(
        &mut self,
        used_states: &mut BTreeSet<u8>,
    ) -> Result<StateSet, XReadParseError> {
        let Some(c) = self.peek_char() else {
            return Err(XReadParseError::new("unexpected eof reading state"));
        };

        if c == '?' || c == '-' {
            self.bump();
            return Ok(StateSet::ALL36);
        }

        if c == '[' {
            self.bump(); /* '[' */
            let mut ss = StateSet::empty();
            loop {
                self.skip_ws();
                let Some(nc) = self.peek_char() else {
                    return Err(XReadParseError::new("unterminated [ ... ] state set"));
                };
                if nc == ']' {
                    self.bump();
                    break;
                }
                let idx = map_state_symbol(nc).ok_or_else(|| {
                    XReadParseError::new(format!("invalid state symbol in [..]: {nc:?}"))
                })?;
                used_states.insert(idx);
                ss.insert(idx);
                self.bump();
            }
            if ss.is_empty() {
                return Err(XReadParseError::new("empty [] state set"));
            }
            return Ok(ss);
        }

        /* singleton */
        let idx = map_state_symbol(c).ok_or_else(|| {
            XReadParseError::new(format!(
                "invalid state symbol: {c:?} (allowed: 0-9 a-z A-Z ? - [..])"
            ))
        })?;
        used_states.insert(idx);
        self.bump();
        Ok(StateSet::singleton(idx))
    }
}

/// 0-9 => 0..9 ; a-z/A-Z => 10..35
fn map_state_symbol(c: char) -> Option<u8> {
    if c.is_ascii_digit() {
        return Some((c as u8) - b'0');
    }
    if c.is_ascii_alphabetic() {
        let lower = c.to_ascii_lowercase() as u8;
        if (b'a'..=b'z').contains(&lower) {
            return Some(10 + (lower - b'a'));
        }
    }
    None
}

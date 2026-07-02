/*!Module: lenient parser that converts lexer tokens into command AST nodes
while preserving enough raw text for block commands such as `xread` and
`tread`.

 */
use std::path::PathBuf;

use crate::{
    ast::*,
    error::{Error, Result},
    lexer::{Token, TokenKind},
    span::Span,
};

/// Lenient command parser: converts token streams into AST Command nodes.
pub struct Parser {
    tokens: Vec<Token>,
    p: usize,
    file: Option<String>,
    raw_text: String,
}

impl Parser {
    /// initializes parser state from tokens, file context, and raw source text
    /// used for fallback spans.
    pub fn new(tokens: Vec<Token>, file: Option<String>, raw_text: String) -> Self {
        Self { tokens, p: 0, file, raw_text }
    }

    /// Parse commands leniently by recording parse errors and resynchronizing at
    /// the next semicolon.
    pub fn parse_script_lenient(mut self) -> (Script, Vec<Error>) {
        let mut commands = Vec::new();
        let mut errors = Vec::new();

        while !self.eof() {
            match self.parse_command() {
                Ok(cmd) => commands.push(cmd),
                Err(e) => {
                    errors.push(e);
                    self.sync_to_next_semi();
                }
            }
        }

        (Script { commands }, errors)
    }
    /// skips malformed input until the next command boundary.

    fn sync_to_next_semi(&mut self) {
        while let Some(tok) = self.peek().cloned() {
            self.bump();
            if matches!(tok.kind, TokenKind::Semi) {
                break;
            }
        }
    }
    /// recognizes every supported legacy command and maps aliases/prefixes into
    /// concrete AST variants.

    fn parse_command(&mut self) -> Result<Spanned<Command>> {
        let start = self.peek_span_start();

        let tok = self.peek().cloned().ok_or_else(|| self.err_parse("unexpected eof", None))?;
        let name = match tok.kind {
            TokenKind::Ident(s) => s,
            _ => return Err(self.err_parse("expected command identifier", Some(tok.span))),
        };
        self.bump();
        /* alias name */
        let name_norm = normalize_cmd_name(&name);

        /* block: xread ... ; */
        if eq_ci(name_norm, "xread") {
            if let Some(tok) = self.peek() {
                if matches!(tok.kind, TokenKind::Semi) {
                    self.bump(); /* consume ';' */
                    let span = Span::new(start, self.last_end());
                    return Ok(Spanned::new(Command::XReadQuery, span));
                }
            }
            let raw = self.collect_raw_until_semi()?;
            let span = Span::new(start, self.last_end());
            return Ok(Spanned::new(Command::XReadRaw(raw), span));
        }

        /* block: xx ... ; */
        if eq_ci(name_norm, "xx") {
            match self.peek().cloned() {
                None => {
                    let span = Span::new(start, self.last_end());
                    return Ok(Spanned::new(Command::XxQuery, span));
                }
                Some(tok) if matches!(tok.kind, TokenKind::Semi) => {
                    self.bump(); /* consume ';' */
                    let span = Span::new(start, self.last_end());
                    return Ok(Spanned::new(Command::XxQuery, span));
                }
                Some(_) => {
                    let _raw = self.collect_raw_until_semi()?;
                    let span = Span::new(start, self.last_end());
                    return Ok(Spanned::new(Command::XxQuery, span));
                }
            }
        }

        let cmd = if eq_ci(name_norm, "log") {
            if self.try_consume_discriminant(TokenKind::Star) {
                self.expect_semi()?;
                Command::LogToggle { enabled: true }
            } else if self.try_consume_discriminant(TokenKind::Minus) {
                self.expect_semi()?;
                Command::LogToggle { enabled: false }
            } else if self.try_consume_discriminant(TokenKind::Slash) {
                self.expect_semi()?;
                Command::Log(LogCmd::Stop)
            } else {
                let path = self.parse_path_like()?;
                self.expect_semi()?;
                Command::Log(LogCmd::Start { path })
            }
        } else if eq_ci(name_norm, "display") {
            let cmd = if self.try_consume_discriminant(TokenKind::Star) {
                Command::Display(DisplayCmd::ToTerminalAlso)
            } else if self.try_consume_discriminant(TokenKind::Minus) {
                Command::Display(DisplayCmd::Disable)
            } else {
                Command::Display(DisplayCmd::ToTerminalAlso)
            };
            self.expect_semi()?;
            cmd
        } else if eq_ci(name_norm, "quote") {
            let s = self.collect_tokens_as_string_until_semi()?;
            Command::Quote(s)
        } else if eq_ci(name_norm, "assist") {
            if self.try_consume_discriminant(TokenKind::Star) {
                self.expect_semi()?;
                Command::Assist(AssistSel::All)
            } else if let Some(tok) = self.peek().cloned() {
                if matches!(tok.kind, TokenKind::Semi) {
                    self.bump();
                    Command::Assist(AssistSel::Default)
                } else {
                    let raw = self.collect_tokens_as_string_until_semi()?;
                    Command::Assist(AssistSel::Filter(raw.trim().to_string()))
                }
            } else {
                return Err(self.err_parse("unexpected eof after assist", None));
            }
        } else if eq_ci(name_norm, "view") {
            if self.try_consume_discriminant(TokenKind::Star) {
                self.expect_semi()?;
                Command::ViewLogToggle { enabled: true }
            } else if self.try_consume_discriminant(TokenKind::Minus) {
                self.expect_semi()?;
                Command::ViewLogToggle { enabled: false }
            } else {
                let path = self.parse_path_like()?;
                self.expect_semi()?;
                Command::View(path)
            }
        } else if eq_ci(name_norm, "procedure") {
            if self.try_consume_discriminant(TokenKind::Star) {
                self.expect_semi()?;
                Command::ProcedureToggle { enabled: true }
            } else if self.try_consume_discriminant(TokenKind::Minus) {
                self.expect_semi()?;
                Command::ProcedureToggle { enabled: false }
            } else if self.peek().map(|t| t.kind == TokenKind::Slash).unwrap_or(false)
                && self.tokens.get(self.p + 1).map(|t| t.kind == TokenKind::Semi).unwrap_or(false)
            {
                self.bump(); /* consume the slash */
                self.expect_semi()?;
                Command::ProcedureClose
            } else {
                let path = self.parse_path_like()?;
                self.expect_semi()?;
                Command::ProcedureOpen(path)
            }
        } else if eq_ci(name_norm, "bytes") {
            self.expect_semi()?;
            Command::Bytes
        } else if eq_ci(name_norm, "batch") {
            /* DOS legacy command in original software; not implemented */
            let raw_args = self.collect_tokens_as_string_until_semi()?;
            Command::Batch(raw_args)
        } else if eq_ci(name_norm, "outgroup") {
            if self.try_consume_discriminant(TokenKind::Eq) {
                let nums = self.parse_outgroup_number_set()?;
                self.expect_semi()?;
                Command::Outgroup(OutgroupCmd::SetByNumbers(nums))
            } else if let Some(tok) = self.peek() {
                if matches!(tok.kind, TokenKind::Semi) {
                    self.bump();
                    Command::Outgroup(OutgroupCmd::Query)
                } else {
                    let names = self.parse_outgroup_name_set()?;
                    self.expect_semi()?;
                    Command::Outgroup(OutgroupCmd::SetByNames(names))
                }
            } else {
                return Err(self.err_parse("unexpected eof after outgroup", None));
            }
        } else if eq_ci(name_norm, "reroot") {
            self.expect_semi()?;
            Command::Reroot
        } else if eq_ci(name_norm, "ccode") {
            if let Some(tok) = self.peek() {
                if matches!(tok.kind, TokenKind::Semi) {
                    self.bump();
                    Command::CCode(CCodeCmd::Query)
                } else {
                    let ops = self.parse_ccode_ops()?; /* This consumes the trailing semicolon. */
                    Command::CCode(CCodeCmd::Apply(ops))
                }
            } else {
                return Err(self.err_parse("unexpected eof after ccode", None));
            }
        } else if eq_ci(name_norm, "ckeep") {
            let slot = self.parse_u8()?;
            self.expect_semi()?;
            Command::CKeep(slot)
        } else if eq_ci(name_norm, "cget") {
            let slot = self.parse_u8()?;
            self.expect_semi()?;
            Command::CGet(slot)
        } else if eq_ci(name_norm, "hennig") {
            let star = self.try_consume_discriminant(TokenKind::Star);
            self.expect_semi()?;
            Command::Hennig { multi: false, star }
        } else if eq_ci(name_norm, "mhennig") {
            let star = self.try_consume_discriminant(TokenKind::Star);
            self.expect_semi()?;
            Command::Hennig { multi: true, star }
        } else if eq_ci(name_norm, "bb") {
            let star = self.try_consume_discriminant(TokenKind::Star);
            self.expect_semi()?;
            Command::Bb { star }
        } else if eq_ci(name_norm, "nelsen") {
            self.expect_semi()?;
            Command::Nelsen
        } else if eq_ci(name_norm, "ie") {
            let cmd = if self.try_consume_discriminant(TokenKind::Star) {
                Command::IeStar /* ie*; */
            } else if self.try_consume_discriminant(TokenKind::Minus) {
                Command::IeDash /* ie-; */
            } else {
                Command::Ie /* ie; */
            };
            self.expect_semi()?;
            cmd
        } else if eq_ci(name_norm, "keep") {
            let slot = self.parse_u8()?;
            self.expect_semi()?;
            Command::Keep(slot)
        } else if eq_ci(name_norm, "get") {
            let slot = self.parse_u8()?;
            self.expect_semi()?;
            Command::Get(slot)
        } else if eq_ci(name_norm, "files") {
            self.expect_semi()?;
            Command::Files
        } else if eq_ci(name_norm, "erase") {
            let slot = self.parse_u8()?;
            self.expect_semi()?;
            Command::Erase(slot)
        } else if eq_ci(name_norm, "tchoose") {
            /*
              tchoose <sel...>;
              Support mixed selectors: single numbers (`0`), ranges (`2.3`),
              and the last tree marker (`/`), e.g. `tchoose 0 2.3 /;`.
            */
            let items = self.parse_tchoose_items()?;
            self.expect_semi()?;
            Command::TChoose(items)
        } else if eq_ci(name_norm, "txascii") {
            let enabled = if self.try_consume_discriminant(TokenKind::Minus) {
                false
            } else if self.try_consume_discriminant(TokenKind::Plus) {
                true
            } else {
                true
            };
            self.expect_semi()?;
            Command::TXAscii(enabled)
        } else if eq_ci(name_norm, "tplot") {
            self.expect_semi()?;
            Command::TPlot
        } else if eq_ci(name_norm, "tread") {
            /* tread 'title' <tree> * <tree> ... ; */
            let title = if let Some(tok) = self.peek().cloned() {
                match tok.kind {
                    TokenKind::StringLit(s) => {
                        self.bump();
                        Some(s)
                    }
                    _ => None,
                }
            } else {
                None
            };

            let raw = self.collect_raw_until_semi()?;
            let trees = split_tread_trees(&raw)
                .map_err(|m| self.err_parse(m, Some(Span::new(start, self.last_end()))))?;
            Command::TRead { title, trees }
        } else if eq_ci(name_norm, "tlist") {
            self.expect_semi()?;
            Command::TList
        } else if eq_ci(name_norm, "tsave") {
            let path = self.parse_path_like()?;
            self.expect_semi()?;
            Command::TSave(path)
        } else if eq_ci(name_norm, "xsteps") {
            /* `xsteps;` defaults to `xsteps l;`. */
            if let Some(tok) = self.peek().cloned() {
                if matches!(tok.kind, TokenKind::Semi) {
                    self.bump(); /* consume ';' */
                    Command::XSteps(vec![XStepsMode::L])
                } else {
                    /*
                      Support compact mode lists like `xsteps lchm;`; the
                      lexer emits `lchm` as a single identifier.
                    */
                    let tok2 = self
                        .peek()
                        .cloned()
                        .ok_or_else(|| self.err_parse("expected xsteps mode(s)", None))?;

                    let s = match tok2.kind {
                        TokenKind::Ident(s) => {
                            self.bump();
                            s
                        }
                        _ => {
                            return Err(self.err_parse(
                                "expected xsteps mode identifier (e.g. l, c, h, m, w, u, or lchm)",
                                Some(tok2.span),
                            ));
                        }
                    };

                    let mut modes: Vec<XStepsMode> = Vec::new();
                    for ch in s.chars() {
                        let m = match ch.to_ascii_lowercase() {
                            'l' => XStepsMode::L,
                            'c' => XStepsMode::C,
                            'h' => XStepsMode::H,
                            'm' => XStepsMode::M,
                            'w' => XStepsMode::W,
                            'u' => XStepsMode::U,
                            _ => {
                                return Err(self.err_parse(
                                    format!("unknown xsteps mode: {ch} (allowed: l c h m w)"),
                                    Some(tok2.span),
                                ));
                            }
                        };
                        modes.push(m);
                    }

                    self.expect_semi()?;
                    if modes.is_empty() {
                        return Err(
                            self.err_parse("xsteps requires at least one mode", Some(tok2.span))
                        );
                    }
                    Command::XSteps(modes)
                }
            } else {
                return Err(self.err_parse("expected ';' or xsteps mode(s)", None));
            }
        } else if eq_ci(name_norm, "steps") {
            self.expect_semi()?;
            Command::Steps
        } else if eq_ci(name_norm, "watch") {
            let on = !self.try_consume_discriminant(TokenKind::Minus);
            self.expect_semi()?;
            Command::Watch(on)
        } else if eq_ci(name_norm, "yama") {
            self.expect_semi()?;
            Command::Yama
        } else {
            let raw_args = self.collect_tokens_as_string_until_semi()?;
            Command::Unknown { name, raw_args }
        };

        let span = Span::new(start, self.last_end());
        Ok(Spanned::new(cmd, span))
    }
    /// parses numeric outgroup selectors and inclusive numeric ranges.

    fn parse_outgroup_number_set(&mut self) -> Result<Vec<usize>> {
        let mut out = Vec::new();

        loop {
            let a = self.parse_u32()? as usize;
            if self.try_consume_discriminant(TokenKind::Dot) {
                let b = self.parse_u32()? as usize;
                if a > b {
                    return Err(self.err_parse(format!("invalid outgroup range {a}.{b}"), None));
                }
                for x in a..=b {
                    out.push(x);
                }
            } else {
                out.push(a);
            }

            let Some(tok) = self.peek() else {
                break;
            };
            if matches!(tok.kind, TokenKind::Semi) {
                break;
            }
            if !matches!(tok.kind, TokenKind::Number(_)) {
                return Err(self.err_parse("expected outgroup number or ';'", Some(tok.span)));
            }
        }

        Ok(out)
    }
    /// parses named outgroup taxa from identifiers or quoted names.

    fn parse_outgroup_name_set(&mut self) -> Result<Vec<String>> {
        let mut names = Vec::new();

        loop {
            let tok = self
                .peek()
                .cloned()
                .ok_or_else(|| self.err_parse("expected outgroup name", None))?;
            match tok.kind {
                TokenKind::Ident(s) => {
                    self.bump();
                    names.push(s);
                }
                TokenKind::StringLit(s) => {
                    self.bump();
                    names.push(s);
                }
                TokenKind::Semi => break,
                _ => return Err(self.err_parse("expected outgroup taxon name", Some(tok.span))),
            }

            let Some(tok2) = self.peek() else {
                break;
            };
            if matches!(tok2.kind, TokenKind::Semi) {
                break;
            }
        }

        if names.is_empty() {
            return Err(self.err_parse("outgroup requires at least one taxon", None));
        }

        Ok(names)
    }
    /// parses `tchoose` selectors, including ranges and the `/` last-tree marker.

    fn parse_tchoose_items(&mut self) -> Result<Vec<crate::ast::TChooseItem>> {
        use crate::ast::TChooseItem;
        use crate::lexer::TokenKind;

        let mut out: Vec<TChooseItem> = Vec::new();

        loop {
            let tok = self
                .peek()
                .cloned()
                .ok_or_else(|| self.err_parse("expected tchoose selector", None))?;

            match tok.kind {
                TokenKind::Semi => break,

                /* tchoose /; */
                TokenKind::Slash => {
                    self.bump();
                    out.push(TChooseItem::Last);
                }

                /* tchoose 0 ...; */
                TokenKind::Number(_) => {
                    let a = self.parse_usize()?;
                    if self.try_consume_discriminant(TokenKind::Dot) {
                        /* tchoose 2.3; */
                        let b = self.parse_usize()?;
                        out.push(TChooseItem::RangeInclusive(a, b));
                    } else {
                        out.push(TChooseItem::Index(a));
                    }
                }

                _ => {
                    return Err(self.err_parse(
                        "expected tchoose selector: number, range a.b, '/', or ';'",
                        Some(tok.span),
                    ));
                }
            }

            /* A following semicolon ends the list; otherwise parse another item. */
            let Some(next) = self.peek() else {
                break;
            };
            if matches!(next.kind, TokenKind::Semi) {
                break;
            }
        }

        if out.is_empty() {
            return Err(
                self.err_parse("tchoose requires at least one selector (e.g. 0, 0.3, /)", None)
            );
        }

        Ok(out)
    }
    /// consumes a non-negative integer as a host-size index.
    fn parse_usize(&mut self) -> Result<usize> {
        let tok = self.peek().cloned().ok_or_else(|| self.err_parse("expected number", None))?;
        match tok.kind {
            TokenKind::Number(s) => {
                self.bump();
                s.parse::<usize>().map_err(|_| self.err_parse("invalid usize", Some(tok.span)))
            }
            _ => Err(self.err_parse("expected number", Some(tok.span))),
        }
    }
    /// accepts identifier or quoted path tokens for file commands.

    fn parse_path_like(&mut self) -> Result<PathBuf> {
        let tok = self.peek().cloned().ok_or_else(|| self.err_parse("expected path", None))?;

        if matches!(tok.kind, TokenKind::StringLit(_)) {
            self.bump();
            return Ok(PathBuf::from(match tok.kind {
                TokenKind::StringLit(s) => s,
                _ => unreachable!(),
            }));
        }

        /* Consume one or more ident segments joined by `/` to support
        directory paths like `test/xxx.ss`. */
        let mut parts = Vec::new();
        loop {
            let tok = self.peek().cloned().ok_or_else(|| self.err_parse("expected path", None))?;
            match tok.kind {
                TokenKind::Ident(s) => {
                    self.bump();
                    parts.push(s);
                }
                TokenKind::Slash => {
                    if parts.is_empty() {
                        parts.push(String::new());
                    }
                    self.bump();
                }
                _ => break,
            }
        }

        if parts.is_empty() {
            let tok = self.peek().cloned().ok_or_else(|| self.err_parse("expected path", None))?;
            return Err(self.err_parse("expected path token", Some(tok.span)));
        }

        Ok(PathBuf::from(parts.join("/")))
    }
    /// consumes a non-negative integer for character selectors and weights.

    fn parse_u32(&mut self) -> Result<u32> {
        let tok = self.peek().cloned().ok_or_else(|| self.err_parse("expected number", None))?;
        match tok.kind {
            TokenKind::Number(s) => {
                self.bump();
                s.parse::<u32>().map_err(|_| self.err_parse("invalid u32", Some(tok.span)))
            }
            _ => Err(self.err_parse("expected number", Some(tok.span))),
        }
    }
    /// consumes a small numeric slot index for tree/code files.

    fn parse_u8(&mut self) -> Result<u8> {
        let tok = self.peek().cloned().ok_or_else(|| self.err_parse("expected number", None))?;
        match tok.kind {
            TokenKind::Number(s) => {
                self.bump();
                s.parse::<u8>().map_err(|_| self.err_parse("invalid u8", Some(tok.span)))
            }
            _ => Err(self.err_parse("expected number", Some(tok.span))),
        }
    }
    /// parses one or more ccode operations separated by `*` and terminated by `;`.

    fn parse_ccode_ops(&mut self) -> Result<Vec<CCodeOp>> {
        let mut ops = Vec::new();

        loop {
            /*
              Multiple operations may be separated with `*`; whitespace has
              already been consumed by the lexer.
            */
            let op = self.parse_ccode_one_op()?;
            ops.push(op);

            /* After one operation, either `*` continues or `;` ends the command. */
            if self.try_consume_discriminant(TokenKind::Star) {
                continue;
            }
            self.expect_semi()?;
            break;
        }

        Ok(ops)
    }
    /// parses a single ccode activation, additive, non-additive, deactivation, or
    /// weight operation.

    fn parse_ccode_one_op(&mut self) -> Result<CCodeOp> {
        /* Every operation starts with `+`, `-`, `[`, `]`, or `/N`. */
        if self.try_consume_discriminant(TokenKind::Plus) {
            let sel = self.parse_char_sel_required()?;
            return Ok(CCodeOp::Additive(sel));
        }
        if self.try_consume_discriminant(TokenKind::Minus) {
            let sel = self.parse_char_sel_required()?;
            return Ok(CCodeOp::NonAdditive(sel));
        }
        if self.try_consume_discriminant(TokenKind::LBracket) {
            let sel = self.parse_char_sel_required()?;
            return Ok(CCodeOp::Active(sel));
        }
        if self.try_consume_discriminant(TokenKind::RBracket) {
            let sel = self.parse_char_sel_required()?;
            return Ok(CCodeOp::Inactive(sel));
        }
        if self.try_consume_discriminant(TokenKind::Slash) {
            /* `/N` and `/.` both mean weight `1`. */
            let w =
                if self.try_consume_discriminant(TokenKind::Dot) { 1 } else { self.parse_u32()? };
            let sel = self.parse_char_sel_required()?;
            return Ok(CCodeOp::Weight(w, sel));
        }

        Err(self.err_parse("ccode op must start with + - [ ] or /N", self.peek().map(|t| t.span)))
    }
    /// parses `.` for all characters or a list of numeric character ranges.

    fn parse_char_sel_required(&mut self) -> Result<CharSel> {
        /* `.` selects every character. */
        if self.try_consume_discriminant(TokenKind::Dot) {
            return Ok(CharSel::All);
        }

        /* Otherwise at least one number or range is required. */
        let mut ranges = Vec::new();
        ranges.push(self.parse_char_range()?);

        /* Additional numbers or ranges continue until `*` or `;`. */
        loop {
            let Some(tok) = self.peek() else {
                break;
            };

            /* Stop at an operation or command separator. */
            if matches!(tok.kind, TokenKind::Star) || matches!(tok.kind, TokenKind::Semi) {
                break;
            }

            /*
              Space-separated numeric selectors already arrive as separate
              number tokens, so parsing can continue directly.
            */
            ranges.push(self.parse_char_range()?);
        }

        Ok(CharSel::List(ranges))
    }
    /// parses a single character number or inclusive `a.b` range.

    fn parse_char_range(&mut self) -> Result<CharRange> {
        let a = self.parse_u32()?;
        if self.try_consume_discriminant(TokenKind::Dot) {
            let b = self.parse_u32()?;
            Ok(CharRange::RangeInclusive(a, b))
        } else {
            Ok(CharRange::Single(a))
        }
    }
    /// preserves token text up to a semicolon for raw block commands.

    fn collect_raw_until_semi(&mut self) -> Result<String> {
        let mut s = String::new();
        while let Some(tok) = self.peek().cloned() {
            if matches!(tok.kind, TokenKind::Semi) {
                self.bump();
                break;
            }
            if !s.is_empty() {
                s.push(' ');
            }
            s.push_str(&token_to_string(&tok));
            self.bump();
        }
        Ok(s)
    }
    /// joins ordinary argument tokens into legacy free-form strings.

    fn collect_tokens_as_string_until_semi(&mut self) -> Result<String> {
        let mut s = String::new();
        while let Some(tok) = self.peek().cloned() {
            if matches!(tok.kind, TokenKind::Semi) {
                self.bump();
                break;
            }
            if !s.is_empty() {
                s.push(' ');
            }
            s.push_str(&token_to_string(&tok));
            self.bump();
        }
        Ok(s)
    }
    /// requires and consumes a command terminator.

    fn expect_semi(&mut self) -> Result<()> {
        self.expect_discriminant(TokenKind::Semi)
    }
    /// consumes a token by variant while ignoring any variant payload.

    fn expect_discriminant(&mut self, kind: TokenKind) -> Result<()> {
        let tok = self.peek().cloned().ok_or_else(|| self.err_parse("unexpected eof", None))?;
        if std::mem::discriminant(&tok.kind) == std::mem::discriminant(&kind) {
            self.bump();
            Ok(())
        } else {
            Err(self.err_parse(format!("expected {:?}, got {:?}", kind, tok.kind), Some(tok.span)))
        }
    }
    /// optionally consumes a token by variant.

    fn try_consume_discriminant(&mut self, kind: TokenKind) -> bool {
        let Some(tok) = self.peek() else {
            return false;
        };
        if std::mem::discriminant(&tok.kind) == std::mem::discriminant(&kind) {
            self.bump();
            true
        } else {
            false
        }
    }
    /// returns the current token without advancing.

    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.p)
    }
    /// advances one token.

    fn bump(&mut self) {
        self.p += 1;
    }
    /// reports whether every token has been consumed.

    fn eof(&self) -> bool {
        self.p >= self.tokens.len()
    }
    /// obtains the current command start byte offset.

    fn peek_span_start(&self) -> usize {
        self.peek().map(|t| t.span.start).unwrap_or(self.raw_text.len())
    }
    /// obtains the previous token end byte offset for spans.

    fn last_end(&self) -> usize {
        if self.p == 0 { 0 } else { self.tokens[self.p - 1].span.end }
    }
    /// builds a parser error using the parser file context.

    fn err_parse(&self, msg: impl Into<String>, span: Option<Span>) -> Error {
        Error::parse(msg, self.file.clone(), span)
    }
}
/// compares command fragments case-insensitively.

fn eq_ci(a: &str, b: &str) -> bool {
    a.eq_ignore_ascii_case(b)
}
/// converts a token back to command text for raw argument preservation.

fn token_to_string(tok: &Token) -> String {
    match &tok.kind {
        TokenKind::Ident(s) => s.clone(),
        TokenKind::Number(s) => s.clone(),
        TokenKind::StringLit(s) => format!("'{}'", s),
        TokenKind::Semi => ";".to_string(),
        TokenKind::Eq => "=".to_string(),
        TokenKind::Slash => "/".to_string(),
        TokenKind::Star => "*".to_string(),
        TokenKind::Plus => "+".to_string(),
        TokenKind::Minus => "-".to_string(),
        TokenKind::Dot => ".".to_string(),
        TokenKind::LBracket => "[".to_string(),
        TokenKind::RBracket => "]".to_string(),
        TokenKind::LParen => "(".to_string(),
        TokenKind::RParen => ")".to_string(),
        TokenKind::Comma => ",".to_string(),
        TokenKind::Question => "?".to_string(),
    }
}
/// splits raw tread text into individual parenthetical tree strings separated
/// by `*`.

fn split_tread_trees(raw: &str) -> std::result::Result<Vec<String>, String> {
    /*
      raw is tokens joined with spaces; we accept that.
      Format expected:
      ( ... )
      * ( ... )
      Multiple '*' separators. Commas may appear.
    */
    let mut trees = Vec::new();

    /*
      Split on '*' token (we inserted spaces around tokens via collect_raw_until_semi)
      So raw may contain " * " or leading/trailing spaces. We treat any '*' char as delimiter.
    */
    let parts: Vec<&str> = raw.split('*').collect();
    for part in parts {
        let t = part.trim();
        if t.is_empty() {
            continue;
        }
        /* allow trailing/leading spaces; ensure no trailing ';' (shouldn't be there) */
        let t = t.trim_end_matches(';').trim();
        if !t.starts_with('(') {
            return Err(format!("tread tree must start with '('; got: {}", preview(t)));
        }
        if !t.ends_with(')') {
            return Err(format!("tread tree must end with ')'; got: {}", preview(t)));
        }
        trees.push(t.to_string());
    }

    if trees.is_empty() {
        return Err("tread contained no trees".to_string());
    }
    Ok(trees)
}
/// shortens long malformed tree snippets in diagnostics.

fn preview(s: &str) -> String {
    const N: usize = 80;
    if s.len() <= N { s.to_string() } else { format!("{}...", &s[..N]) }
}
/// maps command abbreviations and quit aliases to their canonical command
/// names.

fn normalize_cmd_name(name: &str) -> &str {
    let n = name.to_ascii_lowercase();
    match n.as_str() {
        "z" | "zz" | "zzz" => return "yama",
        _ => {}
    }
    /* (canonical, min_prefix_len) */
    const CMDS: &[(&str, usize)] = &[
        ("assist", 1),    /* a... */
        ("batch", 5),     /* batch */
        ("bb", 2),        /* bb */
        ("bytes", 2),     /* by.. */
        ("ccode", 2),     /* cc... */
        ("ckeep", 2),     /* ck... */
        ("cget", 2),      /* cg... */
        ("display", 1),   /* d... */
        ("erase", 1),     /* e... */
        ("files", 1),     /* f... */
        ("get", 1),       /* g... */
        ("hennig", 1),    /* h... */
        ("ie", 1),        /* i... */
        ("keep", 1),      /* k... */
        ("log", 1),       /* l... */
        ("mhennig", 1),   /* m... */
        ("nelsen", 1),    /* n... */
        ("outgroup", 1),  /* o... */
        ("procedure", 1), /* p... */
        ("quote", 1),     /* q... */
        ("reroot", 1),    /* r... */
        ("steps", 1),     /* s... */
        ("tchoose", 2),   /* tc... */
        ("tlist", 2),     /* tl... */
        ("tplot", 2),     /* tp... */
        ("tread", 2),     /* tr... */
        ("tsave", 2),     /* ts... */
        ("txascii", 2),   /* tx... */
        ("view", 1),      /* v... */
        ("watch", 1),     /* w... */
        ("xread", 2),     /* xr... */
        ("xsteps", 2),    /* xs... */
        ("xx", 2),        /* xx */
        ("yama", 1),      /* y... */
    ];

    for &(canon, min_len) in CMDS {
        if n.len() >= min_len && canon.starts_with(&n) {
            return canon;
        }
    }

    name
}

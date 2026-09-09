/*!Module: tokenizer for the legacy command language accepted by scripts,
procedure files, and the interactive REPL.

 */
use crate::{
    error::{Error, Result},
    span::Span,
};

/// Every token variant produced by the legacy-command lexer.
#[derive(Debug, Clone, PartialEq)]
pub enum TokenKind {
    Ident(String),
    Number(String),
    StringLit(String),
    Semi,     /* ; */
    Eq,       /* = */
    Slash,    /* / */
    Star,     /* * */
    Amp,      /* & */
    Plus,     /* + */
    Minus,    /* - */
    Dot,      /* . */
    LBracket, /* [ */
    RBracket, /* ] */

    LParen,   /* ( */
    RParen,   /* ) */
    Comma,    /* , */
    Question, /* ? */
}

/// A single lexed token: kind, source text, and byte span.
#[derive(Debug, Clone)]
pub struct Token {
    pub kind: TokenKind,
    pub span: Span,
}

/// Character-level tokenizer that converts raw script text into Token streams.
pub struct Lexer<'a> {
    input: &'a str,
    i: usize,
    file: Option<String>,
}

impl<'a> Lexer<'a> {
    /// creates a tokenizer over raw source text with optional file context for
    /// diagnostics.
    pub fn new(input: &'a str, file: Option<String>) -> Self {
        Self { input, i: 0, file }
    }
    /// scans the entire source into command tokens and returns a lexical error on
    /// the first unsupported character or unterminated string.

    pub fn tokenize(mut self) -> Result<Vec<Token>> {
        let mut out = Vec::new();
        while self.skip_ws_and_comments() {
            if self.i >= self.input.len() {
                break;
            }
            let c = self.peek().unwrap();
            let start = self.i;

            let kind = match c {
                ';' => {
                    self.i += 1;
                    TokenKind::Semi
                }
                '=' => {
                    self.i += 1;
                    TokenKind::Eq
                }
                '/' => {
                    self.i += 1;
                    TokenKind::Slash
                }
                '*' => {
                    self.i += 1;
                    TokenKind::Star
                }
                '&' => {
                    self.i += 1;
                    TokenKind::Amp
                }
                '+' => {
                    self.i += 1;
                    TokenKind::Plus
                }
                '-' => {
                    self.i += 1;
                    TokenKind::Minus
                }
                '.' => {
                    self.i += 1;
                    TokenKind::Dot
                }
                '[' => {
                    self.i += 1;
                    TokenKind::LBracket
                }
                ']' => {
                    self.i += 1;
                    TokenKind::RBracket
                }
                '(' => {
                    self.i += 1;
                    TokenKind::LParen
                }
                ')' => {
                    self.i += 1;
                    TokenKind::RParen
                }
                ',' => {
                    self.i += 1;
                    TokenKind::Comma
                }
                '?' => {
                    self.i += 1;
                    TokenKind::Question
                }
                '\'' => self.lex_string_lit()?,
                c if is_ident_start(c) => self.lex_ident(),
                c if c.is_ascii_digit() => self.lex_number(),
                other => {
                    return Err(Error::lex(
                        format!("unexpected character: {other:?}"),
                        self.file.clone(),
                        Some(Span::new(start, start + 1)),
                    ));
                }
            };

            let end = self.i;
            out.push(Token { kind, span: Span::new(start, end) });
        }
        Ok(out)
    }
    /// reads the current character without advancing the byte cursor.

    fn peek(&self) -> Option<char> {
        self.input[self.i..].chars().next()
    }
    /// consumes one UTF-8 character and advances the byte cursor.

    fn bump(&mut self) -> Option<char> {
        let c = self.peek()?;
        self.i += c.len_utf8();
        Some(c)
    }
    /// skips whitespace and simple `{...}` comments before token recognition.

    fn skip_ws_and_comments(&mut self) -> bool {
        while self.i < self.input.len() {
            let rest = &self.input[self.i..];

            /* Support simple non-nested `{ ... }` comments. */
            if rest.starts_with('{') {
                if let Some(end) = rest.find('}') {
                    self.i += end + 1;
                    continue;
                } else {
                    /* Leave unterminated comments to the tokenizer's normal
                    error path instead of adding a second error branch here. */
                    return true;
                }
            }

            let c = self.peek().unwrap();
            if c.is_whitespace() {
                self.bump();
                continue;
            }
            break;
        }
        true
    }
    /// consumes command names and file-like identifiers.

    fn lex_ident(&mut self) -> TokenKind {
        let mut s = String::new();
        while let Some(c) = self.peek() {
            if is_ident_continue(c) {
                s.push(c);
                self.bump();
            } else {
                break;
            }
        }
        TokenKind::Ident(s)
    }
    /// consumes ASCII decimal integer tokens.

    fn lex_number(&mut self) -> TokenKind {
        let mut s = String::new();
        while let Some(c) = self.peek() {
            if c.is_ascii_digit() {
                s.push(c);
                self.bump();
            } else {
                break;
            }
        }
        TokenKind::Number(s)
    }
    /// consumes single-quoted strings and reports missing closing quotes.

    fn lex_string_lit(&mut self) -> Result<TokenKind> {
        let quote_start = self.i;
        self.bump(); /* consume ' */

        let mut s = String::new();
        while let Some(c) = self.peek() {
            if c == '\'' {
                self.bump();
                return Ok(TokenKind::StringLit(s));
            } else {
                s.push(c);
                self.bump();
            }
        }

        Err(Error::lex(
            "unterminated string literal".to_string(),
            self.file.clone(),
            Some(Span::new(quote_start, self.i)),
        ))
    }
}
/// decides whether a character can begin an identifier token.

fn is_ident_start(c: char) -> bool {
    c.is_ascii_alphabetic() || c == '_'
}
/// decides whether a character can continue an identifier, including dots used
/// by DOS-like file names.

fn is_ident_continue(c: char) -> bool {
    /* Allow file-like identifiers such as `CZ.SS`. */
    c.is_ascii_alphanumeric() || c == '_' || c == '.'
}

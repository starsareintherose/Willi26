/*!
Module: abstract syntax tree for parsed Willi26 commands and command
argument variants.

 */
use std::path::PathBuf;

use crate::span::Span;

/// A parsed command or AST node tagged with its source byte range for error
/// reporting.
#[derive(Debug, Clone)]
pub struct Spanned<T> {
    pub node: T,
    pub span: Span,
}

impl<T> Spanned<T> {
    /// combines a parsed AST node with its source byte span so later runtime
    /// errors can point back to the command text.
    pub fn new(node: T, span: Span) -> Self {
        Self { node, span }
    }
}

/// A complete parsed script: sequence of commands and any parse errors.
#[derive(Debug, Clone)]
pub struct Script {
    pub commands: Vec<Spanned<Command>>,
}

/// Every supported legacy Hennig86 command, including aliases, toggles, and
/// block commands.
#[derive(Debug, Clone)]
pub enum Command {
    /* logging / UI */
    Log(LogCmd),
    LogToggle { enabled: bool }, /* log*; / log-; */
    Display(DisplayCmd),
    Quote(String),
    Assist(AssistSel),

    /* file & include */
    View(PathBuf),
    ViewLogToggle { enabled: bool }, /* view*; / view-; */
    ProcedureOpen(PathBuf),
    ProcedureClose,                    /* procedure/; */
    ProcedureToggle { enabled: bool }, /* procedure*; / procedure-; */

    /* data */
    XReadQuery,
    XReadRaw(String),

    /* outgroup */
    Outgroup(OutgroupCmd),
    Reroot,

    /* characters */
    CCode(CCodeCmd),
    CKeep(u8),
    CGet(u8),

    /* tree search */
    Hennig { multi: bool, star: bool },
    Bb { star: bool },
    Nelsen,

    Ie,     /* ie; */
    IeStar, /* ie*; */
    IeDash, /* ie-; */

    /* internal tree files */
    Keep(u8),
    Get(u8),
    Files,
    Erase(u8),

    /* reporting / tree info */
    TXAscii(bool),
    TPlot(Vec<TreeSelector>),
    TRead { title: Option<String>, trees: Vec<String> },
    TList(Vec<TreeSelector>),
    TSave(PathBuf),
    Tsvg { tree: TreeSelector, path: PathBuf },
    Apo { tree: TreeSelector, path: PathBuf },

    TChoose(Vec<TreeSelector>),

    XSteps(Vec<XStepsMode>),
    Steps,

    Watch(bool),
    Bytes, /* bytes; */

    /* quit */
    Yama,

    /* xx tree editor */
    XxQuery, /* xx; => interactive editor */

    Batch(String),

    Unknown { name: String, raw_args: String },
}

/// Help selection: default display, all commands, or filter by prefix.
#[derive(Debug, Clone)]
pub enum AssistSel {
    Default,        /* assist; */
    All,            /* assist*; */
    Filter(String), /* assist xx; */
}

/// Log-file control: start, stop, toggle.
#[derive(Debug, Clone)]
pub enum LogCmd {
    Start { path: PathBuf },
    Stop,
}

/// Display-mode control: short, none, or all listings.
#[derive(Debug, Clone)]
pub enum DisplayCmd {
    ToTerminalAlso,
    Disable,
}

/// Outgroup selection: query, set by name/number, or toggle.
#[derive(Debug, Clone)]
pub enum OutgroupCmd {
    SetByNames(Vec<String>),
    SetByNumbers(Vec<usize>),
    Query,
}

/// Character-coding control: query, read file, write file.
#[derive(Debug, Clone)]
pub enum CCodeCmd {
    Query,
    Apply(Vec<CCodeOp>),
}

/// Single character-coding operation: additive, weight, activate, deactivate.
#[derive(Debug, Clone)]
pub enum CCodeOp {
    Additive(CharSel),
    NonAdditive(CharSel),
    Active(CharSel),
    Inactive(CharSel),
    Weight(u32, CharSel),
}

/// Character selection: all characters or a list of ranges.
#[derive(Debug, Clone)]
pub enum CharSel {
    All,
    List(Vec<CharRange>),
}

/// Single character index or an inclusive range a..b.
#[derive(Debug, Clone, Copy)]
pub enum CharRange {
    Single(u32),
    RangeInclusive(u32, u32),
}

/// Tree selection used by commands that operate on a tree subset.
#[derive(Debug, Clone)]
pub enum TreeSelector {
    /// Select one zero-based tree index.
    Index(usize),
    /// Select every tree in an inclusive zero-based range.
    RangeInclusive(usize, usize),
    /// Select the last tree in the current working tree set.
    Last,
}

/// xsteps diagnostic mode: histories, character fits, best/worst, lengths,
/// etc.
#[derive(Debug, Clone)]
pub enum XStepsMode {
    L,
    C,
    H,
    M,
    W,
    U,
}

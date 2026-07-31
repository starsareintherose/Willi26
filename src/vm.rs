/*!Module: command virtual machine that executes parsed Willi26 commands,
maintains runtime state, handles the REPL, coordinates engine calls, and
formats command output.

 */
use std::{
    fs::OpenOptions,
    io::{self, Write},
    path::Path,
};

use crate::{
    ast::*,
    engines::{Engine, trees::TreeSet, util, xx},
    error::{Error, Result},
    io::{load_source, resolve_included_path},
    lexer::Lexer,
    parser::Parser,
    span::Span,
    state::State,
};

const ASSIST_CMD_LIST: &str = "\
apo       assist    batch     bb        bytes     ccode     cget\n\
ckeep     display   erase     files     get       hennig    ie\n\
keep      log       mhennig   nelsen    outgroup  procedure quote\n\
reroot    steps     tchoose   tlist     tplot     tread     tsave\n\
tsvg      txascii   view      watch     xread     xsteps    xx\n\
yama";

/// Command virtual machine: executes parsed commands, manages state, runs
/// REPL.
pub struct Vm {
    pub state: State,
    engine: Engine,
    call_stack: Vec<String>, /* Detect recursive procedure inclusion. */
}

impl Vm {
    /// creates a VM with default state, a stateless engine facade, and an empty
    /// procedure call stack.
    pub fn new() -> Self {
        Self { state: State::new(), engine: Engine::new(), call_stack: Vec::new() }
    }

    /// Run a script file.
    /// IO and lexer failures stop execution; parser and runtime errors are
    /// collected and returned after lenient execution continues.
    pub fn run_file(&mut self, path: &Path) -> Result<Vec<Error>> {
        let src = load_source(path)?;
        let file_id = std::fs::canonicalize(path)
            .unwrap_or_else(|_| path.to_path_buf())
            .display()
            .to_string();

        if self.call_stack.contains(&file_id) {
            return Err(Error::runtime("procedure recursion detected", src.file, None));
        }
        self.call_stack.push(file_id);

        /* Lexer failures stop the run immediately. */
        let tokens = Lexer::new(&src.text, src.file.clone()).tokenize()?;

        /* Parse leniently so later commands can still run after syntax errors. */
        let (script, mut errors) =
            Parser::new(tokens, src.file.clone(), src.text.clone()).parse_script_lenient();

        for cmd in script.commands {
            if let Some(e) = self.exec_lenient(cmd, &src.base_dir) {
                errors.push(e);
            }
        }

        self.call_stack.pop();
        Ok(errors)
    }

    /// Run an in-memory script, primarily for the REPL.
    /// Lexer errors stop execution; parse and runtime errors are collected.
    pub fn run_string(&mut self, text: &str, virtual_file: Option<String>) -> Result<Vec<Error>> {
        let tokens = Lexer::new(text, virtual_file.clone()).tokenize()?;
        let (script, mut errors) =
            Parser::new(tokens, virtual_file.clone(), text.to_string()).parse_script_lenient();

        /* Resolve procedure paths from the current working directory in REPL mode. */
        let base_dir = std::env::current_dir().map_err(|e| Error::io(e, virtual_file.clone()))?;

        for cmd in script.commands {
            if let Some(e) = self.exec_lenient(cmd, &base_dir) {
                errors.push(e);
            }
        }
        Ok(errors)
    }

    /// Banner text for the REPL and `version` command, including the current
    /// crate version.
    fn banner() -> String {
        format!(
            "       Willi26 Version {} Copyright (c) Guoyi Zhang 2026\n                     All rights reserved.\n         This copy produced for the exclusive use of\n                       Every  Cladist.",
            env!("CARGO_PKG_VERSION")
        )
    }

    /// Run the interactive command loop until `yama;`, `quit`, or EOF.
    pub fn repl(&mut self) -> Result<()> {
        /* REPL output is visible on the terminal by default. */
        self.state.display_to_terminal = true;

        self.state.write_line(&Self::banner())?;

        let mut buf = String::new();

        loop {
            if buf.trim().is_empty() {
                print!("*> ");
            } else {
                print!(".> ");
            }
            io::stdout().flush().map_err(Error::from)?;

            let mut line = String::new();
            let n = io::stdin().read_line(&mut line).map_err(Error::from)?;
            if n == 0 {
                self.state.write_line("eof").map_err(Error::from)?;
                break;
            }

            let trimmed = line.trim_end_matches(['\n', '\r']).to_string();
            let t = trimmed.trim();

            /* Top-level `quit` exits the REPL without requiring a semicolon. */
            if buf.trim().is_empty() && t.eq_ignore_ascii_case("quit") {
                break;
            }

            /* Accumulate normal command input. */
            buf.push_str(&trimmed);
            buf.push('\n');

            /*
              Execute only after a semicolon, except for the interactive `xx`
              editor shortcut.
            */
            if !buf.contains(';') {
                let trimmed_buf = buf.trim();
                if !trimmed_buf.eq_ignore_ascii_case("xx") {
                    continue;
                }
            }

            /* Reset the quit flag before executing the buffered command. */
            self.state.should_quit = false;
            /* Execute the current command buffer. */
            match self.run_string(&buf, Some("<repl>".to_string())) {
                Ok(errors) => {
                    for e in errors {
                        eprintln!("{e}");
                    }
                }
                Err(e) => {
                    /* Lexer failures stop the REPL run immediately. */
                    eprintln!("{e}");
                    break;
                }
            }
            if self.state.should_quit {
                break;
            }

            if contains_yama_command(&buf) {
                break;
            }

            buf.clear();
        }

        Ok(())
    }
    /// executes one parsed command and converts strict runtime failures into
    /// collectable errors.

    fn exec_lenient(&mut self, cmd: Spanned<Command>, base_dir: &Path) -> Option<Error> {
        let t0 = if self.state.watch_enabled { Some(std::time::Instant::now()) } else { None };

        let out = if let Err(e) = self.exec_strict(cmd, base_dir) { Some(e) } else { None };

        if let Some(t0) = t0 {
            let d = t0.elapsed();
            let ms = d.as_millis();
            let msg = if ms >= 1000 {
                format!("watch {:.2} s", d.as_secs_f64())
            } else {
                format!("watch {ms} ms")
            };
            let _ = self.state.write_line(&msg);
        }

        out
    }
    /// runs implicit enumeration through the engine and stores the resulting best
    /// tree set in VM state.

    fn run_ie_and_store(
        &mut self,
        tag: &str, /* "[ie]" / "[ie*]" / "[ie-]" */
        mode: crate::engines::ie::IeMode,
    ) -> crate::error::Result<()> {
        ensure_tree_context(&self.state)
            .map_err(|m| crate::error::Error::runtime(m, None, None))?;

        let ds = self.state.dataset.as_ref().unwrap();
        let cfg = self.state.char_config.as_ref().unwrap();
        let outgroup = require_outgroup_index(&self.state)?;

        let outcome = self.engine.ie(ds, cfg, outgroup, mode)?;

        /* TL/CI/RI */
        let min_len = crate::engines::search::minsteps_sum(ds, cfg);
        let max_len = crate::engines::search::maxsteps_sum(ds, cfg);
        let ci = crate::engines::search::calc_ci(min_len, outcome.best_len);
        let ri = crate::engines::search::calc_ri(min_len, max_len, outcome.best_len);

        let ntrees = outcome.trees.len();

        self.write_result_line(tag, outcome.best_len, ci, ri, &format!(" trees={}", ntrees))?;

        /* Store exact-enumeration results in the working tree set. */
        let title = format!("set of {} trees", ntrees);
        self.state.append_to_working_tree_set(title, outcome.trees);

        Ok(())
    }
    /// executes a nested procedure file while honoring the procedure activation
    /// switch and recursion guard.

    fn run_procedure_file(&mut self, path: &std::path::Path) -> Result<Vec<Error>> {
        let src = load_source(path)?;
        let file_id = std::fs::canonicalize(path)
            .unwrap_or_else(|_| path.to_path_buf())
            .display()
            .to_string();

        if self.call_stack.contains(&file_id) {
            return Err(Error::runtime("procedure recursion detected", src.file, None));
        }
        self.call_stack.push(file_id);

        let tokens = Lexer::new(&src.text, src.file.clone()).tokenize()?;
        let (script, mut errors) =
            Parser::new(tokens, src.file.clone(), src.text.clone()).parse_script_lenient();

        /* Resolve nested paths relative to the procedure file. */
        let base_dir = src.base_dir.clone();

        for cmd in script.commands {
            match &cmd.node {
                Command::ProcedureClose => {
                    /* A procedure file closes normally only via `procedure/;`. */
                    self.call_stack.pop();
                    self.state.write_line("procedure --").map_err(Error::from)?;
                    return Ok(errors);
                }
                Command::ProcedureToggle { enabled } => {
                    self.state.procedure_exec_enabled = *enabled;
                    self.state
                        .write_line(&format!(
                            "procedure {} {}",
                            if *enabled { "*" } else { "-" },
                            if *enabled { "enabled" } else { "disabled" }
                        ))
                        .map_err(Error::from)?;
                }
                _ => {
                    if self.state.procedure_exec_enabled {
                        if let Some(e) = self.exec_lenient(cmd, &base_dir) {
                            errors.push(e);
                        }
                    } else {
                        /*
                          Skip ordinary commands while disabled, but keep
                          reading so toggle and close commands can take effect.
                        */
                    }
                }
            }
        }

        self.call_stack.pop();

        /* Reaching EOF without `procedure/;` is a runtime error. */
        Err(Error::runtime(
            format!(
                "procedure not closed: missing procedure/; in {}",
                src.file.clone().unwrap_or_else(|| path.display().to_string())
            ),
            src.file,
            None,
        ))
    }
    /// dispatches every concrete AST command to the corresponding state, parser,
    /// tree, search, report, or file operation.

    fn exec_strict(&mut self, cmd: Spanned<Command>, base_dir: &Path) -> Result<()> {
        let span = cmd.span;

        match cmd.node {
            Command::Log(lc) => match lc {
                LogCmd::Start { path } => {
                    let full = resolve_included_path(base_dir, &path);
                    let f = std::fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(&full)
                        .map_err(|e| Error::io(e, Some(full.display().to_string())))?;

                    self.state.log_file = Some(f);
                    self.state.log_path = Some(full.clone());
                    self.state.log_enabled = true;

                    self.state
                        .write_line(&format!("log {} *", full.display()))
                        .map_err(Error::from)?;
                }
                LogCmd::Stop => {
                    self.state.log_enabled = false;
                    self.state.log_file = None;
                    self.state.log_path = None;

                    self.state.write_line("log --").map_err(Error::from)?;
                }
            },

            Command::LogToggle { enabled } => {
                if enabled {
                    let path = self.state.log_path.clone().ok_or_else(|| {
                        Error::runtime(
                            "log* no log file path set (use: log <path>; first)",
                            None,
                            Some(span),
                        )
                    })?;

                    let f = std::fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(&path)
                        .map_err(|e| Error::io(e, Some(path.display().to_string())))?;

                    self.state.log_file = Some(f);
                    self.state.log_enabled = true;

                    self.state
                        .write_line(&format!("log {} *", path.display()))
                        .map_err(Error::from)?;
                } else {
                    self.state.log_enabled = false;
                    self.state.log_file = None; /* Keep log_path so `log*;` can reopen it. */

                    let lp = self.state.log_path.clone().unwrap_or_default();
                    self.state
                        .write_line(&format!("log {} -", lp.display()))
                        .map_err(Error::from)?;
                }
            }

            Command::ViewLogToggle { enabled } => {
                /* `enabled=true` implements the `view*;` behavior. */
                if !enabled {
                    self.state.write_line("view -").map_err(Error::from)?;
                } else {
                    let path = self.state.log_path.clone().ok_or_else(|| {
                        Error::runtime("view*: no log file path set", None, Some(span))
                    })?;

                    /* Ensure the viewed log is visible on the terminal. */
                    self.state.display_to_terminal = true;

                    let text = std::fs::read_to_string(&path)
                        .map_err(|e| Error::io(e, Some(path.display().to_string())))?;

                    self.state
                        .write_line(&format!("view {}", path.display()))
                        .map_err(Error::from)?;

                    for line in text.lines() {
                        self.state.write_line(line).map_err(Error::from)?;
                    }
                }
            }

            Command::Display(dc) => match dc {
                DisplayCmd::ToTerminalAlso => {
                    self.state.display_to_terminal = true;
                    self.state.write_line("display * terminal + log").map_err(Error::from)?;
                }
                DisplayCmd::Disable => {
                    self.state.display_to_terminal = false;
                    self.state.write_line("display --").map_err(Error::from)?;
                }
            },

            Command::Quote(s) => {
                /* consider push space */
                let out = s.replace(". ,", ";");
                self.state.write_line(&format!("quote {out}")).map_err(Error::from)?;
            }

            Command::Batch(raw) => {
                self.state
                    .write_line(&format!(
                        "DOS legacy command not implemented in this version. args='{}'",
                        raw
                    ))
                    .map_err(Error::from)?;
            }

            Command::Assist(sel) => {
                match sel {
                    AssistSel::Default => {
                        self.state.write_line(ASSIST_CMD_LIST).map_err(Error::from)?;
                    }
                    AssistSel::All => {
                        /* Print full help text for every topic. */
                        for t in crate::help::all_topics() {
                            self.state.write_line("").map_err(Error::from)?;
                            for &line in t.details {
                                self.state.write_line(line).map_err(Error::from)?;
                            }
                        }
                    }
                    AssistSel::Filter(prefix) => {
                        let hits = crate::help::find_topics_by_prefix(&prefix);
                        if hits.is_empty() {
                            self.state.write_line("a;").map_err(Error::from)?;
                            self.state.write_line(&format!(" ? {prefix}")).map_err(Error::from)?;
                        } else {
                            for t in hits {
                                self.state.write_line("").map_err(Error::from)?;
                                for &line in t.details {
                                    self.state.write_line(line).map_err(Error::from)?;
                                }
                            }
                        }
                    }
                }
            }

            Command::View(path) => {
                let full = resolve_included_path(base_dir, &path);
                let text = std::fs::read_to_string(&full)
                    .map_err(|e| Error::io(e, Some(full.display().to_string())))?;
                self.state.write_line(&format!("view {}", full.display())).map_err(Error::from)?;
                for line in text.lines() {
                    self.state.write_line(line).map_err(Error::from)?;
                }
            }

            Command::ProcedureOpen(path) => {
                let full = resolve_included_path(base_dir, &path);
                self.state
                    .write_line(&format!("procedure {} *", full.display()))
                    .map_err(Error::from)?;

                /* Entering a procedure enables command execution by default. */
                self.state.procedure_exec_enabled = true;

                let include_errors = self.run_procedure_file(&full)?;
                if !include_errors.is_empty() {
                    self.state
                        .write_line(&format!(
                            "procedure included file had {} error(s)",
                            include_errors.len()
                        ))
                        .map_err(Error::from)?;
                    for e in include_errors {
                        self.state.write_line(&format!("{e}")).map_err(Error::from)?;
                    }
                }
            }

            Command::ProcedureClose => {
                return Err(Error::runtime(
                    "procedure/; is only valid inside a procedure file",
                    None,
                    Some(span),
                ));
            }

            Command::ProcedureToggle { enabled: _ } => {
                /*
                  Procedure toggles are only meaningful while a procedure file
                  is being interpreted; at top level they are command errors.
                */
                return Err(Error::runtime(
                    "procedure*; / procedure-; are only valid inside a procedure file",
                    None,
                    Some(span),
                ));
            }

            Command::Bytes => match read_mem_available_bytes_linux() {
                Ok(Some(bytes)) => {
                    self.state
                        .write_line(&format!("bytes {}", format_bytes_human(bytes)))
                        .map_err(Error::from)?;
                }
                Ok(None) => {
                    self.state
                        .write_line("bytes This platform is not supported")
                        .map_err(Error::from)?;
                }
                Err(m) => {
                    return Err(Error::runtime(format!("bytes failed: {m}"), None, Some(span)));
                }
            },

            Command::XReadRaw(raw) => match crate::engines::dataset::parse_xread_block(&raw) {
                Ok(ds) => {
                    let nchar = ds.nchar;
                    let title = ds.title.clone().unwrap_or_default();
                    self.state.dataset = Some(ds);
                    self.state.char_config = Some(crate::engines::ccode::CharConfig::new(nchar));
                    self.state.clear_all_trees();
                    self.state.outgroup = Some(crate::state::OutgroupState::new(vec![0]));
                    self.state.write_line("xread").map_err(Error::from)?;
                    self.state.write_line(&title).map_err(Error::from)?;
                }
                Err(e) => {
                    return Err(Error::runtime(
                        format!("xread parse error: {}", e.message),
                        None,
                        Some(span),
                    ));
                }
            },

            Command::XReadQuery => {
                let dump = {
                    let Some(ds) = self.state.dataset.as_ref() else {
                        return Err(Error::runtime("xread: no dataset loaded", None, Some(span)));
                    };

                    let mut out = String::new();
                    out.push_str("xread\n");
                    out.push_str(&format!("'{}'\n", ds.title.clone().unwrap_or_default()));
                    out.push_str(&format!("{} {}\n", ds.nchar, ds.ntax));

                    for t in 0..ds.ntax {
                        let name = &ds.taxa[t];
                        out.push_str(name);
                        out.push(' ');

                        for c in 0..ds.nchar {
                            let ss = ds.matrix[t * ds.nchar + c];
                            out.push_str(&format_stateset(ss));
                        }
                        out.push('\n');
                    }

                    out.push_str(";\n");
                    out
                }; /* End the immutable dataset borrow before writing output. */

                /* It is now safe to mutably borrow state for output. */
                for line in dump.lines() {
                    self.state.write_line(line).map_err(Error::from)?;
                }
            }

            Command::Outgroup(oc) => match oc {
                OutgroupCmd::SetByNames(names) => {
                    self.state.ensure_dataset().map_err(|m| Error::runtime(m, None, Some(span)))?;
                    let ds = self.state.dataset.as_ref().unwrap();

                    let mut taxa = Vec::new();
                    for name in names {
                        let idx = ds.taxa.iter().position(|t| t == &name).ok_or_else(|| {
                            Error::runtime(format!("outgroup ? '{name}'"), None, Some(span))
                        })?;
                        taxa.push(idx);
                    }

                    let resolved =
                        resolve_outgroup_request_from_taxa(&taxa, self.state.working_tree_set());
                    self.state.outgroup = Some(crate::state::OutgroupState::new(resolved.clone()));

                    let desc = self
                        .state
                        .outgroup
                        .as_ref()
                        .unwrap()
                        .describe_with_dataset(self.state.dataset.as_ref());
                    self.state.write_line(&format!("outgroup {desc}")).map_err(Error::from)?;
                }
                OutgroupCmd::SetByNumbers(nums) => {
                    self.state.ensure_dataset().map_err(|m| Error::runtime(m, None, Some(span)))?;
                    let ds = self.state.dataset.as_ref().unwrap();

                    for &t in &nums {
                        if t >= ds.ntax {
                            return Err(Error::runtime(
                                format!("outgroup {t} out of range (ntax={})", ds.ntax),
                                None,
                                Some(span),
                            ));
                        }
                    }

                    let resolved =
                        resolve_outgroup_request_from_taxa(&nums, self.state.working_tree_set());
                    self.state.outgroup = Some(crate::state::OutgroupState::new(resolved.clone()));

                    let desc = self
                        .state
                        .outgroup
                        .as_ref()
                        .unwrap()
                        .describe_with_dataset(self.state.dataset.as_ref());
                    self.state.write_line(&format!("outgroup {desc}")).map_err(Error::from)?;
                }
                OutgroupCmd::Query => {
                    let desc = self
                        .state
                        .outgroup
                        .as_ref()
                        .map(|og| og.describe_with_dataset(self.state.dataset.as_ref()))
                        .unwrap_or_else(|| "none".to_string());

                    self.state.write_line(&format!("outgroup {desc}")).map_err(Error::from)?;
                }
            },

            Command::Reroot => {
                self.state.ensure_dataset().map_err(|m| Error::runtime(m, None, Some(span)))?;

                let outgroup =
                    self.state.outgroup.as_ref().cloned().ok_or_else(|| {
                        Error::runtime("reroot no outgroup set", None, Some(span))
                    })?;

                let ts = self.state.working_tree_set().cloned().ok_or_else(|| {
                    Error::runtime("reroot no trees in working slot 0", None, Some(span))
                })?;

                let mut rerooted = Vec::with_capacity(ts.trees.len());
                for tr in &ts.trees {
                    let rr = tr.reroot_by_outgroup_set(&outgroup.taxa).map_err(|m| {
                        Error::runtime(format!("reroot failed {m}"), None, Some(span))
                    })?;
                    rerooted.push(rr);
                }

                self.state.set_working_tree_set(TreeSet::new(ts.title, rerooted));
                self.state.write_line("reroot").map_err(Error::from)?;
            }

            Command::CCode(cc) => match cc {
                CCodeCmd::Query => {
                    let dump = {
                        let Some(cfg) = self.state.char_config.as_ref() else {
                            return Err(Error::runtime(
                                "ccode no dataset loaded",
                                None,
                                Some(span),
                            ));
                        };

                        let mut out = String::new();
                        out.push_str("Ccode\n");
                        let per_row: usize = 5;
                        let mut col = 0usize;
                        for (i, ch) in cfg.chars.iter().enumerate() {
                            if col == 0 {
                                if i > 0 {
                                    out.push('\n');
                                }
                                out.push_str("   ");
                            }
                            let flag = if ch.additive { "+" } else { "-" };
                            let bracket = if ch.active { "[" } else { "]" };
                            let s = format!("{flag}{bracket}/{}  {i}", ch.weight);
                            out.push_str(&format!("{:16}", s));
                            col += 1;
                            if col == per_row {
                                col = 0;
                            }
                        }
                        if col != 0 {
                            out.push_str("   ;");
                        }
                        out.push('\n');
                        out
                    };

                    for line in dump.lines() {
                        self.state.write_line(line).map_err(Error::from)?;
                    }
                }
                CCodeCmd::Apply(ops) => {
                    let Some(cfg) = self.state.char_config.as_mut() else {
                        return Err(Error::runtime("ccode no dataset loaded", None, Some(span)));
                    };
                    if let Err(msg) = apply_ccode_ops(cfg, &ops) {
                        return Err(Error::runtime(msg, None, Some(span)));
                    }
                    self.state.write_line("ccode updated").map_err(Error::from)?;
                }
            },

            Command::CKeep(slot) => {
                let Some(cfg) = self.state.char_config.as_ref() else {
                    return Err(Error::runtime("ckeep no dataset loaded", None, Some(span)));
                };
                self.state.ccode_slots.insert(slot, cfg.clone());
                self.state
                    .write_line(&format!("ckeep saved ccode to {slot}"))
                    .map_err(Error::from)?;
            }

            Command::CGet(slot) => {
                let Some(saved) = self.state.ccode_slots.get(&slot).cloned() else {
                    return Err(Error::runtime(
                        format!("cget failed {slot} not found"),
                        None,
                        Some(span),
                    ));
                };
                /* Saved coding must match the current dataset character count. */
                if let Some(cur) = self.state.char_config.as_ref() {
                    if cur.len() != saved.len() {
                        return Err(Error::runtime(
                            format!(
                                "cget failed {slot} has nchar={}, current nchar={}",
                                saved.len(),
                                cur.len()
                            ),
                            None,
                            Some(span),
                        ));
                    }
                }
                self.state.char_config = Some(saved);
                self.state
                    .write_line(&format!("cget loaded ccode from {slot}"))
                    .map_err(Error::from)?;
            }

            Command::Hennig { multi, star } => {
                self.state.ensure_dataset().map_err(|m| Error::runtime(m, None, Some(span)))?;

                let ds = self.state.dataset.as_ref().unwrap();
                let cfg = self.state.char_config.as_ref().unwrap();

                let outgroup = 0usize;

                let outcome = self.engine.hennig(ds, cfg, multi, star, outgroup)?;

                let title = if multi {
                    if star { "MHennig* trees" } else { "MHennig trees" }
                } else {
                    if star { "Hennig* trees" } else { "Hennig trees" }
                }
                .to_string();

                let ntrees = outcome.trees.len();
                let best_len = outcome.best_len;
                let ci = outcome.ci;
                let ri = outcome.ri;
                let trees = outcome.trees;

                self.state.append_to_working_tree_set(title, trees);

                if !multi {
                    self.write_result_line(
                        &format!("hennig{}", if star { "*" } else { "" }),
                        best_len,
                        ci,
                        ri,
                        "",
                    )?;
                } else {
                    self.write_result_line(
                        &format!("mhennig{}", if star { "*" } else { "" }),
                        best_len,
                        ci,
                        ri,
                        &format!(" trees={}", ntrees),
                    )?;
                }
            }

            Command::Bb { star } => {
                self.state.ensure_dataset().map_err(|m| Error::runtime(m, None, Some(span)))?;

                let ds = self.state.dataset.as_ref().unwrap();
                let cfg = self.state.char_config.as_ref().unwrap();

                let ts0 = self
                    .state
                    .working_tree_set()
                    .cloned()
                    .ok_or_else(|| Error::runtime("bb no trees in 0", None, Some(span)))?;

                if ts0.trees.is_empty() {
                    return Err(Error::runtime("bb 0 is empty", None, Some(span)));
                }

                let outgroup = 0usize;
                let mut ws = crate::engines::search::ScoreWorkspace::new();

                /*
                  Score all working trees and keep near-best starting trees
                  within a slack to give TBR closure more diversity.
                */
                const BB_START_SLACK: u64 = 5;
                let mut start_best_len = u64::MAX;
                let mut start_bucket: Vec<crate::engines::trees::Tree> = Vec::new();

                for tr in &ts0.trees {
                    let norm = tr.reroot_by_outgroup(outgroup).unwrap_or_else(|_| tr.clone());
                    let len = crate::engines::search::score_tree(ds, cfg, &norm, &mut ws);

                    if len < start_best_len {
                        start_best_len = len;
                        start_bucket.clear();
                        start_bucket.push(norm);
                    } else if len <= start_best_len + BB_START_SLACK {
                        start_bucket.push(norm);
                    }
                }

                /* Canonicalize and deduplicate the starting bucket. */
                let mut seen_start = std::collections::BTreeSet::<String>::new();
                let mut start_trees: Vec<crate::engines::trees::Tree> = Vec::new();
                for mut tr in start_bucket {
                    tr.canonicalize();
                    let key = tr.to_storage_string();
                    if seen_start.insert(key) {
                        start_trees.push(tr);
                    }
                }

                if start_trees.is_empty() {
                    return Err(Error::runtime("bb no start trees after uniq", None, Some(span)));
                }

                /* Run closure-style branch breaking for `bb` or `bb*`. */
                let (mut best_trees, best_len) =
                    self.engine.bb(ds, cfg, &start_trees, outgroup, star)?;

                /*
                  Canonicalize and deduplicate once more before replacing the
                  working slot.
                */
                let mut seen_out = std::collections::BTreeSet::<String>::new();
                let mut out_trees = Vec::new();
                for mut tr in best_trees.drain(..) {
                    tr.canonicalize();
                    let key = tr.to_storage_string();
                    if seen_out.insert(key) {
                        out_trees.push(tr);
                    }
                }

                let n = out_trees.len();

                self.state.set_working_tree_set(TreeSet::new(
                    if star { "bb* trees" } else { "bb trees" },
                    out_trees,
                ));

                self.state
                    .write_line(&format!(
                        "bb{} start_trees={} start_best_TL={} best_TL={} trees={}",
                        if star { "*" } else { "" },
                        start_trees.len(),
                        start_best_len,
                        best_len,
                        n
                    ))
                    .map_err(Error::from)?;
            }

            Command::Nelsen => {
                ensure_tree_context(&self.state)
                    .map_err(|m| Error::runtime(m, None, Some(span)))?;

                let ds = self.state.dataset.as_ref().unwrap();
                let cfg = self.state.char_config.as_ref().unwrap();

                let Some(ts) = self.state.current_tree_set() else {
                    return Err(Error::runtime("nelsen: no current tree set", None, Some(span)));
                };

                if ts.trees.is_empty() {
                    return Err(Error::runtime(
                        "nelsen: current tree set is empty",
                        None,
                        Some(span),
                    ));
                }

                let outgroups = self
                    .state
                    .outgroup
                    .as_ref()
                    .map(|og| og.taxa.clone())
                    .ok_or_else(|| Error::runtime("outgroup not set", None, Some(span)))?;

                if outgroups.is_empty() {
                    return Err(Error::runtime("nelsen outgroup set is empty", None, Some(span)));
                }

                let outcome = self.engine.nelsen_consensus(ds, cfg, &ts.trees, &outgroups)?;

                /*
                  Append the consensus tree to the current slot rather than
                  replacing the source trees.
                */
                self.state.append_to_current_tree_set(
                    format!("strict consensus of {} trees", outcome.uniq_trees),
                    vec![outcome.tree],
                );

                self.write_result_line(
                    "nelsen",
                    outcome.best_len,
                    outcome.ci,
                    outcome.ri,
                    &format!(" input={} uniq={} trees=1", outcome.input_trees, outcome.uniq_trees),
                )?;
            }

            Command::Ie => {
                self.run_ie_and_store("ie", crate::engines::ie::IeMode::KeepUpTo(100))?;
            }

            Command::IeStar => {
                self.run_ie_and_store("ie*", crate::engines::ie::IeMode::KeepAll)?;
            }

            Command::IeDash => {
                self.run_ie_and_store("ie-", crate::engines::ie::IeMode::KeepOne)?;
            }

            Command::TRead { title, trees } => {
                if let Err(msg) = ensure_tree_context(&self.state) {
                    return Err(Error::runtime(msg, None, Some(span)));
                }

                let ds = self.state.dataset.as_ref().unwrap();

                let normalized = normalize_and_validate_trees(ds, &trees)
                    .map_err(|m| Error::runtime(m, None, Some(span)))?;

                let mut parsed = Vec::with_capacity(normalized.len());
                for s in normalized {
                    let tr = crate::engines::trees::Tree::from_storage_string(&s, ds.ntax)
                        .map_err(|m| {
                            Error::runtime(format!("tread parse failed: {m}"), None, Some(span))
                        })?;
                    parsed.push(tr);
                }

                let final_title: String =
                    match title.as_ref().map(|s| s.trim()).filter(|s| !s.is_empty()) {
                        Some(non_empty) => non_empty.to_string(),
                        None => format!("set of {} trees", parsed.len()),
                    };

                self.state.append_to_working_tree_set(final_title.clone(), parsed);

                let current_len =
                    self.state.working_tree_set().map(|ts| ts.trees.len()).unwrap_or(0);

                self.state
                    .write_line(&format!(
                        "tread\ncurrent treeset '{}' now has {} tree(s)",
                        final_title, current_len
                    ))
                    .map_err(Error::from)?;
            }

            Command::Keep(slot) => {
                if let Err(msg) = ensure_tree_context(&self.state) {
                    return Err(Error::runtime(msg, None, Some(span)));
                }

                let ts =
                    self.state.tree_slots.remove(&State::WORKING_TREE_SLOT).ok_or_else(|| {
                        Error::runtime("keep failed: no trees in working set 0", None, Some(span))
                    })?;

                self.state.tree_slots.insert(slot, ts);
                /*
                  `keep` moves the working slot into a vault slot and makes it
                  current.
                */
                self.state.current_tree_slot = Some(slot);

                self.state
                    .write_line(&format!("keep moved working tree set into {slot}"))
                    .map_err(Error::from)?;
            }

            Command::Get(slot) => {
                if let Err(msg) = ensure_tree_context(&self.state) {
                    return Err(Error::runtime(msg, None, Some(span)));
                }

                let ts = self.state.tree_slots.get(&slot).cloned().ok_or_else(|| {
                    Error::runtime(
                        format!("get failed: tree set {slot} not found"),
                        None,
                        Some(span),
                    )
                })?;

                /* Append the selected vault slot into working slot 0. */
                self.state.append_to_working_tree_set(ts.title.clone(), ts.trees.clone());

                /* Refresh an auto-generated working title with the real count. */
                self.state.refresh_working_title_auto();

                let n = self.state.working_tree_set().map(|x| x.trees.len()).unwrap_or(0);
                self.state
                    .write_line(&format!("get appended {slot} into working 0 (now {n} tree(s))"))
                    .map_err(Error::from)?;
            }

            Command::Files => {
                if let Err(msg) = ensure_tree_context(&self.state) {
                    return Err(Error::runtime(msg, None, Some(span)));
                }

                let mut keys: Vec<u8> = self.state.tree_slots.keys().copied().collect();
                keys.sort();

                if keys.is_empty() {
                    self.state.write_line("files\nno tree").map_err(Error::from)?;
                } else {
                    self.state.write_line("files tree:").map_err(Error::from)?;
                    for k in keys {
                        let ts = &self.state.tree_slots[&k];
                        let mut tags = Vec::new();
                        if k == State::WORKING_TREE_SLOT {
                            tags.push("working");
                        }
                        if self.state.current_tree_slot == Some(k) {
                            tags.push("current");
                        }
                        let tag_s = if tags.is_empty() {
                            "".to_string()
                        } else {
                            format!(" [{}]", tags.join(","))
                        };

                        self.state
                            .write_line(&format!(
                                "  tree file {k}: '{}' ({} tree(s)){tag_s}",
                                ts.title,
                                ts.trees.len()
                            ))
                            .map_err(Error::from)?;
                    }
                }
            }

            Command::Erase(slot) => {
                if let Err(msg) = ensure_tree_context(&self.state) {
                    return Err(Error::runtime(msg, None, Some(span)));
                }
                self.state.tree_slots.remove(&slot);
                if self.state.current_tree_slot == Some(slot) {
                    self.state.current_tree_slot = None;
                }
                self.state
                    .write_line(&format!("erase removed tree set {slot}"))
                    .map_err(Error::from)?;
            }

            Command::TChoose(items) => {
                ensure_tree_context(&self.state)
                    .map_err(|m| Error::runtime(m, None, Some(span)))?;

                let n = self.state.working_tree_set_len();
                if n == 0 {
                    return Err(Error::runtime(
                        "tchoose no trees in working tree set 0",
                        None,
                        Some(span),
                    ));
                }

                let picked = expand_tree_selectors("tchoose", &items, n, span)?;

                /* Clone selected trees first to avoid conflicting state borrows. */
                let new_ts = {
                    let ts = self.state.working_tree_set().unwrap();
                    let mut new_trees = Vec::with_capacity(picked.len());
                    for &i in &picked {
                        new_trees.push(ts.trees[i].clone());
                    }
                    crate::engines::trees::TreeSet::new(ts.title.clone(), new_trees)
                };

                self.state.set_working_tree_set(new_ts);

                let show = picked.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" ");
                self.state
                    .write_line(&format!("tchoose kept {}/{} tree(s): {}", picked.len(), n, show))
                    .map_err(Error::from)?;
            }

            Command::TXAscii(enabled) => {
                self.state.tree_plot_style = if enabled {
                    crate::engines::tplot::TPlotStyle::Unicode
                } else {
                    crate::engines::tplot::TPlotStyle::Plain
                };

                self.state.write_line(&format!("txascii *")).map_err(Error::from)?;
            }

            Command::TPlot(trees) => {
                if let Err(msg) = ensure_tree_context(&self.state) {
                    return Err(Error::runtime(msg, None, Some(span)));
                }

                let rendered_all = {
                    let ds = self.state.dataset.as_ref().ok_or_else(|| {
                        Error::runtime("tplot: dataset not loaded", None, Some(span))
                    })?;

                    let ts = self.state.working_tree_set().ok_or_else(|| {
                        Error::runtime("tplot: no trees in working tree set 0", None, Some(span))
                    })?;

                    if ts.trees.is_empty() {
                        return Err(Error::runtime(
                            "tplot: working tree set is empty",
                            None,
                            Some(span),
                        ));
                    }

                    let outgroups = self.state.outgroup.as_ref().map(|og| og.taxa.as_slice());
                    let style = self.state.tree_plot_style;
                    let picked = if trees.is_empty() {
                        (0..ts.trees.len()).collect::<Vec<_>>()
                    } else {
                        expand_tree_selectors("tplot", &trees, ts.trees.len(), span)?
                    };

                    let rendered_all: std::result::Result<Vec<(usize, Vec<String>)>, String> =
                        picked
                            .into_iter()
                            .map(|i| {
                                crate::engines::tplot::render_tree(
                                    &ts.trees[i],
                                    &ds.taxa,
                                    outgroups,
                                    style,
                                )
                                .map(|lines| (i, lines))
                            })
                            .collect();

                    rendered_all.map_err(|m| {
                        Error::runtime(format!("tplot failed: {m}"), None, Some(span))
                    })?
                };

                for (i, lines) in rendered_all {
                    self.state.write_line(&format!("TREE {}", i)).map_err(Error::from)?;
                    for line in lines {
                        self.state.write_line(&line).map_err(Error::from)?;
                    }
                }
            }

            Command::Tsvg { tree, path } => {
                if let Err(msg) = ensure_tree_context(&self.state) {
                    return Err(Error::runtime(msg, None, Some(span)));
                }

                let text = {
                    let ds = self.state.dataset.as_ref().ok_or_else(|| {
                        Error::runtime("tsvg: dataset not loaded", None, Some(span))
                    })?;

                    let ts = self.state.working_tree_set().ok_or_else(|| {
                        Error::runtime("tsvg: no trees in working tree set 0", None, Some(span))
                    })?;

                    let picked = expand_tree_selectors("tsvg", &[tree], ts.trees.len(), span)?;
                    let i = picked[0];
                    let outgroups = self.state.outgroup.as_ref().map(|og| og.taxa.as_slice());
                    let style = self.state.tree_plot_style;

                    crate::engines::tplot::render_tree_svg(&ts.trees[i], &ds.taxa, outgroups, style)
                        .map_err(|m| {
                            Error::runtime(format!("tsvg failed: {m}"), None, Some(span))
                        })?
                };

                let mut f = OpenOptions::new()
                    .create(true)
                    .truncate(true)
                    .write(true)
                    .open(&path)
                    .map_err(|e| Error::io(e, Some(path.display().to_string())))?;

                write!(f, "{text}").map_err(|e| Error::io(e, Some(path.display().to_string())))?;

                self.state
                    .write_line(&format!("tsvg wrote {}", path.display()))
                    .map_err(Error::from)?;
            }

            Command::Apo { tree, path } => {
                if let Err(msg) = ensure_tree_context(&self.state) {
                    return Err(Error::runtime(msg, None, Some(span)));
                }

                let text = {
                    let ds = self.state.dataset.as_ref().ok_or_else(|| {
                        Error::runtime("apo: dataset not loaded", None, Some(span))
                    })?;

                    let cfg = self.state.char_config.as_ref().ok_or_else(|| {
                        Error::runtime("apo: ccode/char config not loaded", None, Some(span))
                    })?;

                    let ts = self.state.working_tree_set().ok_or_else(|| {
                        Error::runtime("apo: no trees in working tree set 0", None, Some(span))
                    })?;

                    let picked = expand_tree_selectors("apo", &[tree], ts.trees.len(), span)?;
                    let i = picked[0];
                    let annotations =
                        crate::engines::apo::detect_apomorphic_changes(ds, cfg, &ts.trees[i])
                            .map_err(|m| {
                                Error::runtime(format!("apo failed: {m}"), None, Some(span))
                            })?;

                    crate::engines::tplot::render_tree_apo_svg(&ts.trees[i], &ds.taxa, &annotations)
                        .map_err(|m| Error::runtime(format!("apo failed: {m}"), None, Some(span)))?
                };

                let mut f = OpenOptions::new()
                    .create(true)
                    .truncate(true)
                    .write(true)
                    .open(&path)
                    .map_err(|e| Error::io(e, Some(path.display().to_string())))?;

                write!(f, "{text}").map_err(|e| Error::io(e, Some(path.display().to_string())))?;

                self.state
                    .write_line(&format!("apo wrote {}", path.display()))
                    .map_err(Error::from)?;
            }

            Command::TList(trees) => {
                if let Err(msg) = ensure_tree_context(&self.state) {
                    return Err(Error::runtime(msg, None, Some(span)));
                }

                let dump = {
                    let Some(ts) = self.state.working_tree_set() else {
                        return Err(Error::runtime("tlist no current tree set", None, Some(span)));
                    };

                    if trees.is_empty() {
                        dump_tread_block_from_treeset(&ts.title, &ts.trees)
                    } else {
                        if ts.trees.is_empty() {
                            return Err(Error::runtime(
                                "tlist no trees in working tree set 0",
                                None,
                                Some(span),
                            ));
                        }
                        let picked = expand_tree_selectors("tlist", &trees, ts.trees.len(), span)?;
                        let selected =
                            picked.into_iter().map(|i| ts.trees[i].clone()).collect::<Vec<_>>();
                        dump_tread_block_from_treeset(&ts.title, &selected)
                    }
                };

                for line in dump.lines() {
                    self.state.write_line(line).map_err(Error::from)?;
                }
            }

            Command::TSave(path) => {
                if let Err(msg) = ensure_tree_context(&self.state) {
                    return Err(Error::runtime(msg, None, Some(span)));
                }
                let dump = {
                    let Some(ts) = self.state.working_tree_set() else {
                        return Err(Error::runtime("tsave no current tree set", None, Some(span)));
                    };
                    dump_tread_block_from_treeset(&ts.title, &ts.trees)
                };

                let mut f = OpenOptions::new()
                    .create(true)
                    .truncate(true)
                    .write(true)
                    .open(&path)
                    .map_err(|e| Error::io(e, Some(path.display().to_string())))?;

                write!(f, "{dump}").map_err(|e| Error::io(e, Some(path.display().to_string())))?;

                self.state
                    .write_line(&format!("tsave wrote {}", path.display()))
                    .map_err(Error::from)?;
            }

            Command::XSteps(modes) => {
                if let Err(msg) = ensure_tree_context(&self.state) {
                    return Err(Error::runtime(msg, None, Some(span)));
                }

                for mode in modes {
                    match mode {
                        crate::ast::XStepsMode::U => {
                            let old_ts =
                                self.state.working_tree_set().cloned().ok_or_else(|| {
                                    Error::runtime(
                                        "xsteps u no trees in working tree set 0",
                                        None,
                                        Some(span),
                                    )
                                })?;

                            if old_ts.trees.is_empty() {
                                return Err(Error::runtime(
                                    "xsteps u working tree set 0 is empty",
                                    None,
                                    Some(span),
                                ));
                            }

                            let outgroup = require_outgroup_index(&self.state)?;
                            let ds = self.state.dataset.as_ref().ok_or_else(|| {
                                Error::runtime("xsteps u dataset not loaded", None, Some(span))
                            })?;
                            let cfg = self.state.char_config.as_ref().ok_or_else(|| {
                                Error::runtime("xsteps u ccode config not loaded", None, Some(span))
                            })?;
                            let mut ws = crate::engines::search::ScoreWorkspace::new();

                            let mut seen = std::collections::HashSet::<u64>::new();
                            let mut uniq = Vec::<crate::engines::trees::Tree>::new();

                            for tr0 in old_ts.trees {
                                let mut tr = tr0.reroot_by_outgroup(outgroup).unwrap_or(tr0);
                                tr.canonicalize();
                                crate::engines::search::collapse_zero_length_branches(
                                    &mut tr, ds, cfg, &mut ws,
                                );
                                let h = crate::engines::search::collapsed_topology_hash(&tr);
                                if seen.insert(h) {
                                    uniq.push(tr);
                                }
                            }

                            let new_n = uniq.len();
                            self.state.set_working_tree_set(TreeSet::new(old_ts.title, uniq));

                            self.state
                                .write_line(&format!(
                                    "xsteps u working treeset 0 uniq => {} tree(s)",
                                    new_n
                                ))
                                .map_err(Error::from)?;
                        }

                        crate::ast::XStepsMode::L => {
                            let ds = self.require_dataset(span)?;
                            let cfg = self.require_char_config(span)?;
                            let trees = self.require_working_trees(span)?;

                            let rows = crate::engines::xsteps::analyze_lengths(ds, cfg, &trees);

                            self.state.write_line("Tree lengths").map_err(Error::from)?;
                            self.state.write_line(" ").map_err(Error::from)?;

                            const COLS: usize = 10;
                            const W: usize = 6;

                            let mut header = "".to_string();
                            for c in 0..COLS {
                                header.push_str(&format!("{:>W$}", c));
                            }
                            self.state
                                .write_line(&format!("{:>W$}{header}", ""))
                                .map_err(Error::from)?;
                            self.state
                                .write_line(&format!("{:>W$}{}", "", " ".repeat(COLS * W)))
                                .map_err(Error::from)?;

                            let mut start = 0usize;
                            while start < rows.len() {
                                let end = (start + COLS).min(rows.len());
                                let mut line = format!("{:>W$}", start);
                                for i in start..end {
                                    line.push_str(&format!("{:>W$}", rows[i].length));
                                }
                                self.state.write_line(&line).map_err(Error::from)?;
                                start = end;
                            }
                        }

                        crate::ast::XStepsMode::C => {
                            let ds = self.require_dataset(span)?;
                            let cfg = self.require_char_config(span)?;
                            let trees = self.require_working_trees(span)?;

                            let reports =
                                crate::engines::xsteps::analyze_character_fits(ds, cfg, &trees)
                                    .map_err(|m| {
                                        Error::runtime(
                                            format!("xsteps c failed: {m}"),
                                            None,
                                            Some(span),
                                        )
                                    })?;

                            self.state
                                .write_line(&format!("xsteps file 0 {} trees", reports.len()))
                                .map_err(Error::from)?;

                            for rep in reports {
                                let ri_s = rep
                                    .ri
                                    .map(|x| format!("{:.2}", x))
                                    .unwrap_or_else(|| "NA".to_string());

                                self.state
                                    .write_line(&format!(
                                        "tree {} length {} ci {:.2} ri {}",
                                        rep.tree_index, rep.length, rep.ci, ri_s
                                    ))
                                    .map_err(Error::from)?;

                                self.state
                                    .write_line("character/steps/ci/ri")
                                    .map_err(Error::from)?;

                                let nchars = rep.rows.len();
                                let mut chunk_start = 0usize;
                                while chunk_start < nchars {
                                    let chunk_end = (chunk_start + 10).min(nchars);

                                    let mut char_line = String::new();
                                    let mut step_line = String::new();
                                    let mut ci_line = String::new();
                                    let mut ri_line = String::new();
                                    for idx in chunk_start..chunk_end {
                                        let r = &rep.rows[idx];
                                        char_line.push_str(&format!("{:>4}", r.character));
                                        step_line.push_str(&format!("{:>4}", r.steps));
                                        ci_line.push_str(&format!("{:>4}", (r.ci * 100.0) as u64));
                                        let ri = r.ri.map(|x| (x * 100.0) as u64);
                                        ri_line.push_str(&match ri {
                                            Some(v) => format!("{:>4}", v),
                                            None => "   0".to_string(),
                                        });
                                    }
                                    self.state.write_line(&char_line).map_err(Error::from)?;
                                    self.state.write_line(&step_line).map_err(Error::from)?;
                                    self.state.write_line(&ci_line).map_err(Error::from)?;
                                    self.state.write_line(&ri_line).map_err(Error::from)?;
                                    self.state.write_line("").map_err(Error::from)?;
                                    chunk_start = chunk_end;
                                }
                            }
                        }

                        crate::ast::XStepsMode::H => {
                            let ds = self.require_dataset(span)?;
                            let cfg = self.require_char_config(span)?;
                            let trees = self.require_working_trees(span)?;

                            let reports =
                                crate::engines::xsteps::analyze_histories(ds, cfg, &trees)
                                    .map_err(|m| {
                                        Error::runtime(
                                            format!("xsteps h failed: {m}"),
                                            None,
                                            Some(span),
                                        )
                                    })?;

                            self.state
                                .write_line(&format!("xsteps file 0 {} trees", reports.len()))
                                .map_err(Error::from)?;

                            for tree_rep in reports {
                                self.state
                                    .write_line(&format!("tree {}", tree_rep.tree_index))
                                    .map_err(Error::from)?;

                                for hist in tree_rep.histories {
                                    self.state
                                        .write_line(&format!("character {}", hist.character))
                                        .map_err(Error::from)?;

                                    let node_ids = hist
                                        .internal_nodes
                                        .iter()
                                        .map(|x| x.to_string())
                                        .collect::<Vec<_>>()
                                        .join(" ");
                                    self.state.write_line(&node_ids).map_err(Error::from)?;

                                    let states = hist.node_states.join(" ");
                                    self.state.write_line(&states).map_err(Error::from)?;
                                }
                            }
                        }

                        crate::ast::XStepsMode::M => {
                            let ds = self.require_dataset(span)?;
                            let cfg = self.require_char_config(span)?;
                            let trees = self.require_working_trees(span)?;

                            let rows =
                                crate::engines::xsteps::analyze_best_worst_fits_across_trees(
                                    ds, cfg, &trees,
                                )
                                .map_err(|m| {
                                    Error::runtime(
                                        format!("xsteps m failed: {m}"),
                                        None,
                                        Some(span),
                                    )
                                })?;

                            self.state
                                .write_line("character best_steps worst_steps best_ci worst_ci best_ri worst_ri")
                                .map_err(Error::from)?;

                            for r in rows {
                                let best_ri_s = r
                                    .best_ri
                                    .map(|x| format!("{:.2}", x))
                                    .unwrap_or_else(|| "NA".to_string());
                                let worst_ri_s = r
                                    .worst_ri
                                    .map(|x| format!("{:.2}", x))
                                    .unwrap_or_else(|| "NA".to_string());

                                self.state
                                    .write_line(&format!(
                                        "{} {} {} {:.2} {:.2} {} {}",
                                        r.character,
                                        r.best_steps,
                                        r.worst_steps,
                                        r.best_ci,
                                        r.worst_ci,
                                        best_ri_s,
                                        worst_ri_s
                                    ))
                                    .map_err(Error::from)?;
                            }
                        }

                        crate::ast::XStepsMode::W => {
                            let (ds_nchar, last_tree) = {
                                let ds = self.state.dataset.as_ref().ok_or_else(|| {
                                    Error::runtime("xsteps w dataset not loaded", None, Some(span))
                                })?;

                                let ts = self.state.working_tree_set().ok_or_else(|| {
                                    Error::runtime(
                                        "xsteps w no trees in working tree set 0",
                                        None,
                                        Some(span),
                                    )
                                })?;

                                let tr = ts.trees.last().cloned().ok_or_else(|| {
                                    Error::runtime(
                                        "xsteps w working tree set 0 is empty",
                                        None,
                                        Some(span),
                                    )
                                })?;

                                (ds.nchar, tr)
                            };

                            let rows = {
                                let ds = self.state.dataset.as_ref().unwrap();
                                crate::engines::xsteps::successive_weights_rc10_from_tree(
                                    ds, &last_tree,
                                )
                                .map_err(|m| {
                                    Error::runtime(
                                        format!("xsteps w failed: {m}"),
                                        None,
                                        Some(span),
                                    )
                                })?
                            };

                            {
                                let cfg = self.state.char_config.as_mut().ok_or_else(|| {
                                    Error::runtime(
                                        "xsteps w char config not loaded",
                                        None,
                                        Some(span),
                                    )
                                })?;

                                crate::engines::xsteps::apply_successive_weights_to_ccode(
                                    cfg, ds_nchar, &rows,
                                )
                                .map_err(|m| {
                                    Error::runtime(
                                        format!("xsteps w writeback failed: {m}"),
                                        None,
                                        Some(span),
                                    )
                                })?;
                            }

                            self.state.write_line("xsteps w *").map_err(Error::from)?;
                        }
                    }
                }
            }

            Command::Steps => {
                /* steps: display min/max steps per character (dataset-only) */
                self.state.ensure_dataset().map_err(|m| Error::runtime(m, None, Some(span)))?;

                let ds = self.state.dataset.as_ref().unwrap();

                let rows = crate::engines::xsteps::character_min_max_steps(ds);

                self.state.write_line("min/max steps per character").map_err(Error::from)?;

                self.state.write_line("character").map_err(Error::from)?;
                let chars =
                    rows.iter().map(|r| r.character.to_string()).collect::<Vec<_>>().join(" ");
                self.state.write_line(&chars).map_err(Error::from)?;

                self.state.write_line("min").map_err(Error::from)?;
                let mins =
                    rows.iter().map(|r| r.min_steps.to_string()).collect::<Vec<_>>().join(" ");
                self.state.write_line(&mins).map_err(Error::from)?;

                self.state.write_line("max").map_err(Error::from)?;
                let maxs =
                    rows.iter().map(|r| r.max_steps.to_string()).collect::<Vec<_>>().join(" ");
                self.state.write_line(&maxs).map_err(Error::from)?;
            }

            Command::Watch(on) => {
                self.state.watch_enabled = on;
                self.state.write_line(&format!("watch *")).map_err(Error::from)?;
            }

            Command::XxQuery => {
                ensure_tree_context(&self.state)
                    .map_err(|m| Error::runtime(m, None, Some(span)))?;

                let ds = self
                    .state
                    .dataset
                    .as_ref()
                    .ok_or_else(|| Error::runtime("xx: dataset not loaded", None, Some(span)))?
                    .clone();

                let cfg = self
                    .state
                    .char_config
                    .as_ref()
                    .ok_or_else(|| Error::runtime("xx: char config not loaded", None, Some(span)))?
                    .clone();

                let ts = self.state.working_tree_set().cloned().ok_or_else(|| {
                    Error::runtime("xx: no trees in working tree set 0", None, Some(span))
                })?;

                let last_tree = ts.trees.last().cloned().ok_or_else(|| {
                    Error::runtime("xx: working tree set 0 is empty", None, Some(span))
                })?;

                let result = xx::run_xx_editor(last_tree, ds, cfg)
                    .map_err(|m| Error::runtime(format!("xx failed: {m}"), None, Some(span)))?;

                match result {
                    xx::XxEditResult::Cancelled => {
                        self.state.write_line("xx cancelled").map_err(Error::from)?;
                    }
                    xx::XxEditResult::Saved { tree, char_config } => {
                        /* Save the edited tree back into the last working tree. */
                        if let Some(ts0) = self.state.tree_slots.get_mut(&State::WORKING_TREE_SLOT)
                        {
                            if let Some(last) = ts0.trees.last_mut() {
                                *last = tree;
                            } else {
                                ts0.trees.push(tree);
                            }
                        } else {
                            self.state
                                .set_working_tree_set(TreeSet::new("xx edited tree", vec![tree]));
                        }

                        self.state.char_config = Some(char_config);
                        self.state.write_line("xx saved").map_err(Error::from)?;
                    }
                }
            }

            Command::Yama => {
                self.state.should_quit = true;
                self.state.write_line("yama").map_err(Error::from)?;
            }

            Command::Unknown { name, raw_args } => {
                self.state.write_line(ASSIST_CMD_LIST).map_err(Error::from)?;
                let tip = if raw_args.trim().is_empty() {
                    format!("{name} ? please use a; or a*; to get help")
                } else {
                    format!("{name} {raw_args} ? please use a; or a*; to get help")
                };
                self.state.write_line(&tip).map_err(Error::from)?;
            }
        }

        Ok(())
    }

    /// Helper: require dataset is loaded.
    fn require_dataset(&self, span: Span) -> Result<&crate::engines::dataset::Dataset> {
        self.state
            .dataset
            .as_ref()
            .ok_or_else(|| Error::runtime("xsteps dataset not loaded", None, Some(span)))
    }

    /// Helper: require char config is loaded.
    fn require_char_config(&self, span: Span) -> Result<&crate::engines::ccode::CharConfig> {
        self.state
            .char_config
            .as_ref()
            .ok_or_else(|| Error::runtime("xsteps ccode/char config not loaded", None, Some(span)))
    }

    /// Helper: clone and return the working tree set.
    fn require_working_trees(&self, span: Span) -> Result<Vec<crate::engines::trees::Tree>> {
        let ts = self.state.working_tree_set().ok_or_else(|| {
            Error::runtime("xsteps no trees in working tree set 0", None, Some(span))
        })?;
        if ts.trees.is_empty() {
            return Err(Error::runtime("xsteps working tree set is empty", None, Some(span)));
        }
        Ok(ts.trees.clone())
    }

    /// Helper: write a result line with tag, TL, CI, optional RI, and an extra suffix.
    fn write_result_line(
        &mut self,
        tag: &str,
        best_len: u64,
        ci: f64,
        ri: Option<f64>,
        extra: &str,
    ) -> Result<()> {
        match ri {
            Some(ri) => {
                self.state.write_line(&format!("{tag} TL={best_len} CI={ci:.2} RI={ri:.2}{extra}"))
            }
            None => self.state.write_line(&format!("{tag} TL={best_len} CI={ci:.2} RI=NA{extra}")),
        }
        .map_err(Error::from)
    }
}
/// detects whether a REPL input buffer contains a quit command that should end
/// interactive execution.

fn contains_yama_command(buf: &str) -> bool {
    /*
      This lightweight check is enough because commands are semicolon
      terminated; full parsing is handled by the normal command path.
    */
    let lower = buf.to_ascii_lowercase();
    lower.contains("yama") && lower.contains("yama;")
}
/// formats a dataset state set for human-readable matrix and diagnostic
/// output.

fn format_stateset(ss: crate::engines::dataset::StateSet) -> String {
    let bits = ss.bits();

    /* Unknown/all-36 state sets print as `?`. */
    if bits == crate::engines::dataset::StateSet::ALL36.bits() {
        return "?".to_string();
    }

    /* Singleton states print as a single symbol. */
    if bits.count_ones() == 1 {
        let idx = bits.trailing_zeros() as u8;
        return util::idx_to_symbol(idx).to_string();
    }

    /* Polymorphic states print as `[a b c]` with spaces between symbols. */
    let mut syms = Vec::new();
    for idx in 0..36u8 {
        if (bits >> idx) & 1 == 1 {
            syms.push(util::idx_to_symbol(idx));
        }
    }
    let inside = syms.into_iter().map(|c| c.to_string()).collect::<Vec<_>>().join(" ");
    format!("[{inside}]")
}
/// applies parsed character-coding operations to a character configuration
/// with range validation.

fn apply_ccode_ops(
    cfg: &mut crate::engines::ccode::CharConfig,
    ops: &[crate::ast::CCodeOp],
) -> std::result::Result<(), String> {
    for op in ops {
        match op {
            crate::ast::CCodeOp::Additive(sel) => {
                for idx in iter_char_sel(cfg.len(), sel)? {
                    cfg.chars[idx].additive = true;
                }
            }
            crate::ast::CCodeOp::NonAdditive(sel) => {
                for idx in iter_char_sel(cfg.len(), sel)? {
                    cfg.chars[idx].additive = false;
                }
            }
            crate::ast::CCodeOp::Active(sel) => {
                for idx in iter_char_sel(cfg.len(), sel)? {
                    cfg.chars[idx].active = true;
                }
            }
            crate::ast::CCodeOp::Inactive(sel) => {
                for idx in iter_char_sel(cfg.len(), sel)? {
                    cfg.chars[idx].active = false;
                }
            }
            crate::ast::CCodeOp::Weight(w, sel) => {
                if *w == 0 {
                    return Err("weight must be >= 1".to_string());
                }
                for idx in iter_char_sel(cfg.len(), sel)? {
                    cfg.chars[idx].weight = *w;
                }
            }
        }
    }
    Ok(())
}
/// verifies that dataset, coding, and current trees are available before tree
/// diagnostics or editing commands run.

fn ensure_tree_context(state: &crate::state::State) -> std::result::Result<(), String> {
    if state.dataset.is_none() {
        Err("dataset not loaded: tree operations require xread first".to_string())
    } else {
        Ok(())
    }
}
/// expands a parsed character selection into zero-based indices.

fn iter_char_sel(
    nchar: usize,
    sel: &crate::ast::CharSel,
) -> std::result::Result<Vec<usize>, String> {
    let mut out = Vec::new();
    match sel {
        crate::ast::CharSel::All => {
            out.extend(0..nchar);
        }
        crate::ast::CharSel::List(ranges) => {
            for r in ranges {
                match *r {
                    crate::ast::CharRange::Single(a) => {
                        let a = a as usize;
                        if a >= nchar {
                            return Err(format!("character {a} out of range (nchar={nchar})"));
                        }
                        out.push(a);
                    }
                    crate::ast::CharRange::RangeInclusive(a, b) => {
                        let a = a as usize;
                        let b = b as usize;
                        if a > b {
                            return Err(format!("invalid range {a}.{b}"));
                        }
                        if b >= nchar {
                            return Err(format!("range end {b} out of range (nchar={nchar})"));
                        }
                        out.extend(a..=b);
                    }
                }
            }
        }
    }
    Ok(out)
}
/// expands tree selectors into zero-based tree indices, preserving command order
/// and removing duplicates.

fn expand_tree_selectors(
    command: &str,
    items: &[crate::ast::TreeSelector],
    n: usize,
    span: Span,
) -> Result<Vec<usize>> {
    if n == 0 {
        return Err(Error::runtime(
            format!("{command} no trees in working tree set 0"),
            None,
            Some(span),
        ));
    }

    let mut picked: Vec<usize> = Vec::new();

    for it in items {
        match *it {
            crate::ast::TreeSelector::Last => picked.push(n - 1),
            crate::ast::TreeSelector::Index(i) => {
                if i >= n {
                    return Err(Error::runtime(
                        format!("{command}: index {i} out of range (n={n})"),
                        None,
                        Some(span),
                    ));
                }
                picked.push(i);
            }
            crate::ast::TreeSelector::RangeInclusive(a, b) => {
                if a > b {
                    return Err(Error::runtime(
                        format!("{command} invalid range (a > b)"),
                        None,
                        Some(span),
                    ));
                }
                if b >= n {
                    return Err(Error::runtime(
                        format!("{command} range end {b} out of range (n={n})"),
                        None,
                        Some(span),
                    ));
                }
                picked.extend(a..=b);
            }
        }
    }

    if picked.is_empty() {
        return Err(Error::runtime(format!("{command} no trees selected"), None, Some(span)));
    }

    let mut seen = std::collections::HashSet::<usize>::new();
    picked.retain(|i| seen.insert(*i));
    Ok(picked)
}
/// parses, validates, and normalizes all tread tree strings against the
/// current dataset.

fn normalize_and_validate_trees(
    ds: &crate::engines::dataset::Dataset,
    trees: &[String],
) -> std::result::Result<Vec<String>, String> {
    trees.iter().map(|t| normalize_and_validate_one_tree(ds, t)).collect()
}
/// tokenizes one parenthetical tree, checks taxa, and converts it into
/// internal storage form.

fn normalize_and_validate_one_tree(
    ds: &crate::engines::dataset::Dataset,
    tree: &str,
) -> std::result::Result<String, String> {
    let chars: Vec<char> = tree.chars().collect();
    let mut i = 0usize;
    let mut out_parts: Vec<String> = Vec::new();

    while i < chars.len() {
        let c = chars[i];

        if c.is_whitespace() {
            i += 1;
            continue;
        }

        if c == '(' || c == ')' || c == ',' {
            out_parts.push(c.to_string());
            i += 1;
            continue;
        }

        /* number token */
        if c.is_ascii_digit() {
            let start = i;
            i += 1;
            while i < chars.len() && chars[i].is_ascii_digit() {
                i += 1;
            }
            let s: String = chars[start..i].iter().collect();
            let idx: usize = s.parse().map_err(|_| format!("invalid taxon number: {s}"))?;
            if idx >= ds.ntax {
                return Err(format!("taxon number {idx} out of range (ntax={})", ds.ntax));
            }

            out_parts.push(idx.to_string());
            continue;
        }

        /* identifier token */
        if c.is_ascii_alphabetic() || c == '_' {
            let start = i;
            i += 1;
            while i < chars.len() && (chars[i].is_ascii_alphanumeric() || chars[i] == '_') {
                i += 1;
            }
            let name: String = chars[start..i].iter().collect();
            let idx = ds
                .taxa
                .iter()
                .position(|t| t == &name)
                .ok_or_else(|| format!("unknown taxon name in tree: {name}"))?;

            out_parts.push(idx.to_string());
            continue;
        }

        return Err(format!("invalid character in tree: {:?}", c));
    }

    let out = out_parts.join(" ");

    if !out.starts_with('(') || !out.ends_with(')') {
        return Err("tree must start with '(' and end with ')'".to_string());
    }

    validate_tree_parentheses(&out)?;
    Ok(out)
}
/// performs a lightweight balance check for raw parenthetical tree text.

fn validate_tree_parentheses(tree: &str) -> std::result::Result<(), String> {
    let mut bal = 0i32;
    for c in tree.chars() {
        match c {
            '(' => bal += 1,
            ')' => {
                bal -= 1;
                if bal < 0 {
                    return Err("unbalanced parentheses in tree".to_string());
                }
            }
            _ => {}
        }
    }
    if bal != 0 {
        return Err("unbalanced parentheses in tree".to_string());
    }
    Ok(())
}
/// serializes an internal tree set back into a `tread` block for file saving
/// and external validation.

fn dump_tread_block_from_treeset(title: &str, trees: &[crate::engines::trees::Tree]) -> String {
    let mut out = String::new();
    out.push_str(&format!("tread '{}'\n", title));

    for (i, t) in trees.iter().enumerate() {
        let pt = t.to_hennig_spine_string();

        if i + 1 == trees.len() {
            out.push_str(&format!("{pt} ;\n"));
        } else {
            out.push_str(&format!("{pt} *\n"));
        }
    }

    out.push_str(";\n");
    out
}
/// maps outgroup names or numeric requests into zero-based taxon indices using
/// the loaded dataset.

fn resolve_outgroup_request_from_taxa(
    requested: &[usize],
    ts: Option<&crate::engines::trees::TreeSet>,
) -> Vec<usize> {
    let mut taxa = requested.to_vec();
    taxa.sort_unstable();
    taxa.dedup();

    if taxa.is_empty() {
        return taxa;
    }

    if taxa.len() == 1 {
        return taxa;
    }

    if let Some(ts) = ts {
        let all_ok = ts.trees.iter().all(|tr| {
            tr.reroot_by_outgroup_set(&taxa).is_ok() && tr.find_monophyletic_node(&taxa).is_some()
        });
        if all_ok {
            return taxa;
        }
    }

    vec![taxa[0]]
}
/// reads Linux memory availability for the legacy `bytes;` command when
/// `/proc/meminfo` is available.

fn read_mem_available_bytes_linux() -> std::result::Result<Option<u64>, String> {
    let meminfo_path = std::path::Path::new("/proc/meminfo");
    if !meminfo_path.exists() {
        return Err("/proc/meminfo not found".to_string());
    }

    let content = std::fs::read_to_string(meminfo_path).map_err(|e| e.to_string())?;

    for line in content.lines() {
        if !line.starts_with("MemAvailable:") {
            continue;
        }

        let value_kb = line
            .split_whitespace()
            .nth(1)
            .ok_or_else(|| "invalid MemAvailable line".to_string())?
            .parse::<u64>()
            .map_err(|_| "invalid MemAvailable number".to_string())?;

        return Ok(Some(value_kb * 1024));
    }

    Err("MemAvailable not found in /proc/meminfo".to_string())
}
/// converts a byte count into compact human-readable text.

fn format_bytes_human(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];

    let mut value = bytes as f64;
    let mut unit_idx = 0usize;

    while value >= 1024.0 && unit_idx + 1 < UNITS.len() {
        value /= 1024.0;
        unit_idx += 1;
    }

    if unit_idx == 0 {
        format!("{} {}", bytes, UNITS[unit_idx])
    } else {
        format!("{:.1} {}", value, UNITS[unit_idx])
    }
}
/// returns the selected outgroup anchor or reports that no usable outgroup is
/// configured.

fn require_outgroup_index(state: &crate::state::State) -> crate::error::Result<usize> {
    state
        .outgroup
        .as_ref()
        .and_then(|og| og.first())
        .ok_or_else(|| crate::error::Error::runtime("outgroup not set", None, None))
}

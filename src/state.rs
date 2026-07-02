/*!Module: mutable VM state for display/logging, loaded datasets, character
coding slots, tree slots, outgroup selection, plot style, stopwatch state,
and process-control flags.

 */
use std::{
    collections::HashMap,
    fs::File,
    io::{self, Write},
};

use crate::engines::tplot::TPlotStyle;

/// Outgroup taxon selection: normalized, sorted, deduplicated list.
#[derive(Debug, Clone)]
pub struct OutgroupState {
    pub taxa: Vec<usize>,
}

impl OutgroupState {
    /// normalizes a user-supplied outgroup taxon list by sorting and deduplicating
    /// it.
    pub fn new(mut taxa: Vec<usize>) -> Self {
        taxa.sort_unstable();
        taxa.dedup();
        Self { taxa }
    }
    /// returns the lowest selected outgroup taxon for commands that need a single
    /// root anchor.

    pub fn first(&self) -> Option<usize> {
        self.taxa.first().copied()
    }
    /// formats selected outgroups as either numeric indices or index:name pairs
    /// when a dataset is loaded.

    pub fn describe_with_dataset(&self, ds: Option<&crate::engines::dataset::Dataset>) -> String {
        if self.taxa.is_empty() {
            return "none".to_string();
        }

        match ds {
            Some(ds) => {
                let parts = self
                    .taxa
                    .iter()
                    .map(|&t| {
                        let name = ds.taxa.get(t).cloned().unwrap_or_else(|| format!("?{t}"));
                        format!("{t}:{name}")
                    })
                    .collect::<Vec<_>>();
                parts.join(" ")
            }
            None => self.taxa.iter().map(|t| t.to_string()).collect::<Vec<_>>().join(" "),
        }
    }
}

/// Mutable VM state: datasets, coding, tree slots, display/logging, stopwatch.
#[derive(Debug)]
pub struct State {
    pub display_to_terminal: bool,
    pub log_file: Option<File>,
    pub dataset: Option<crate::engines::dataset::Dataset>,
    /* outgroup */
    pub outgroup: Option<OutgroupState>,

    /* character config slots */
    pub char_config: Option<crate::engines::ccode::CharConfig>,
    pub ccode_slots: std::collections::HashMap<u8, crate::engines::ccode::CharConfig>,

    /* trees: slot -> TreeSet(title + trees) */
    pub current_tree_slot: Option<u8>,
    pub tree_slots: HashMap<u8, crate::engines::trees::TreeSet>,

    pub tree_plot_style: TPlotStyle,

    pub log_path: Option<std::path::PathBuf>,
    pub log_enabled: bool,

    pub procedure_exec_enabled: bool,
    /* watch */
    pub watch_enabled: bool,
    pub watch_t0: Option<std::time::Instant>,

    pub should_quit: bool,
}

impl State {
    /// creates the default command-runtime state with no dataset, no tree files,
    /// logging disabled, and Unicode plotting enabled.
    pub fn new() -> Self {
        Self {
            display_to_terminal: false,
            log_path: None,
            log_file: None,
            log_enabled: false,
            procedure_exec_enabled: true,
            should_quit: false,
            dataset: None,
            outgroup: None,
            char_config: None,
            ccode_slots: HashMap::new(),
            current_tree_slot: None,
            tree_slots: HashMap::new(),
            tree_plot_style: TPlotStyle::Unicode,
            watch_enabled: false,
            watch_t0: None,
        }
    }
    /// writes one output line to the terminal and/or active log file according to
    /// the current display switches.

    pub fn write_line(&mut self, s: &str) -> io::Result<()> {
        if self.display_to_terminal {
            println!("{s}");
        }
        if self.log_enabled {
            if let Some(f) = &mut self.log_file {
                writeln!(f, "{s}")?;
            }
        }
        Ok(())
    }
    /// emits stopwatch timing information and advances the watch checkpoint when
    /// `watch;` is active.

    pub fn watch_note(&mut self, label: &str) -> io::Result<()> {
        if !self.watch_enabled {
            return Ok(());
        }
        let now = std::time::Instant::now();
        let dt = self.watch_t0.map(|t0| now.duration_since(t0));
        self.watch_t0 = Some(now);

        if let Some(dt) = dt {
            self.write_line(&format!("[watch] {label}: {:.3}s", dt.as_secs_f64()))?;
        } else {
            self.write_line(&format!("[watch] {label}: start"))?;
        }
        Ok(())
    }
    /// validates that character data has been loaded before commands that require
    /// a dataset run.

    pub fn ensure_dataset(&self) -> Result<(), String> {
        if self.dataset.is_some() {
            Ok(())
        } else {
            Err("dataset not loaded (use xread ... ; or procedure <file>; first)".to_string())
        }
    }

    /// Default tree slot index (0) used by search and tree-manipulation commands.
    pub const WORKING_TREE_SLOT: u8 = 0;
    /// returns the implicit working tree slot used by search and tree-editing
    /// commands.

    pub fn working_tree_set(&self) -> Option<&crate::engines::trees::TreeSet> {
        self.tree_slots.get(&Self::WORKING_TREE_SLOT)
    }
    /// appends newly produced trees to the working slot while preserving or
    /// initializing its title.

    pub fn append_to_working_tree_set(
        &mut self,
        title: String,
        mut trees: Vec<crate::engines::trees::Tree>,
    ) {
        let slot = Self::WORKING_TREE_SLOT;

        if let Some(ts) = self.tree_slots.get_mut(&slot) {
            ts.trees.append(&mut trees);
            if ts.title.trim().is_empty() {
                ts.title = title;
            }
        } else {
            self.tree_slots.insert(slot, crate::engines::trees::TreeSet::new(title, trees));
        }

        self.current_tree_slot = Some(slot);
        self.refresh_working_title_auto();
    }
    /// replaces the working tree slot and makes it current.

    pub fn set_working_tree_set(&mut self, ts: crate::engines::trees::TreeSet) {
        let slot = Self::WORKING_TREE_SLOT;
        self.tree_slots.insert(slot, ts);
        self.current_tree_slot = Some(slot);
    }
    /// updates auto-generated working titles so they match the current number of
    /// stored trees.

    pub fn refresh_working_title_auto(&mut self) {
        let slot = Self::WORKING_TREE_SLOT;
        let Some(ts) = self.tree_slots.get_mut(&slot) else {
            return;
        };

        let n = ts.trees.len();

        if ts.title.trim().to_ascii_lowercase().starts_with("set of ") {
            ts.title = format!("set of {} trees", n);
        }
    }
    /// returns the number of trees in the working slot without exposing the slot
    /// contents.

    pub fn working_tree_set_len(&self) -> usize {
        self.working_tree_set().map(|ts| ts.trees.len()).unwrap_or(0)
    }
    /// returns the tree set selected by `get`, search, or tree-loading commands.

    pub fn current_tree_set(&self) -> Option<&crate::engines::trees::TreeSet> {
        self.current_tree_slot.and_then(|slot| self.tree_slots.get(&slot))
    }
    /// removes every internal tree slot.

    pub fn clear_all_trees(&mut self) {
        self.current_tree_slot = None;
        self.tree_slots.clear();
    }
    /// appends trees to whichever internal slot is current, creating slot 0 when
    /// none exists.

    pub fn append_to_current_tree_set(
        &mut self,
        title: String,
        mut trees: Vec<crate::engines::trees::Tree>,
    ) {
        let slot = self.current_tree_slot.unwrap_or(0);

        if let Some(ts) = self.tree_slots.get_mut(&slot) {
            ts.trees.append(&mut trees);
        } else {
            self.tree_slots.insert(slot, crate::engines::trees::TreeSet::new(title, trees));
            self.current_tree_slot = Some(slot);
        }
    }
}

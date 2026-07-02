/*!Module: high-level engine facade that coordinates search, branch swapping,
implicit enumeration, and Nelson consensus operations for the VM.

 */
/// Fast TBR branch-swapping (=branch-breaking) search operators.
pub mod branchswap;
/// Character-coding configuration (additive, weight, active/inactive).
pub mod ccode;
/// Strict consensus (nelsen) tree construction.
pub mod con;
/// Character matrix parsing, storage, and state-set utilities.
pub mod dataset;
/// Implicit enumeration (exact branch-and-bound search).
pub mod ie;
/// Core parsimony search: Wagner trees, scoring, topology hashing, branch
/// collapse.
pub mod search;
/// Tree plotting and ASCII diagram generation.
pub mod tplot;
/// Tree data structures: Node, Tree, TreeSet, serialization, rerooting.
pub mod trees;
/// Shared utilities: state formatting, plot ordering, additive range math.
pub mod util;
/// Tree diagnostics: character fit, ancestral states, successive weighting.
pub mod xsteps;
/// Interactive tree editor for character-state inspection and topology
/// editing.
pub mod xx;

use crate::error::Result;

const TREE_BUFFER_LIMIT: usize = 100;
const MHENNIG_REPS_NONSTAR: usize = 500;
const MHENNIG_REPS_STAR: usize = 1000;

/// Stateless engine facade implementing legacy command dispatch for search,
/// consensus, and enumeration.
pub struct Engine;

impl Engine {
    /// constructs the stateless engine facade.
    pub fn new() -> Self {
        Self
    }
    /// validates dataset/config/outgroup inputs, runs single or multiple Wagner
    /// search, and optionally improves results with TBR swapping.

    pub fn hennig(
        &self,
        ds: &dataset::Dataset,
        cfg: &ccode::CharConfig,
        multi: bool,
        star: bool,
        outgroup: usize,
    ) -> Result<search::SearchOutcome> {
        if ds.ntax < 3 {
            return Err(crate::error::Error::runtime("hennig: ntax must be >= 3", None, None));
        }

        if cfg.len() != ds.nchar {
            return Err(crate::error::Error::runtime(
                format!("hennig: char_config nchar={} != dataset nchar={}", cfg.len(), ds.nchar),
                None,
                None,
            ));
        }

        if outgroup >= ds.ntax {
            return Err(crate::error::Error::runtime(
                format!("hennig: outgroup {outgroup} out of range for ntax={}", ds.ntax),
                None,
                None,
            ));
        }

        let initial = if !multi {
            search::hennig(ds, cfg, 0x1234_5678_9abc_def0, outgroup)
        } else {
            let reps = if star { MHENNIG_REPS_STAR } else { MHENNIG_REPS_NONSTAR };
            let keep_limit = if star { None } else { Some(TREE_BUFFER_LIMIT) };
            search::mhennig_limited(ds, cfg, reps, 0x3141_5926_5358_9793, outgroup, keep_limit)
        };

        if !star {
            return Ok(initial);
        }

        /*
          Star mode: run full branch-breaking closure (bb*) on the Wagner
          trees, including TBR pre-optimisation and collapse-aware dedup.
          This makes h* = h + bb* and mh* = mh + bb*.
        */
        let keep_limit = None; /* bb*: unlimited tree buffer */
        let (trees, best_len) =
            branchswap::branch_break_closure(ds, cfg, &initial.trees, outgroup, keep_limit)
                .map_err(|m| crate::error::Error::runtime(m, None, None))?;

        let min_len = search::minsteps_sum(ds, cfg);
        let max_len = search::maxsteps_sum(ds, cfg);
        let ci = search::calc_ci(min_len, best_len);
        let ri = search::calc_ri(min_len, max_len, best_len);

        Ok(search::SearchOutcome { trees, best_len, ci, ri })
    }
    /// runs branch-breaking closure from the current tree set using TBR and the
    /// requested tree-buffer policy.

    pub fn bb(
        &self,
        ds: &dataset::Dataset,
        cfg: &ccode::CharConfig,
        start_trees: &[crate::engines::trees::Tree],
        outgroup: usize,
        star: bool,
    ) -> Result<(Vec<crate::engines::trees::Tree>, u64)> {
        let keep_limit = if star { None } else { Some(TREE_BUFFER_LIMIT) };

        let (trees, best_len) =
            branchswap::branch_break_closure(ds, cfg, start_trees, outgroup, keep_limit)
                .map_err(|m| crate::error::Error::runtime(m, None, None))?;

        Ok((trees, best_len))
    }
    /// delegates exact implicit enumeration and maps string failures into runtime
    /// errors.

    pub fn ie(
        &self,
        ds: &dataset::Dataset,
        cfg: &ccode::CharConfig,
        outgroup: usize,
        mode: ie::IeMode,
    ) -> Result<ie::IeOutcome> {
        let out = ie::ie_exact(ds, cfg, outgroup, mode)
            .map_err(|m| crate::error::Error::runtime(m, None, None))?;
        Ok(out)
    }
    /// delegates strict Nelson consensus construction and scoring for the current
    /// tree set.

    pub fn nelsen_consensus(
        &self,
        ds: &dataset::Dataset,
        cfg: &ccode::CharConfig,
        trees: &[crate::engines::trees::Tree],
        outgroups: &[usize],
    ) -> Result<con::ConsensusOutcome> {
        crate::engines::con::nelsen_consensus(ds, cfg, trees, outgroups)
            .map_err(|m| crate::error::Error::runtime(m, None, None))
    }
}

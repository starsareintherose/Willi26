/*!Module: high-level engine facade that coordinates search, branch swapping,
implicit enumeration, and Nelson consensus operations for the VM.

 */
/// Apomorphy/homoplasy mapping for annotated SVG output.
pub mod apo;
/// Fast TBR branch-swapping (=branch-breaking) search operators.
pub mod branchswap;
/// Character-coding configuration (additive, weight, active/inactive).
pub mod ccode;
/// Strict consensus (nelsen) tree construction.
pub mod consensus;
/// Character matrix parsing, storage, and state-set utilities.
pub mod dataset;
/// Implicit enumeration (exact branch-and-bound search).
pub mod ie;
/// Bootstrap/jackknife-style resampling support calculations.
pub mod resample;
/// Deterministic random-number mixing.
pub(crate) mod rng;
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
/* Number of random-addition replicates used by multiple Hennig searches. */
const MHENNIG_REPS: usize = 10;

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

        /* Interactive Hennig searches use the reproducible seed sequence rooted
        at 42. Resampling paths may still pass their own explicit seeds. */
        let seed = search::entropy_seed();
        if !star {
            return Ok(if !multi {
                search::hennig(ds, cfg, seed, outgroup)
            } else {
                search::mhennig_limited(
                    ds,
                    cfg,
                    MHENNIG_REPS,
                    seed,
                    outgroup,
                    Some(TREE_BUFFER_LIMIT),
                )
            });
        }

        /* Star Hennig swaps every random Wagner start independently, retaining
        one local optimum from each start before selecting the shortest
        results. */
        let reps = if multi { MHENNIG_REPS } else { 1 };
        let starts: Vec<crate::engines::trees::Tree> = (0..reps)
            .map(|rep| {
                let rep_seed = if rep == 0 { seed } else { search::entropy_seed() };
                search::hennig(ds, cfg, rep_seed, outgroup).trees.remove(0)
            })
            .collect();

        /* Independent starts can be swapped concurrently without cloning the
        immutable dataset or character configuration. */
        let swapped: std::result::Result<Vec<_>, String> = std::thread::scope(|scope| {
            let handles: Vec<_> = starts
                .iter()
                .map(|start| {
                    scope.spawn(move || branchswap::branch_swap_one(ds, cfg, start, outgroup))
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().map_err(|_| "mhennig*: TBR worker panicked".to_string())?)
                .collect()
        });
        let swapped = swapped.map_err(|m| crate::error::Error::runtime(m, None, None))?;

        let mut trees = Vec::new();
        let mut seen = std::collections::BTreeSet::new();
        let mut best_len = u64::MAX;
        for (mut tree, len) in swapped {
            tree.canonicalize();
            if len < best_len {
                best_len = len;
                trees.clear();
                seen.clear();
            }
            if len == best_len {
                let key = tree.to_storage_string();
                if seen.insert(key) {
                    trees.push(tree);
                }
            }
        }

        Ok(search::search_outcome(ds, cfg, trees, best_len))
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
    ) -> Result<consensus::ConsensusOutcome> {
        crate::engines::consensus::nelsen_consensus(ds, cfg, trees, outgroups)
            .map_err(|m| crate::error::Error::runtime(m, None, None))
    }
}

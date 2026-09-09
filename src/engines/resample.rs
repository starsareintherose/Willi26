/*!Module: bootstrap resampling support for target-tree clades.

The implementation estimates support for clades already present in a target tree.
Each pseudoreplicate is searched, summarized with a strict consensus, and then
matched back to target-tree internal nodes by descendant taxon set.
 */
use std::collections::{BTreeSet, HashMap};

use crate::ast::ResampleSearchStep;
use crate::engines::{
    branchswap, ccode::CharConfig, consensus, dataset::Dataset, ie, search, trees::Tree,
};

const BOOT_SEARCH_REPS: usize = 16;
const BOOT_TREE_LIMIT: usize = 100;

#[derive(Debug, Clone)]
pub struct BootstrapSupportResult {
    pub target_tree: Tree,
    pub labels: HashMap<usize, String>,
}

#[derive(Debug, Clone, Copy)]
pub enum ResampleMethod {
    Bootstrap,
    Jackknife { delete_percent: f64 },
    Symmetric { delete_percent: f64 },
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct TaxonSet {
    bits: Vec<u64>,
}

impl TaxonSet {
    fn new(ntax: usize) -> Self {
        Self { bits: vec![0; ntax.div_ceil(64)] }
    }

    fn singleton(ntax: usize, taxon: usize) -> Self {
        let mut out = Self::new(ntax);
        out.insert(taxon);
        out
    }

    fn insert(&mut self, taxon: usize) {
        self.bits[taxon / 64] |= 1u64 << (taxon % 64);
    }

    fn union_inplace(&mut self, other: &Self) {
        for (a, b) in self.bits.iter_mut().zip(&other.bits) {
            *a |= *b;
        }
    }

    fn contains(&self, taxon: usize) -> bool {
        ((self.bits[taxon / 64] >> (taxon % 64)) & 1) == 1
    }

    fn count(&self) -> usize {
        self.bits.iter().map(|x| x.count_ones() as usize).sum()
    }
}

pub fn bootstrap_support(
    ds: &Dataset,
    cfg: &CharConfig,
    target_tree: &Tree,
    outgroups: &[usize],
    replications: usize,
    method: ResampleMethod,
    search_steps: &[ResampleSearchStep],
) -> Result<BootstrapSupportResult, String> {
    if replications == 0 {
        return Err("resample boot requires replications > 0".to_string());
    }
    if ds.nchar == 0 {
        return Err("resample boot requires at least one character".to_string());
    }
    if target_tree.ntax != ds.ntax {
        return Err("resample boot target tree ntax does not match dataset".to_string());
    }
    if cfg.len() != ds.nchar {
        return Err(format!(
            "resample boot ccode nchar={} != dataset nchar={}",
            cfg.len(),
            ds.nchar
        ));
    }

    let target_clades = target_node_clades(target_tree, outgroups);
    let mut counts = vec![0usize; target_clades.len()];
    let outgroup = outgroups.first().copied().unwrap_or(0);
    let search_steps = if search_steps.is_empty() {
        &[ResampleSearchStep::MHennig { star: false }, ResampleSearchStep::Bb { star: false }][..]
    } else {
        search_steps
    };

    for rep in 0..replications {
        let rep_cfg = resampled_config(cfg, ds.nchar, method, 0xa11c_e5eed ^ mix64(rep as u64));
        let searched = run_resample_search(
            ds,
            &rep_cfg,
            outgroup,
            search_steps,
            0x5eed_0000_0000_0001 ^ mix64(rep as u64),
        )?;
        let searched = keep_replicate_optimal_trees(ds, &rep_cfg, searched);

        if searched.is_empty() {
            continue;
        }

        let consensus_tree = consensus::strict_consensus_with_outgroups(&searched, outgroups)
            .map_err(|m| format!("resample boot replicate {} consensus failed: {m}", rep + 1))?;
        let rep_clades = tree_clade_set(&consensus_tree, outgroups);

        for (i, (_, clade)) in target_clades.iter().enumerate() {
            if rep_clades.contains(clade) {
                counts[i] += 1;
            }
        }
    }

    let mut labels = HashMap::new();
    for ((node, _), count) in target_clades.iter().zip(counts) {
        let pct = (count * 100 + replications / 2) / replications;
        labels.insert(*node, pct.to_string());
    }

    Ok(BootstrapSupportResult { target_tree: target_tree.clone(), labels })
}

fn resampled_config(
    cfg: &CharConfig,
    nchar: usize,
    method: ResampleMethod,
    seed: u64,
) -> CharConfig {
    match method {
        ResampleMethod::Bootstrap => bootstrap_config(cfg, nchar, seed),
        ResampleMethod::Jackknife { delete_percent } => {
            jackknife_config(cfg, nchar, seed, delete_percent)
        }
        ResampleMethod::Symmetric { delete_percent } => {
            symmetric_config(cfg, nchar, seed, delete_percent)
        }
    }
}

fn bootstrap_config(cfg: &CharConfig, nchar: usize, seed: u64) -> CharConfig {
    let mut rng = Rng::new(seed);
    let mut counts = vec![0u32; nchar];
    let active: Vec<usize> = cfg
        .chars
        .iter()
        .enumerate()
        .filter_map(|(i, setting)| setting.active.then_some(i))
        .collect();
    for _ in 0..active.len() {
        let c = active[rng.next_usize(active.len())];
        counts[c] = counts[c].saturating_add(1);
    }

    let mut out = cfg.clone();
    for (i, setting) in out.chars.iter_mut().enumerate() {
        if !cfg.chars[i].active || counts[i] == 0 {
            setting.active = false;
            setting.weight = 1;
        } else {
            setting.active = true;
            setting.weight = cfg.chars[i].weight.saturating_mul(counts[i]);
        }
    }
    out
}

/// A replicate contributes the strict consensus of its MPTs, not near-optimal
/// search survivors. `mhennig_limited` intentionally retains a slack buffer,
/// so normalize its output before consensus support is counted.
fn keep_replicate_optimal_trees(ds: &Dataset, cfg: &CharConfig, trees: Vec<Tree>) -> Vec<Tree> {
    let mut ws = search::ScoreWorkspace::new();
    let mut best = u64::MAX;
    let mut out = Vec::new();

    for tree in trees {
        let len = search::score_tree(ds, cfg, &tree, &mut ws);
        if len < best {
            best = len;
            out.clear();
            out.push(tree);
        } else if len == best {
            out.push(tree);
        }
    }

    out
}

fn symmetric_config(cfg: &CharConfig, nchar: usize, seed: u64, delete_percent: f64) -> CharConfig {
    let mut rng = Rng::new(seed);
    let mut out = cfg.clone();
    for i in 0..nchar {
        if !cfg.chars[i].active {
            out.chars[i].active = false;
            out.chars[i].weight = 1;
            continue;
        }

        let roll = rng.next_percent();
        if roll < delete_percent {
            out.chars[i].active = false;
            out.chars[i].weight = 1;
        } else if roll < delete_percent * 2.0 {
            out.chars[i] = cfg.chars[i];
            out.chars[i].weight = cfg.chars[i].weight.saturating_mul(2);
        } else {
            out.chars[i] = cfg.chars[i];
        }
    }
    out
}

fn jackknife_config(cfg: &CharConfig, nchar: usize, seed: u64, delete_percent: f64) -> CharConfig {
    let mut rng = Rng::new(seed);
    let mut out = cfg.clone();
    for i in 0..nchar {
        if !cfg.chars[i].active || rng.next_percent() < delete_percent {
            out.chars[i].active = false;
            out.chars[i].weight = 1;
        } else {
            out.chars[i] = cfg.chars[i];
        }
    }
    out
}

fn run_resample_search(
    ds: &Dataset,
    cfg: &CharConfig,
    outgroup: usize,
    steps: &[ResampleSearchStep],
    seed: u64,
) -> Result<Vec<Tree>, String> {
    let mut trees = Vec::<Tree>::new();

    for (idx, step) in steps.iter().enumerate() {
        match step {
            ResampleSearchStep::Hennig { star } => {
                trees = search::hennig(ds, cfg, seed ^ mix64(idx as u64), outgroup).trees;
                if *star {
                    trees = branchswap::branch_break_closure(ds, cfg, &trees, outgroup, None)
                        .map(|(trees, _)| trees)
                        .unwrap_or(trees);
                }
            }
            ResampleSearchStep::MHennig { star } => {
                trees = search::mhennig_limited(
                    ds,
                    cfg,
                    BOOT_SEARCH_REPS,
                    seed ^ mix64(idx as u64),
                    outgroup,
                    if *star { None } else { Some(BOOT_TREE_LIMIT) },
                )
                .trees;
                if *star {
                    trees = branchswap::branch_break_closure(ds, cfg, &trees, outgroup, None)
                        .map(|(trees, _)| trees)
                        .unwrap_or(trees);
                }
            }
            ResampleSearchStep::Bb { star } => {
                if trees.is_empty() {
                    trees = search::mhennig_limited(
                        ds,
                        cfg,
                        BOOT_SEARCH_REPS,
                        seed ^ mix64(idx as u64),
                        outgroup,
                        Some(BOOT_TREE_LIMIT),
                    )
                    .trees;
                }
                trees = branchswap::branch_break_closure(
                    ds,
                    cfg,
                    &trees,
                    outgroup,
                    if *star { None } else { Some(BOOT_TREE_LIMIT) },
                )
                .map(|(trees, _)| trees)
                .unwrap_or(trees);
            }
            ResampleSearchStep::Ie { star, dash } => {
                let mode = if *star {
                    ie::IeMode::KeepAll
                } else if *dash {
                    ie::IeMode::KeepOne
                } else {
                    ie::IeMode::KeepUpTo(100)
                };
                trees = ie::ie_exact(ds, cfg, outgroup, mode)?.trees;
            }
        }
    }

    Ok(trees)
}

fn target_node_clades(tree: &Tree, outgroups: &[usize]) -> Vec<(usize, TaxonSet)> {
    let mut out = Vec::new();
    collect_node_clades_rec(tree, tree.root, outgroups, &mut out);
    out
}

fn tree_clade_set(tree: &Tree, outgroups: &[usize]) -> BTreeSet<TaxonSet> {
    target_node_clades(tree, outgroups).into_iter().map(|(_, c)| c).collect()
}

fn collect_node_clades_rec(
    tree: &Tree,
    node: usize,
    outgroups: &[usize],
    out: &mut Vec<(usize, TaxonSet)>,
) -> TaxonSet {
    let n = &tree.nodes[node];
    if let Some(taxon) = n.taxon {
        return TaxonSet::singleton(tree.ntax, taxon);
    }

    let mut acc = TaxonSet::new(tree.ntax);
    for &child in &n.children {
        let child_set = collect_node_clades_rec(tree, child, outgroups, out);
        acc.union_inplace(&child_set);
    }

    let count = acc.count();
    let contains_outgroup = outgroups.iter().any(|&og| og < tree.ntax && acc.contains(og));

    if node != tree.root && count > 1 && count < tree.ntax && !contains_outgroup {
        out.push((node, acc.clone()));
    }

    acc
}

struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Self { state: seed | 1 }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        mix64(self.state)
    }

    fn next_usize(&mut self, n: usize) -> usize {
        if n <= 1 { 0 } else { (self.next_u64() as usize) % n }
    }

    fn next_percent(&mut self) -> f64 {
        (self.next_u64() as f64 / u64::MAX as f64) * 100.0
    }
}

fn mix64(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

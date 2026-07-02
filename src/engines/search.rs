/*!Module: core parsimony search and scoring. It contains Wagner/Hennig tree
construction, repeated multi-start search, compressed character-pattern
scoring, additive ordered-state scoring, incremental candidate evaluation,
topology hashing, CI/RI summaries, and zero-length branch collapse.

 */
use std::collections::{BTreeSet, HashMap};

use crate::engines::{ccode::CharConfig, dataset::Dataset, trees::Tree, util};

/*Shared constants */
/// Maximum number of ordered character states (0..=35).
pub const NSTATES: usize = 36;

/// Sentinel infinity value for additive DP cost arrays.
pub const INF: u64 = u64::MAX / 8;

/* Fitch merge */
/// Non-additive (Fitch) merge of two state sets.
///
/// Returns `(merged_bits, step_cost)`. If the intersection is non-empty it is
/// kept with zero cost; otherwise the union is taken with cost 1.
#[inline(always)]
pub fn merge_bits(a: u64, b: u64) -> (u64, u64) {
    let inter = a & b;
    if inter != 0 { (inter, 0) } else { (a | b, 1) }
}

/// Result of a Wagner or branch-breaking search: best length, tree set, and
/// statistics.
#[derive(Debug, Clone)]
pub struct SearchOutcome {
    pub trees: Vec<Tree>,
    pub best_len: u64,
    pub ci: f64,
    pub ri: Option<f64>,
}

#[derive(Debug, Clone)]
struct CharGroup {
    chars: Vec<usize>,
    total_weight: u64,
    additive: bool,
    all_contiguous: bool,
}
/// compresses active characters with identical taxon-state vectors and scoring
/// mode into weighted groups.

fn build_char_groups(ds: &Dataset, cfg: &CharConfig) -> Vec<CharGroup> {
    type PatternKey = (Vec<u64>, bool);
    let nchar = ds.nchar;
    let ntax = ds.ntax;

    let mut groups: HashMap<PatternKey, CharGroup> = HashMap::new();

    for c in 0..nchar {
        let cs = &cfg.chars[c];
        if !cs.active {
            continue;
        }

        let bits_vec: Vec<u64> = (0..ntax).map(|t| ds.matrix[t * nchar + c].bits()).collect();

        let additive = cs.additive;
        let all_contiguous =
            if additive { bits_vec.iter().all(|&b| util::is_contiguous(b)) } else { false };
        let key = (bits_vec, additive);
        groups
            .entry(key)
            .and_modify(|g: &mut CharGroup| {
                g.chars.push(c);
                g.total_weight += cs.weight as u64;
            })
            .or_insert_with(|| CharGroup {
                chars: vec![c],
                total_weight: cs.weight as u64,
                additive,
                all_contiguous,
            });
    }

    let mut result: Vec<CharGroup> = groups.into_values().collect();

    result.sort_by(|a, b| b.total_weight.cmp(&a.total_weight));

    result
}

/// Compute a dedup hash that collapses zero-length branches first,
/// then hashes the collapsed topology via clade sets.
///
/// This treats different binary resolutions of the same collapsed tree as
/// identical, matching IE's behaviour for all search and branch-breaking
/// commands (h, h*, mh, mh*, bb, bb*).
pub fn collapsed_dedup_hash(
    tree: &Tree,
    ds: &Dataset,
    cfg: &CharConfig,
    ws: &mut ScoreWorkspace,
) -> u64 {
    ws.scratch_nodes.clear();
    ws.scratch_nodes.extend_from_slice(&tree.nodes);
    ws.scratch_root = tree.root;
    ws.scratch_ntax = tree.ntax;
    let mut scratch = Tree {
        nodes: std::mem::take(&mut ws.scratch_nodes),
        root: ws.scratch_root,
        ntax: ws.scratch_ntax,
    };
    collapse_zero_length_branches(&mut scratch, ds, cfg, ws);
    let h = collapsed_topology_hash(&scratch);
    ws.scratch_nodes = scratch.nodes;
    ws.scratch_root = scratch.root;
    ws.scratch_ntax = scratch.ntax;
    h
}

#[derive(Debug, Clone)]
pub(crate) struct ScoreWorkspace {
    post: Vec<usize>,
    st: Vec<u64>,
    rg: Vec<util::AdditiveRange>,
    char_steps: Vec<u64>,

    /* cache leaf states for current dataset shape */
    leaf_cache_ntax: usize,
    leaf_cache_nchar: usize,
    leaf_bits: Vec<u64>,
    leaf_rg: Vec<util::AdditiveRange>,

    /* character groups for pattern-based acceleration */
    char_groups: Vec<CharGroup>,
    groups_initialized: bool,

    /* reusable Sankoff cost array (nnode * NSTATES) */
    sankoff_costs: Vec<u64>,

    /* reusable buffers for compute_rule1_collapsible_edges */
    coll_down: Vec<u64>,
    coll_child_transition: Vec<u64>,
    coll_up: Vec<u64>,
    coll_post: Vec<usize>,

    /* scratch tree for clone-free hashing */
    pub(crate) scratch_nodes: Vec<crate::engines::trees::Node>,
    pub(crate) scratch_root: usize,
    pub(crate) scratch_ntax: usize,
}

impl ScoreWorkspace {
    /// allocates reusable scoring buffers and caches.
    pub(crate) fn new() -> Self {
        Self {
            post: Vec::new(),
            st: Vec::new(),
            rg: Vec::new(),
            char_steps: Vec::new(),

            leaf_cache_ntax: 0,
            leaf_cache_nchar: 0,
            leaf_bits: Vec::new(),
            leaf_rg: Vec::new(),

            char_groups: Vec::new(),
            groups_initialized: false,

            sankoff_costs: Vec::new(),

            coll_down: Vec::new(),
            coll_child_transition: Vec::new(),
            coll_up: Vec::new(),
            coll_post: Vec::new(),

            scratch_nodes: Vec::new(),
            scratch_root: 0,
            scratch_ntax: 0,
        }
    }
    /// lazily initializes compressed character groups for the current
    /// dataset/configuration.

    fn ensure_groups(&mut self, ds: &Dataset, cfg: &CharConfig) {
        if !self.groups_initialized {
            self.char_groups = build_char_groups(ds, cfg);
            self.groups_initialized = true;
        }
    }
    /// refreshes postorder traversal and node buffers for a tree shape.

    fn prepare(&mut self, tree: &Tree, nchar: usize) {
        let n = tree.nodes.len();

        self.post.clear();
        self.post.reserve(n);
        postorder_fill(tree, tree.root, &mut self.post);

        if self.st.len() < n {
            self.st.resize(n, 0);
        }
        if self.rg.len() < n {
            self.rg.resize(n, util::AdditiveRange::default());
        }
        if self.char_steps.len() < nchar {
            self.char_steps.resize(nchar, 0);
        }
        let needed_sankoff = n * NSTATES;
        if self.sankoff_costs.len() < needed_sankoff {
            self.sankoff_costs.resize(needed_sankoff, INF);
        }
    }
    /// refreshes cached taxon-state bitsets and ranges for a dataset shape.

    fn ensure_leaf_cache(&mut self, ds: &Dataset) {
        if self.leaf_cache_ntax == ds.ntax
            && self.leaf_cache_nchar == ds.nchar
            && self.leaf_bits.len() == ds.ntax * ds.nchar
            && self.leaf_rg.len() == ds.ntax * ds.nchar
        {
            return;
        }

        let n = ds.ntax * ds.nchar;
        self.leaf_cache_ntax = ds.ntax;
        self.leaf_cache_nchar = ds.nchar;
        self.leaf_bits.resize(n, 0);
        self.leaf_rg.resize(n, util::AdditiveRange::default());

        for t in 0..ds.ntax {
            for c in 0..ds.nchar {
                let idx = t * ds.nchar + c;
                let bits = ds.matrix[idx].bits();
                self.leaf_bits[idx] = bits;
                self.leaf_rg[idx] = util::AdditiveRange::from_bits(bits);
            }
        }
    }
}

/* public API */
/// builds one greedy Wagner tree and returns its length/statistics.

pub fn hennig(ds: &Dataset, cfg: &CharConfig, seed: u64, outgroup: usize) -> SearchOutcome {
    let mut tree = hennig_once_tree(ds, cfg, seed, outgroup);
    let mut ws = ScoreWorkspace::new();

    tree.canonicalize();

    let best_len = score_tree(ds, cfg, &tree, &mut ws);
    let min_len = minsteps_sum(ds, cfg);
    let max_len = maxsteps_sum(ds, cfg);
    let ci = calc_ci(min_len, best_len);
    let ri = calc_ri(min_len, max_len, best_len);

    SearchOutcome { trees: vec![tree], best_len, ci, ri }
}
/// runs repeated Wagner starts with optional best-tree buffer limiting for
/// legacy non-star behavior.

pub fn mhennig_limited(
    ds: &Dataset,
    cfg: &CharConfig,
    reps: usize,
    seed0: u64,
    outgroup: usize,
    keep_limit: Option<usize>,
) -> SearchOutcome {
    let mut best_len = u64::MAX;
    let mut best_keys: BTreeSet<u64> = BTreeSet::new();
    let mut best_trees: Vec<Tree> = Vec::new();
    let mut ws = ScoreWorkspace::new();

    const MHENNIG_SLACK: u64 = 2;

    for i in 0..reps {
        let seed = seed0 ^ mix64(i as u64);
        let mut tr = hennig_once_tree(ds, cfg, seed, outgroup);
        let len = score_tree(ds, cfg, &tr, &mut ws);

        tr.canonicalize();
        let key = collapsed_dedup_hash(&tr, ds, cfg, &mut ws);

        if len < best_len {
            best_len = len;
            best_keys.clear();
            best_trees.clear();
            best_keys.insert(key);
            best_trees.push(tr);
        } else if len <= best_len + MHENNIG_SLACK {
            if best_keys.insert(key) {
                let can_keep = keep_limit.map(|lim| best_trees.len() < lim).unwrap_or(true);
                if can_keep {
                    best_trees.push(tr);
                }
            }
        }
    }

    let min_len = minsteps_sum(ds, cfg);
    let max_len = maxsteps_sum(ds, cfg);
    let ci = calc_ci(min_len, best_len);
    let ri = calc_ri(min_len, max_len, best_len);

    SearchOutcome { trees: best_trees, best_len, ci, ri }
}

#[derive(Clone, Copy, Debug, Default)]
struct LocalEval {
    st: u64,
    rg: util::AdditiveRange,
    step: u64,
}

#[derive(Debug, Clone)]
struct WagnerCharCache {
    nnode: usize,
    nchar: usize,

    st: Vec<u64>,
    rg: Vec<util::AdditiveRange>,

    node_steps_nonadd: Vec<u64>,
    node_steps_add: Vec<u64>,

    char_steps: Vec<u64>,
    total: u64,

    subtree_min_taxon: Vec<usize>,

    leaf_bits: Vec<u64>,
    leaf_rg: Vec<util::AdditiveRange>,

    char_groups: Vec<CharGroup>,
    groups_initialized: bool,
}

impl WagnerCharCache {
    /// lazily initializes grouped patterns for insertion scoring.
    fn ensure_groups(&mut self, ds: &Dataset, cfg: &CharConfig) {
        if !self.groups_initialized {
            self.char_groups = build_char_groups(ds, cfg);
            self.groups_initialized = true;
        }
    }
    /// builds leaf caches and node caches for one partial Wagner tree.

    fn new(ds: &Dataset, cfg: &CharConfig, tree: &Tree) -> Self {
        let nnode = tree.nodes.len().max(ds.ntax.saturating_mul(2).saturating_add(8));

        let nchar = ds.nchar;

        let mut out = Self {
            nnode,
            nchar,
            st: vec![0; nnode * nchar],
            rg: vec![util::AdditiveRange::default(); nnode * nchar],
            node_steps_nonadd: vec![0; nnode * nchar],
            node_steps_add: vec![0; nnode * nchar],
            char_steps: vec![0; nchar],
            total: 0,
            subtree_min_taxon: vec![usize::MAX; nnode],
            leaf_bits: vec![0; ds.ntax * nchar],
            leaf_rg: vec![util::AdditiveRange::default(); ds.ntax * nchar],
            char_groups: Vec::new(),
            groups_initialized: false,
        };

        out.rebuild_leaf_cache(ds);
        out.rebuild(ds, cfg, tree);
        out
    }
    /// maps taxon/character coordinates into leaf-cache indices.

    fn leaf_idx(&self, taxon: usize, c: usize) -> usize {
        taxon * self.nchar + c
    }
    /// fills leaf bitset/range caches from the dataset matrix.

    fn rebuild_leaf_cache(&mut self, ds: &Dataset) {
        let need = ds.ntax * ds.nchar;

        if self.leaf_bits.len() != need {
            self.leaf_bits.resize(need, 0);
        }
        if self.leaf_rg.len() != need {
            self.leaf_rg.resize(need, util::AdditiveRange::default());
        }

        for t in 0..ds.ntax {
            for c in 0..ds.nchar {
                let idx = t * ds.nchar + c;
                let bits = ds.matrix[idx].bits();
                self.leaf_bits[idx] = bits;
                self.leaf_rg[idx] = util::AdditiveRange::from_bits(bits);
            }
        }
    }
    /// grows cache arrays when a partial tree gains internal nodes.

    fn ensure_shape(&mut self, tree: &Tree, nchar: usize) {
        let needed_stride = tree.nodes.len().max(self.nnode);

        if self.nchar != nchar {
            self.nchar = nchar;
            self.nnode = needed_stride;

            let flat = self.nnode * self.nchar;
            self.st = vec![0; flat];
            self.rg = vec![util::AdditiveRange::default(); flat];
            self.node_steps_nonadd = vec![0; flat];
            self.node_steps_add = vec![0; flat];
            self.char_steps = vec![0; self.nchar];
            self.subtree_min_taxon = vec![usize::MAX; self.nnode];

            self.total = 0;
            return;
        }

        if needed_stride <= self.nnode {
            return;
        }

        let old_stride = self.nnode;
        let new_stride = needed_stride;

        self.st = grow_flat_copy(&self.st, old_stride, new_stride, self.nchar);
        self.rg = grow_flat_copy(&self.rg, old_stride, new_stride, self.nchar);
        self.node_steps_nonadd =
            grow_flat_copy(&self.node_steps_nonadd, old_stride, new_stride, self.nchar);
        self.node_steps_add =
            grow_flat_copy(&self.node_steps_add, old_stride, new_stride, self.nchar);

        let mut new_min = vec![usize::MAX; new_stride];
        new_min[..old_stride].copy_from_slice(&self.subtree_min_taxon[..old_stride]);
        self.subtree_min_taxon = new_min;

        self.nnode = new_stride;
    }
    /// maps character/node coordinates into flat cache arrays.

    fn idx(&self, c: usize, v: usize) -> usize {
        c * self.nnode + v
    }
    /// recomputes all cached node states, steps, and subtree minimum taxa for the
    /// current tree.

    fn rebuild(&mut self, ds: &Dataset, cfg: &CharConfig, tree: &Tree) {
        self.ensure_shape(tree, ds.nchar);
        self.ensure_groups(ds, cfg);

        self.total = 0;
        self.char_steps.fill(0);

        for x in &mut self.subtree_min_taxon {
            *x = usize::MAX;
        }

        let mut post = Vec::with_capacity(tree.nodes.len());
        postorder_fill(tree, tree.root, &mut post);

        for &v in &post {
            if let Some(t) = tree.nodes[v].taxon {
                self.subtree_min_taxon[v] = t;
            } else {
                let mut m = usize::MAX;
                for &ch in &tree.nodes[v].children {
                    m = m.min(self.subtree_min_taxon[ch]);
                }
                self.subtree_min_taxon[v] = m;
            }
        }

        for group in &self.char_groups {
            let rep = group.chars[0];
            let rep_w = cfg.chars[rep].weight as u64;

            if !group.additive {
                let mut unweighted_total = 0u64;

                for (i, node) in tree.nodes.iter().enumerate() {
                    let fi = self.idx(rep, i);
                    if let Some(t) = node.taxon {
                        self.st[fi] = self.leaf_bits[self.leaf_idx(t, rep)];
                        self.node_steps_nonadd[fi] = 0;
                    }
                }

                for &v in &post {
                    let node = &tree.nodes[v];
                    if node.taxon.is_some() {
                        continue;
                    }

                    let ev = eval_node_nonadd_from_children(
                        tree,
                        rep,
                        self.nnode,
                        &self.st,
                        &node.children,
                    );
                    let fi = self.idx(rep, v);
                    self.st[fi] = ev.st;
                    /* Store WEIGHTED for rep char */
                    self.node_steps_nonadd[fi] = ev.step * rep_w;
                    unweighted_total += ev.step;
                }

                for &c in &group.chars {
                    let w = cfg.chars[c].weight as u64;
                    let weighted = unweighted_total * w;
                    self.char_steps[c] = weighted;
                    if c != rep {
                        for v in 0..self.nnode {
                            let fi_rep = self.idx(rep, v);
                            let fi_c = self.idx(c, v);
                            self.st[fi_c] = self.st[fi_rep];
                            /* Scale from rep's weighted value */
                            self.node_steps_nonadd[fi_c] =
                                self.node_steps_nonadd[fi_rep] / rep_w * w;
                        }
                    }
                }
                self.total += unweighted_total * group.total_weight;
            } else {
                let mut unweighted_total = 0u64;

                for (i, node) in tree.nodes.iter().enumerate() {
                    let fi = self.idx(rep, i);
                    if let Some(t) = node.taxon {
                        self.rg[fi] = self.leaf_rg[self.leaf_idx(t, rep)];
                        self.node_steps_add[fi] = 0;
                    }
                }

                for &v in &post {
                    let node = &tree.nodes[v];
                    if node.taxon.is_some() {
                        continue;
                    }

                    let ev = eval_node_add_from_children(
                        tree,
                        rep,
                        self.nnode,
                        &self.rg,
                        &node.children,
                    );
                    let fi = self.idx(rep, v);
                    self.rg[fi] = ev.rg;
                    /* Store WEIGHTED for rep char */
                    self.node_steps_add[fi] = ev.step * rep_w;
                    unweighted_total += ev.step;
                }

                for &c in &group.chars {
                    let w = cfg.chars[c].weight as u64;
                    let weighted = unweighted_total * w;
                    self.char_steps[c] = weighted;
                    if c != rep {
                        for v in 0..self.nnode {
                            let fi_rep = self.idx(rep, v);
                            let fi_c = self.idx(c, v);
                            self.rg[fi_c] = self.rg[fi_rep];
                            self.node_steps_add[fi_c] = self.node_steps_add[fi_rep] / rep_w * w;
                        }
                    }
                }
                self.total += unweighted_total * group.total_weight;
            }
        }
    }
    /// repairs subtree-minimum taxon keys along paths affected by an insertion.

    fn update_subtree_min_after_insert(
        &mut self,
        tree: &Tree,
        child: usize,
        clip_taxon: usize,
        new_internal: usize,
    ) {
        /*
          Make sure all leaf entries have valid representatives.
          This is cheap compared with scoring and prevents usize::MAX leaks
          when subtree_min_taxon was resized or not initialized for a taxon.
        */
        for t in 0..tree.ntax {
            if t < self.subtree_min_taxon.len() {
                self.subtree_min_taxon[t] = t;
            }
        }

        let child_min = self
            .subtree_min_taxon
            .get(child)
            .copied()
            .filter(|&x| x != usize::MAX)
            .unwrap_or(child);

        self.subtree_min_taxon[new_internal] = child_min.min(clip_taxon);

        let mut cur = Some(new_internal);
        while let Some(v) = cur {
            if let Some(t) = tree.nodes[v].taxon {
                self.subtree_min_taxon[v] = t;
            } else {
                let mut m = usize::MAX;
                for &ch in &tree.nodes[v].children {
                    let cm = self
                        .subtree_min_taxon
                        .get(ch)
                        .copied()
                        .filter(|&x| x != usize::MAX)
                        .unwrap_or(ch);
                    m = m.min(cm);
                }
                self.subtree_min_taxon[v] = m;
            }

            cur = tree.nodes[v].parent;
        }
    }
    /// updates caches after an insertion is accepted into the growing Wagner tree.

    fn apply_committed_insert(
        &mut self,
        ds: &Dataset,
        cfg: &CharConfig,
        tree: &Tree,
        parent: usize,
        child: usize,
        clip_taxon: usize,
        new_internal: usize,
        path_buf: &mut Vec<usize>,
    ) {
        self.ensure_shape(tree, ds.nchar);

        self.update_subtree_min_after_insert(tree, child, clip_taxon, new_internal);

        let stride = self.nnode;

        fill_path_to_root(tree, new_internal, path_buf);

        for group in &self.char_groups {
            let rep = group.chars[0];
            let rep_w = cfg.chars[rep].weight as u64;

            if !group.additive {
                let base = self.idx(rep, 0);
                let clip_bits = self.leaf_bits[self.leaf_idx(clip_taxon, rep)];
                let child_bits = self.st[base + child];

                let (below_state, below_step) = merge_bits(clip_bits, child_bits);

                let ni_fi = base + new_internal;
                self.st[ni_fi] = below_state;
                self.node_steps_nonadd[ni_fi] = below_step * rep_w;

                for &v in path_buf.iter() {
                    if v == new_internal {
                        continue;
                    }

                    let fi = base + v;
                    let ev = eval_node_nonadd_committed(tree, rep, stride, &self.st, v);

                    self.st[fi] = ev.st;
                    self.node_steps_nonadd[fi] = ev.step * rep_w;
                }

                let mut total_delta = 0i128;
                for &c in &group.chars {
                    let w = cfg.chars[c].weight as u64;
                    let base_c = self.idx(c, 0);

                    let old_char_total = self.char_steps[c] as i128;
                    let mut new_char_total = old_char_total;

                    new_char_total += (below_step * w) as i128;
                    self.node_steps_nonadd[base_c + new_internal] = below_step * w;
                    if c != rep {
                        self.st[base_c + new_internal] = below_state;
                    }

                    for &v in path_buf.iter() {
                        if v == new_internal {
                            continue;
                        }
                        let fi = base_c + v;
                        let old_step = self.node_steps_nonadd[fi];
                        let ev = eval_node_nonadd_committed(tree, c, stride, &self.st, v);
                        let new_step = ev.step * w;
                        self.node_steps_nonadd[fi] = new_step;
                        if c != rep {
                            self.st[fi] = ev.st;
                        }
                        new_char_total = new_char_total
                            .saturating_add((new_step as i128).saturating_sub(old_step as i128));
                    }

                    self.char_steps[c] = new_char_total as u64;
                    total_delta += new_char_total - old_char_total;
                }

                if total_delta >= 0 {
                    self.total = self.total.saturating_add(total_delta as u64);
                } else {
                    self.total = self.total.saturating_sub((-total_delta) as u64);
                }
            } else {
                let base = self.idx(rep, 0);
                let clip_rg = self.leaf_rg[self.leaf_idx(clip_taxon, rep)];
                let child_rg = self.rg[base + child];

                let (below_rg, below_step) = clip_rg.merge(child_rg);

                let ni_fi = base + new_internal;
                self.rg[ni_fi] = below_rg;
                self.node_steps_add[ni_fi] = below_step * rep_w;

                for &v in path_buf.iter() {
                    if v == new_internal {
                        continue;
                    }

                    let fi = base + v;
                    let ev = eval_node_add_committed(tree, rep, stride, &self.rg, v);

                    self.rg[fi] = ev.rg;
                    self.node_steps_add[fi] = ev.step * rep_w;
                }

                let mut total_delta = 0i128;
                for &c in &group.chars {
                    let w = cfg.chars[c].weight as u64;
                    let base_c = self.idx(c, 0);

                    let old_char_total = self.char_steps[c] as i128;
                    let mut new_char_total = old_char_total;

                    new_char_total += (below_step * w) as i128;
                    self.node_steps_add[base_c + new_internal] = below_step * w;
                    if c != rep {
                        self.rg[base_c + new_internal] = below_rg;
                    }

                    for &v in path_buf.iter() {
                        if v == new_internal {
                            continue;
                        }
                        let fi = base_c + v;
                        let old_step = self.node_steps_add[fi];
                        let ev = eval_node_add_committed(tree, c, stride, &self.rg, v);
                        let new_step = ev.step * w;
                        self.node_steps_add[fi] = new_step;
                        if c != rep {
                            self.rg[fi] = ev.rg;
                        }
                        new_char_total = new_char_total
                            .saturating_add((new_step as i128).saturating_sub(old_step as i128));
                    }

                    self.char_steps[c] = new_char_total as u64;
                    total_delta += new_char_total - old_char_total;
                }

                if total_delta >= 0 {
                    self.total = self.total.saturating_add(total_delta as u64);
                } else {
                    self.total = self.total.saturating_sub((-total_delta) as u64);
                }
            }
        }

        let _ = parent;
    }
}

fn grow_flat_copy<T: Copy + Default>(
    old: &[T],
    old_stride: usize,
    new_stride: usize,
    nchar: usize,
) -> Vec<T> {
    let mut new_vec = vec![T::default(); new_stride * nchar];

    for c in 0..nchar {
        let old_start = c * old_stride;
        let new_start = c * new_stride;

        new_vec[new_start..new_start + old_stride]
            .copy_from_slice(&old[old_start..old_start + old_stride]);
    }

    new_vec
}
/// writes a node-to-root parent chain into a reusable vector.

fn fill_path_to_root(tree: &Tree, mut v: usize, out: &mut Vec<usize>) {
    out.clear();
    out.push(v);
    while let Some(p) = tree.nodes[v].parent {
        out.push(p);
        v = p;
    }
}
/// evaluates one non-additive internal node from candidate child states.

fn eval_node_nonadd_from_children(
    tree: &Tree,
    c: usize,
    nnode: usize,
    st: &[u64],
    children: &[usize],
) -> LocalEval {
    let _ = tree;

    let mut iter = children.iter().copied();
    let first = iter.next().expect("internal node must have children");
    let base = c * nnode;

    let mut acc = st[base + first];
    let mut steps = 0u64;

    for ch in iter {
        let b = st[base + ch];
        let (merged, step) = merge_bits(acc, b);
        acc = merged;
        steps += step;
    }

    LocalEval { st: acc, rg: util::AdditiveRange::default(), step: steps }
}
/// evaluates one additive internal node from candidate child intervals.

fn eval_node_add_from_children(
    tree: &Tree,
    c: usize,
    nnode: usize,
    rg: &[util::AdditiveRange],
    children: &[usize],
) -> LocalEval {
    let _ = tree;

    let mut iter = children.iter().copied();
    let first = iter.next().expect("internal node must have children");
    let base = c * nnode;

    let mut acc = rg[base + first];
    let mut steps = 0u64;

    for ch in iter {
        let b = rg[base + ch];
        let (merged, gap) = acc.merge(b);
        acc = merged;
        steps += gap;
    }

    LocalEval { st: 0, rg: acc, step: steps }
}
/// evaluates one committed non-additive node in the Wagner cache.

fn eval_node_nonadd_committed(
    tree: &Tree,
    c: usize,
    nnode: usize,
    st: &[u64],
    v: usize,
) -> LocalEval {
    let children = &tree.nodes[v].children;
    eval_node_nonadd_from_children(tree, c, nnode, st, children)
}
/// evaluates one committed additive node in the Wagner cache.

fn eval_node_add_committed(
    tree: &Tree,
    c: usize,
    nnode: usize,
    rg: &[util::AdditiveRange],
    v: usize,
) -> LocalEval {
    let children = &tree.nodes[v].children;
    eval_node_add_from_children(tree, c, nnode, rg, children)
}
/// incrementally scores inserting one taxon on one edge of the partial Wagner
/// tree.

fn score_insert_candidate_wagner_with_path(
    ds: &Dataset,
    cfg: &CharConfig,
    tree: &Tree,
    cache: &WagnerCharCache,
    parent: usize,
    child: usize,
    clip_taxon: usize,
    path: &[usize],
) -> u64 {
    let stride = cache.nnode;

    debug_assert_eq!(path.first().copied(), Some(parent));

    let mut total = cache.total;

    for group in &cache.char_groups {
        let rep = group.chars[0];
        let rep_w = cfg.chars[rep].weight as u64;
        let base = rep * stride;

        if !group.additive {
            let clip_bits = ds.matrix[clip_taxon * ds.nchar + rep].bits();
            let child_bits = cache.st[base + child];

            let inter = clip_bits & child_bits;
            let (mut below_state, below_step) =
                if inter != 0 { (inter, 0u64) } else { (clip_bits | child_bits, 1u64) };

            let mut unweighted_delta: i128 = below_step as i128;

            for i in 0..path.len() {
                let v = path[i];
                let child_on_path = if i == 0 { child } else { path[i - 1] };

                let old_fi = base + v;
                unweighted_delta -= (cache.node_steps_nonadd[old_fi] / rep_w) as i128;

                let children = &tree.nodes[v].children;
                let mut local_step = 0u64;

                let mut iter = children.iter().copied();
                let first = iter.next().expect("internal node must have children");

                let mut acc =
                    if first == child_on_path { below_state } else { cache.st[base + first] };

                for ch in iter {
                    let b = if ch == child_on_path { below_state } else { cache.st[base + ch] };

                    let (merged, step) = merge_bits(acc, b);
                    acc = merged;
                    local_step += step;
                }

                below_state = acc;
                unweighted_delta += local_step as i128;
            }

            let group_delta = unweighted_delta * group.total_weight as i128;
            if group_delta < 0 {
                total = total.saturating_sub((-group_delta) as u64);
            } else {
                total = total.saturating_add(group_delta as u64);
            }
        } else {
            let clip_rg = cache.leaf_rg[cache.leaf_idx(clip_taxon, rep)];
            let child_rg = cache.rg[base + child];

            let (mut below_rg, below_step) = clip_rg.merge(child_rg);

            let mut unweighted_delta: i128 = below_step as i128;

            for i in 0..path.len() {
                let v = path[i];
                let child_on_path = if i == 0 { child } else { path[i - 1] };

                let old_fi = base + v;
                unweighted_delta -= (cache.node_steps_add[old_fi] / rep_w) as i128;

                let children = &tree.nodes[v].children;
                let mut local_step = 0u64;

                let mut iter = children.iter().copied();
                let first = iter.next().expect("internal node must have children");

                let mut acc =
                    if first == child_on_path { below_rg } else { cache.rg[base + first] };

                for ch in iter {
                    let b = if ch == child_on_path { below_rg } else { cache.rg[base + ch] };

                    let (merged, gap) = acc.merge(b);
                    acc = merged;
                    local_step += gap;
                }

                below_rg = acc;
                unweighted_delta += local_step as i128;
            }

            let group_delta = unweighted_delta * group.total_weight as i128;
            if group_delta < 0 {
                total = total.saturating_sub((-group_delta) as u64);
            } else {
                total = total.saturating_add(group_delta as u64);
            }
        }
    }

    total
}

/* hennig core */

const WAGNER_MIN_EXACT_EDGES: usize = 24;
const WAGNER_CHEAP_WINDOW: u64 = 2;
const WAGNER_VALIDATE_CACHE: bool = false;
/// constructs one complete Wagner tree from a deterministic taxon order and
/// best insertion edges.

fn hennig_once_tree(ds: &Dataset, cfg: &CharConfig, seed: u64, outgroup: usize) -> Tree {
    assert!(outgroup < ds.ntax, "outgroup out of range");

    let mut rng = Rng::new(seed ^ 0xa076_1d64_78bd_642f);
    let ingroup = random_order_ingroup_taxa(ds, seed, outgroup);

    let a = ingroup[0];
    let b = ingroup[1];
    let mut tr = Tree::new_outgrouped(ds.ntax, outgroup, a, b);

    let mut cache = WagnerCharCache::new(ds, cfg, &tr);
    let mut path_buf = Vec::<usize>::with_capacity(ds.ntax * 2 + 8);
    let mut commit_path_buf = Vec::<usize>::with_capacity(ds.ntax * 2 + 8);

    for &clip in &ingroup[2..] {
        let edges = tr.rooted_edges();

        let mut ranked_edges: Vec<((usize, usize), u64)> = edges
            .into_iter()
            .filter(|&(_, child)| child != outgroup)
            .map(|(parent, child)| {
                let cheap = candidate_edge_cheap_score(ds, cfg, &cache, clip, child);
                ((parent, child), cheap)
            })
            .collect();

        ranked_edges.sort_by(|a, b| a.1.cmp(&b.1));

        let cheap_best = ranked_edges.first().map(|x| x.1).unwrap_or(0);

        let mut best_score = u64::MAX;
        let mut best_edges = Vec::<(usize, usize)>::new();

        for (i, &((parent, child), cheap)) in ranked_edges.iter().enumerate() {
            if i >= WAGNER_MIN_EXACT_EDGES
                && cheap > cheap_best.saturating_add(WAGNER_CHEAP_WINDOW)
                && !best_edges.is_empty()
            {
                break;
            }

            fill_path_to_root(&tr, parent, &mut path_buf);

            let s = score_insert_candidate_wagner_with_path(
                ds, cfg, &tr, &cache, parent, child, clip, &path_buf,
            );

            if s < best_score {
                best_score = s;
                best_edges.clear();
                best_edges.push((parent, child));
            } else if s == best_score {
                best_edges.push((parent, child));
            }
        }

        let pick = if best_edges.len() <= 1 { 0 } else { rng.next_usize(best_edges.len()) };

        let (parent, child) = best_edges[pick];
        let new_internal = tr.insert_taxon_on_edge(parent, child, clip);

        cache.apply_committed_insert(
            ds,
            cfg,
            &tr,
            parent,
            child,
            clip,
            new_internal,
            &mut commit_path_buf,
        );

        if WAGNER_VALIDATE_CACHE {
            let mut ws_check = ScoreWorkspace::new();
            let exact = score_tree(ds, cfg, &tr, &mut ws_check);
            assert_eq!(
                cache.total, exact,
                "Wagner cache mismatch after inserting taxon {clip}: cache={} exact={}",
                cache.total, exact
            );
        }
    }

    tr
}
/// builds a shuffled ingroup order for alternate starts.

fn random_order_ingroup_taxa(ds: &Dataset, seed: u64, outgroup: usize) -> Vec<usize> {
    let mut taxa: Vec<usize> = (0..ds.ntax).filter(|&t| t != outgroup).collect();

    let mut rng = Rng::new(seed ^ 0x9e37_79b9_7f4a_7c15);
    rng.shuffle(&mut taxa);

    taxa
}
/// computes a lightweight edge heuristic before full insertion scoring.

fn candidate_edge_cheap_score(
    ds: &Dataset,
    _cfg: &CharConfig,
    cache: &WagnerCharCache,
    clip_taxon: usize,
    child: usize,
) -> u64 {
    let rep_taxon = cache
        .subtree_min_taxon
        .get(child)
        .copied()
        .filter(|&t| t < ds.ntax)
        .unwrap_or_else(|| if child < ds.ntax { child } else { 0 });

    debug_assert!(rep_taxon < ds.ntax, "invalid rep_taxon={} for child={}", rep_taxon, child);

    let mut score = 0u64;

    for group in &cache.char_groups {
        let rep = group.chars[0];
        let w = group.total_weight;

        if !group.additive {
            let a = cache.leaf_bits[cache.leaf_idx(clip_taxon, rep)];
            let b = cache.leaf_bits[cache.leaf_idx(rep_taxon, rep)];

            if a == crate::engines::dataset::StateSet::ALL36.bits()
                || b == crate::engines::dataset::StateSet::ALL36.bits()
            {
                continue;
            }

            if (a & b) == 0 {
                score += w;
            }
        } else {
            let a_bits = cache.leaf_bits[cache.leaf_idx(clip_taxon, rep)];
            let b_bits = cache.leaf_bits[cache.leaf_idx(rep_taxon, rep)];

            if a_bits == crate::engines::dataset::StateSet::ALL36.bits()
                || b_bits == crate::engines::dataset::StateSet::ALL36.bits()
            {
                continue;
            }

            let a = cache.leaf_rg[cache.leaf_idx(clip_taxon, rep)];
            let b = cache.leaf_rg[cache.leaf_idx(rep_taxon, rep)];
            let (_, gap) = a.merge(b);
            score += gap * w;
        }
    }

    score
}

/*  scoring  */
/// computes total parsimony length for a full tree using reusable workspace
/// buffers.

pub(crate) fn score_tree(
    ds: &Dataset,
    cfg: &CharConfig,
    tree: &Tree,
    ws: &mut ScoreWorkspace,
) -> u64 {
    ws.prepare(tree, ds.nchar);
    ws.ensure_leaf_cache(ds);
    ws.ensure_groups(ds, cfg);

    let nchar = ds.nchar;
    let mut total: u64 = 0;

    let n_groups = ws.char_groups.len();
    for gi in 0..n_groups {
        let group_additive = ws.char_groups[gi].additive;
        let rep = ws.char_groups[gi].chars[0];
        let total_weight = ws.char_groups[gi].total_weight;

        if !group_additive {
            let mut steps: u64 = 0;

            for (i, node) in tree.nodes.iter().enumerate() {
                if let Some(t) = node.taxon {
                    ws.st[i] = ws.leaf_bits[t * nchar + rep];
                }
            }

            for &v in &ws.post {
                let node = &tree.nodes[v];
                if node.taxon.is_some() {
                    continue;
                }

                let mut iter = node.children.iter().copied();
                let first = iter.next().expect("internal node must have children");
                let mut acc = ws.st[first];

                for ch in iter {
                    let b = ws.st[ch];
                    let (merged, step) = merge_bits(acc, b);
                    acc = merged;
                    steps += step;
                }

                ws.st[v] = acc;
            }

            let weighted = steps * total_weight;
            ws.char_steps[rep] = weighted;
            total += weighted;
        } else if ws.char_groups[gi].all_contiguous {
            /* Fast-range interval scoring for contiguous additive characters. */
            let mut steps = 0u64;
            for (i, node) in tree.nodes.iter().enumerate() {
                if let Some(t) = node.taxon {
                    ws.rg[i] = ws.leaf_rg[t * nchar + rep];
                }
            }
            for &v in &ws.post {
                let node = &tree.nodes[v];
                if node.taxon.is_some() {
                    continue;
                }
                let mut iter = node.children.iter().copied();
                let first = iter.next().expect("internal node must have children");
                let mut acc = ws.rg[first];
                for ch in iter {
                    let (merged, gap) = acc.merge(ws.rg[ch]);
                    acc = merged;
                    steps += gap;
                }
                ws.rg[v] = acc;
            }
            let weighted = steps * total_weight;
            ws.char_steps[rep] = weighted;
            total += weighted;
        } else {
            let steps =
                sankoff_ordered_char_steps_fast(ds, tree, rep, &ws.post, &mut ws.sankoff_costs);
            let weighted = steps * total_weight;
            ws.char_steps[rep] = weighted;
            total += weighted;
        }
    }

    total
}

pub(crate) fn score_tree_bounded(
    ds: &Dataset,
    cfg: &CharConfig,
    tree: &Tree,
    ws: &mut ScoreWorkspace,
    limit: u64,
) -> Option<u64> {
    ws.prepare(tree, ds.nchar);
    ws.ensure_leaf_cache(ds);
    ws.ensure_groups(ds, cfg);

    let nchar = ds.nchar;
    let mut total: u64 = 0;

    let n_groups = ws.char_groups.len();
    for gi in 0..n_groups {
        let group_additive = ws.char_groups[gi].additive;
        let rep = ws.char_groups[gi].chars[0];
        let total_weight = ws.char_groups[gi].total_weight;

        if !group_additive {
            let mut steps: u64 = 0;

            for (i, node) in tree.nodes.iter().enumerate() {
                if let Some(t) = node.taxon {
                    ws.st[i] = ws.leaf_bits[t * nchar + rep];
                }
            }

            for &v in &ws.post {
                let node = &tree.nodes[v];
                if node.taxon.is_some() {
                    continue;
                }

                let mut iter = node.children.iter().copied();
                let first = iter.next().expect("internal node must have children");
                let mut acc = ws.st[first];

                for ch in iter {
                    let b = ws.st[ch];
                    let (merged, step) = merge_bits(acc, b);
                    acc = merged;
                    steps += step;
                }

                ws.st[v] = acc;
            }

            let weighted = steps * total_weight;
            ws.char_steps[rep] = weighted;
            total = total.saturating_add(weighted);
        } else if ws.char_groups[gi].all_contiguous {
            let mut steps = 0u64;
            for (i, node) in tree.nodes.iter().enumerate() {
                if let Some(t) = node.taxon {
                    ws.rg[i] = ws.leaf_rg[t * nchar + rep];
                }
            }
            for &v in &ws.post {
                let node = &tree.nodes[v];
                if node.taxon.is_some() {
                    continue;
                }
                let mut iter = node.children.iter().copied();
                let first = iter.next().expect("internal node must have children");
                let mut acc = ws.rg[first];
                for ch in iter {
                    let (merged, gap) = acc.merge(ws.rg[ch]);
                    acc = merged;
                    steps += gap;
                }
                ws.rg[v] = acc;
            }
            let weighted = steps * total_weight;
            ws.char_steps[rep] = weighted;
            total = total.saturating_add(weighted);
        } else {
            let steps =
                sankoff_ordered_char_steps_fast(ds, tree, rep, &ws.post, &mut ws.sankoff_costs);
            let weighted = steps * total_weight;
            ws.char_steps[rep] = weighted;
            total = total.saturating_add(weighted);
        }

        if total > limit {
            return None;
        }
    }

    Some(total)
}
/// applies ordered-state transition costs to a cost vector for Sankoff
/// scoring.
#[inline(always)]
pub(crate) fn min_cost_with_transition(costs: &[u64]) -> [u64; NSTATES] {
    let mut forward = [INF; NSTATES];
    forward[0] = costs[0];
    for s in 1..NSTATES {
        forward[s] = costs[s].min(forward[s - 1].saturating_add(1));
    }

    let mut backward = [INF; NSTATES];
    backward[NSTATES - 1] = costs[NSTATES - 1];
    for s in (0..NSTATES - 1).rev() {
        backward[s] = costs[s].min(backward[s + 1].saturating_add(1));
    }

    let mut result = [INF; NSTATES];
    for s in 0..NSTATES {
        result[s] = forward[s].min(backward[s]);
    }
    result
}
/// computes exact ordered-state steps for one additive character with dynamic
/// programming.

fn sankoff_ordered_char_steps_fast(
    ds: &Dataset,
    tree: &Tree,
    ch: usize,
    postorder: &[usize],
    costs: &mut [u64],
) -> u64 {
    costs.fill(INF);

    for (i, node) in tree.nodes.iter().enumerate() {
        if let Some(t) = node.taxon {
            let bits = ds.matrix[t * ds.nchar + ch].bits();
            let base = i * NSTATES;
            for s in 0..NSTATES {
                if ((bits >> s) & 1) == 1 {
                    costs[base + s] = 0;
                }
            }
        }
    }

    for &v in postorder {
        let node = &tree.nodes[v];
        if node.taxon.is_some() {
            continue;
        }

        let base = v * NSTATES;
        let mut iter = node.children.iter();
        let first = *iter.next().expect("internal node must have children");
        let mut acc = min_cost_with_transition(&costs[first * NSTATES..][..NSTATES]);
        for &ch in iter {
            let cm = min_cost_with_transition(&costs[ch * NSTATES..][..NSTATES]);
            for s in 0..NSTATES {
                acc[s] = acc[s].saturating_add(cm[s]);
            }
        }
        costs[base..base + NSTATES].copy_from_slice(&acc);
    }

    let root_base = tree.root * NSTATES;
    costs[root_base..root_base + NSTATES].iter().copied().min().unwrap_or(0)
}
/// records a postorder traversal of a rooted tree.

fn postorder_fill(tree: &Tree, v: usize, out: &mut Vec<usize>) {
    let n = &tree.nodes[v];
    if n.taxon.is_some() {
        out.push(v);
        return;
    }
    for &c in &n.children {
        postorder_fill(tree, c, out);
    }
    out.push(v);
}

/*  rooted canonical hash  */
/// hashes a rooted topology deterministically after canonical binary-view
/// traversal.

pub(crate) fn rooted_topology_hash(tree: &Tree) -> u64 {
    fn mix(a: u64, b: u64) -> u64 {
        let mut x = a ^ b.rotate_left(17) ^ 0x9e37_79b9_7f4a_7c15;
        x ^= x >> 30;
        x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x ^= x >> 27;
        x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
        x ^ (x >> 31)
    }

    let bv = match BinaryView::from_tree(tree) {
        Ok(v) => v,
        Err(_) => {
            return 0;
        }
    };
    /// local recursive topology hash helper.

    fn rec(tree: &Tree, bv: &BinaryView, v: usize) -> u64 {
        let node = &tree.nodes[v];
        if let Some(t) = node.taxon {
            return mix(0x1234_5678_9abc_def0, t as u64);
        }

        let a = bv.left[v].unwrap();
        let b = bv.right[v].unwrap();

        let ha = rec(tree, bv, a);
        let hb = rec(tree, bv, b);

        let (x, y) = if ha <= hb { (ha, hb) } else { (hb, ha) };

        let mut acc = 0xfeed_face_cafe_babeu64;
        acc = mix(acc, x);
        acc = mix(acc, y);
        acc
    }

    rec(tree, &bv, bv.root)
}

/*CI / RI  */
/// computes consistency index from minimum and observed tree length.

pub(crate) fn calc_ci(min_len: u64, tree_len: u64) -> f64 {
    if tree_len == 0 { 0.0 } else { min_len as f64 / tree_len as f64 }
}
/// computes retention index when the denominator is defined.

pub(crate) fn calc_ri(min_len: u64, max_len: u64, tree_len: u64) -> Option<f64> {
    if max_len == min_len {
        None
    } else {
        Some((max_len - tree_len) as f64 / (max_len - min_len) as f64)
    }
}
/// computes the weighted dataset-wide minimum possible steps.

pub(crate) fn minsteps_sum(ds: &Dataset, cfg: &CharConfig) -> u64 {
    let mut total = 0u64;

    for c in 0..ds.nchar {
        let cs = cfg.chars[c];
        if !cs.active {
            continue;
        }
        let w = cs.weight as u64;

        if !cs.additive {
            let mut union_bits = 0u64;
            for t in 0..ds.ntax {
                let bits = ds.matrix[t * ds.nchar + c].bits();
                if bits != crate::engines::dataset::StateSet::ALL36.bits() {
                    union_bits |= bits;
                }
            }
            let k = union_bits.count_ones() as u64;
            if k > 0 {
                total += (k - 1) * w;
            }
        } else {
            let mut lo: Option<u8> = None;
            let mut hi: u8 = 0;
            for t in 0..ds.ntax {
                let bits = ds.matrix[t * ds.nchar + c].bits();
                if bits == crate::engines::dataset::StateSet::ALL36.bits() {
                    continue;
                }
                let r = util::AdditiveRange::from_bits(bits);
                lo = Some(lo.map(|x| x.min(r.lo)).unwrap_or(r.lo));
                hi = hi.max(r.hi);
            }
            if let Some(lo) = lo {
                total += (hi - lo) as u64 * w;
            }
        }
    }

    total
}
/// computes the weighted dataset-wide maximum possible steps.

pub(crate) fn maxsteps_sum(ds: &Dataset, cfg: &CharConfig) -> u64 {
    let mut total = 0u64;

    for c in 0..ds.nchar {
        let cs = cfg.chars[c];
        if !cs.active {
            continue;
        }
        let w = cs.weight as u64;

        if !cs.additive {
            let mut informative_taxa = 0u64;
            let mut freq = [0u64; 36];
            let mut union_bits = 0u64;

            for t in 0..ds.ntax {
                let bits = ds.matrix[t * ds.nchar + c].bits();
                if bits == crate::engines::dataset::StateSet::ALL36.bits() {
                    continue;
                }

                informative_taxa += 1;
                union_bits |= bits;

                if bits.count_ones() == 1 {
                    let idx = bits.trailing_zeros() as usize;
                    freq[idx] += 1;
                }
            }

            if informative_taxa > 0 {
                let maxf = freq.into_iter().max().unwrap_or(0);
                let themin_char = {
                    let k = union_bits.count_ones() as u64;
                    if k > 0 { k - 1 } else { 0 }
                };

                let mut themax_char = informative_taxa.saturating_sub(maxf);
                if themax_char < themin_char {
                    themax_char = themin_char;
                }
                total += themax_char * w;
            }
        } else {
            let mut informative_taxa = 0u64;
            let mut lo: Option<u8> = None;
            let mut hi: u8 = 0;

            for t in 0..ds.ntax {
                let bits = ds.matrix[t * ds.nchar + c].bits();
                if bits == crate::engines::dataset::StateSet::ALL36.bits() {
                    continue;
                }
                informative_taxa += 1;
                let r = util::AdditiveRange::from_bits(bits);
                lo = Some(lo.map(|x| x.min(r.lo)).unwrap_or(r.lo));
                hi = hi.max(r.hi);
            }

            if let Some(lo) = lo {
                let span = (hi - lo) as u64;
                let themin_char = span;
                let mut themax_char = informative_taxa.saturating_sub(1) * span;
                if themax_char < themin_char {
                    themax_char = themin_char;
                }
                total += themax_char * w;
            }
        }
    }

    total
}

/* Random */

struct Rng {
    s: u64,
}

impl Rng {
    /// initializes the deterministic pseudo-random generator.
    fn new(seed: u64) -> Self {
        Self { s: if seed == 0 { 0x9e37_79b9_7f4a_7c15 } else { seed } }
    }
    /// advances and returns the next random word.

    fn next_u64(&mut self) -> u64 {
        let mut x = self.s;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.s = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    /// maps the next random word into a bounded index.

    fn next_usize(&mut self, upper: usize) -> usize {
        if upper <= 1 {
            return 0;
        }
        (self.next_u64() % (upper as u64)) as usize
    }

    fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = self.next_usize(i + 1);
            v.swap(i, j);
        }
    }
}
/// scrambles seeds and iteration counters into reproducible start seeds.

fn mix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

#[derive(Debug, Clone)]
pub(crate) struct BinaryView {
    pub root: usize,
    pub left: Vec<Option<usize>>,
    pub right: Vec<Option<usize>>,
}

impl BinaryView {
    /// converts a possibly-polytomous rooted tree into the binary/scoring view
    /// used by candidate caches.
    pub(crate) fn from_tree(tree: &Tree) -> Result<Self, String> {
        let n = tree.nodes.len();
        let mut left = vec![None; n];
        let mut right = vec![None; n];
        /// local BinaryView builder that recursively copies nodes and child links.

        fn dfs(
            tree: &Tree,
            v: usize,
            left_arr: &mut [Option<usize>],
            right_arr: &mut [Option<usize>],
        ) -> Result<(), String> {
            let node = &tree.nodes[v];
            if node.taxon.is_some() {
                return Ok(());
            }

            if node.children.len() != 2 {
                return Err(format!(
                    "BinaryView requires strictly binary tree; node {v} has {} children",
                    node.children.len()
                ));
            }

            let a = node.children[0];
            let b = node.children[1];

            left_arr[v] = Some(a);
            right_arr[v] = Some(b);

            dfs(tree, a, left_arr, right_arr)?;
            dfs(tree, b, left_arr, right_arr)?;
            Ok(())
        }

        dfs(tree, tree.root, &mut left, &mut right)?;

        Ok(Self { root: tree.root, left, right })
    }
}

/*  Zero-length branch collapse */

/// collapse of zero-length branches.
///
/// A branch is zero-length (unsupported) if for EVERY character, there exists a
/// most parsimonious reconstruction where the branch has length 0.
/// Equivalently: the MP-sets (possible state sets) of the two endpoint nodes
/// intersect for every character.
///
/// This function marks every internal branch whose minimum length is zero for
/// every active character group, then contracts all marked branches in one pass.
/// This matches minimum branch length collapse: the temporarily
/// collapsed tree may be longer than the original binary resolution.
pub fn collapse_zero_length_branches(
    tree: &mut Tree,
    ds: &Dataset,
    cfg: &CharConfig,
    ws: &mut ScoreWorkspace,
) -> bool {
    ws.ensure_groups(ds, cfg);

    if ws.char_groups.is_empty() {
        return false;
    }

    let mut collapse_edge = ws.compute_rule1_collapsible_edges(tree, ds);
    let mut any_collapsed = false;

    for (ch, node) in tree.nodes.iter().enumerate() {
        let Some(par) = node.parent else {
            collapse_edge[ch] = false;
            continue;
        };
        if par == tree.root || ch == tree.root {
            collapse_edge[ch] = false;
            continue;
        }
        if tree.nodes[par].taxon.is_some() || node.taxon.is_some() {
            collapse_edge[ch] = false;
            continue;
        }
        if collapse_edge[ch] {
            any_collapsed = true;
        } else {
            collapse_edge[ch] = false;
        }
    }

    if !any_collapsed {
        return false;
    }

    contract_marked_edges(tree, &collapse_edge);
    tree.canonicalize();
    true
}

impl ScoreWorkspace {
    /// determines which parent-child edges can be collapsed without increasing
    /// parsimony length.
    fn compute_rule1_collapsible_edges(&mut self, tree: &Tree, ds: &Dataset) -> Vec<bool> {
        let nnode = tree.nodes.len();
        let nchar = ds.nchar;
        let mut collapsible = vec![true; nnode];

        let groups = &self.char_groups;
        let needed = nnode * NSTATES;
        if self.coll_down.len() < needed {
            self.coll_down.resize(needed, INF);
            self.coll_child_transition.resize(needed, INF);
            self.coll_up.resize(needed, INF);
        }
        self.coll_post.clear();
        postorder_fill(tree, tree.root, &mut self.coll_post);
        let down = &mut self.coll_down[..needed];
        let child_transition = &mut self.coll_child_transition[..needed];
        let up = &mut self.coll_up[..needed];
        let post = &self.coll_post;
        for group in groups {
            let rep = group.chars[0];
            down.fill(INF);

            for (i, node) in tree.nodes.iter().enumerate() {
                if let Some(t) = node.taxon {
                    let leaf_bits = ds.matrix[t * nchar + rep].bits();
                    let base = i * NSTATES;
                    for s in 0..NSTATES {
                        if ((leaf_bits >> s) & 1) == 1 {
                            down[base + s] = 0;
                        }
                    }
                }
            }

            for &v in post {
                let node = &tree.nodes[v];
                if node.taxon.is_some() {
                    continue;
                }

                let base = v * NSTATES;
                let mut iter = node.children.iter();
                let first = *iter.next().expect("internal node must have children");
                let mut acc =
                    transition_min_costs(&down[first * NSTATES..][..NSTATES], group.additive);
                for &ch in iter {
                    let cm = transition_min_costs(&down[ch * NSTATES..][..NSTATES], group.additive);
                    for s in 0..NSTATES {
                        acc[s] = acc[s].saturating_add(cm[s]);
                    }
                }
                down[base..base + NSTATES].copy_from_slice(&acc);
            }

            for v in 0..nnode {
                let base = v * NSTATES;
                let mins = transition_min_costs(&down[base..base + NSTATES], group.additive);
                child_transition[base..base + NSTATES].copy_from_slice(&mins);
            }

            let root_base = tree.root * NSTATES;
            let tree_opt = down[root_base..root_base + NSTATES].iter().copied().min().unwrap_or(0);

            up.fill(INF);
            for s in 0..NSTATES {
                up[root_base + s] = 0;
            }

            for &p in post.iter().rev() {
                let node = &tree.nodes[p];
                if node.taxon.is_some() {
                    continue;
                }

                let pbase = p * NSTATES;
                let mut total_at_parent = [0u64; NSTATES];
                total_at_parent.copy_from_slice(&up[pbase..pbase + NSTATES]);

                for &child in &node.children {
                    let cbase = child * NSTATES;
                    for s in 0..NSTATES {
                        total_at_parent[s] =
                            total_at_parent[s].saturating_add(child_transition[cbase + s]);
                    }
                }

                for &child in &node.children {
                    let cbase = child * NSTATES;
                    let mut outside_child = [INF; NSTATES];
                    for s in 0..NSTATES {
                        outside_child[s] =
                            total_at_parent[s].saturating_sub(child_transition[cbase + s]);
                    }

                    if tree.nodes[child].taxon.is_none() && p != tree.root && child != tree.root {
                        let mut zero_len_possible = false;
                        for s in 0..NSTATES {
                            if outside_child[s].saturating_add(down[cbase + s]) == tree_opt {
                                zero_len_possible = true;
                                break;
                            }
                        }
                        if !zero_len_possible {
                            collapsible[child] = false;
                        }
                    }

                    let child_up = transition_min_costs(&outside_child, group.additive);
                    up[cbase..cbase + NSTATES].copy_from_slice(&child_up);
                }
            }
        }

        collapsible
    }
}
/// applies additive or non-additive transition rules to a state-cost vector.

#[inline(always)]
fn transition_min_costs(costs: &[u64], additive: bool) -> [u64; NSTATES] {
    if additive {
        return min_cost_with_transition(costs);
    }

    let min_any = costs.iter().copied().min().unwrap_or(INF);
    let min_change = min_any.saturating_add(1);
    let mut out = [INF; NSTATES];
    for s in 0..NSTATES {
        out[s] = costs[s].min(min_change);
    }
    out
}
/// rebuilds a tree after expanding children through marked collapsible edges.

fn contract_marked_edges(tree: &mut Tree, collapse_edge: &[bool]) {
    fn expanded_children(tree: &Tree, collapse_edge: &[bool], child: usize, out: &mut Vec<usize>) {
        if collapse_edge.get(child).copied().unwrap_or(false) {
            for &grandchild in &tree.nodes[child].children {
                expanded_children(tree, collapse_edge, grandchild, out);
            }
        } else {
            out.push(child);
        }
    }

    let old_children: Vec<Vec<usize>> = tree.nodes.iter().map(|n| n.children.clone()).collect();
    let mut replacements = vec![Vec::<usize>::new(); tree.nodes.len()];

    for v in 0..tree.nodes.len() {
        if tree.nodes[v].taxon.is_some() || collapse_edge.get(v).copied().unwrap_or(false) {
            continue;
        }
        for &child in &old_children[v] {
            expanded_children(tree, collapse_edge, child, &mut replacements[v]);
        }
    }

    for v in 0..tree.nodes.len() {
        if collapse_edge.get(v).copied().unwrap_or(false) {
            tree.nodes[v].children.clear();
            tree.nodes[v].parent = None;
            continue;
        }
        if tree.nodes[v].taxon.is_some() {
            continue;
        }
        let children = std::mem::take(&mut replacements[v]);
        tree.nodes[v].children = children.clone();
        for child in children {
            tree.nodes[child].parent = Some(v);
        }
    }
    tree.nodes[tree.root].parent = None;
}

/// Compute a clade-set based topology hash for deduplication purposes.
/// Unlike `rooted_topology_hash` (which requires a strictly binary tree),
/// this works for any tree and gives the same hash for trees that share
/// the same set of collapsed clades. Polytomy child order is ignored.
pub fn collapsed_topology_hash(tree: &Tree) -> u64 {
    fn mix(a: u64, b: u64) -> u64 {
        let mut x = a ^ b.rotate_left(17) ^ 0x9e37_79b9_7f4a_7c15;
        x ^= x >> 30;
        x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x ^= x >> 27;
        x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
        x ^ (x >> 31)
    }
    /// local helper counting bits across a clade bit-vector.

    fn bit_count(bits: &[u64]) -> u32 {
        bits.iter().map(|w| w.count_ones()).sum()
    }
    /// local helper hashing a clade bit-vector.

    fn hash_bits(bits: &[u64]) -> u64 {
        let mut h = mix(0x6a09_e667_f3bc_c909, bits.len() as u64);
        for &word in bits {
            h = mix(h, word);
        }
        h
    }

    fn rec(tree: &Tree, v: usize, words: usize, clades: &mut Vec<Vec<u64>>) -> Vec<u64> {
        let mut bits = vec![0u64; words];
        if let Some(t) = tree.nodes[v].taxon {
            bits[t / 64] |= 1u64 << (t % 64);
            return bits;
        }

        for &child in &tree.nodes[v].children {
            let child_bits = rec(tree, child, words, clades);
            for (dst, src) in bits.iter_mut().zip(child_bits) {
                *dst |= src;
            }
        }

        let count = bit_count(&bits);
        if v != tree.root && count > 1 && count < tree.ntax as u32 {
            clades.push(bits.clone());
        }

        bits
    }

    let words = tree.ntax.div_ceil(64);
    let mut clades = Vec::<Vec<u64>>::new();
    rec(tree, tree.root, words, &mut clades);
    clades.sort();

    let mut h = mix(0xfeed_face_cafe_babe, tree.ntax as u64);
    h = mix(h, clades.len() as u64);
    for clade in clades {
        h = mix(h, hash_bits(&clade));
    }
    h
}

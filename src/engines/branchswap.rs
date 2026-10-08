/*!Module: branch swapping and branch-breaking search operators. It implements
fast TBR side-message scoring, candidate ranking,
deduplication, and closure search over improved trees.

 */
use std::collections::{HashMap, HashSet};
use std::{collections::VecDeque, rc::Rc};

use crate::engines::{
    ccode::CharConfig,
    dataset::Dataset,
    search::{
        BinaryView, INF, NSTATES, ScoreWorkspace, collapsed_dedup_hash, merge_bits,
        min_cost_with_transition, rooted_topology_hash, score_tree,
    },
    trees::Tree,
    util::{self, is_contiguous},
};

/* Small fast-pass limits.  They affect traversal order only: before a tree is
declared locally optimal, an uncapped pass verifies the whole TBR
neighbourhood. */
const FAST_RECONNECTS_PER_SIDE: usize = 8;
const FAST_CANDIDATES_PER_CUT: usize = 16;
/*  Pattern – detect & compress identical character vectors  */

#[derive(Clone, Copy, Debug)]
enum AdditiveStrategy {
    FastRange,
    ExactCosts,
}

#[derive(Clone, Debug)]
struct PatternInfo {
    total_weight: u64,
    strategy: AdditiveStrategy,
    leaf_bits: Vec<u64>,                   /* [ntax] */
    leaf_ranges: Vec<util::AdditiveRange>, /* [ntax] */
}

#[derive(Clone, Copy)]
struct ActivePattern {
    pattern_id: usize,
    total_weight: u64,
    additive: bool,
}

/// Return whether an additive character can be scored as ordinary Fitch data.
/// This is safe only when all observed states occupy at most two adjacent
/// states; then every additive change has the same cost as a non-additive
/// change for the observed data.
fn can_downgrade_to_nonadditive(leaf_bits: &[u64]) -> bool {
    let mut global_lo = 36u8;
    let mut global_hi = 0u8;
    for &bits in leaf_bits {
        if bits == 0 || bits == crate::engines::dataset::StateSet::ALL36.bits() {
            continue;
        }
        if !is_contiguous(bits) {
            return false;
        }
        let lo = bits.trailing_zeros() as u8; /* low */
        let hi = 63 - bits.leading_zeros() as u8; /* high */
        global_lo = global_lo.min(lo);
        global_hi = global_hi.max(hi);
    }
    /* Only a two-state contiguous span preserves Fitch-equivalent costs. */
    global_hi.saturating_sub(global_lo) <= 1
}
/// groups active characters by identical taxon-state vectors and effective
/// additive mode for faster candidate scoring.

fn build_active_patterns(ds: &Dataset, cfg: &CharConfig) -> (Vec<ActivePattern>, Vec<PatternInfo>) {
    let nchar = ds.nchar;
    let ntax = ds.ntax;

    /* Group by the exact leaf-state vector and the effective scoring mode. */
    type PatternKey = (Vec<u64>, bool);
    let mut pattern_map: HashMap<PatternKey, (usize, u64)> = HashMap::new();
    let mut patterns: Vec<PatternInfo> = Vec::new();

    for c in 0..nchar {
        let cs = cfg.chars[c];
        if !cs.active {
            continue;
        }
        let bits_vec: Vec<u64> = (0..ntax).map(|t| ds.matrix[t * nchar + c].bits()).collect();

        let effective_additive = if cs.additive {
            /*
              Downgradable additive characters are cheaper to score as Fitch
              characters and produce identical costs under the condition above.
            */
            !can_downgrade_to_nonadditive(&bits_vec)
        } else {
            false
        };

        let key = (bits_vec.clone(), effective_additive);
        let entry = pattern_map.entry(key).or_insert_with(|| {
            let pid = patterns.len();
            patterns.push(PatternInfo {
                total_weight: 0,
                strategy: AdditiveStrategy::FastRange,
                leaf_bits: vec![0; ntax],
                leaf_ranges: vec![util::AdditiveRange::default(); ntax],
            });
            (pid, 0u64)
        });
        entry.1 += cs.weight as u64;
    }

    /* Populate the concrete leaf caches for each compressed pattern. */
    let mut active_patterns = Vec::new();
    for ((bits_vec, additive), (pid, total_weight)) in pattern_map.into_iter() {
        /*
          Choose exact costs only when non-contiguous additive states make
          interval merging insufficient.
        */
        let all_contig = bits_vec.iter().all(|&b| is_contiguous(b));
        let strategy = if additive && !all_contig {
            AdditiveStrategy::ExactCosts
        } else {
            AdditiveStrategy::FastRange
        };

        /* Fill leaf caches used by side-message scoring. */
        for t in 0..ntax {
            patterns[pid].leaf_bits[t] = bits_vec[t];
            patterns[pid].leaf_ranges[t] = util::AdditiveRange::from_bits(bits_vec[t]);
        }
        patterns[pid].total_weight = total_weight;
        patterns[pid].strategy = strategy;

        active_patterns.push(ActivePattern { pattern_id: pid, total_weight, additive });
    }

    /* High-weight patterns first improve limit-based pruning. */
    active_patterns.sort_by(|a, b| {
        b.total_weight
            .cmp(&a.total_weight)
            .then_with(|| b.additive.cmp(&a.additive))
            .then_with(|| a.pattern_id.cmp(&b.pattern_id))
    });

    (active_patterns, patterns)
}

/*  Side message – per pattern, unweighted  */

#[derive(Debug, Clone)]
struct SidePatternInfo {
    raw_len: u64,
    st: u64,
    rg: util::AdditiveRange,
    costs: Option<Rc<[u64; NSTATES]>>,
}

#[derive(Debug, Clone)]
struct SideMessage {
    infos: Vec<SidePatternInfo>,
    weighted_len: u64,
}

impl SideMessage {
    fn compute_weighted_len(&self, active: &[ActivePattern], patterns: &[PatternInfo]) -> u64 {
        self.infos
            .iter()
            .enumerate()
            .map(|(i, info)| info.raw_len * patterns[active[i].pattern_id].total_weight)
            .sum()
    }
}

type SideMessageRef = Rc<SideMessage>;

struct SideMessageMemo {
    nnode: usize,
    slots: Vec<Option<SideMessageRef>>,
}

impl SideMessageMemo {
    /// allocates a dense directed-edge memo table for side messages.
    fn new(nnode: usize) -> Self {
        Self { nnode, slots: vec![None; nnode.saturating_mul(nnode)] }
    }
    /// maps a directed `(from, blocked)` edge to a memo slot.

    #[inline(always)]
    fn idx(&self, from: usize, blocked: usize) -> usize {
        from * self.nnode + blocked
    }
    /// retrieves a cached side message by directed edge.

    fn get(&self, from: usize, blocked: usize) -> Option<SideMessageRef> {
        self.slots[self.idx(from, blocked)].as_ref().map(Rc::clone)
    }
    /// stores a computed side message in the memo table.

    fn put(&mut self, from: usize, blocked: usize, msg: SideMessageRef) {
        let idx = self.idx(from, blocked);
        self.slots[idx] = Some(msg);
    }

    /// Drops cut-specific values while retaining the allocated dense table.
    fn clear(&mut self) {
        self.slots.fill(None);
    }
}
/// creates a side message for a terminal taxon.

fn leaf_side_message(
    taxon: usize,
    active: &[ActivePattern],
    patterns: &[PatternInfo],
) -> SideMessage {
    let mut infos = Vec::with_capacity(active.len());
    for ap in active {
        let pat = &patterns[ap.pattern_id];
        let info = match pat.strategy {
            AdditiveStrategy::FastRange if ap.additive => {
                let rg = pat.leaf_ranges[taxon];
                SidePatternInfo { raw_len: 0, st: 0, rg, costs: None }
            }
            AdditiveStrategy::ExactCosts => {
                let bits = pat.leaf_bits[taxon];
                let mut arr = [INF; NSTATES];
                for s in 0..NSTATES {
                    if (bits >> s) & 1 == 1 {
                        arr[s] = 0;
                    }
                }
                SidePatternInfo {
                    raw_len: 0,
                    st: 0,
                    rg: util::AdditiveRange::default(),
                    costs: Some(Rc::new(arr)),
                }
            }
            _ => {
                let st = pat.leaf_bits[taxon];
                SidePatternInfo { raw_len: 0, st, rg: util::AdditiveRange::default(), costs: None }
            }
        };
        infos.push(info);
    }
    SideMessage { infos, weighted_len: 0 }
}
/// combines two exact ordered-state Sankoff arrays and returns the minimum
/// merge cost.

fn merge_two_cost_arrays(a: &[u64; NSTATES], b: &[u64; NSTATES]) -> ([u64; NSTATES], u64) {
    /*
      Each child contributes min_k (child[k] + |s-k|) independently.
    */
    let a_min = min_cost_with_transition(a.as_slice());
    let b_min = min_cost_with_transition(b.as_slice());

    let mut merged = [INF; NSTATES];
    let mut min_total = INF;
    for s in 0..NSTATES {
        let v = a_min[s].saturating_add(b_min[s]);
        merged[s] = v;
        min_total = min_total.min(v);
    }
    (merged, min_total)
}
/// merges child side messages into a parent side message over every active
/// compressed pattern.

fn merge_side_messages(
    active: &[ActivePattern],
    patterns: &[PatternInfo],
    children: &[&SideMessage],
) -> SideMessage {
    let n = active.len();
    let mut merged_infos = Vec::with_capacity(n);

    for (i, ap) in active.iter().enumerate() {
        let pat = &patterns[ap.pattern_id];
        let first = &children[0].infos[i];
        let mut info = first.clone();

        for child in &children[1..] {
            let ci = &child.infos[i];
            match pat.strategy {
                AdditiveStrategy::FastRange if ap.additive => {
                    let (new_rg, gap) = info.rg.merge(ci.rg);
                    info.rg = new_rg;
                    info.raw_len += ci.raw_len + gap;
                }
                AdditiveStrategy::ExactCosts => {
                    let my_costs = info.costs.as_ref().unwrap();
                    let other_costs = ci.costs.as_ref().unwrap();
                    let (merged_costs, min) = merge_two_cost_arrays(my_costs, other_costs);
                    info.raw_len = min;
                    info.costs = Some(Rc::new(merged_costs));
                }
                _ => {
                    /* Non-additive Fitch merge. */
                    let (merged_bits, step) = merge_bits(info.st, ci.st);
                    info.st = merged_bits;
                    info.raw_len += ci.raw_len + step;
                }
            }
        }
        merged_infos.push(info);
    }

    let mut msg = SideMessage { infos: merged_infos, weighted_len: 0 };
    msg.weighted_len = msg.compute_weighted_len(active, patterns);
    msg
}
/// recursively computes and memoizes the message for one side of a blocked
/// edge.

fn compute_side_message_memo(
    tree: &Tree,
    active: &[ActivePattern],
    patterns: &[PatternInfo],
    neigh: &[Vec<usize>],
    from: usize,
    blocked: usize,
    cut_edge: (usize, usize),
    memo: &mut SideMessageMemo,
) -> SideMessageRef {
    if let Some(cached) = memo.get(from, blocked) {
        return cached;
    }

    if let Some(t) = tree.nodes[from].taxon {
        let msg = leaf_side_message(t, active, patterns);
        let rc = Rc::new(msg);
        memo.put(from, blocked, Rc::clone(&rc));
        return rc;
    }

    let mut children_refs = Vec::new();
    for &nb in &neigh[from] {
        if nb == blocked || is_same_edge(from, nb, cut_edge) {
            continue;
        }
        children_refs.push(compute_side_message_memo(
            tree, active, patterns, neigh, nb, from, cut_edge, memo,
        ));
    }

    let msg = if children_refs.is_empty() {
        /*
          Isolated internal nodes accept every state at zero cost for Fitch and
          interval scoring; exact-cost scoring uses INF to mark invalid state
          arrays.
        */
        let mut infos = Vec::with_capacity(active.len());
        for ap in active {
            let pat = &patterns[ap.pattern_id];
            let info = match pat.strategy {
                AdditiveStrategy::FastRange if ap.additive => SidePatternInfo {
                    raw_len: 0,
                    st: 0,
                    rg: util::AdditiveRange { lo: 0, hi: 35 },
                    costs: None,
                },
                AdditiveStrategy::ExactCosts => {
                    let costs = Rc::new([INF; NSTATES]);
                    SidePatternInfo {
                        raw_len: INF, /* Unreachable cost marks this option invalid. */
                        st: 0,
                        rg: util::AdditiveRange::default(),
                        costs: Some(costs),
                    }
                }
                _ => SidePatternInfo {
                    raw_len: 0,
                    st: crate::engines::dataset::StateSet::ALL36.bits(),
                    rg: util::AdditiveRange::default(),
                    costs: None,
                },
            };
            infos.push(info);
        }
        let mut msg = SideMessage { infos, weighted_len: 0 };
        msg.weighted_len = msg.compute_weighted_len(active, patterns);
        msg
    } else {
        let refs: Vec<&SideMessage> = children_refs.iter().map(|rc| rc.as_ref()).collect();
        merge_side_messages(active, patterns, &refs)
    };

    let rc = Rc::new(msg);
    memo.put(from, blocked, Rc::clone(&rc));
    rc
}

/* TBR option structures */

#[derive(Debug, Clone)]
struct TbrInsertOption {
    /* None denotes a one-node (normally terminal) component, which reconnects
    directly instead of subdividing an edge. */
    edge: Option<(usize, usize)>,
    reconnect: (usize, usize),
    base_cost: u64,
    a_to_b: SideMessageRef,
    b_to_a: SideMessageRef,
}

impl TbrInsertOption {
    /// returns the local side-message length used to rank component insertion
    /// options.
    #[inline(always)]
    fn base_len(&self) -> u64 {
        self.base_cost
    }
}

#[derive(Debug, Clone)]
struct TbrComponentInfo {
    left_edges: Vec<(usize, usize)>,
    right_edges: Vec<(usize, usize)>,
    left_options: Vec<TbrInsertOption>,
    right_options: Vec<TbrInsertOption>,
}

#[derive(Debug, Clone, Copy)]
struct TbrRankedCandidate {
    delta_len: u64,
    left_idx: usize,
    right_idx: usize,
}

/* TBR search context */

struct TbrSearchContext {
    neigh: Vec<Vec<usize>>,
    edges: Vec<(usize, usize)>,
    queue: VecDeque<usize>,
    mark: Vec<u8>,
    side_memo: SideMessageMemo,
    build_edges: Vec<(usize, usize)>,
}

impl TbrSearchContext {
    /// prepares adjacency, active patterns, and side-message memoization for a
    /// source tree.
    fn new(tree: &Tree) -> Self {
        Self {
            neigh: undirected_neighbors(tree),
            edges: tree.undirected_edges(),
            queue: VecDeque::with_capacity(tree.nodes.len()),
            mark: vec![0u8; tree.nodes.len()],
            side_memo: SideMessageMemo::new(tree.nodes.len()),
            build_edges: Vec::with_capacity(tree.nodes.len() + 4),
        }
    }
    /// lists undirected edges inside one component after removing the cut edge.

    fn component_edges_after_cut(
        &mut self,
        cut_edge: (usize, usize),
    ) -> (Vec<(usize, usize)>, Vec<(usize, usize)>) {
        self.mark.fill(0);
        bfs_mark_component(
            &self.neigh,
            cut_edge.0,
            cut_edge.0,
            cut_edge.1,
            1,
            &mut self.mark,
            &mut self.queue,
        );
        bfs_mark_component(
            &self.neigh,
            cut_edge.1,
            cut_edge.0,
            cut_edge.1,
            2,
            &mut self.mark,
            &mut self.queue,
        );

        let mut left_edges = Vec::new();
        let mut right_edges = Vec::new();
        for &e in &self.edges {
            if e == cut_edge {
                continue;
            }
            let (a, b) = e;
            if self.mark[a] == 1 && self.mark[b] == 1 {
                left_edges.push(e);
            } else if self.mark[a] == 2 && self.mark[b] == 2 {
                right_edges.push(e);
            }
        }
        (left_edges, right_edges)
    }
}
/// marks one connected component while ignoring a blocked cut edge.

fn bfs_mark_component(
    neigh: &[Vec<usize>],
    start: usize,
    cut_a: usize,
    cut_b: usize,
    value: u8,
    mark: &mut [u8],
    q: &mut VecDeque<usize>,
) {
    q.clear();
    mark[start] = value;
    q.push_back(start);
    while let Some(u) = q.pop_front() {
        for &v in &neigh[u] {
            if is_same_edge(u, v, (cut_a, cut_b)) {
                continue;
            }
            if mark[v] == 0 {
                mark[v] = value;
                q.push_back(v);
            }
        }
    }
}

/*  Component options construction  */
/// builds component anchors and insertion options for one TBR cut.

fn make_tbr_component_info_for_cut(
    tree: &Tree,
    active: &[ActivePattern],
    patterns: &[PatternInfo],
    outgroup: usize,
    ctx: &mut TbrSearchContext,
    cut_edge: (usize, usize),
) -> Option<TbrComponentInfo> {
    if cut_edge.0 == outgroup || cut_edge.1 == outgroup {
        return None;
    }
    let (left_edges, right_edges) = ctx.component_edges_after_cut(cut_edge);
    ctx.side_memo.clear();
    let left_options = if left_edges.is_empty() {
        vec![make_singleton_option(
            tree,
            active,
            patterns,
            &ctx.neigh,
            &mut ctx.side_memo,
            cut_edge.0,
            cut_edge.1,
            cut_edge,
        )]
    } else {
        make_insert_options_for_component(
            tree,
            active,
            patterns,
            &ctx.neigh,
            &mut ctx.side_memo,
            &left_edges,
            cut_edge,
        )
    };
    let right_options = if right_edges.is_empty() {
        vec![make_singleton_option(
            tree,
            active,
            patterns,
            &ctx.neigh,
            &mut ctx.side_memo,
            cut_edge.1,
            cut_edge.0,
            cut_edge,
        )]
    } else {
        make_insert_options_for_component(
            tree,
            active,
            patterns,
            &ctx.neigh,
            &mut ctx.side_memo,
            &right_edges,
            cut_edge,
        )
    };

    Some(TbrComponentInfo { left_edges, right_edges, left_options, right_options })
}
/// turns component edges into ranked places where a reconnected edge can be
/// inserted.

fn make_insert_options_for_component(
    tree: &Tree,
    active: &[ActivePattern],
    patterns: &[PatternInfo],
    neigh: &[Vec<usize>],
    memo: &mut SideMessageMemo,
    component_edges: &[(usize, usize)],
    cut_edge: (usize, usize),
) -> Vec<TbrInsertOption> {
    let mut out = Vec::with_capacity(component_edges.len());

    for &(a, b) in component_edges {
        let edge = sorted_edge(a, b);
        let a_to_b = compute_side_message_memo(tree, active, patterns, neigh, a, b, cut_edge, memo);
        let b_to_a = compute_side_message_memo(tree, active, patterns, neigh, b, a, cut_edge, memo);
        let base_cost = a_to_b.weighted_len + b_to_a.weighted_len;

        out.push(TbrInsertOption { edge: Some(edge), reconnect: edge, base_cost, a_to_b, b_to_a });
    }
    out
}

/// Builds a reconnection option for a single-node TBR component.
fn make_singleton_option(
    tree: &Tree,
    active: &[ActivePattern],
    patterns: &[PatternInfo],
    neigh: &[Vec<usize>],
    memo: &mut SideMessageMemo,
    node: usize,
    blocked: usize,
    cut_edge: (usize, usize),
) -> TbrInsertOption {
    let real =
        compute_side_message_memo(tree, active, patterns, neigh, node, blocked, cut_edge, memo);
    let neutral = Rc::new(neutral_side_message(active, patterns));
    TbrInsertOption {
        edge: None,
        reconnect: (node, node),
        base_cost: real.weighted_len,
        a_to_b: real,
        b_to_a: neutral,
    }
}

fn neutral_side_message(active: &[ActivePattern], patterns: &[PatternInfo]) -> SideMessage {
    let infos = active
        .iter()
        .map(|ap| match patterns[ap.pattern_id].strategy {
            AdditiveStrategy::FastRange if ap.additive => SidePatternInfo {
                raw_len: 0,
                st: 0,
                rg: util::AdditiveRange { lo: 0, hi: 35 },
                costs: None,
            },
            AdditiveStrategy::ExactCosts => SidePatternInfo {
                raw_len: 0,
                st: 0,
                rg: util::AdditiveRange::default(),
                costs: Some(Rc::new([0; NSTATES])),
            },
            _ => SidePatternInfo {
                raw_len: 0,
                st: crate::engines::dataset::StateSet::ALL36.bits(),
                rg: util::AdditiveRange::default(),
                costs: None,
            },
        })
        .collect();
    SideMessage { infos, weighted_len: 0 }
}

/*  Candidate scoring */
/// Scores a paired TBR reconnection up to `limit`.

fn score_tbr_candidate_limited(
    active: &[ActivePattern],
    patterns: &[PatternInfo],
    left: &TbrInsertOption,
    right: &TbrInsertOption,
    limit: u64,
) -> Option<u64> {
    let mut total = 0u64;

    for (i, ap) in active.iter().enumerate() {
        let pat = &patterns[ap.pattern_id];
        let w = pat.total_weight;

        let base = left.a_to_b.infos[i].raw_len
            + left.b_to_a.infos[i].raw_len
            + right.a_to_b.infos[i].raw_len
            + right.b_to_a.infos[i].raw_len;

        let raw = match pat.strategy {
            AdditiveStrategy::FastRange if ap.additive => {
                let (left_rg, left_gap) = left.a_to_b.infos[i].rg.merge(left.b_to_a.infos[i].rg);
                let (right_rg, right_gap) =
                    right.a_to_b.infos[i].rg.merge(right.b_to_a.infos[i].rg);
                let (_root_rg, bridge_gap) = left_rg.merge(right_rg);
                base + left_gap + right_gap + bridge_gap
            }
            AdditiveStrategy::ExactCosts => {
                let left_costs = merge_two_cost_arrays(
                    left.a_to_b.infos[i].costs.as_ref().unwrap(),
                    left.b_to_a.infos[i].costs.as_ref().unwrap(),
                )
                .0;
                let right_costs = merge_two_cost_arrays(
                    right.a_to_b.infos[i].costs.as_ref().unwrap(),
                    right.b_to_a.infos[i].costs.as_ref().unwrap(),
                )
                .0;
                let (_root_costs, root_min) = merge_two_cost_arrays(&left_costs, &right_costs);
                root_min
            }
            _ => {
                let (left_st, left_step) =
                    merge_bits(left.a_to_b.infos[i].st, left.b_to_a.infos[i].st);
                let (right_st, right_step) =
                    merge_bits(right.a_to_b.infos[i].st, right.b_to_a.infos[i].st);
                let (_root_st, bridge_step) = merge_bits(left_st, right_st);
                base + left_step + right_step + bridge_step
            }
        };

        total = total.saturating_add(raw.saturating_mul(w));
        if total > limit {
            return None;
        }
    }
    Some(total)
}
/// Scores all TBR reconnections that can meet `score_limit`.
fn exhaustive_candidates_for_cut(
    active: &[ActivePattern],
    patterns: &[PatternInfo],
    info: &TbrComponentInfo,
    score_limit: u64,
) -> Vec<TbrRankedCandidate> {
    let left_indices = sorted_option_indices(&info.left_options, None);
    let right_indices = sorted_option_indices(&info.right_options, None);

    let mut candidates =
        Vec::with_capacity(info.left_options.len().saturating_mul(info.right_options.len()));
    for li in left_indices {
        let left = &info.left_options[li];
        if left.base_len() > score_limit {
            break;
        }
        for &ri in &right_indices {
            let right = &info.right_options[ri];
            if left.base_len().saturating_add(right.base_len()) > score_limit {
                break;
            }
            if let Some(delta_len) =
                score_tbr_candidate_limited(active, patterns, left, right, score_limit)
            {
                candidates.push(TbrRankedCandidate { delta_len, left_idx: li, right_idx: ri });
            }
        }
    }
    candidates.sort_unstable_by_key(|c| (c.delta_len, c.left_idx, c.right_idx));
    candidates
}

/// Returns ranked TBR candidates for the first-improvement pass.
fn fast_candidates_for_cut(
    active: &[ActivePattern],
    patterns: &[PatternInfo],
    info: &TbrComponentInfo,
    score_limit: u64,
) -> Vec<TbrRankedCandidate> {
    let left = sorted_option_indices(&info.left_options, Some(FAST_RECONNECTS_PER_SIDE));
    let right = sorted_option_indices(&info.right_options, Some(FAST_RECONNECTS_PER_SIDE));

    let mut out =
        Vec::with_capacity(left.len().saturating_mul(right.len()).min(FAST_CANDIDATES_PER_CUT));
    for li in left {
        for &ri in &right {
            let l = &info.left_options[li];
            let r = &info.right_options[ri];
            if l.base_len().saturating_add(r.base_len()) > score_limit {
                continue;
            }
            if let Some(delta_len) =
                score_tbr_candidate_limited(active, patterns, l, r, score_limit)
            {
                out.push(TbrRankedCandidate { delta_len, left_idx: li, right_idx: ri });
            }
        }
    }
    out.sort_unstable_by_key(|c| (c.delta_len, c.left_idx, c.right_idx));
    out.truncate(FAST_CANDIDATES_PER_CUT);
    out
}

fn sorted_option_indices(options: &[TbrInsertOption], limit: Option<usize>) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..options.len()).collect();
    indices.sort_unstable_by_key(|&index| options[index].base_len());
    if let Some(limit) = limit {
        indices.truncate(limit);
    }
    indices
}
/* Public entry points */
/// repeatedly explores equal/better TBR neighbors until closure.

pub fn branch_break_closure(
    ds: &Dataset,
    cfg: &CharConfig,
    start_trees: &[Tree],
    outgroup: usize,
    keep_limit: Option<usize>,
) -> Result<(Vec<Tree>, u64), String> {
    /*
      Greedy TBR pre-optimisation on every start tree.
      `branch_swap_best_tree` uses TBR , and lifting each
      Wagner-built tree to its local optimum gives the subsequent closure
      phase a much richer frontier to explore.
    */
    let valid_starts: Vec<&Tree> = start_trees.iter().filter(|tr| is_strictly_binary(tr)).collect();
    let improved: Result<Vec<(Tree, u64)>, String> = std::thread::scope(|scope| {
        let handles: Vec<_> = valid_starts
            .iter()
            .map(|tr| scope.spawn(move || branch_swap_best_tree_tbr_fast(ds, cfg, tr, outgroup)))
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().map_err(|_| "branch breaking worker panicked".to_string())?)
            .collect()
    });

    let mut pre_opt: Vec<Tree> = Vec::new();
    let mut seen_pre = HashSet::<u64>::new();
    for (improved, _) in improved? {
        let h = rooted_topology_hash(&improved);
        if seen_pre.insert(h) {
            pre_opt.push(improved);
        }
    }

    branch_break_closure_tbr_fast(ds, cfg, &pre_opt, outgroup, keep_limit, false)
}

/// Improves one start tree by TBR and retains one local optimum.
pub fn branch_swap_one(
    ds: &Dataset,
    cfg: &CharConfig,
    start: &Tree,
    outgroup: usize,
) -> Result<(Tree, u64), String> {
    branch_swap_best_tree_tbr_fast(ds, cfg, start, outgroup)
}

/*  Fast TBR best-tree search  */
/// evaluates a ranked fast-TBR neighborhood and returns the best candidate.

fn branch_swap_best_tree_tbr_fast(
    ds: &Dataset,
    cfg: &CharConfig,
    start: &Tree,
    outgroup: usize,
) -> Result<(Tree, u64), String> {
    let (active_patterns, patterns) = build_active_patterns(ds, cfg);
    let mut score_ws = ScoreWorkspace::new();

    let mut current = normalize_tree(start, outgroup);
    let mut current_len = score_tree(ds, cfg, &current, &mut score_ws);

    loop {
        let mut improved = false;

        /* Two passes: a small first-improvement pass, followed only when
        necessary by an exhaustive proof that no improving reconnect was
        omitted by the fast ranking. */
        for exhaustive in [false, true] {
            let mut ctx = TbrSearchContext::new(&current);
            let cut_edges = ctx.edges.clone();
            for cut_edge in cut_edges {
                let Some(info) = make_tbr_component_info_for_cut(
                    &current,
                    &active_patterns,
                    &patterns,
                    outgroup,
                    &mut ctx,
                    cut_edge,
                ) else {
                    continue;
                };

                let limit = current_len.saturating_sub(1);
                let ranked = if exhaustive {
                    exhaustive_candidates_for_cut(&active_patterns, &patterns, &info, limit)
                } else {
                    fast_candidates_for_cut(&active_patterns, &patterns, &info, limit)
                };
                if ranked.is_empty() {
                    continue;
                }

                for best in &ranked {
                    let left_opt = &info.left_options[best.left_idx];
                    let right_opt = &info.right_options[best.right_idx];
                    let cand_tree = match build_tbr_tree_from_options(
                        &current,
                        outgroup,
                        &info,
                        left_opt,
                        right_opt,
                        &mut ctx.build_edges,
                    ) {
                        Ok(t) if is_strictly_binary(&t) => t,
                        _ => continue,
                    };
                    let len = best.delta_len;
                    #[cfg(debug_assertions)]
                    {
                        let checked = score_tree(ds, cfg, &cand_tree, &mut score_ws);
                        debug_assert_eq!(len, checked, "indirect TBR score mismatch");
                    }
                    if len < current_len {
                        current = canonicalize_outgrouped_candidate(cand_tree);
                        current_len = len;
                        improved = true;
                        break;
                    }
                }
                if improved {
                    break;
                }
            }
            if improved {
                break;
            }
        }

        if !improved {
            return Ok((current, current_len));
        }
    }
}

/*  Fast TBR closure search  */
/// repeatedly expands best TBR candidates while maintaining a deduplicated
/// best-tree frontier.

fn branch_break_closure_tbr_fast(
    ds: &Dataset,
    cfg: &CharConfig,
    start_trees: &[Tree],
    outgroup: usize,
    keep_limit: Option<usize>,
    collapsed_frontier: bool,
) -> Result<(Vec<Tree>, u64), String> {
    let (active_patterns, patterns) = build_active_patterns(ds, cfg);
    let mut score_ws = ScoreWorkspace::new();

    let mut best_len = u64::MAX;
    let mut best_seen = HashSet::<u64>::new();
    let mut best_trees = Vec::<Tree>::new();
    /* Binary resolutions with the same collapsed report can have different
    improving TBR paths, so frontier identity must remain rooted/binary. */
    let mut frontier_seen = HashSet::<u64>::new();
    let mut frontier = Vec::<Tree>::new();
    let mut rooted_seen = HashSet::<u64>::new();
    for tr in start_trees {
        if !is_strictly_binary(tr) {
            continue;
        }
        let norm = normalize_tree(tr, outgroup);
        let h_coll = collapsed_dedup_hash(&norm, ds, cfg, &mut score_ws);
        let h_root = rooted_topology_hash(&norm);
        if !frontier_seen.insert(h_root) {
            continue;
        }
        let len = score_tree(ds, cfg, &norm, &mut score_ws);

        if len < best_len {
            best_len = len;
            best_seen.clear();
            best_trees.clear();
            frontier.clear();
            frontier_seen.clear();
            rooted_seen.clear();

            best_seen.insert(h_coll);
            frontier_seen.insert(h_root);
            best_trees.push(norm.clone());
            frontier.push(norm);

            if best_buffer_full(best_trees.len(), keep_limit) {
                return Ok((best_trees, best_len));
            }
        } else if len == best_len {
            if best_seen.insert(h_coll) {
                best_trees.push(norm.clone());
                if best_buffer_full(best_trees.len(), keep_limit) {
                    return Ok((best_trees, best_len));
                }
            }
            frontier.push(norm);
        }
    }

    if best_trees.is_empty() {
        return Err("branch_break_closure: no valid binary start trees".to_string());
    }

    loop {
        let mut next_frontier = Vec::new();
        let mut new_found = false;
        let mut buffer_full = false;

        for tr in &frontier {
            if buffer_full {
                break;
            }
            let mut ctx = TbrSearchContext::new(tr);
            let cut_edges = ctx.edges.clone();

            for cut_edge in cut_edges {
                if buffer_full {
                    break;
                }
                let Some(info) = make_tbr_component_info_for_cut(
                    tr,
                    &active_patterns,
                    &patterns,
                    outgroup,
                    &mut ctx,
                    cut_edge,
                ) else {
                    continue;
                };

                /* Both commands traverse the complete equal-tree closure.
                `bb` stops only after its 100-tree buffer fills; `bb*`
                removes that storage limit. */
                let ranked =
                    exhaustive_candidates_for_cut(&active_patterns, &patterns, &info, best_len);

                for rc in &ranked {
                    if buffer_full {
                        break;
                    }
                    let delta = rc.delta_len;
                    if delta > best_len {
                        continue;
                    }
                    let left_opt = &info.left_options[rc.left_idx];
                    let right_opt = &info.right_options[rc.right_idx];
                    let cand_tree = match build_tbr_tree_from_options(
                        tr,
                        outgroup,
                        &info,
                        left_opt,
                        right_opt,
                        &mut ctx.build_edges,
                    ) {
                        Ok(t) if is_strictly_binary(&t) => t,
                        _ => continue,
                    };
                    /* Deduplicate and score before canonical child sorting. */
                    let h_root = rooted_topology_hash(&cand_tree);

                    /*
                      Cheap rooted-hash prefilter: skip if this binary topology
                      has already been seen in the current frontier iteration.
                    */
                    if !rooted_seen.insert(h_root) {
                        continue;
                    }

                    /* The side-message calculation is the exact indirect TBR
                    length: fixed component lengths plus reconnect cost.
                    Debug builds cross-check it with a full down-pass; release
                    builds avoid that O(T*C) duplicate work. */
                    let exact_len = rc.delta_len;
                    #[cfg(debug_assertions)]
                    {
                        let checked = score_tree(ds, cfg, &cand_tree, &mut score_ws);
                        debug_assert_eq!(exact_len, checked, "indirect TBR score mismatch");
                    }

                    if exact_len > best_len {
                        continue;
                    }

                    let norm = canonicalize_outgrouped_candidate(cand_tree);
                    let h_coll = collapsed_dedup_hash(&norm, ds, cfg, &mut score_ws);

                    if exact_len < best_len {
                        best_len = exact_len;
                        best_seen.clear();
                        best_trees.clear();
                        next_frontier.clear();
                        frontier_seen.clear();
                        rooted_seen.clear();
                        rooted_seen.insert(h_root);

                        best_seen.insert(h_coll);
                        best_trees.push(norm.clone());
                        frontier_seen.insert(h_root);
                        next_frontier.push(norm);
                        new_found = true;

                        if best_buffer_full(best_trees.len(), keep_limit) {
                            buffer_full = true;
                        }
                        continue;
                    }

                    /* equal length */
                    let mut useful = false;
                    if best_seen.insert(h_coll) {
                        best_trees.push(norm.clone());
                        useful = true;
                        if best_buffer_full(best_trees.len(), keep_limit) {
                            buffer_full = true;
                        }
                    }
                    if frontier_seen.insert(h_root) {
                        next_frontier.push(norm);
                        useful = true;
                    }
                    if useful {
                        new_found = true;
                    }
                }
            }
        }

        if buffer_full {
            return Ok((best_trees, best_len));
        }
        if !new_found {
            break;
        }
        frontier = next_frontier;
        if frontier.is_empty() {
            break;
        }
    }

    /* First pass preserves every rooted binary route that may improve TL.
    Once it has finished, bounded bb may use a collapsed-MPT frontier only
    to collect additional equal-length reported trees up to its 100 slots.
    This cannot discard a shorter route already considered by the first
    pass; if a shorter tree is nevertheless found, it is adopted. */
    if !collapsed_frontier
        && keep_limit.is_some()
        && !best_buffer_full(best_trees.len(), keep_limit)
    {
        return branch_break_closure_tbr_fast(ds, cfg, &best_trees, outgroup, keep_limit, true);
    }

    Ok((best_trees, best_len))
}
/// constructs a TBR candidate from two ranked component insertion options.

fn build_tbr_tree_from_options(
    tree: &Tree,
    outgroup: usize,
    info: &TbrComponentInfo,
    left_opt: &TbrInsertOption,
    right_opt: &TbrInsertOption,
    edges: &mut Vec<(usize, usize)>,
) -> Result<Tree, String> {
    edges.clear();
    edges.reserve(
        (info.left_edges.len() + info.right_edges.len() + 5).saturating_sub(edges.capacity()),
    );
    for &e in &info.left_edges {
        if Some(e) != left_opt.edge {
            edges.push(e);
        }
    }
    for &e in &info.right_edges {
        if Some(e) != right_opt.edge {
            edges.push(e);
        }
    }
    let mut next = next_internal_node_id(edges, tree.ntax);
    let left_anchor = if left_opt.edge.is_some() {
        let id = next;
        next += 1;
        let (a, b) = left_opt.reconnect;
        edges.push(sorted_edge(a, id));
        edges.push(sorted_edge(b, id));
        id
    } else {
        left_opt.reconnect.0
    };
    let right_anchor = if right_opt.edge.is_some() {
        let id = next;
        let (a, b) = right_opt.reconnect;
        edges.push(sorted_edge(a, id));
        edges.push(sorted_edge(b, id));
        id
    } else {
        right_opt.reconnect.0
    };
    edges.push(sorted_edge(left_anchor, right_anchor));
    edges.sort_unstable();
    edges.dedup();
    Tree::from_undirected_edges_with_outgroup(tree.ntax, edges, outgroup)
}

/* Helpers  */
/// reroots and canonicalizes a candidate by the outgroup.

fn normalize_tree(tree: &Tree, outgroup: usize) -> Tree {
    let mut norm = tree.reroot_by_outgroup(outgroup).unwrap_or_else(|_| tree.clone());
    norm.canonicalize();
    norm
}

/// TBR candidate construction already roots on the chosen outgroup; avoid a
/// redundant full reroot before deterministic child ordering.
fn canonicalize_outgrouped_candidate(mut tree: Tree) -> Tree {
    tree.canonicalize();
    tree
}
/// tests whether every internal node has exactly two children.

fn is_strictly_binary(tree: &Tree) -> bool {
    BinaryView::from_tree(tree).is_ok()
}
/// checks whether the best-tree buffer reached its optional storage limit.

fn best_buffer_full(len: usize, keep_limit: Option<usize>) -> bool {
    keep_limit.map(|limit| len >= limit).unwrap_or(false)
}
/// builds an adjacency list from a rooted tree.

fn undirected_neighbors(tree: &Tree) -> Vec<Vec<usize>> {
    let mut g = vec![Vec::new(); tree.nodes.len()];
    for (u, node) in tree.nodes.iter().enumerate() {
        for &v in &node.children {
            g[u].push(v);
            g[v].push(u);
        }
    }
    g
}
/// compares an edge to two endpoints ignoring direction.

#[inline(always)]
fn is_same_edge(a: usize, b: usize, e: (usize, usize)) -> bool {
    (a == e.0 && b == e.1) || (a == e.1 && b == e.0)
}
/// returns the next unused internal node id in an edge list.

fn next_internal_node_id(edges: &[(usize, usize)], ntax: usize) -> usize {
    let max_id = edges.iter().flat_map(|&(a, b)| [a, b]).max().unwrap_or(ntax.saturating_sub(1));
    max_id + 1
}
/// normalizes an undirected edge endpoint order.

fn sorted_edge(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
}

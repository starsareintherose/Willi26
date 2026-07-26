/*!Module: Nelson/strict consensus support, including clade extraction,
outgroup normalization, duplicate collapsed-topology removal, and consensus
tree scoring.

 */
use std::collections::{BTreeSet, HashSet};

use crate::engines::{
    ccode::CharConfig,
    dataset::Dataset,
    search::{
        ScoreWorkspace, calc_ci, calc_ri, collapse_zero_length_branches, collapsed_topology_hash,
        maxsteps_sum, minsteps_sum, score_tree,
    },
    trees::Tree,
};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct TaxonSet {
    bits: Vec<u64>,
}

impl TaxonSet {
    /// creates an empty fixed-width bitset for a dataset taxon count.
    fn new(ntax: usize) -> Self {
        let words = ntax.div_ceil(64);
        Self { bits: vec![0; words] }
    }
    /// creates a one-taxon set.

    fn singleton(ntax: usize, t: usize) -> Self {
        let mut s = Self::new(ntax);
        s.insert(t);
        s
    }
    /// adds a taxon to the set.

    fn insert(&mut self, t: usize) {
        let w = t / 64;
        let b = t % 64;
        self.bits[w] |= 1u64 << b;
    }
    /// unions another taxon set into this one.

    fn union_inplace(&mut self, other: &Self) {
        for (a, b) in self.bits.iter_mut().zip(&other.bits) {
            *a |= *b;
        }
    }
    /// tests membership for one taxon.

    fn contains(&self, t: usize) -> bool {
        let w = t / 64;
        let b = t % 64;
        ((self.bits[w] >> b) & 1) == 1
    }
    /// tests whether every bit in this set is present in another set.

    fn is_subset_of(&self, other: &Self) -> bool {
        self.bits.iter().zip(&other.bits).all(|(a, b)| (*a & !*b) == 0)
    }
    /// compares two taxon sets exactly.

    fn is_equal(&self, other: &Self) -> bool {
        self.bits == other.bits
    }
    /// tests whether two clades share any taxon.

    fn intersects(&self, other: &Self) -> bool {
        self.bits.iter().zip(&other.bits).any(|(a, b)| (*a & *b) != 0)
    }
    /// checks clade compatibility by nesting or disjointness.

    fn is_compatible_with(&self, other: &Self) -> bool {
        self.is_subset_of(other) || other.is_subset_of(self) || !self.intersects(other)
    }
    /// counts selected taxa.

    fn count(&self) -> usize {
        self.bits.iter().map(|w| w.count_ones() as usize).sum()
    }
    /// expands the bitset into sorted taxon indices.

    fn taxa(&self) -> Vec<usize> {
        let mut out = Vec::new();
        for (wi, &word) in self.bits.iter().enumerate() {
            let mut w = word;
            while w != 0 {
                let b = w.trailing_zeros() as usize;
                out.push(wi * 64 + b);
                w &= w - 1;
            }
        }
        out
    }
}
/// recursively gathers rooted clades from a tree.

fn collect_clades_rec(tree: &Tree, v: usize, ntax: usize, out: &mut Vec<TaxonSet>) -> TaxonSet {
    let node = &tree.nodes[v];
    if let Some(t) = node.taxon {
        return TaxonSet::singleton(ntax, t);
    }

    let mut acc = TaxonSet::new(ntax);
    for &c in &node.children {
        let sub = collect_clades_rec(tree, c, ntax, out);
        acc.union_inplace(&sub);
    }

    out.push(acc.clone());
    acc
}
/// collects informative ingroup clades after excluding outgroups, whole-tree
/// clades, and singleton clades.

fn rooted_clades_excluding_outgroups(tree: &Tree, outgroups: &[usize]) -> Vec<TaxonSet> {
    let mut out = Vec::new();
    let whole = collect_clades_rec(tree, tree.root, tree.ntax, &mut out);
    /* Compute the full ingroup set so we can exclude it later. */
    let mut ingroup = TaxonSet::new(tree.ntax);
    for t in 0..tree.ntax {
        if !outgroups.iter().any(|&og| og == t) {
            ingroup.insert(t);
        }
    }

    out.into_iter()
        .filter(|c| {
            let k = c.count();
            k > 1
                && !c.is_equal(&whole)
                && !c.is_equal(&ingroup)
                && !outgroups.iter().any(|&og| c.contains(og))
        })
        .collect()
}
/// reroots and canonicalizes a tree according to the selected outgroup set.

fn normalize_by_outgroups(tree: &Tree, outgroups: &[usize]) -> Tree {
    let mut norm = tree
        .reroot_by_outgroup_set(outgroups)
        .or_else(|_| tree.reroot_by_outgroup(outgroups[0]))
        .unwrap_or_else(|_| tree.clone());
    norm.canonicalize();
    norm
}
/// removes duplicated trees after outgroup rooting and zero-length branch
/// collapse.

pub fn uniq_rooted_by_outgroups(
    trees: &[Tree],
    outgroups: &[usize],
    ds: &Dataset,
    cfg: &CharConfig,
    ws: &mut ScoreWorkspace,
) -> Vec<Tree> {
    let mut seen = HashSet::<u64>::new();
    let mut out = Vec::new();

    for tr in trees {
        let mut norm = normalize_by_outgroups(tr, outgroups);
        /*
          Collapse zero-length branches so trees differing only in
          unsupported resolution are treated as identical.
        */
        collapse_zero_length_branches(&mut norm, ds, cfg, ws);
        let h = collapsed_topology_hash(&norm);
        if seen.insert(h) {
            out.push(norm);
        }
    }

    out
}

/// Result of consensus tree construction: the consensus tree and a summary
/// message.
#[derive(Debug, Clone)]
pub struct ConsensusOutcome {
    pub tree: Tree,
    pub input_trees: usize,
    pub uniq_trees: usize,
    pub best_len: u64,
    pub ci: f64,
    pub ri: Option<f64>,
}

#[derive(Debug, Clone)]
struct ClusterNode {
    set: TaxonSet,
    parent: Option<usize>,
    children: Vec<usize>,
}
/// constructs a rooted consensus tree from compatible ingroup clades and
/// explicit outgroup children.

fn build_tree_from_rooted_clades_with_outgroups(
    ntax: usize,
    outgroups: &[usize],
    clades: &[TaxonSet],
) -> Result<Tree, String> {
    let mut ogset = TaxonSet::new(ntax);
    for &t in outgroups {
        ogset.insert(t);
    }

    let mut ingroup_set = TaxonSet::new(ntax);
    for t in 0..ntax {
        if !ogset.contains(t) {
            ingroup_set.insert(t);
        }
    }

    let mut clusters = Vec::<ClusterNode>::new();
    clusters.push(ClusterNode { set: ingroup_set.clone(), parent: None, children: Vec::new() });

    for c in clades {
        if outgroups.iter().any(|&og| c.contains(og)) {
            continue;
        }
        if c.count() <= 1 {
            continue;
        }
        if c.is_equal(&ingroup_set) {
            continue;
        }
        clusters.push(ClusterNode { set: c.clone(), parent: None, children: Vec::new() });
    }

    for i in 1..clusters.len() {
        let mut best_parent: Option<usize> = None;
        let mut best_size = usize::MAX;

        for j in 0..clusters.len() {
            if i == j {
                continue;
            }
            if clusters[i].set.is_subset_of(&clusters[j].set)
                && !clusters[i].set.is_equal(&clusters[j].set)
            {
                let sz = clusters[j].set.count();
                if sz < best_size {
                    best_size = sz;
                    best_parent = Some(j);
                }
            }
        }

        let p = best_parent.unwrap_or(0); /* fallback to full ingroup */
        clusters[i].parent = Some(p);
        clusters[p].children.push(i);
    }

    let mut tree = Tree::new(ntax);
    /// recursively turns one cluster node into tree nodes while adding uncovered
    /// terminal taxa.

    fn build_cluster(idx: usize, clusters: &[ClusterNode], tree: &mut Tree, ntax: usize) -> usize {
        let cluster = &clusters[idx];
        let mut members = cluster.set.taxa();
        members.retain(|&t| t < ntax);

        let mut covered = TaxonSet::new(ntax);
        for &ch in &cluster.children {
            covered.union_inplace(&clusters[ch].set);
        }

        let mut child_nodes = Vec::<usize>::new();

        for &ch in &cluster.children {
            let nid = build_cluster(ch, clusters, tree, ntax);
            child_nodes.push(nid);
        }

        for t in members {
            if !covered.contains(t) {
                child_nodes.push(t);
            }
        }

        child_nodes.sort_unstable();

        if child_nodes.len() == 1 { child_nodes[0] } else { tree.add_internal_node(child_nodes) }
    }

    let ingroup_root = build_cluster(0, &clusters, &mut tree, ntax);

    let mut root_children = Vec::<usize>::new();
    root_children.extend_from_slice(outgroups);
    root_children.push(ingroup_root);
    root_children.sort_unstable();

    let root = if root_children.len() == 1 {
        root_children[0]
    } else {
        tree.add_internal_node(root_children)
    };
    tree.root = root;

    Ok(tree)
}

/// Build a strict consensus tree by retaining only clades that occur in every
/// normalized input tree. Unsupported branches naturally collapse into
/// polytomies because no non-common clade is passed to the tree builder.
pub fn strict_consensus_with_outgroups(
    trees: &[Tree],
    outgroups: &[usize],
) -> Result<Tree, String> {
    if trees.is_empty() {
        return Err("nelsen: no trees available".to_string());
    }
    if outgroups.is_empty() {
        return Err("nelsen: outgroup set is empty".to_string());
    }

    let ntax = trees[0].ntax;
    for (i, t) in trees.iter().enumerate() {
        if t.ntax != ntax {
            return Err(format!("nelsen: tree {i} has different ntax"));
        }
    }

    let normalized: Vec<Tree> =
        trees.iter().map(|t| normalize_by_outgroups(t, outgroups)).collect();

    let first_clades = rooted_clades_excluding_outgroups(&normalized[0], outgroups);
    let mut common = BTreeSet::<TaxonSet>::new();
    for c in first_clades {
        common.insert(c);
    }

    for tr in &normalized[1..] {
        let set: BTreeSet<TaxonSet> =
            rooted_clades_excluding_outgroups(tr, outgroups).into_iter().collect();
        common = common.intersection(&set).cloned().collect();
    }

    let common_vec: Vec<TaxonSet> = common.into_iter().collect();
    for i in 0..common_vec.len() {
        for j in i + 1..common_vec.len() {
            if !common_vec[i].is_compatible_with(&common_vec[j]) {
                return Err("nelsen: incompatible clades survived strict intersection".to_string());
            }
        }
    }

    build_tree_from_rooted_clades_with_outgroups(ntax, outgroups, &common_vec)
}
/// deduplicates input trees, builds strict consensus, scores it, and returns
/// summary statistics.

pub fn nelsen_consensus(
    ds: &Dataset,
    cfg: &CharConfig,
    trees: &[Tree],
    outgroups: &[usize],
) -> Result<ConsensusOutcome, String> {
    if trees.is_empty() {
        return Err("nelsen: no trees available".to_string());
    }

    let mut ws = ScoreWorkspace::new();
    let uniq = uniq_rooted_by_outgroups(trees, outgroups, ds, cfg, &mut ws);
    let consensus = strict_consensus_with_outgroups(&uniq, outgroups)?;

    let best_len = score_tree(ds, cfg, &consensus, &mut ws);
    let min_len = minsteps_sum(ds, cfg);
    let max_len = maxsteps_sum(ds, cfg);
    let ci = calc_ci(min_len, best_len);
    let ri = calc_ri(min_len, max_len, best_len);

    Ok(ConsensusOutcome {
        tree: consensus,
        input_trees: trees.len(),
        uniq_trees: uniq.len(),
        best_len,
        ci,
        ri,
    })
}

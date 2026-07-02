/*!Module: exact implicit-enumeration search over all legal taxon insertion
sequences, with branch-and-bound pruning and topology deduplication.

 */
use std::collections::HashSet;

use crate::engines::{
    ccode::CharConfig,
    dataset::Dataset,
    search::{
        ScoreWorkspace, collapse_zero_length_branches, collapsed_topology_hash,
        rooted_topology_hash, score_tree,
    },
    trees::Tree,
};

/// Implicit enumeration mode: keep one tree, keep all, or unlimited.
#[derive(Debug, Clone, Copy)]
pub enum IeMode {
    KeepUpTo(usize), /* ie; */
    KeepAll,         /* ie*; */
    KeepOne,         /* ie-; */
}

/// Result of implicit enumeration: best length, tree set, and search
/// statistics.
#[derive(Debug, Clone)]
pub struct IeOutcome {
    pub trees: Vec<Tree>,
    pub best_len: u64,
}
/// validates inputs, builds an initial bound, enumerates all rooted insertion
/// trees, and returns the optimal trees requested by mode.

pub fn ie_exact(
    ds: &Dataset,
    cfg: &CharConfig,
    outgroup: usize,
    mode: IeMode,
) -> Result<IeOutcome, String> {
    if ds.ntax < 3 {
        return Err("ie: ntax must be >= 3".to_string());
    }
    if outgroup >= ds.ntax {
        return Err(format!("ie: outgroup {outgroup} out of range"));
    }
    if cfg.len() != ds.nchar {
        return Err(format!("ie: char_config nchar={} != dataset nchar={}", cfg.len(), ds.nchar));
    }

    /*
      Build a quick feasible tree first so the exact enumerator starts with a
      finite pruning bound.
    */
    let mut ws = ScoreWorkspace::new();
    let init = greedy_build_initial_bound(ds, cfg, outgroup, &mut ws)?;
    let bound = score_tree(ds, cfg, &init, &mut ws);

    /* Start enumeration from the canonical three-taxon tree */
    let start = initial_three_taxon_tree(ds.ntax, outgroup)?;

    let mut ctx =
        Ctx { ds, cfg, outgroup, mode, ws, bound, best: Vec::new(), best_seen: HashSet::new() };

    /* Keep a deterministic ingroup order for exact enumeration. */
    let ingroup: Vec<usize> = (0..ds.ntax).filter(|&t| t != outgroup).collect();

    /*
      The start tree already contains the outgroup and the first two ingroup
      taxa; recurse over all remaining taxa.
    */
    let fixed = start.fixed_taxa.clone();
    let mut remaining: Vec<usize> = ingroup.into_iter().filter(|t| !fixed.contains(t)).collect();

    /* Recursively enumerate every legal insertion sequence. */
    dfs_enum(&mut ctx, &start.tree, &mut remaining)?;

    /* Return trees in a stable topology-hash order for deterministic output. */
    ctx.best.sort_by_key(|t| rooted_topology_hash(t));
    Ok(IeOutcome { trees: ctx.best, best_len: ctx.bound })
}

struct Ctx<'a> {
    ds: &'a Dataset,
    cfg: &'a CharConfig,
    outgroup: usize,
    mode: IeMode,
    ws: ScoreWorkspace,
    bound: u64,
    best: Vec<Tree>,
    best_seen: HashSet<u64>,
}

/// Store the canonical three-taxon start tree and the taxa already fixed in it
/// so enumeration can continue with the fourth taxon.
struct ThreeStart {
    tree: Tree,
    fixed_taxa: Vec<usize>,
}

/// Build the start tree `(outgroup (a b))`, where `a` and `b` are the smallest
/// two ingroup taxon indices.
fn initial_three_taxon_tree(ntax: usize, outgroup: usize) -> Result<ThreeStart, String> {
    let mut ingroup: Vec<usize> = (0..ntax).filter(|&t| t != outgroup).collect();
    if ingroup.len() < 2 {
        return Err("ie: need at least 2 ingroup taxa".to_string());
    }
    ingroup.sort_unstable();
    let a = ingroup[0];
    let b = ingroup[1];

    let tr = Tree::new_outgrouped(ntax, outgroup, a, b);
    Ok(ThreeStart { tree: tr, fixed_taxa: vec![a, b] })
}

/// Construct a feasible tree by greedy insertion and use it as the initial
/// branch-and-bound limit.
fn greedy_build_initial_bound(
    ds: &Dataset,
    cfg: &CharConfig,
    outgroup: usize,
    ws: &mut ScoreWorkspace,
) -> Result<Tree, String> {
    let mut ingroup: Vec<usize> = (0..ds.ntax).filter(|&x| x != outgroup).collect();
    ingroup.sort_unstable();

    let a = ingroup[0];
    let b = ingroup[1];
    let mut tr = Tree::new_outgrouped(ds.ntax, outgroup, a, b);

    for &clip in &ingroup[2..] {
        let edges = tr.rooted_edges();

        let mut best_len = u64::MAX;
        let mut best_edge = None::<(usize, usize)>;

        for &(p, c) in &edges {
            if c == outgroup {
                continue;
            }
            let mut tmp = tr.clone();
            tmp.insert_taxon_on_edge(p, c, clip);
            let len = score_tree(ds, cfg, &tmp, ws);
            if len < best_len {
                best_len = len;
                best_edge = Some((p, c));
            }
        }

        let (p, c) = best_edge.ok_or("ie: no insertion edge found")?;
        tr.insert_taxon_on_edge(p, c, clip);
    }

    Ok(tr)
}
/// recursively expands remaining taxa on every legal edge, prunes partial
/// trees above the bound, and records best collapsed topologies.

fn dfs_enum(ctx: &mut Ctx, cur: &Tree, remaining: &mut Vec<usize>) -> Result<(), String> {
    if remaining.is_empty() {
        /* complete tree: evaluate */
        let mut tmp = cur.clone();
        tmp.canonicalize();
        let len = score_tree(ctx.ds, ctx.cfg, &tmp, &mut ctx.ws);

        if len < ctx.bound {
            ctx.bound = len;
            ctx.best.clear();
            ctx.best_seen.clear();
        }
        if len == ctx.bound {
            /*
              Normalize by outgroup, then collapse before
              computing dedup hash.  This ensures trees that differ only
              in unsupported resolution are treated as identical
            */
            let h = {
                let og = [ctx.outgroup];
                let mut norm = tmp
                    .reroot_by_outgroup_set(&og)
                    .or_else(|_| tmp.reroot_by_outgroup(ctx.outgroup))
                    .unwrap_or_else(|_| tmp.clone());
                norm.canonicalize();
                collapse_zero_length_branches(&mut norm, ctx.ds, ctx.cfg, &mut ctx.ws);
                collapsed_topology_hash(&norm)
            };
            if ctx.best_seen.insert(h) {
                match ctx.mode {
                    IeMode::KeepOne => {
                        ctx.best.clear();
                        ctx.best.push(tmp);
                        /*
                          KeepOne still cannot stop at the first tree matching
                          the current bound, because a lower bound may still be
                          found later without a stronger subtree lower bound.
                        */
                    }
                    IeMode::KeepAll => ctx.best.push(tmp),
                    IeMode::KeepUpTo(lim) => {
                        if ctx.best.len() < lim {
                            ctx.best.push(tmp);
                        }
                    }
                }
            }
        }
        return Ok(());
    }

    /*
      Pick the next taxon in fixed order.
    */
    let next_tax = remaining.pop().unwrap();

    /*
      Try every directed insertion edge returned as `(parent, child)` by
      `Tree::rooted_edges()`.
    */
    let edges = cur.rooted_edges();

    for &(p, c) in &edges {
        /*
          Avoid inserting below the outgroup edge; it adds redundant rooted
          placements for this enumeration.
        */
        if c == ctx.outgroup {
            continue;
        }

        let mut t2 = cur.clone();
        t2.insert_taxon_on_edge(p, c, next_tax);

        /*
          Branch-and-bound check: the current partial-tree score is a lower
          bound on every completion of this partial topology.
        */
        let partial_len = score_tree(ctx.ds, ctx.cfg, &t2, &mut ctx.ws);

        if partial_len > ctx.bound {
            continue; /* prune */
        }

        dfs_enum(ctx, &t2, remaining)?;
    }

    /* restore */
    remaining.push(next_tax);
    Ok(())
}

/*!Module: branch mapping for apomorphy/homoplasy SVG output.

*/
use std::collections::HashMap;

use crate::engines::{
    ccode::CharConfig,
    dataset::{Dataset, StateSet},
    search::{INF, NSTATES},
    trees::Tree,
    util,
};

/// Type of mapped character-state change for SVG apomorphy output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApoChangeKind {
    Apomorphy,
    Homoplasy,
}

/// One character-state change mapped to a directed parent-child branch.
#[derive(Debug, Clone)]
pub struct ApoBranchChange {
    pub parent: usize,
    pub child: usize,
    pub character: usize,
    pub state: String,
    pub kind: ApoChangeKind,
}

/// Maps active character-state changes to tree branches and classifies each
/// change as apomorphy or homoplasy. Root state is left unconstrained so basal
/// ambiguity is preserved; additive and non-additive characters differ only in
/// transition costs and marker admissibility.
pub fn detect_apomorphic_changes(
    ds: &Dataset,
    cfg: &CharConfig,
    tr: &Tree,
) -> Result<Vec<ApoBranchChange>, String> {
    if cfg.len() != ds.nchar {
        return Err(format!("char_config nchar={} != dataset nchar={}", cfg.len(), ds.nchar));
    }
    if tr.ntax != ds.ntax {
        return Err(format!("tree ntax {} != dataset ntax {}", tr.ntax, ds.ntax));
    }

    let descendant_masks = descendant_taxon_masks(tr, ds.ntax)?;
    let mut raw = Vec::<RawApoChange>::new();
    let mut counts = HashMap::<(usize, u64), usize>::new();

    for ch in 0..ds.nchar {
        if !cfg.chars[ch].active {
            continue;
        }

        for change in detect_character_changes(ds, tr, ch, cfg.chars[ch].additive)? {
            *counts.entry((change.character, change.state_bits)).or_insert(0) += 1;
            raw.push(change);
        }
    }

    let mut out = Vec::with_capacity(raw.len());
    for change in raw {
        let count = counts.get(&(change.character, change.state_bits)).copied().unwrap_or(0);
        let kind = if count > 1
            || state_exists_outside_subtree(
                ds,
                &descendant_masks[change.child],
                change.character,
                change.state_bits,
            )? {
            ApoChangeKind::Homoplasy
        } else {
            ApoChangeKind::Apomorphy
        };

        out.push(ApoBranchChange {
            parent: change.parent,
            child: change.child,
            character: change.character,
            state: format_state_bits(change.state_bits),
            kind,
        });
    }

    out.sort_by_key(|x| (x.parent, x.child, x.character, x.state.clone()));
    Ok(out)
}

#[derive(Debug, Clone)]
struct RawApoChange {
    parent: usize,
    child: usize,
    character: usize,
    state_bits: u64,
}

/// Computes descendant-taxon membership for every node, used to decide whether
/// a state also appears outside the branch being annotated.
fn descendant_taxon_masks(tr: &Tree, ntax: usize) -> Result<Vec<Vec<bool>>, String> {
    let mut masks = vec![vec![false; ntax]; tr.nodes.len()];
    let mut postorder = Vec::with_capacity(tr.nodes.len());
    postorder_fill(tr, tr.root, &mut postorder);

    for &v in &postorder {
        if let Some(taxon) = tr.nodes[v].taxon {
            if taxon >= ntax {
                return Err(format!("taxon {taxon} out of range (ntax={ntax})"));
            }
            masks[v][taxon] = true;
        } else {
            let mut mask = vec![false; ntax];
            for &child in &tr.nodes[v].children {
                for taxon in 0..ntax {
                    mask[taxon] |= masks[child][taxon];
                }
            }
            masks[v] = mask;
        }
    }

    Ok(masks)
}

/// Returns true when the same state is present in terminal taxa outside the
/// child subtree, making a unique branch change homoplasious.
fn state_exists_outside_subtree(
    ds: &Dataset,
    subtree: &[bool],
    ch: usize,
    state_bits: u64,
) -> Result<bool, String> {
    for (taxon, in_subtree) in subtree.iter().copied().enumerate() {
        if in_subtree {
            continue;
        }
        let bits = dataset_state_at(ds, taxon, ch)?.bits();
        if bits != StateSet::ALL36.bits() && bits & state_bits != 0 {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Optimizes one character over the whole tree and returns raw branch changes
/// before apomorphy/homoplasy classification.
fn detect_character_changes(
    ds: &Dataset,
    tr: &Tree,
    ch: usize,
    additive: bool,
) -> Result<Vec<RawApoChange>, String> {
    let n_nodes = tr.nodes.len();
    let mut down = vec![INF; n_nodes * NSTATES];
    let mut up = vec![INF; n_nodes * NSTATES];
    let mut postorder = Vec::with_capacity(n_nodes);
    postorder_fill(tr, tr.root, &mut postorder);

    for (node_id, node) in tr.nodes.iter().enumerate() {
        if let Some(taxon) = node.taxon {
            let bits = normalize_sankoff_leaf_bits(dataset_state_at(ds, taxon, ch)?.bits());
            for state in 0..NSTATES {
                if (bits >> state) & 1 == 1 {
                    down[node_id * NSTATES + state] = 0;
                }
            }
        }
    }

    for &node_id in &postorder {
        let node = &tr.nodes[node_id];
        if node.taxon.is_some() {
            continue;
        }
        if node.children.is_empty() {
            return Err(format!("internal node {node_id} has no children"));
        }

        for parent_state in 0..NSTATES {
            let mut total = 0u64;
            for &child in &node.children {
                total =
                    add_cost(total, sankoff_min_child_cost(&down, child, parent_state, additive));
            }
            down[node_id * NSTATES + parent_state] = total;
        }
    }

    let root_base = tr.root * NSTATES;
    let best = down[root_base..root_base + NSTATES].iter().copied().min().unwrap_or(INF);
    if best >= INF {
        return Err(format!("character {ch} has no finite reconstruction"));
    }

    for state in 0..NSTATES {
        if down[root_base + state] == best {
            up[root_base + state] = 0;
        }
    }

    for &parent in postorder.iter().rev() {
        let node = &tr.nodes[parent];
        if node.taxon.is_some() {
            continue;
        }

        let mut child_costs = Vec::<(usize, [u64; NSTATES])>::with_capacity(node.children.len());
        for &child in &node.children {
            let mut costs = [INF; NSTATES];
            for (parent_state, slot) in costs.iter_mut().enumerate() {
                *slot = sankoff_min_child_cost(&down, child, parent_state, additive);
            }
            child_costs.push((child, costs));
        }

        for &(child, _) in &child_costs {
            for child_state in 0..NSTATES {
                let mut best_up = INF;
                for parent_state in 0..NSTATES {
                    let mut total = up[parent * NSTATES + parent_state];
                    for &(sibling, sibling_costs) in &child_costs {
                        if sibling != child {
                            total = add_cost(total, sibling_costs[parent_state]);
                        }
                    }
                    total = add_cost(
                        total,
                        sankoff_transition_cost(parent_state, child_state, additive),
                    );
                    best_up = best_up.min(total);
                }
                up[child * NSTATES + child_state] = best_up;
            }
        }
    }

    let mut final_bits = vec![0u64; n_nodes];
    for node_id in 0..n_nodes {
        let base = node_id * NSTATES;
        for state in 0..NSTATES {
            if add_cost(down[base + state], up[base + state]) == best {
                final_bits[node_id] |= 1u64 << state;
            }
        }
    }

    let mut out = Vec::new();
    for (parent, child) in tr.rooted_edges() {
        let parent_bits = final_bits[parent];
        let child_bits = final_bits[child];
        if !is_marker_change(parent_bits, child_bits, additive) {
            continue;
        }
        out.push(RawApoChange { parent, child, character: ch, state_bits: child_bits });
    }

    Ok(out)
}

/// Decides whether optimized endpoint state sets are specific enough to draw a
/// marker. Ordered characters allow interval-to-singleton changes outside the
/// parent interval; unordered characters require singleton endpoints.
fn is_marker_change(parent_bits: u64, child_bits: u64, additive: bool) -> bool {
    if additive {
        /* Ordered characters mark a move to a singleton child state when that
        state lies outside the parent's optimal interval. */
        parent_bits != 0 && child_bits.count_ones() == 1 && parent_bits & child_bits == 0
    } else {
        /* Unordered characters are only unambiguous when both endpoint
        optimizations are singleton states. */
        parent_bits.count_ones() == 1 && child_bits.count_ones() == 1 && parent_bits != child_bits
    }
}

/// Treats missing/unknown states as all possible states for optimization.
fn normalize_sankoff_leaf_bits(bits: u64) -> u64 {
    if bits == 0 || bits == StateSet::ALL36.bits() { StateSet::ALL36.bits() } else { bits }
}

/// Computes the minimum cost of connecting a parent state to one child subtree.
fn sankoff_min_child_cost(down: &[u64], child: usize, parent_state: usize, additive: bool) -> u64 {
    let base = child * NSTATES;
    let mut best = INF;
    for child_state in 0..NSTATES {
        let cost = add_cost(
            down[base + child_state],
            sankoff_transition_cost(parent_state, child_state, additive),
        );
        best = best.min(cost);
    }
    best
}

/// Transition cost between concrete states: ordered distance for additive
/// characters, unit change cost for non-additive characters.
fn sankoff_transition_cost(parent_state: usize, child_state: usize, additive: bool) -> u64 {
    if additive {
        parent_state.abs_diff(child_state) as u64
    } else if parent_state == child_state {
        0
    } else {
        1
    }
}

/// Adds dynamic-programming costs while preserving the INF sentinel.
fn add_cost(a: u64, b: u64) -> u64 {
    if a >= INF || b >= INF { INF } else { a.saturating_add(b).min(INF) }
}

/// Safely retrieves one taxon-character state from the row-major matrix.
fn dataset_state_at(ds: &Dataset, taxon: usize, ch: usize) -> Result<StateSet, String> {
    let idx = taxon
        .checked_mul(ds.nchar)
        .and_then(|x| x.checked_add(ch))
        .ok_or_else(|| "dataset index overflow".to_string())?;

    ds.matrix
        .get(idx)
        .copied()
        .ok_or_else(|| format!("dataset state index out of range: taxon={taxon} ch={ch}"))
}

/// Collects node ids in postorder for bottom-up dynamic programming.
fn postorder_fill(tr: &Tree, v: usize, out: &mut Vec<usize>) {
    for &child in &tr.nodes[v].children {
        postorder_fill(tr, child, out);
    }
    out.push(v);
}

/// Formats a low-36-bit state set as Hennig/TNT-style state symbols.
fn format_state_bits(bits: u64) -> String {
    if bits == 0 {
        return "?".to_string();
    }

    let mut s = String::new();
    for i in 0..36usize {
        if ((bits >> i) & 1) == 1 {
            s.push(util::idx_to_symbol(i as u8));
        }
    }
    s
}

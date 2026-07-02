/*!Module: tree diagnostics for xsteps commands, including tree lengths,
per-character fits, hypothetical ancestor states, successive weighting, and
best/worst summaries across trees.

 */
use crate::engines::{
    ccode::{CharConfig, CharSetting},
    dataset::{Dataset, StateSet},
    search::{INF, NSTATES, calc_ci, calc_ri, maxsteps_sum, minsteps_sum, score_tree},
    trees::Tree,
    util,
};

/// One row in the tree-length report: total weighted length.
#[derive(Debug, Clone)]
pub struct TreeLengthRow {
    pub length: u64,
}

/// Per-character fit statistics: steps, CI, RI for one tree.
#[derive(Debug, Clone)]
pub struct CharacterFitRow {
    pub character: usize,
    pub steps: u64,
    pub ci: f64,
    pub ri: Option<f64>,
}

/// Character-by-character report for one tree with state assignments.
#[derive(Debug, Clone)]
pub struct TreeCharacterReport {
    pub tree_index: usize,
    pub length: u64,
    pub ci: f64,
    pub ri: Option<f64>,
    pub rows: Vec<CharacterFitRow>,
}

/// Ancestral state and synapomorphy report for one character across trees.
#[derive(Debug, Clone)]
pub struct CharacterHistoryReport {
    pub character: usize,
    pub internal_nodes: Vec<usize>,
    pub node_states: Vec<String>,
}

/// Complete ancestral-state reconstruction for all active characters on one
/// tree.
#[derive(Debug, Clone)]
pub struct TreeHistoryReport {
    pub tree_index: usize,
    pub histories: Vec<CharacterHistoryReport>,
}

/// Full analysis of one character on one tree: steps, states, and per-node
/// assignments.
#[derive(Debug, Clone)]
pub struct SingleCharAnalysis {
    pub steps: u64,
    pub node_states: Vec<u64>, /* per node: low 36 bits used */
}
/// scores each supplied tree and returns one tree-length row per tree.

pub fn analyze_lengths(ds: &Dataset, cfg: &CharConfig, trees: &[Tree]) -> Vec<TreeLengthRow> {
    let mut rows = Vec::with_capacity(trees.len());
    let mut ws = crate::engines::search::ScoreWorkspace::new();

    for tr in trees.iter() {
        let len = score_tree(ds, cfg, tr, &mut ws);
        rows.push(TreeLengthRow { length: len });
    }

    rows
}
/// computes per-character fit rows and whole-tree CI/RI summaries for each
/// tree.

pub fn analyze_character_fits(
    ds: &Dataset,
    cfg: &CharConfig,
    trees: &[Tree],
) -> Result<Vec<TreeCharacterReport>, String> {
    let mut out = Vec::with_capacity(trees.len());

    for (tree_index, tr) in trees.iter().enumerate() {
        let rows = analyze_tree_character_fits(ds, cfg, tr)?;
        let length: u64 = rows.iter().map(|r| r.steps).sum();

        let min_len = minsteps_sum(ds, cfg);
        let max_len = maxsteps_sum(ds, cfg);
        let ci = calc_ci(min_len, length);
        let ri = calc_ri(min_len, max_len, length);

        out.push(TreeCharacterReport { tree_index, length, ci, ri, rows });
    }

    Ok(out)
}
/// computes possible ancestral state strings for every internal node and
/// character.

pub fn analyze_histories(
    ds: &Dataset,
    _cfg: &CharConfig,
    trees: &[Tree],
) -> Result<Vec<TreeHistoryReport>, String> {
    let mut out = Vec::with_capacity(trees.len());

    for (tree_index, tr) in trees.iter().enumerate() {
        let mut histories = Vec::with_capacity(ds.nchar);

        let internal_nodes: Vec<usize> = tr
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(i, n)| if n.taxon.is_none() { Some(i) } else { None })
            .collect();

        for ch in 0..ds.nchar {
            let analyzed = analyze_single_character(ds, tr, ch)?;
            let mut node_states = Vec::with_capacity(internal_nodes.len());

            for &nid in &internal_nodes {
                node_states.push(format_state_bits(analyzed.node_states[nid]));
            }

            histories.push(CharacterHistoryReport {
                character: ch,
                internal_nodes: internal_nodes.clone(),
                node_states,
            });
        }

        out.push(TreeHistoryReport { tree_index, histories });
    }

    Ok(out)
}
/// evaluates every character on one tree using the active character
/// configuration.

fn analyze_tree_character_fits(
    ds: &Dataset,
    cfg: &CharConfig,
    tr: &Tree,
) -> Result<Vec<CharacterFitRow>, String> {
    let mut rows = Vec::with_capacity(ds.nchar);

    for ch in 0..ds.nchar {
        let cs = cfg.chars[ch];
        let analyzed = analyze_single_character_with_setting(ds, tr, ch, cs)?;
        let obs = analyzed.steps;
        let min_s = minsteps_char_with_setting(ds, ch, cs);
        let max_s = maxsteps_char_with_setting(ds, ch, cs);

        let ci = if obs == 0 { 1.0 } else { (min_s as f64) / (obs as f64) };

        let ri = if max_s > min_s {
            Some((max_s.saturating_sub(obs)) as f64 / (max_s - min_s) as f64)
        } else {
            None
        };

        rows.push(CharacterFitRow { character: ch, steps: obs, ci, ri });
    }

    Ok(rows)
}
/// diagnoses one character as active, non-additive, and unit-weighted for
/// display-only history output.

pub fn analyze_single_character(
    ds: &Dataset,
    tr: &Tree,
    ch: usize,
) -> Result<SingleCharAnalysis, String> {
    analyze_single_character_with_setting(
        ds,
        tr,
        ch,
        CharSetting { active: true, additive: false, weight: 1 },
    )
}
/// dispatches inactive, additive, and non-additive single-character scoring
/// paths.

fn analyze_single_character_with_setting(
    ds: &Dataset,
    tr: &Tree,
    ch: usize,
    cs: CharSetting,
) -> Result<SingleCharAnalysis, String> {
    if ch >= ds.nchar {
        return Err(format!("character {ch} out of range"));
    }

    if tr.ntax != ds.ntax {
        return Err(format!("tree ntax {} != dataset ntax {}", tr.ntax, ds.ntax));
    }

    if !cs.active {
        return Ok(SingleCharAnalysis { steps: 0, node_states: vec![0u64; tr.nodes.len()] });
    }

    if cs.additive {
        return analyze_single_additive_character(ds, tr, ch, cs.weight as u64);
    }

    let mut node_states = vec![0u64; tr.nodes.len()];
    let mut steps = 0u64;
    /// local Fitch recursion that propagates non-additive state sets upward and
    /// counts union events.

    fn rec(
        ds: &Dataset,
        tr: &Tree,
        v: usize,
        ch: usize,
        node_states: &mut [u64],
        steps: &mut u64,
    ) -> Result<u64, String> {
        let node = &tr.nodes[v];

        if let Some(taxon) = node.taxon {
            let ss = dataset_state_at(ds, taxon, ch)?;
            let bits = ss.bits();
            node_states[v] = bits;
            return Ok(bits);
        }

        if node.children.is_empty() {
            return Err(format!("internal node {v} has no children"));
        }

        let mut child_bits = Vec::with_capacity(node.children.len());
        for &c in &node.children {
            child_bits.push(rec(ds, tr, c, ch, node_states, steps)?);
        }

        let mut inter = child_bits[0];
        for &b in &child_bits[1..] {
            inter &= b;
        }

        let here = if inter != 0 {
            inter
        } else {
            let mut uni = 0u64;
            for &b in &child_bits {
                uni |= b;
            }
            *steps += 1;
            uni
        };

        node_states[v] = here;
        Ok(here)
    }

    rec(ds, tr, tr.root, ch, &mut node_states, &mut steps)?;

    Ok(SingleCharAnalysis { steps: steps * cs.weight as u64, node_states })
}

#[derive(Clone, Copy, Debug, Default)]
struct OrderedRange {
    lo: u8,
    hi: u8,
}
/// applies ordered-state Sankoff scoring to one additive character and records
/// optimal root-state sets.

fn analyze_single_additive_character(
    ds: &Dataset,
    tr: &Tree,
    ch: usize,
    weight: u64,
) -> Result<SingleCharAnalysis, String> {
    let mut postorder = Vec::with_capacity(tr.nodes.len());
    postorder_fill(tr, tr.root, &mut postorder);

    let mut costs = vec![INF; tr.nodes.len() * NSTATES];
    let mut node_states = vec![0u64; tr.nodes.len()];

    for (i, node) in tr.nodes.iter().enumerate() {
        if let Some(taxon) = node.taxon {
            let bits = dataset_state_at(ds, taxon, ch)?.bits();
            let base = i * NSTATES;
            for s in 0..NSTATES {
                if ((bits >> s) & 1) == 1 {
                    costs[base + s] = 0;
                }
            }
            node_states[i] = bits;
        }
    }

    for &v in &postorder {
        let node = &tr.nodes[v];
        if node.taxon.is_some() {
            continue;
        }

        if node.children.is_empty() {
            return Err(format!("internal node {v} has no children"));
        }

        let base = v * NSTATES;
        for s in 0..NSTATES {
            let mut total = 0u64;
            for &child in &node.children {
                let child_base = child * NSTATES;
                let mut best = INF;
                for k in 0..NSTATES {
                    let edge = s.abs_diff(k) as u64;
                    best = best.min(costs[child_base + k].saturating_add(edge));
                }
                total = total.saturating_add(best);
            }
            costs[base + s] = total;
        }

        let best = costs[base..base + NSTATES].iter().copied().min().unwrap_or(0);
        let mut bits = 0u64;
        for s in 0..NSTATES {
            if costs[base + s] == best {
                bits |= 1u64 << s;
            }
        }
        node_states[v] = bits;
    }

    let root_base = tr.root * NSTATES;
    let steps = costs[root_base..root_base + NSTATES].iter().copied().min().unwrap_or(0);

    Ok(SingleCharAnalysis { steps: steps * weight, node_states })
}
/// collects a postorder node traversal for dynamic programs.

fn postorder_fill(tr: &Tree, v: usize, out: &mut Vec<usize>) {
    for &child in &tr.nodes[v].children {
        postorder_fill(tr, child, out);
    }
    out.push(v);
}
/// converts a bitset into the lowest/highest observed ordered state interval.

fn range_from_bits(bits: u64) -> OrderedRange {
    if bits == StateSet::ALL36.bits() {
        return OrderedRange { lo: 0, hi: 35 };
    }

    let mut lo: Option<u8> = None;
    let mut hi = 0u8;
    for i in 0..36u8 {
        if ((bits >> i) & 1) == 1 {
            if lo.is_none() {
                lo = Some(i);
            }
            hi = i;
        }
    }

    OrderedRange { lo: lo.unwrap_or(0), hi }
}
/// safely retrieves one taxon-character state set from the row-major matrix.

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
/// computes default non-additive minimum steps for one character.

fn minsteps_char(ds: &Dataset, ch: usize) -> u64 {
    minsteps_char_with_setting(ds, ch, CharSetting { active: true, additive: false, weight: 1 })
}
/// computes weighted minimum possible steps for one character under
/// active/additive settings.

fn minsteps_char_with_setting(ds: &Dataset, ch: usize, cs: CharSetting) -> u64 {
    if !cs.active {
        return 0;
    }

    let weight = cs.weight as u64;

    if cs.additive {
        let mut lo: Option<u8> = None;
        let mut hi = 0u8;

        for taxon in 0..ds.ntax {
            let idx = taxon * ds.nchar + ch;
            let bits = ds.matrix[idx].bits();
            if bits == StateSet::ALL36.bits() {
                continue;
            }
            let rg = range_from_bits(bits);
            lo = Some(lo.map(|x| x.min(rg.lo)).unwrap_or(rg.lo));
            hi = hi.max(rg.hi);
        }

        return lo.map(|x| (hi - x) as u64 * weight).unwrap_or(0);
    }

    let mut seen = 0u64;

    for taxon in 0..ds.ntax {
        let idx = taxon * ds.nchar + ch;
        let bits = ds.matrix[idx].bits();
        if bits != StateSet::ALL36.bits() {
            seen |= bits;
        }
    }

    let k = seen.count_ones() as u64;
    if k == 0 { 0 } else { (k - 1) * weight }
}
/// computes default non-additive maximum steps for one character.

fn maxsteps_char(ds: &Dataset, ch: usize) -> u64 {
    maxsteps_char_with_setting(ds, ch, CharSetting { active: true, additive: false, weight: 1 })
}
/// computes weighted maximum possible steps for one character under
/// active/additive settings.

fn maxsteps_char_with_setting(ds: &Dataset, ch: usize, cs: CharSetting) -> u64 {
    if !cs.active {
        return 0;
    }

    let weight = cs.weight as u64;

    if cs.additive {
        let mut vals = Vec::<u8>::new();
        for taxon in 0..ds.ntax {
            let idx = taxon * ds.nchar + ch;
            let bits = ds.matrix[idx].bits();
            if bits == StateSet::ALL36.bits() || bits.count_ones() != 1 {
                continue;
            }
            vals.push(bits.trailing_zeros() as u8);
        }

        if vals.is_empty() {
            return 0;
        }

        let mut best = 0u64;
        for &candidate in &vals {
            let sum = vals.iter().map(|&v| v.abs_diff(candidate) as u64).sum::<u64>();
            best = best.max(sum);
        }
        return best * weight;
    }

    use std::collections::HashMap;

    let mut counts = HashMap::<u64, usize>::new();
    let mut total = 0usize;

    for taxon in 0..ds.ntax {
        let idx = taxon * ds.nchar + ch;
        let bits = ds.matrix[idx].bits();

        /*
          Only singleton states are treated as fixed here; polymorphic and
          missing states do not refine the maximum-step statistic.
        */
        if bits != 0 && bits.count_ones() == 1 {
            *counts.entry(bits).or_insert(0) += 1;
            total += 1;
        }
    }

    if total == 0 {
        return 0;
    }

    let max_count = counts.values().copied().max().unwrap_or(0);
    (total.saturating_sub(max_count)) as u64 * weight
}
/// formats low-36-bit state masks for xsteps output.

pub fn format_state_bits(bits: u64) -> String {
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

/// Minimum and maximum possible steps for one character from the dataset
/// alone.
#[derive(Debug, Clone)]
pub struct CharMinMaxRow {
    pub character: usize,
    pub min_steps: u64,
    pub max_steps: u64,
}

/// For `steps;`: show min/max steps per character (dataset-only, independent of trees).
pub fn character_min_max_steps(ds: &Dataset) -> Vec<CharMinMaxRow> {
    let mut out = Vec::with_capacity(ds.nchar);
    for ch in 0..ds.nchar {
        out.push(CharMinMaxRow {
            character: ch,
            min_steps: minsteps_char(ds, ch),
            max_steps: maxsteps_char(ds, ch),
        });
    }
    out
}

/// Successive-weighting result: character, RC value, and rc10 bucket.
#[derive(Debug, Clone)]
pub struct SaWeightRow {
    pub character: usize,
    pub rc10: u32, /* 0/1/3/10 (Hennig86-style bucket) */
}
/// computes Hennig86-style RC and rc10 weights from one diagnosed tree.

pub fn successive_weights_rc10_from_tree(
    ds: &crate::engines::dataset::Dataset,
    tr: &crate::engines::trees::Tree,
) -> Result<Vec<SaWeightRow>, String> {
    if ds.nchar == 0 {
        return Err("xsteps w: dataset has 0 characters".to_string());
    }

    let mut rows = Vec::with_capacity(ds.nchar);

    for ch in 0..ds.nchar {
        let analyzed = analyze_single_character(ds, tr, ch)?;
        let steps = analyzed.steps;

        let min_s = minsteps_char(ds, ch);
        let max_s = maxsteps_char(ds, ch);

        /* CI: min/obs (obs==0 => treat as perfect) */
        let ci = if steps == 0 { 1.0 } else { (min_s as f64) / (steps as f64) };

        /* RI: (max-obs)/(max-min), if defined; NA if max==min */
        let ri = if max_s > min_s {
            Some((max_s.saturating_sub(steps)) as f64 / (max_s - min_s) as f64)
        } else {
            None
        };

        /*
          Hennig86 behavior guess:
          RI=NA => treat as 1.0 (full credit)
        */
        let ri_eff = ri.unwrap_or(1.0);

        let rc = ci * ri_eff;
        let rc10 = rc_to_rc10_floor(rc);

        rows.push(SaWeightRow { character: ch, rc10 });
    }

    Ok(rows)
}
/// converts a floating RC value to the integer weight bucket.

fn rc_to_rc10_floor(rc: f64) -> u32 {
    if !rc.is_finite() || rc <= 0.0 {
        0
    } else if rc >= 1.0 {
        10
    } else {
        (rc * 10.0).floor() as u32
    }
}
/// writes successive-weight  back to the mutable character
/// configuration.

pub fn apply_successive_weights_to_ccode(
    cfg: &mut crate::engines::ccode::CharConfig,
    nchar: usize,
    rows: &[SaWeightRow],
) -> Result<(), String> {
    if cfg.len() != nchar {
        return Err(format!("ccode length mismatch: cfg.len()={} nchar={}", cfg.len(), nchar));
    }

    for r in rows {
        if r.character >= nchar {
            return Err(format!("character {} out of range (nchar={})", r.character, nchar));
        }

        /* Write Hennig86-style rc10   */
        cfg.chars[r.character].weight = r.rc10;
    }

    Ok(())
}

/// Best and worst fit statistics for one character across a tree set.
#[derive(Debug, Clone)]
pub struct XStepsMBestWorstRow {
    pub character: usize,

    pub best_steps: u64,
    pub worst_steps: u64,

    pub best_ci: f64,
    pub worst_ci: f64,

    pub best_ri: Option<f64>,
    pub worst_ri: Option<f64>,
}
/// summarizes best and worst per-character fits across all current trees.

pub fn analyze_best_worst_fits_across_trees(
    ds: &Dataset,
    cfg: &CharConfig,
    trees: &[Tree],
) -> Result<Vec<XStepsMBestWorstRow>, String> {
    if trees.is_empty() {
        return Err("xsteps m: no trees provided".to_string());
    }
    if ds.nchar == 0 {
        return Err("xsteps m: dataset has 0 characters".to_string());
    }
    if cfg.len() != ds.nchar {
        return Err(format!(
            "xsteps m: char_config nchar={} != dataset nchar={}",
            cfg.len(),
            ds.nchar
        ));
    }

    let reports = analyze_character_fits(ds, cfg, trees)?;

    let mut best_steps = vec![u64::MAX; ds.nchar];
    let mut worst_steps = vec![0u64; ds.nchar];

    for rep in &reports {
        for row in &rep.rows {
            let ch = row.character;
            let s = row.steps;
            if s < best_steps[ch] {
                best_steps[ch] = s;
            }
            if s > worst_steps[ch] {
                worst_steps[ch] = s;
            }
        }
    }

    let mut out = Vec::with_capacity(ds.nchar);
    for ch in 0..ds.nchar {
        let min_s = minsteps_char(ds, ch);
        let max_s = maxsteps_char(ds, ch);

        let best_obs = best_steps[ch];
        let worst_obs = worst_steps[ch];

        let best_ci = if best_obs == 0 { 1.0 } else { (min_s as f64) / (best_obs as f64) };
        let worst_ci = if worst_obs == 0 { 1.0 } else { (min_s as f64) / (worst_obs as f64) };

        let best_ri = if max_s > min_s {
            Some((max_s.saturating_sub(best_obs)) as f64 / (max_s - min_s) as f64)
        } else {
            None
        };
        let worst_ri = if max_s > min_s {
            Some((max_s.saturating_sub(worst_obs)) as f64 / (max_s - min_s) as f64)
        } else {
            None
        };

        out.push(XStepsMBestWorstRow {
            character: ch,
            best_steps: best_obs,
            worst_steps: worst_obs,
            best_ci,
            worst_ci,
            best_ri,
            worst_ri,
        });
    }

    Ok(out)
}

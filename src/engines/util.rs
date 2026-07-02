/*!Module: shared utility functions used across engine modules for formatting,
plot ordering, additive range logic, and character-set conversion.
*/

use crate::engines::{dataset::StateSet, trees::Tree};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct PlotOrderKey {
    leaf_count: usize,
    min_taxon: usize,
    max_taxon: usize,
}

pub(crate) fn order_children_for_plot(tree: &mut Tree) {
    if tree.nodes.is_empty() || tree.root >= tree.nodes.len() {
        return;
    }
    let _ = order_children_for_plot_rec(tree, tree.root);
}

fn order_children_for_plot_rec(tree: &mut Tree, node: usize) -> PlotOrderKey {
    if let Some(taxon) = tree.nodes[node].taxon {
        return PlotOrderKey { leaf_count: 1, min_taxon: taxon, max_taxon: taxon };
    }

    let children = tree.nodes[node].children.clone();
    let mut keyed_children = Vec::with_capacity(children.len());

    for child in children {
        let key = order_children_for_plot_rec(tree, child);
        keyed_children.push((key, child));
    }

    keyed_children.sort_by(|(ka, a), (kb, b)| ka.cmp(kb).then_with(|| a.cmp(b)));

    tree.nodes[node].children = keyed_children.iter().map(|(_, child)| *child).collect();

    let leaf_count = keyed_children.iter().map(|(key, _)| key.leaf_count).sum();
    let min_taxon = keyed_children.iter().map(|(key, _)| key.min_taxon).min().unwrap_or(usize::MAX);
    let max_taxon = keyed_children.iter().map(|(key, _)| key.max_taxon).max().unwrap_or(0);

    PlotOrderKey { leaf_count, min_taxon, max_taxon }
}

/// Converts a numeric state index (0..35) to its display character (0-9, a-z,
/// A-Z).
pub fn idx_to_symbol(idx: u8) -> char {
    match idx {
        0..=9 => (b'0' + idx) as char,
        10..=35 => (b'a' + (idx - 10)) as char,
        _ => '?',
    }
}

pub(crate) fn apply_charset(lines: Vec<String>, unicode: bool) -> Vec<String> {
    if unicode {
        return lines;
    }
    lines
        .into_iter()
        .map(|line| {
            line.chars()
                .map(|ch| match ch {
                    '┌' => ',',
                    '├' => '+',
                    '└' => '`',
                    '│' => '|',
                    '─' => '-',
                    _ => ch,
                })
                .collect()
        })
        .collect()
}

/* AdditiveRange: shared interval type for ordered character scoring */

#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct AdditiveRange {
    pub lo: u8,
    pub hi: u8,
}

impl AdditiveRange {
    /// Converts a state bitset into an AdditiveRange (lo, hi, contiguous).
    pub fn from_bits(bits: u64) -> Self {
        if bits == StateSet::ALL36.bits() {
            return Self { lo: 0, hi: 35 };
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
        Self { lo: lo.unwrap_or(0), hi }
    }

    /* Merge two intervals using the Farris inner-interval rule.
    Returns the tightest interval containing the optimal parent state
    and the step cost (gap) between them. */
    /// Merges two AdditiveRanges and returns the edge cost to connect them.
    pub fn merge(self, other: Self) -> (Self, u64) {
        let ix_lo = self.lo.max(other.lo);
        let ix_hi = self.hi.min(other.hi);

        if ix_lo <= ix_hi {
            return (Self { lo: ix_lo, hi: ix_hi }, 0);
        }

        let (lo, hi, gap) = if self.hi < other.lo {
            (self.hi, other.lo, (other.lo - self.hi) as u64)
        } else {
            (other.hi, self.lo, (self.lo - other.hi) as u64)
        };
        (Self { lo, hi }, gap)
    }
}

/* Whether a state bitset represents one contiguous ordered-state interval. */
pub(crate) fn is_contiguous(bits: u64) -> bool {
    if bits == 0 || bits == StateSet::ALL36.bits() {
        return true;
    }
    let lo = bits.trailing_zeros() as u8;
    let hi = 63 - bits.leading_zeros() as u8;
    let mask = ((1u64 << (hi - lo + 1)) - 1) << lo;
    (bits & mask) == mask
}

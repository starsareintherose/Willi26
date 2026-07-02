/*!Module: text tree plotting for `tplot`, supporting Unicode and plain ASCII
connector styles.

 */
use crate::engines::{trees::Tree, util};

/// Tree plot character style: normal ASCII or extended line-drawing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TPlotStyle {
    Plain,
    Unicode,
}
/// orders a cloned tree for stable display, builds its layout, and returns
/// rendered lines in the requested character set.

pub fn render_tree(
    tree: &Tree,
    taxon_names: &[String],
    outgroups: Option<&[usize]>,
    style: TPlotStyle,
) -> Result<Vec<String>, String> {
    let mut rooted = tree.clone();
    let _ = outgroups;

    util::order_children_for_plot(&mut rooted);

    let layout = Layout::build(&rooted, rooted.root, taxon_names)?;
    let width = layout.width();

    let mut out = Vec::with_capacity(layout.lines.len());
    for row in &layout.lines {
        let mut line = String::with_capacity(width);
        for cell in row {
            line.push(*cell);
        }
        while line.ends_with(' ') {
            line.pop();
        }
        out.push(line);
    }
    Ok(util::apply_charset(out, style == TPlotStyle::Unicode))
}

#[derive(Debug, Clone)]
struct Layout {
    lines: Vec<Vec<char>>,
    anchor_row: usize,
}

impl Layout {
    /// recursively lays out a subtree with child connector rows and leaf labels.
    fn build(tree: &Tree, node: usize, taxon_names: &[String]) -> Result<Self, String> {
        let n = &tree.nodes[node];

        if let Some(t) = n.taxon {
            let label = taxon_names.get(t).cloned().unwrap_or_else(|| format!("{t}"));
            return Ok(Self { lines: vec![label.chars().collect()], anchor_row: 0 });
        }

        if n.children.is_empty() {
            return Ok(Self { lines: vec![format!("#{node}").chars().collect()], anchor_row: 0 });
        }

        let mut child_blocks = Vec::with_capacity(n.children.len());
        for &ch in &n.children {
            child_blocks.push(Self::build(tree, ch, taxon_names)?);
        }

        let gap = 1usize;
        let mut total_rows = 0usize;
        let mut child_anchor_rows = Vec::with_capacity(child_blocks.len());

        for (i, block) in child_blocks.iter().enumerate() {
            if i > 0 {
                total_rows += gap;
            }
            child_anchor_rows.push(total_rows + block.anchor_row);
            total_rows += block.lines.len();
        }

        let anchor_row = match (child_anchor_rows.first(), child_anchor_rows.last()) {
            (Some(&a), Some(&b)) => (a + b) / 2,
            _ => 0,
        };

        let child_width = child_blocks.iter().map(|b| b.width()).max().unwrap_or(0);
        let left_width = 3usize; /* connector area */
        let total_width = left_width + child_width;
        let mut lines = vec![vec![' '; total_width]; total_rows];

        for (block_idx, block) in child_blocks.iter().enumerate() {
            let top =
                if block_idx == 0 { 0 } else { child_anchor_rows[block_idx] - block.anchor_row };

            for (r, src_line) in block.lines.iter().enumerate() {
                let dst_r = top + r;
                for (c, ch) in src_line.iter().enumerate() {
                    lines[dst_r][left_width + c] = *ch;
                }
            }
        }

        for r in 0..total_rows {
            if r == anchor_row {
                continue;
            }
            let intersects_child_anchor = child_anchor_rows.contains(&r);
            if (r > anchor_row
                && child_anchor_rows.iter().any(|&x| x < r)
                && child_anchor_rows.iter().any(|&x| x > r))
                || (r < anchor_row
                    && child_anchor_rows.iter().any(|&x| x < r)
                    && child_anchor_rows.iter().any(|&x| x > r))
                || (r > anchor_row
                    && child_anchor_rows.iter().any(|&x| x >= r)
                    && child_anchor_rows[0] <= r)
                || (r < anchor_row
                    && child_anchor_rows.last().copied().unwrap_or(0) >= r
                    && child_anchor_rows.iter().any(|&x| x <= r))
            {
                if !intersects_child_anchor {
                    lines[r][0] = '│';
                }
            }
        }

        for (idx, &r) in child_anchor_rows.iter().enumerate() {
            lines[r][0] = if idx == 0 {
                '┌'
            } else if idx + 1 == child_anchor_rows.len() {
                '└'
            } else {
                '├'
            };
            lines[r][1] = '─';
            lines[r][2] = '─';
        }

        if !child_anchor_rows.is_empty() {
            let min_r = *child_anchor_rows.first().unwrap();
            let max_r = *child_anchor_rows.last().unwrap();
            for r in min_r..=max_r {
                if r != anchor_row && !child_anchor_rows.contains(&r) {
                    lines[r][0] = '│';
                }
            }
            if lines[anchor_row][0] == ' ' {
                lines[anchor_row][0] = '│';
            }
        }

        /* internal-node label hook reserved here; currently hidden by default. */
        Ok(Self { lines, anchor_row })
    }
    /// returns the maximum rendered row width for padding.

    fn width(&self) -> usize {
        self.lines.iter().map(|x| x.len()).max().unwrap_or(0)
    }
}

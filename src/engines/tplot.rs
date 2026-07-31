/*!Module: text tree plotting for `tplot`, supporting Unicode and plain ASCII
connector styles.

 */
use std::collections::HashMap;

use crate::engines::{
    apo::{ApoBranchChange, ApoChangeKind},
    trees::Tree,
    util,
};

type ApoEdgeMap<'a> = HashMap<(usize, usize), Vec<&'a ApoBranchChange>>;

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

/// renders a tree into a standalone SVG document using the same character-grid
/// layout as `render_tree`.
pub fn render_tree_svg(
    tree: &Tree,
    taxon_names: &[String],
    outgroups: Option<&[usize]>,
    style: TPlotStyle,
) -> Result<String, String> {
    let mut rooted = tree.clone();
    let _ = outgroups;
    let _ = style;

    util::order_children_for_plot(&mut rooted);

    let layout = Layout::build(&rooted, rooted.root, taxon_names)?;
    let cols = layout.width().max(1);
    let rows = layout.lines.len().max(1);

    let margin = 12.0;
    let cell_w = 9.0;
    let cell_h = 18.0;
    let label_gap = 4.0;
    let font_size = 14.0;
    let width = margin * 2.0 + label_gap + cols as f64 * cell_w;
    let height = margin * 2.0 + rows as f64 * cell_h;

    let mut svg = String::new();
    svg.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.0}\" height=\"{:.0}\" viewBox=\"0 0 {:.0} {:.0}\">\n",
        width, height, width, height
    ));
    svg.push_str("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");
    svg.push_str(
        "<g fill=\"none\" stroke=\"black\" stroke-width=\"3\" stroke-linecap=\"square\">\n",
    );
    for branch_path in connector_paths(&layout, margin, cell_w, cell_h) {
        svg.push_str("<path d=\"");
        svg.push_str(&branch_path);
        svg.push_str("\"/>\n");
    }
    svg.push_str("</g>\n");
    svg.push_str(&format!(
        "<g font-family=\"Arial, DejaVu Sans Mono, Consolas, monospace\" font-size=\"{font_size}\" font-style=\"italic\" font-weight=\"bold\" fill=\"black\">\n"
    ));
    for (r, row) in layout.lines.iter().enumerate() {
        let mut c = 0usize;
        while c < row.len() {
            if row[c] == ' ' || is_connector(row[c]) {
                c += 1;
                continue;
            }

            let start = c;
            let mut label = String::new();
            while c < row.len() && row[c] != ' ' && !is_connector(row[c]) {
                label.push(row[c]);
                c += 1;
            }

            let x = margin + start as f64 * cell_w + label_gap;
            let y = margin + r as f64 * cell_h + cell_h / 2.0;
            svg.push_str(&format!("<text x=\"{x:.1}\" y=\"{y:.1}\" dominant-baseline=\"middle\">"));
            push_xml_escaped(&mut svg, &label);
            svg.push_str("</text>\n");
        }
    }
    svg.push_str("</g>\n</svg>\n");
    Ok(svg)
}

/// Renders one tree as SVG with apomorphy/homoplasy markers drawn on horizontal
/// branches. Black circles are apomorphies and white circles are homoplasies.
pub fn render_tree_apo_svg(
    tree: &Tree,
    taxon_names: &[String],
    annotations: &[ApoBranchChange],
) -> Result<String, String> {
    let mut rooted = tree.clone();
    util::order_children_for_plot(&mut rooted);

    let mut by_edge = ApoEdgeMap::new();
    for a in annotations {
        by_edge.entry((a.parent, a.child)).or_default().push(a);
    }
    for changes in by_edge.values_mut() {
        changes.sort_by_key(|x| x.character);
    }

    let n = rooted.nodes.len();
    let leaf_count = rooted.nodes.iter().filter(|x| x.taxon.is_some()).count().max(1);
    let margin_x = 36.0;
    let margin_y = 58.0;
    let root_stem = 21.0;
    let leaf_gap = 72.0;
    let branch_base = 72.0;
    /* Extra branch length per marker; branch_base supplies most marker spacing. */
    let marker_gap = 6.0;
    let label_gap = 14.0;
    let font_size = 14.0;
    let mark_font_size = 11.0;
    let radius = 6.0;

    let mut y = vec![0.0; n];
    let mut next_leaf = 0usize;
    assign_apo_y(&rooted, rooted.root, margin_y, leaf_gap, &mut next_leaf, &mut y)?;

    let mut x = vec![0.0; n];
    x[rooted.root] = margin_x + root_stem;
    assign_apo_x(&rooted, rooted.root, &by_edge, branch_base, marker_gap, &mut x)?;

    let max_label_chars = taxon_names.iter().map(|s| s.chars().count()).max().unwrap_or(1);
    let max_x = x.iter().copied().fold(margin_x, f64::max);
    let width = max_x + label_gap + max_label_chars as f64 * 9.0 + margin_x;
    let height = margin_y * 2.0 + (leaf_count.saturating_sub(1)) as f64 * leaf_gap + 22.0;

    let mut svg = String::new();
    svg.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.0}\" height=\"{:.0}\" viewBox=\"0 0 {:.0} {:.0}\">\n",
        width, height, width, height
    ));
    svg.push_str("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");

    svg.push_str(
        "<g fill=\"none\" stroke=\"black\" stroke-width=\"3\" stroke-linecap=\"square\">\n",
    );
    svg.push_str(&format!(
        "<path d=\"M {:.1} {:.1} H {:.1}\"/>\n",
        margin_x, y[rooted.root], x[rooted.root]
    ));
    for (parent, node) in rooted.nodes.iter().enumerate() {
        if node.children.is_empty() {
            continue;
        }

        let ymin = node.children.iter().map(|&c| y[c]).fold(f64::INFINITY, f64::min);
        let ymax = node.children.iter().map(|&c| y[c]).fold(f64::NEG_INFINITY, f64::max);
        if ymin.is_finite() && ymax.is_finite() && ymin < ymax {
            svg.push_str(&format!("<path d=\"M {:.1} {:.1} V {:.1}\"/>\n", x[parent], ymin, ymax));
        }

        for &child in &node.children {
            svg.push_str(&format!(
                "<path d=\"M {:.1} {:.1} H {:.1}\"/>\n",
                x[parent], y[child], x[child]
            ));
        }
    }
    svg.push_str("</g>\n");

    svg.push_str(&format!(
        "<g font-family=\"Arial, DejaVu Sans, sans-serif\" font-size=\"{mark_font_size}\" text-anchor=\"middle\" fill=\"black\">\n"
    ));
    for ((parent, child), changes) in &by_edge {
        if changes.is_empty() {
            continue;
        }
        let start = x[*parent];
        let end = x[*child];
        let step = (end - start) / (changes.len() + 1) as f64;
        for (idx, change) in changes.iter().enumerate() {
            let cx = start + step * (idx + 1) as f64;
            let cy = y[*child];
            svg.push_str(&format!(
                "<text x=\"{cx:.1}\" y=\"{:.1}\">{}</text>\n",
                cy - 13.0,
                change.character
            ));
            match change.kind {
                ApoChangeKind::Apomorphy => {
                    svg.push_str(&format!(
                        "<circle cx=\"{cx:.1}\" cy=\"{cy:.1}\" r=\"{radius:.1}\" fill=\"black\" stroke=\"black\"/>\n"
                    ));
                }
                ApoChangeKind::Homoplasy => {
                    svg.push_str(&format!(
                        "<circle cx=\"{cx:.1}\" cy=\"{cy:.1}\" r=\"{radius:.1}\" fill=\"white\" stroke=\"black\" stroke-width=\"2\"/>\n"
                    ));
                }
            }
            svg.push_str(&format!("<text x=\"{cx:.1}\" y=\"{:.1}\">", cy + 23.0));
            push_xml_escaped(&mut svg, &change.state);
            svg.push_str("</text>\n");
        }
    }
    svg.push_str("</g>\n");

    svg.push_str(&format!(
        "<g font-family=\"Arial, DejaVu Sans, sans-serif\" font-size=\"{font_size}\" font-style=\"italic\" font-weight=\"bold\" fill=\"black\">\n"
    ));
    for (node_id, node) in rooted.nodes.iter().enumerate() {
        let Some(taxon) = node.taxon else {
            continue;
        };
        let label = taxon_names.get(taxon).map(String::as_str).unwrap_or("?");
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" dominant-baseline=\"middle\">",
            x[node_id] + label_gap,
            y[node_id]
        ));
        push_xml_escaped(&mut svg, label);
        svg.push_str("</text>\n");
    }
    svg.push_str("</g>\n</svg>\n");
    Ok(svg)
}

/// Assigns y positions from leaf order and places internal nodes midway between
/// their first and last descendant leaves.
fn assign_apo_y(
    tree: &Tree,
    node: usize,
    margin_y: f64,
    leaf_gap: f64,
    next_leaf: &mut usize,
    y: &mut [f64],
) -> Result<f64, String> {
    if tree.nodes[node].taxon.is_some() {
        let here = margin_y + *next_leaf as f64 * leaf_gap;
        *next_leaf += 1;
        y[node] = here;
        return Ok(here);
    }

    if tree.nodes[node].children.is_empty() {
        return Err(format!("internal node {node} has no children"));
    }

    let mut sum = 0.0;
    for &child in &tree.nodes[node].children {
        sum += assign_apo_y(tree, child, margin_y, leaf_gap, next_leaf, y)?;
    }
    let here = sum / tree.nodes[node].children.len() as f64;
    y[node] = here;
    Ok(here)
}

/// Assigns x positions recursively, lengthening only branches that carry
/// apomorphy/homoplasy markers.
fn assign_apo_x(
    tree: &Tree,
    node: usize,
    by_edge: &ApoEdgeMap<'_>,
    branch_base: f64,
    marker_gap: f64,
    x: &mut [f64],
) -> Result<(), String> {
    for &child in &tree.nodes[node].children {
        let n_marks = by_edge.get(&(node, child)).map(|x| x.len()).unwrap_or(0);
        let branch_len =
            if n_marks == 0 { branch_base } else { branch_base + n_marks as f64 * marker_gap };
        x[child] = x[node] + branch_len;
        assign_apo_x(tree, child, by_edge, branch_base, marker_gap, x)?;
    }
    Ok(())
}

/// converts connector characters in a text layout directly into editable SVG
/// path segments.
fn connector_paths(layout: &Layout, margin: f64, cell_w: f64, cell_h: f64) -> Vec<String> {
    let mut out = Vec::new();

    for (r, row) in layout.lines.iter().enumerate() {
        let mut c = 0usize;
        while c < row.len() {
            let Some((start, mut end)) = horizontal_part(row, c) else {
                c += 1;
                continue;
            };

            c += 1;
            while c < row.len() && row[c] == '─' {
                end = horizontal_end(row, c);
                c += 1;
            }

            let y = svg_y(r as i32 * 2 + 1, margin, cell_h);
            let x1 = svg_x(start, margin, cell_w);
            let x2 = svg_x(end, margin, cell_w);
            out.push(format!("M {x1:.1} {y:.1} H {x2:.1}"));
        }
    }

    for c in 0..layout.width() {
        let x = c as i32 * 2 + 1;
        let mut start = None::<i32>;
        let mut end = 0i32;

        for (r, row) in layout.lines.iter().enumerate() {
            let Some((part_start, part_end, cut)) =
                vertical_part(row.get(c).copied().unwrap_or(' '), r)
            else {
                finish_vertical(&mut out, x, &mut start, end, margin, cell_w, cell_h);
                continue;
            };

            if start.is_none() || part_start > end {
                finish_vertical(&mut out, x, &mut start, end, margin, cell_w, cell_h);
                start = Some(part_start);
            }
            end = end.max(part_end);

            if cut {
                let y = r as i32 * 2 + 1;
                if let Some(s) = start {
                    push_vertical(&mut out, x, s, y, margin, cell_w, cell_h);
                }
                start = Some(y);
            }
        }

        finish_vertical(&mut out, x, &mut start, end, margin, cell_w, cell_h);
    }

    out
}

/// returns the horizontal grid span contributed by one connector cell.
fn horizontal_part(row: &[char], c: usize) -> Option<(i32, i32)> {
    let x = c as i32 * 2 + 1;
    match row[c] {
        '─' => Some((horizontal_start(row, c), horizontal_end(row, c))),
        '┌' | '└' | '├' => Some((x, horizontal_end(row, c))),
        _ => None,
    }
}

/// computes the left edge of a horizontal connector, extending into an adjacent
/// vertical connector when needed to avoid a visible gap.
fn horizontal_start(row: &[char], c: usize) -> i32 {
    if c > 0 && is_vertical_connector(row[c - 1]) { c as i32 * 2 - 1 } else { c as i32 * 2 }
}

/// computes the right edge of a horizontal connector, extending into an adjacent
/// vertical connector when needed to avoid a visible gap.
fn horizontal_end(row: &[char], c: usize) -> i32 {
    if row.get(c + 1).copied().map(is_vertical_connector).unwrap_or(false) {
        c as i32 * 2 + 3
    } else {
        c as i32 * 2 + 2
    }
}

/// returns the vertical grid span contributed by one connector cell and whether
/// that cell is an editable branch intersection.
fn vertical_part(ch: char, r: usize) -> Option<(i32, i32, bool)> {
    let y = r as i32 * 2 + 1;
    match ch {
        '│' => Some((y - 1, y + 1, false)),
        '┌' => Some((y, y + 1, true)),
        '└' => Some((y - 1, y, true)),
        '├' => Some((y - 1, y + 1, true)),
        _ => None,
    }
}

/// emits the active vertical path segment, if any, and clears the active start.
fn finish_vertical(
    out: &mut Vec<String>,
    x: i32,
    start: &mut Option<i32>,
    end: i32,
    margin: f64,
    cell_w: f64,
    cell_h: f64,
) {
    if let Some(s) = start.take() {
        push_vertical(out, x, s, end, margin, cell_w, cell_h);
    }
}

/// appends one non-empty vertical SVG path segment.
fn push_vertical(
    out: &mut Vec<String>,
    x: i32,
    start: i32,
    end: i32,
    margin: f64,
    cell_w: f64,
    cell_h: f64,
) {
    if start >= end {
        return;
    }
    let x = svg_x(x, margin, cell_w);
    let y1 = svg_y(start, margin, cell_h);
    let y2 = svg_y(end, margin, cell_h);
    out.push(format!("M {x:.1} {y1:.1} V {y2:.1}"));
}

/// maps an integer connector-grid x coordinate into an SVG x coordinate.
fn svg_x(grid_x: i32, margin: f64, cell_w: f64) -> f64 {
    margin + grid_x as f64 * cell_w / 2.0
}

/// maps an integer connector-grid y coordinate into an SVG y coordinate.
fn svg_y(grid_y: i32, margin: f64, cell_h: f64) -> f64 {
    margin + grid_y as f64 * cell_h / 2.0
}

/// reports whether a character belongs to the text-tree connector alphabet.
fn is_connector(ch: char) -> bool {
    matches!(ch, '─' | '│' | '┌' | '└' | '├')
}

/// reports whether a connector character carries a vertical branch segment.
fn is_vertical_connector(ch: char) -> bool {
    matches!(ch, '│' | '┌' | '└' | '├')
}

/// appends text with XML special characters escaped for SVG output.
fn push_xml_escaped(out: &mut String, text: &str) {
    for ch in text.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(ch),
        }
    }
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

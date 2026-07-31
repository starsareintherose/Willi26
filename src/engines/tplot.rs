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

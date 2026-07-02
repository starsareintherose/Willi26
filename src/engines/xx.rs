/*!Module: interactive `xx` tree editor and tree/character diagnostic renderer.

*/
use std::io::{self, Write};

use crate::engines::{
    ccode::CharConfig,
    dataset::{Dataset, StateSet},
    search::ScoreWorkspace,
    tplot::TPlotStyle,
    trees::Tree,
    util,
    xsteps::analyze_single_character,
};

/// Result of one xx editor command: continue, saved, or quit.
#[derive(Debug, Clone)]
pub enum XxEditResult {
    Cancelled,
    Saved { tree: Tree, char_config: CharConfig },
}
/// drives the edit loop, renders the current tree view, reads commands, and
/// returns either cancellation or saved tree/config changes.

pub fn run_xx_editor(
    mut tree: Tree,
    ds: Dataset,
    mut cfg: CharConfig,
) -> Result<XxEditResult, String> {
    if ds.nchar == 0 {
        return Err("xx dataset has 0 characters".to_string());
    }
    if cfg.len() != ds.nchar {
        return Err(format!("xx: char config nchar={} != dataset nchar={}", cfg.len(), ds.nchar));
    }

    let mut current_char: usize = 0;

    loop {
        render_xx_view(&tree, &ds, &cfg, current_char)?;

        print!("xx> ");
        io::stdout().flush().map_err(|e| e.to_string())?;

        let mut line = String::new();
        let n = io::stdin().read_line(&mut line).map_err(|e| e.to_string())?;
        if n == 0 {
            return Ok(XxEditResult::Cancelled);
        }

        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let outcome = apply_xx_line(line, &mut tree, current_char, &mut cfg)?;
        match outcome {
            XxLineOutcome::Stay { new_char } => {
                current_char = new_char;
            }
            XxLineOutcome::Cancel => {
                return Ok(XxEditResult::Cancelled);
            }
            XxLineOutcome::Save => {
                return Ok(XxEditResult::Saved { tree, char_config: cfg });
            }
        }
    }
}

enum XxLineOutcome {
    Stay { new_char: usize },
    Cancel,
    Save,
}
/// interprets one editor command line for character switches, weight changes,
/// subtree moves, subtree deletion, save, or cancel.

fn apply_xx_line(
    line: &str,
    tree: &mut Tree,
    current_char: usize,
    cfg: &mut CharConfig,
) -> Result<XxLineOutcome, String> {
    let tokens = tokenize_editor_line(line)?;
    if tokens.is_empty() {
        return Ok(XxLineOutcome::Stay { new_char: current_char });
    }

    let mut i = 0usize;
    let mut cur = current_char;

    while i < tokens.len() {
        match tokens[i].as_str() {
            ";" => return Ok(XxLineOutcome::Cancel),
            "=" => return Ok(XxLineOutcome::Save),

            "+" => {
                cfg.chars[cur].additive = true;
                i += 1;
            }
            "-" => {
                cfg.chars[cur].additive = false;
                i += 1;
            }
            "[" => {
                cfg.chars[cur].active = true;
                i += 1;
            }
            "]" => {
                cfg.chars[cur].active = false;
                i += 1;
            }

            tok if tok.starts_with('/') && tok.len() > 1 => {
                let w: u32 =
                    tok[1..].parse().map_err(|_| format!("xx invalid weight token: {tok}"))?;
                if w == 0 {
                    return Err("xx weight must be >= 1".to_string());
                }
                cfg.chars[cur].weight = w;
                i += 1;
            }

            "\\" => {
                if i + 2 >= tokens.len() {
                    return Err("xx: '\\ a b' requires two node numbers".to_string());
                }
                let from: usize = tokens[i + 1]
                    .parse()
                    .map_err(|_| format!("xx invalid source node: {}", tokens[i + 1]))?;
                let to: usize = tokens[i + 2]
                    .parse()
                    .map_err(|_| format!("xx invalid destination node: {}", tokens[i + 2]))?;

                move_subtree(tree, from, to)?;
                i += 3;
            }

            "\\\\" => {
                if i + 1 >= tokens.len() {
                    return Err("xx: '\\\\ n' requires a node number".to_string());
                }
                let node: usize = tokens[i + 1]
                    .parse()
                    .map_err(|_| format!("xx invalid delete node: {}", tokens[i + 1]))?;

                delete_subtree(tree, node)?;
                i += 2;
            }

            tok if tok.chars().all(|c| c.is_ascii_digit()) => {
                let ch: usize =
                    tok.parse().map_err(|_| format!("xx: invalid character number: {tok}"))?;
                if ch >= cfg.len() {
                    return Err(format!("xx character {} out of range (nchar={})", ch, cfg.len()));
                }
                cur = ch;
                i += 1;
            }

            other => {
                return Err(format!("xx unsupported command token: {other}"));
            }
        }
    }

    Ok(XxLineOutcome::Stay { new_char: cur })
}
/// tokenizes the compact xx command syntax into editor tokens.

fn tokenize_editor_line(line: &str) -> Result<Vec<String>, String> {
    let mut out = Vec::<String>::new();
    let chars: Vec<char> = line.chars().collect();
    let mut i = 0usize;

    while i < chars.len() {
        let c = chars[i];

        if c.is_whitespace() {
            i += 1;
            continue;
        }

        match c {
            ';' | '=' | '+' | '-' | '[' | ']' => {
                out.push(c.to_string());
                i += 1;
            }
            '/' => {
                i += 1;
                let start = i;
                while i < chars.len() && chars[i].is_ascii_digit() {
                    i += 1;
                }
                if start == i {
                    return Err("xx '/' must be followed by digits".to_string());
                }
                let num: String = chars[start..i].iter().collect();
                out.push(format!("/{num}"));
            }
            '\\' => {
                if i + 1 < chars.len() && chars[i + 1] == '\\' {
                    out.push("\\\\".to_string());
                    i += 2;
                } else {
                    out.push("\\".to_string());
                    i += 1;
                }
            }
            d if d.is_ascii_digit() => {
                let start = i;
                i += 1;
                while i < chars.len() && chars[i].is_ascii_digit() {
                    i += 1;
                }
                let num: String = chars[start..i].iter().collect();
                out.push(num);
            }
            _ => {
                return Err(format!("xx invalid character in command line: {:?}", c));
            }
        }
    }

    Ok(out)
}
/// prints the current annotated tree, total length, and active
/// character-coding switches.

fn render_xx_view(
    tree: &Tree,
    ds: &Dataset,
    cfg: &CharConfig,
    current_char: usize,
) -> Result<(), String> {
    let lines = render_tree_with_character_labels(tree, ds, current_char, TPlotStyle::Unicode)?;

    let tl = {
        let mut ws = ScoreWorkspace::new();
        crate::engines::search::score_tree(ds, cfg, tree, &mut ws)
    };

    println!();
    for line in lines {
        println!("{line}");
    }

    let chcfg = &cfg.chars[current_char];
    let active_s = if chcfg.active { "[" } else { "]" };
    let additive_s = if chcfg.additive { "+" } else { "-" };

    println!();
    println!("=sv /w +a -na [i ]ni \\mv \\\\rm");
    println!("TL={} char={} weight={} {} {}", tl, current_char, chcfg.weight, active_s, additive_s);
    println!(
        "commands: ; quit-no-save | = save | N switch-char | /W set-weight | + - [ ] | \\ A B | \\\\ N"
    );

    Ok(())
}
/// combines single-character analysis with a deterministic tree layout.

fn render_tree_with_character_labels(
    tree: &Tree,
    ds: &Dataset,
    ch: usize,
    style: TPlotStyle,
) -> Result<Vec<String>, String> {
    let analyzed = analyze_single_character(ds, tree, ch)?;

    let mut rooted = tree.clone();
    util::order_children_for_plot(&mut rooted);

    let layout = build_xx_layout(&rooted, rooted.root, ds, &analyzed.node_states, ch)?;
    Ok(util::apply_charset(layout.lines_to_strings(), style == TPlotStyle::Unicode))
}

#[derive(Debug, Clone)]
struct XxLayout {
    lines: Vec<Vec<char>>,
    anchor_row: usize,
}

impl XxLayout {
    /// returns the maximum width of the mutable layout rows.
    fn width(&self) -> usize {
        self.lines.iter().map(|x| x.len()).max().unwrap_or(0)
    }
    /// converts padded layout rows into trimmed output strings.

    fn lines_to_strings(self) -> Vec<String> {
        let width = self.width();
        let mut out = Vec::with_capacity(self.lines.len());
        for row in self.lines {
            let mut line = String::with_capacity(width);
            for ch in row {
                line.push(ch);
            }
            while line.ends_with(' ') {
                line.pop();
            }
            out.push(line);
        }
        out
    }
}
/// recursively builds an annotated tree plot containing leaf states and
/// internal node state sets.

fn build_xx_layout(
    tree: &Tree,
    node: usize,
    ds: &Dataset,
    node_states: &[u64],
    ch: usize,
) -> Result<XxLayout, String> {
    let n = &tree.nodes[node];

    if let Some(t) = n.taxon {
        let name = ds.taxa.get(t).cloned().unwrap_or_else(|| format!("{t}"));
        let ss = ds.matrix[t * ds.nchar + ch];
        let label = format!("{t}{name}={}", format_stateset_inline(ss));
        return Ok(XxLayout { lines: vec![label.chars().collect()], anchor_row: 0 });
    }

    if n.children.is_empty() {
        let label = format!("{node}={}", format_state_bits_inline(node_states[node]));
        return Ok(XxLayout { lines: vec![label.chars().collect()], anchor_row: 0 });
    }

    let mut child_blocks = Vec::with_capacity(n.children.len());
    for &chd in &n.children {
        child_blocks.push(build_xx_layout(tree, chd, ds, node_states, ch)?);
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

    let mut anchor_row = match (child_anchor_rows.first(), child_anchor_rows.last()) {
        (Some(&a), Some(&b)) => (a + b) / 2,
        _ => 0,
    };

    /*
      If this is a polytomy (more than 2 children), shift the label down by 1 row
      for better readability, but keep it in-bounds.
    */
    if child_anchor_rows.len() > 2 && anchor_row + 1 < total_rows {
        anchor_row += 1;
    }

    let child_width = child_blocks.iter().map(|b| b.width()).max().unwrap_or(0);

    let my_label = format!("{node}={}", format_state_bits_inline(node_states[node]));
    let my_label_chars: Vec<char> = my_label.chars().collect();

    let left_width = 3usize;
    let total_width = left_width + child_width.max(my_label_chars.len());
    let mut lines = vec![vec![' '; total_width]; total_rows];

    for (block_idx, block) in child_blocks.iter().enumerate() {
        let top = if block_idx == 0 { 0 } else { child_anchor_rows[block_idx] - block.anchor_row };

        for (r, src_line) in block.lines.iter().enumerate() {
            let dst_r = top + r;
            for (c, ch) in src_line.iter().enumerate() {
                lines[dst_r][left_width + c] = *ch;
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

    let label_col0 = if left_width >= 2 { left_width - 2 } else { left_width };
    for (i, ch) in my_label_chars.iter().enumerate() {
        let col = label_col0 + i;
        if col >= lines[anchor_row].len() {
            break;
        }

        /* Avoid overwriting the leftmost connector column (0). */
        if col == 0 {
            continue;
        }

        /*
          overwrite the horizontal dashes (cols 1..2) for internal node labels
        */
        lines[anchor_row][col] = *ch;
    }

    Ok(XxLayout { lines, anchor_row })
}
/// renders a dataset state set as a compact leaf label.

fn format_stateset_inline(ss: StateSet) -> String {
    let bits = ss.bits();
    if bits == StateSet::ALL36.bits() {
        return "?".to_string();
    }
    if bits.count_ones() == 1 {
        let idx = bits.trailing_zeros() as u8;
        return util::idx_to_symbol(idx).to_string();
    }

    let mut syms = Vec::new();
    for idx in 0..36u8 {
        if (bits >> idx) & 1 == 1 {
            syms.push(util::idx_to_symbol(idx));
        }
    }
    let inside = syms.into_iter().map(|c| c.to_string()).collect::<Vec<_>>().join("");
    format!("[{inside}]")
}
/// renders an inferred internal-node state mask.

fn format_state_bits_inline(bits: u64) -> String {
    if bits == 0 {
        return "?".to_string();
    }
    if bits == StateSet::ALL36.bits() {
        return "?".to_string();
    }

    let mut syms = Vec::new();
    for idx in 0..36u8 {
        if (bits >> idx) & 1 == 1 {
            syms.push(util::idx_to_symbol(idx));
        }
    }
    if syms.len() == 1 {
        syms[0].to_string()
    } else {
        format!("[{}]", syms.into_iter().map(|c| c.to_string()).collect::<Vec<_>>().join(""))
    }
}
/// detaches one subtree and reattaches it below or beside another node while
/// avoiding cycles.

fn move_subtree(tree: &mut Tree, from: usize, to: usize) -> Result<(), String> {
    if from >= tree.nodes.len() {
        return Err(format!("xx source node {} out of range", from));
    }
    if to >= tree.nodes.len() {
        return Err(format!("xx destination node {} out of range", to));
    }
    if from == tree.root {
        return Err("xx cannot move root".to_string());
    }
    if from == to {
        return Err("xx cannot move a node under itself".to_string());
    }
    if is_descendant(tree, to, from) {
        return Err("xx cannot move a node under its own descendant".to_string());
    }

    let old_parent =
        tree.nodes[from].parent.ok_or_else(|| "xx: source node has no parent".to_string())?;

    /* First detach the moved subtree from its old parent. */
    {
        let siblings = &mut tree.nodes[old_parent].children;
        if let Some(pos) = siblings.iter().position(|&x| x == from) {
            siblings.remove(pos);
        } else {
            return Err("xx source node not found in its parent children".to_string());
        }
    }

    /* Then attach it according to the destination node type. */
    if tree.nodes[to].taxon.is_some() {
        /*
          Make `from` the sister of the destination tip.
          old: parent(..., to, ...)
          new: parent(..., new_internal, ...)
          new_internal -> [to, from]
        */
        let tip_parent =
            tree.nodes[to].parent.ok_or_else(|| "xx: destination tip has no parent".to_string())?;

        /*
          Reject moves that would make `from` an ancestor of its destination
          parent and therefore create a cycle.
        */
        if is_descendant(tree, tip_parent, from) {
            return Err(
                "xx cannot move a node to become sister of one of its own ancestors' child"
                    .to_string(),
            );
        }

        let new_internal_id = tree.nodes.len();
        tree.nodes.push(crate::engines::trees::Node {
            parent: Some(tip_parent),
            children: vec![to, from],
            taxon: None,
        });

        /*
          Replace the destination tip with the new internal node under the old
          tip parent.
        */
        {
            let ch = &mut tree.nodes[tip_parent].children;
            if let Some(pos) = ch.iter().position(|&x| x == to) {
                ch[pos] = new_internal_id;
            } else {
                return Err("xx destination tip not found in parent children".to_string());
            }
        }

        tree.nodes[to].parent = Some(new_internal_id);
        tree.nodes[from].parent = Some(new_internal_id);
    } else {
        /* Internal destinations can receive the moved subtree directly. */
        tree.nodes[to].children.push(from);
        tree.nodes[from].parent = Some(to);
    }

    normalize_after_edit(tree)?;
    Ok(())
}
/// removes a node/subtree from its parent and normalizes the edited tree.

fn delete_subtree(tree: &mut Tree, node: usize) -> Result<(), String> {
    if node >= tree.nodes.len() {
        return Err(format!("xx node {} out of range", node));
    }
    if node == tree.root {
        return Err("xx cannot delete root".to_string());
    }

    let parent = tree.nodes[node].parent.ok_or_else(|| "xx node has no parent".to_string())?;

    {
        let siblings = &mut tree.nodes[parent].children;
        if let Some(pos) = siblings.iter().position(|&x| x == node) {
            siblings.remove(pos);
        } else {
            return Err("xx node not found in parent children".to_string());
        }
    }

    normalize_after_edit(tree)?;
    Ok(())
}
/// tests ancestry to reject cycle-creating edits.

fn is_descendant(tree: &Tree, candidate: usize, ancestor: usize) -> bool {
    let mut cur = Some(candidate);
    while let Some(v) = cur {
        if v == ancestor {
            return true;
        }
        cur = tree.nodes[v].parent;
    }
    false
}
/// removes unreachable nodes, collapses unary internals, and repairs
/// root/parent indices after editing.

fn normalize_after_edit(tree: &mut Tree) -> Result<(), String> {
    if tree.nodes.is_empty() || tree.root >= tree.nodes.len() {
        return Err("xx tree is empty or root invalid".to_string());
    }

    let mut reachable = vec![false; tree.nodes.len()];
    mark_reachable(tree, tree.root, &mut reachable);

    let mut map = vec![usize::MAX; tree.nodes.len()];
    let mut new_nodes = Vec::new();

    for (old_id, node) in tree.nodes.iter().enumerate() {
        if reachable[old_id] {
            map[old_id] = new_nodes.len();
            new_nodes.push(node.clone());
        }
    }

    for node in &mut new_nodes {
        node.parent = node.parent.and_then(|p| {
            let np = map[p];
            if np == usize::MAX { None } else { Some(np) }
        });

        node.children = node
            .children
            .iter()
            .copied()
            .filter_map(|c| {
                let nc = map[c];
                if nc == usize::MAX { None } else { Some(nc) }
            })
            .collect();
    }

    let new_root = map[tree.root];
    if new_root == usize::MAX {
        return Err("xx root became unreachable".to_string());
    }

    tree.nodes = new_nodes;
    tree.root = new_root;

    loop {
        let changed = collapse_unary_nodes(tree)?;
        if !changed {
            break;
        }
    }

    collapse_root_if_unary(tree)?;
    Ok(())
}
/// marks every node reachable from the current root.

fn mark_reachable(tree: &Tree, node: usize, reachable: &mut [bool]) {
    if reachable[node] {
        return;
    }
    reachable[node] = true;
    for &ch in &tree.nodes[node].children {
        mark_reachable(tree, ch, reachable);
    }
}
/// collapses one non-root internal node with a single child and reports
/// whether a change was made.

fn collapse_unary_nodes(tree: &mut Tree) -> Result<bool, String> {
    for v in 0..tree.nodes.len() {
        if v == tree.root {
            continue;
        }
        if tree.nodes[v].taxon.is_some() {
            continue;
        }
        if tree.nodes[v].children.len() != 1 {
            continue;
        }

        let child = tree.nodes[v].children[0];
        let parent = tree.nodes[v]
            .parent
            .ok_or_else(|| "xx unary internal node has no parent".to_string())?;

        {
            let pch = &mut tree.nodes[parent].children;
            if let Some(pos) = pch.iter().position(|&x| x == v) {
                pch[pos] = child;
            } else {
                return Err("xx parent missing unary child".to_string());
            }
        }

        tree.nodes[child].parent = Some(parent);
        tree.nodes[v].children.clear();

        return Ok(true);
    }

    Ok(false)
}
/// repeatedly promotes the only child of a unary root.

fn collapse_root_if_unary(tree: &mut Tree) -> Result<(), String> {
    loop {
        if tree.root >= tree.nodes.len() {
            return Err("xx invalid root".to_string());
        }
        if tree.nodes[tree.root].taxon.is_some() {
            break;
        }
        if tree.nodes[tree.root].children.len() != 1 {
            break;
        }

        let child = tree.nodes[tree.root].children[0];
        tree.nodes[child].parent = None;
        tree.root = child;
    }

    Ok(())
}

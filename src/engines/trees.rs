/*!Module: rooted tree data structures, tree-set storage, parent/child editing,
rerooting, component extraction, and internal storage-string conversion.

 */
/// Named collection of trees with an auto-generated or user-supplied title.
#[derive(Debug, Clone)]
pub struct TreeSet {
    pub title: String,
    pub trees: Vec<Tree>,
}

impl TreeSet {
    /// packages a title and list of trees for VM tree slots.
    pub fn new(title: impl Into<String>, trees: Vec<Tree>) -> Self {
        Self { title: title.into(), trees }
    }
}

/// Single node in a phylogenetic tree: children, parent, taxon binding.
#[derive(Debug, Clone)]
pub struct Node {
    pub parent: Option<usize>,
    pub children: Vec<usize>,
    pub taxon: Option<usize>,
}

impl Node {
    /// creates a terminal node bound to one taxon index.
    pub fn new_leaf(taxon: usize) -> Self {
        Self { parent: None, children: Vec::new(), taxon: Some(taxon) }
    }
    /// creates an unlabeled internal node with no children yet.

    pub fn new_internal() -> Self {
        Self { parent: None, children: Vec::new(), taxon: None }
    }
    /// reports whether a node represents a terminal taxon.

    pub fn is_leaf(&self) -> bool {
        self.taxon.is_some()
    }
}

/// Rooted phylogenetic tree with an arena-based node storage.
#[derive(Debug, Clone)]
pub struct Tree {
    pub nodes: Vec<Node>,
    pub root: usize,
    pub ntax: usize,
}

impl Tree {
    /// creates a tree containing one detached leaf node per taxon.
    pub fn new(ntax: usize) -> Self {
        let mut nodes = Vec::with_capacity(ntax);
        for i in 0..ntax {
            nodes.push(Node::new_leaf(i));
        }
        Self { nodes, root: 0, ntax }
    }
    /// reports whether a node index points to a terminal node.

    pub fn is_leaf(&self, v: usize) -> bool {
        self.nodes[v].is_leaf()
    }
    /// counts descendant terminal taxa below a node.

    fn leaf_count(&self, v: usize) -> usize {
        if self.nodes[v].taxon.is_some() {
            1
        } else {
            self.nodes[v].children.iter().map(|&c| self.leaf_count(c)).sum()
        }
    }
    /// appends an internal node, assigns its children, and updates each child's
    /// parent pointer.

    pub fn add_internal_node(&mut self, children: Vec<usize>) -> usize {
        let id = self.nodes.len();
        let mut n = Node::new_internal();
        n.children = children.clone();
        self.nodes.push(n);
        for &c in &children {
            self.nodes[c].parent = Some(id);
        }
        id
    }
    /// parses the internal whitespace-separated tree storage format into a rooted
    /// tree.

    pub fn from_storage_string(s: &str, ntax: usize) -> std::result::Result<Self, String> {
        let tokens: Vec<&str> = s.split_whitespace().collect();
        if tokens.is_empty() {
            return Err("empty tree".to_string());
        }

        let mut tree = Tree::new(ntax);
        let mut p = 0usize;
        /// local recursive parser for internal tree storage strings.

        fn parse_rec(
            tokens: &[&str],
            p: &mut usize,
            tree: &mut Tree,
            ntax: usize,
        ) -> std::result::Result<usize, String> {
            if *p >= tokens.len() {
                return Err("unexpected eof while parsing tree".to_string());
            }

            let tok = tokens[*p];
            if tok == "(" {
                *p += 1;
                let mut children = Vec::new();

                while *p < tokens.len() && tokens[*p] != ")" {
                    let child = parse_rec(tokens, p, tree, ntax)?;
                    children.push(child);
                }

                if *p >= tokens.len() || tokens[*p] != ")" {
                    return Err("expected ')'".to_string());
                }
                *p += 1;

                if children.len() < 2 {
                    return Err("internal node must have at least 2 children".to_string());
                }

                let id = tree.add_internal_node(children);
                Ok(id)
            } else {
                let tax =
                    tok.parse::<usize>().map_err(|_| format!("invalid taxon token: {tok}"))?;
                if tax >= ntax {
                    return Err(format!("taxon {tax} out of range (ntax={ntax})"));
                }
                *p += 1;
                Ok(tax)
            }
        }

        let root = parse_rec(&tokens, &mut p, &mut tree, ntax)?;
        if p != tokens.len() {
            return Err("extra tokens after tree".to_string());
        }
        tree.root = root;
        tree.nodes[root].parent = None;
        Ok(tree)
    }
    /// serializes a rooted tree to the internal storage format used by tree files
    /// and hashes.

    pub fn to_storage_string(&self) -> String {
        let mut out = String::new();
        self.write_storage_rec(self.root, &mut out);
        out.trim_end().to_string()
    }
    /// recursively writes storage-format tokens.

    fn write_storage_rec(&self, v: usize, out: &mut String) {
        let n = &self.nodes[v];
        if let Some(t) = n.taxon {
            out.push_str(&t.to_string());
            out.push(' ');
            return;
        }

        out.push_str("( ");
        for &c in &n.children {
            self.write_storage_rec(c, out);
        }
        out.push_str(") ");
    }
    /// sorts child order recursively for deterministic output and hashing.

    pub fn canonicalize(&mut self) {
        self.canonicalize_rec(self.root);
    }
    /// canonicalizes one subtree and returns its ordering key.

    fn canonicalize_rec(&mut self, v: usize) -> String {
        if self.is_leaf(v) {
            return format!("{} ", self.nodes[v].taxon.unwrap());
        }

        let children = self.nodes[v].children.clone();
        let mut child_keys: Vec<(String, usize)> =
            children.into_iter().map(|c| (self.canonicalize_rec(c), c)).collect();

        child_keys.sort_by(|a, b| a.0.cmp(&b.0));
        self.nodes[v].children = child_keys.iter().map(|(_, c)| *c).collect();

        let mut s = String::from("( ");
        for (k, _) in &child_keys {
            s.push_str(k);
        }
        s.push_str(") ");
        s
    }
    /// reroots a tree so a single outgroup leaf is sister to the ingroup root.

    pub fn reroot_by_outgroup(&self, outgroup: usize) -> std::result::Result<Tree, String> {
        if outgroup >= self.ntax {
            return Err(format!("outgroup {outgroup} out of range"));
        }

        if !self.nodes[outgroup].is_leaf() {
            return Err(format!("outgroup {outgroup} is not a leaf"));
        }

        let neigh = self.undirected_neighbors();
        if neigh[outgroup].len() != 1 {
            return Err("outgroup leaf does not have degree 1".to_string());
        }

        let attach = neigh[outgroup][0];
        let mut new_tree = Tree::new(self.ntax);
        /// local reroot helper that rebuilds rooted child lists from an undirected
        /// traversal.

        fn build_rooted(
            old: &Tree,
            neigh: &[Vec<usize>],
            cur: usize,
            parent: usize,
            new_tree: &mut Tree,
        ) -> std::result::Result<usize, String> {
            if let Some(t) = old.nodes[cur].taxon {
                return Ok(t);
            }

            let children_old: Vec<usize> =
                neigh[cur].iter().copied().filter(|&x| x != parent).collect();

            if children_old.is_empty() {
                return Err("reroot produced empty internal continuation".to_string());
            }

            let mut new_children = Vec::with_capacity(children_old.len());
            for ch in children_old {
                let nid = build_rooted(old, neigh, ch, cur, new_tree)?;
                new_children.push(nid);
            }

            if new_children.len() == 1 {
                Ok(new_children[0])
            } else {
                Ok(new_tree.add_internal_node(new_children))
            }
        }

        let ing_root = build_rooted(self, &neigh, attach, outgroup, &mut new_tree)?;
        let root = new_tree.add_internal_node(vec![outgroup, ing_root]);
        new_tree.root = root;
        Ok(new_tree)
    }
    /// reroots by a monophyletic outgroup set when possible, otherwise falls back
    /// to the first outgroup.
    pub fn reroot_by_outgroup_set(&self, outgroups: &[usize]) -> std::result::Result<Tree, String> {
        if outgroups.is_empty() {
            return Err("empty outgroup set".to_string());
        }

        let mut og = outgroups.to_vec();
        og.sort_unstable();
        og.dedup();

        for &t in &og {
            if t >= self.ntax {
                return Err(format!("outgroup {t} out of range"));
            }
            if !self.nodes[t].is_leaf() {
                return Err(format!("outgroup {t} is not a leaf"));
            }
        }

        if og.len() == 1 {
            return self.reroot_by_outgroup(og[0]);
        }

        if self.is_monophyletic_taxa_set(&og) {
            self.reroot_by_monophyletic_outgroup_set(&og)
        } else {
            self.reroot_by_outgroup(og[0])
        }
    }
    /// formats a stable Hennig-style parenthetical string that prefers smaller
    /// subtrees on the left.

    pub fn to_hennig_spine_string(&self) -> String {
        fn min_taxon(tr: &Tree, v: usize) -> usize {
            if let Some(t) = tr.nodes[v].taxon {
                t
            } else {
                tr.nodes[v].children.iter().map(|&c| min_taxon(tr, c)).min().unwrap_or(usize::MAX)
            }
        }
        /// local formatter used by `to_hennig_spine_string`.

        fn rec(tr: &Tree, v: usize) -> String {
            if let Some(t) = tr.nodes[v].taxon {
                return t.to_string();
            }

            let ch = &tr.nodes[v].children;
            if ch.len() != 2 {
                return tr.to_hennig_string_from(v);
            }

            let a = ch[0];
            let b = ch[1];
            let ca = tr.leaf_count(a);
            let cb = tr.leaf_count(b);

            let (small, big) = if ca < cb {
                (a, b)
            } else if cb < ca {
                (b, a)
            } else {
                let ma = min_taxon(tr, a);
                let mb = min_taxon(tr, b);
                if ma <= mb { (a, b) } else { (b, a) }
            };

            format!("({} {})", rec(tr, small), rec(tr, big))
        }

        rec(self, self.root)
    }
    /// formats a subtree without the spine ordering heuristic.

    fn to_hennig_string_from(&self, v: usize) -> String {
        let n = &self.nodes[v];
        if let Some(t) = n.taxon {
            return t.to_string();
        }
        let parts = n.children.iter().map(|&c| self.to_hennig_string_from(c)).collect::<Vec<_>>();
        format!("({})", parts.join(" "))
    }
    /// reroots at the split between a monophyletic outgroup clade and the
    /// remaining ingroup.

    fn reroot_by_monophyletic_outgroup_set(
        &self,
        outgroups: &[usize],
    ) -> std::result::Result<Tree, String> {
        let target = self
            .find_monophyletic_node(outgroups)
            .ok_or_else(|| "outgroup set is not monophyletic".to_string())?;

        if target == self.root {
            return Ok(self.clone());
        }

        let parent = self.nodes[target]
            .parent
            .ok_or_else(|| "monophyletic outgroup node has no parent".to_string())?;

        let neigh = self.undirected_neighbors();
        let mut new_tree = Tree::new(self.ntax);

        fn build_rooted(
            old: &Tree,
            neigh: &[Vec<usize>],
            cur: usize,
            parent: usize,
            new_tree: &mut Tree,
        ) -> std::result::Result<usize, String> {
            if let Some(t) = old.nodes[cur].taxon {
                return Ok(t);
            }

            let children_old: Vec<usize> =
                neigh[cur].iter().copied().filter(|&x| x != parent).collect();

            if children_old.is_empty() {
                return Err("reroot produced empty internal continuation".to_string());
            }

            let mut new_children = Vec::with_capacity(children_old.len());
            for ch in children_old {
                let nid = build_rooted(old, neigh, ch, cur, new_tree)?;
                new_children.push(nid);
            }

            if new_children.len() == 1 {
                Ok(new_children[0])
            } else {
                Ok(new_tree.add_internal_node(new_children))
            }
        }

        let out_root = build_rooted(self, &neigh, target, parent, &mut new_tree)?;
        let in_root = build_rooted(self, &neigh, parent, target, &mut new_tree)?;
        let root = new_tree.add_internal_node(vec![out_root, in_root]);
        new_tree.root = root;
        Ok(new_tree)
    }
    /// tests whether the selected taxa form a clade.

    fn is_monophyletic_taxa_set(&self, taxa: &[usize]) -> bool {
        self.find_monophyletic_node(taxa).is_some()
    }
    /// finds the first internal node whose descendant taxa exactly cover a
    /// requested set.

    pub fn find_monophyletic_node(&self, taxa: &[usize]) -> Option<usize> {
        if taxa.is_empty() {
            return None;
        }

        let mut target = vec![false; self.ntax];
        for &t in taxa {
            if t >= self.ntax {
                return None;
            }
            target[t] = true;
        }
        let need = taxa.len();
        /// local clade-count traversal used by monophyly detection.

        fn dfs(
            tree: &Tree,
            v: usize,
            target: &[bool],
            need: usize,
            found: &mut Option<usize>,
        ) -> usize {
            if let Some(t) = tree.nodes[v].taxon {
                return usize::from(target[t]);
            }

            let mut cnt = 0usize;
            for &c in &tree.nodes[v].children {
                cnt += dfs(tree, c, target, need, found);
            }

            if cnt == need && found.is_none() {
                found.replace(v);
            }
            cnt
        }

        let mut found = None;
        dfs(self, self.root, &target, need, &mut found);
        found
    }
    /// returns sorted, deduplicated edges regardless of parent direction.

    pub fn undirected_edges(&self) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        for (u, n) in self.nodes.iter().enumerate() {
            for &v in &n.children {
                out.push(sorted_edge(u, v));
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }
    /// rebuilds a rooted tree from an undirected edge set using an outgroup leaf
    /// as the root anchor.

    pub fn from_undirected_edges_with_outgroup(
        ntax: usize,
        edges: &[(usize, usize)],
        outgroup: usize,
    ) -> Result<Tree, String> {
        if outgroup >= ntax {
            return Err(format!("outgroup {outgroup} out of range"));
        }

        let max_node =
            edges.iter().flat_map(|(a, b)| [*a, *b]).max().unwrap_or(ntax.saturating_sub(1));

        let mut neigh = vec![Vec::<usize>::new(); max_node + 1];
        for &(a, b) in edges {
            neigh[a].push(b);
            neigh[b].push(a);
        }

        if neigh[outgroup].len() != 1 {
            return Err("cannot root tree: outgroup must have degree 1".to_string());
        }

        let attach = neigh[outgroup][0];
        let mut new_tree = Tree::new(ntax);
        /// local helper that reconstructs rooted subtrees from undirected edges.

        fn build(
            cur: usize,
            parent: usize,
            ntax: usize,
            neigh: &[Vec<usize>],
            new_tree: &mut Tree,
        ) -> Result<usize, String> {
            if cur < ntax {
                return Ok(cur);
            }

            let nexts: Vec<usize> = neigh[cur].iter().copied().filter(|&x| x != parent).collect();

            if nexts.is_empty() {
                return Err("empty continuation while rebuilding rooted tree".to_string());
            }

            let mut kids = Vec::new();
            for ch in nexts {
                kids.push(build(ch, cur, ntax, neigh, new_tree)?);
            }

            if kids.len() == 1 { Ok(kids[0]) } else { Ok(new_tree.add_internal_node(kids)) }
        }

        let ing_root = build(attach, outgroup, ntax, &neigh, &mut new_tree)?;
        let root = new_tree.add_internal_node(vec![outgroup, ing_root]);
        new_tree.root = root;
        Ok(new_tree)
    }
    /// builds an adjacency list from parent/child edges.
    fn undirected_neighbors(&self) -> Vec<Vec<usize>> {
        let mut g = vec![Vec::<usize>::new(); self.nodes.len()];
        for (i, n) in self.nodes.iter().enumerate() {
            for &c in &n.children {
                g[i].push(c);
                g[c].push(i);
            }
        }
        g
    }
    /// creates the canonical three-taxon rooted start tree.

    pub fn new_outgrouped(ntax: usize, outgroup: usize, a: usize, b: usize) -> Self {
        assert!(outgroup < ntax);
        assert!(a < ntax && b < ntax);
        assert!(outgroup != a && outgroup != b && a != b);

        let mut tree = Tree::new(ntax);
        let ing = tree.add_internal_node(vec![a, b]);
        let root = tree.add_internal_node(vec![outgroup, ing]);
        tree.root = root;
        tree
    }
    /// subdivides a rooted edge and inserts a new taxon as sister to the old
    /// child.

    pub fn insert_taxon_on_edge(&mut self, parent: usize, child: usize, taxon: usize) -> usize {
        let pos = self.nodes[parent]
            .children
            .iter()
            .position(|&x| x == child)
            .expect("parent must reference child");

        let new_internal = self.add_internal_node(vec![child, taxon]);
        self.nodes[parent].children[pos] = new_internal;
        self.nodes[new_internal].parent = Some(parent);
        new_internal
    }
    /// returns every directed parent-child edge.

    pub fn rooted_edges(&self) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        self.collect_rooted_edges(self.root, &mut out);
        out
    }
    /// recursively appends directed edges below one node.

    fn collect_rooted_edges(&self, v: usize, out: &mut Vec<(usize, usize)>) {
        for &c in &self.nodes[v].children {
            out.push((v, c));
            self.collect_rooted_edges(c, out);
        }
    }
}
/// returns an undirected edge with endpoints in ascending order.

fn sorted_edge(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
}

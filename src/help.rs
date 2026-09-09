/*!Module: static help text for the legacy command vocabulary.

*/
#[derive(Debug, Clone, Copy)]
pub struct HelpTopic {
    pub name: &'static str,
    pub summary: &'static str,
    pub details: &'static [&'static str],
}
/// returns the full compile-time help table in display order.

pub fn all_topics() -> &'static [HelpTopic] {
    debug_assert!(TOPICS.iter().all(|t| !t.summary.is_empty()));
    TOPICS
}
/// performs case-insensitive prefix matching for `assist xx;` style filtered
/// help.

pub fn find_topics_by_prefix(prefix: &str) -> Vec<&'static HelpTopic> {
    let p = prefix.trim().to_ascii_lowercase();
    if p.is_empty() {
        return vec![];
    }

    TOPICS.iter().filter(|t| t.name.starts_with(&p)).collect()
}

static TOPICS: &[HelpTopic] = &[
    HelpTopic {
        name: "apo",
        summary: "show or save tree with apomorphy/homoplasy markers",
        details: &[
            "apo n;           show tree n with apomorphy/homoplasy character labels",
            "apo n f;         save tree n as SVG file f with apomorphy/homoplasy markers",
            "apo+;            same with apo, but fast optimization / ACCTRAN",
            "apo-;            same with apo, slow optimization / DELTRAN",
            "                 black dots are apomorphies; white dots are homoplasies",
            "                 if n is /, last tree is chosen",
        ],
    },
    HelpTopic {
        name: "assist",
        summary: "list available commands",
        details: &[
            "assist;         list available commands",
            "assist *;       show full help for all commands",
            "assist xx;      show help for commands starting with 'xx'",
        ],
    },
    HelpTopic {
        name: "batch",
        summary: "turn on batch switch (DOS legacy, not implemented)",
        details: &[
            "batch;          turn on batch switch (DOS legacy; not implemented)",
            "batch -;        turn it off (DOS legacy, not implemented)",
        ],
    },
    HelpTopic {
        name: "bb",
        summary: "produce multiple trees by branch breaking",
        details: &[
            "bb;             produce multiple trees by branch breaking",
            "bb *;           use all available tree space",
        ],
    },
    HelpTopic {
        name: "bytes",
        summary: "display bytes of free ram",
        details: &["bytes;          display bytes of free ram"],
    },
    HelpTopic {
        name: "ccode",
        summary: "control character coding",
        details: &[
            "ccode [opt];    control character coding",
            "                [opt] / set weight [ activate ] deactivate",
            "                [opt] + additive - nonadditive",
            "                [opt] * connect different opts",
            "ccode ;         display codings",
        ],
    },
    HelpTopic {
        name: "cget",
        summary: "set coding from code file x",
        details: &["cget x;         set coding from code file x"],
    },
    HelpTopic {
        name: "ckeep",
        summary: "save current coding in code file x",
        details: &["ckeep x;        save current coding in code file x"],
    },
    HelpTopic {
        name: "display",
        summary: "short listings to display",
        details: &[
            "display;        short listings to display",
            "display -;      no listing to display",
            "display *;      all listing to display",
        ],
    },
    HelpTopic {
        name: "erase",
        summary: "delete tree files in scope s",
        details: &["erase  s;       delete tree files in scope s"],
    },
    HelpTopic {
        name: "files",
        summary: "display directory of treefiles",
        details: &["files;          display directory of treefiles"],
    },
    HelpTopic {
        name: "get",
        summary: "make tree file x current",
        details: &["get x;          make tree file x current"],
    },
    HelpTopic {
        name: "hennig",
        summary: "calculate single tree",
        details: &["hennig;         calculate single tree", "hennig *;       use branch breaking"],
    },
    HelpTopic {
        name: "ie",
        summary: "find trees by implicit enumeration",
        details: &[
            "ie;             find trees by implicit enumeration",
            "ie -;           find just one tree",
            "ie *;           use all available tree space",
        ],
    },
    HelpTopic {
        name: "keep",
        summary: "save current tree file as tree file x",
        details: &["keep x;         save current tree file as tree file x"],
    },
    HelpTopic {
        name: "log",
        summary: "control log file",
        details: &[
            "log  n;         open dos file n as new log file",
            "log  -;         deactivate log file",
            "log  *;         activate it",
            "log  /;         close it",
        ],
    },
    HelpTopic {
        name: "mhennig",
        summary: "calculate multiple trees",
        details: &[
            "mhennig;        calculate multiple trees",
            "mhennig *;      use branch breaking",
        ],
    },
    HelpTopic {
        name: "nelsen",
        summary: "calculate nelson consensus tree",
        details: &["nelsen;         calculate nelson consensus tree"],
    },
    HelpTopic {
        name: "naked",
        summary: "show or hide tree node labels",
        details: &[
            "naked;         show node-label display status",
            "naked =;       hide tree node labels (default)",
            "naked -;       show tree node labels in tplot/ttags SVG",
        ],
    },
    HelpTopic {
        name: "outgroup",
        summary: "control outgroup",
        details: &[
            "outgroup xx;    control outgroup",
            "outgroup = x;   set outgroup to list",
            "outgroup ;      list outgroup alphabetically",
        ],
    },
    HelpTopic {
        name: "procedure",
        summary: "open dos file n as procedure file",
        details: &[
            "procedure n;    open dos file n as procedure file",
            "procedure -;    deactivate procedure file",
            "procedure *;    activate it",
            "procedure /;    close it",
        ],
    },
    HelpTopic {
        name: "quote",
        summary: "copy message",
        details: &["quote;          copy message"],
    },
    HelpTopic {
        name: "reroot",
        summary: "reroot current treefile according to current outgroup",
        details: &["reroot;         current treefile according to current outgroup"],
    },
    HelpTopic {
        name: "resample",
        summary: "resampling support on a target tree",
        details: &[
            "resample boot n; estimate bootstrap support for target-tree nodes",
            "resample boot n from t; use tree t as target (default: from /)",
            "resample jak n; independent-deletion jackknife support",
            "resample jak n delete 36.8; set deletion percentage",
            "resample sym n; symmetric resampling (default delete/up 33%)",
            "resample sym n delete 33; set equal down/up percentage",
            "resample boot n [mh; bb;]; set replicate search strategy",
            "resample boot n [mh*; bb*;]; use unlimited-tree star search steps",
            "resample boot n [ie;], [ie*;], or [ie-;]; use implicit enumeration",
            "                 support labels are stored in ttags for SVG output",
        ],
    },
    HelpTopic {
        name: "steps",
        summary: "display max/min steps per character",
        details: &["steps;          display max/min steps per character"],
    },
    HelpTopic {
        name: "tchoose",
        summary: "select trees in scope s from current tree file",
        details: &[
            "tchoose  s;     select trees in scope s from current tree file",
            "                s accepts n, a.b, and / for the last tree",
        ],
    },
    HelpTopic {
        name: "tlist",
        summary: "display trees in parenthetical notation",
        details: &[
            "tlist;          display trees in parenthetical notation",
            "tlist s;        display trees in scope s",
            "                s accepts n, a.b, and / for the last tree",
        ],
    },
    HelpTopic {
        name: "tplot",
        summary: "produce tree diagrams",
        details: &[
            "tplot;          produce tree diagrams",
            "tplot s;        produce tree diagrams in scope s",
            "                s accepts n, a.b, and / for the last tree",
        ],
    },
    HelpTopic { name: "tread", summary: "read trees", details: &["tread;          read trees"] },
    HelpTopic {
        name: "tsave",
        summary: "save current tree file on dos file n",
        details: &["tsave  n;       save current tree file on dos file n"],
    },
    HelpTopic {
        name: "tsvg",
        summary: "save one tree diagram as SVG file",
        details: &[
            "tsvg n f;       save tree n as SVG file f",
            "                if n is /, last tree is chosen",
        ],
    },
    HelpTopic {
        name: "ttags",
        summary: "store tree node labels and export SVG",
        details: &[
            "ttags;         show ttags status",
            "ttags =;       enable labels and select the last tree if available",
            "ttags -;       clear labels and target tree",
            "ttags +N txt;  write txt to node N of the target tree",
            "ttags & f;     write target tree with labels to SVG file f",
        ],
    },
    HelpTopic {
        name: "txascii",
        summary: "use extended ascii characters in tree plots",
        details: &[
            "txascii;        use extended ascii characters in tree plots",
            "txascii -;      don't",
        ],
    },
    HelpTopic {
        name: "view",
        summary: "inspect dos file n",
        details: &[
            "view  n;        inspect dos file n",
            "view  *;        close and inspect current log file",
        ],
    },
    HelpTopic {
        name: "watch",
        summary: "turn on stopwatch",
        details: &["watch;          turn on stopwatch", "watch -;        turn it off"],
    },
    HelpTopic {
        name: "xread",
        summary: "read character data",
        details: &["xread;          read character data"],
    },
    HelpTopic {
        name: "xsteps",
        summary: "diagnose trees in current tree file",
        details: &[
            "xsteps;         diagnose trees in current tree file",
            "xsteps h;       list possible states for hypothetical ancestors",
            "xsteps c;       list character fits",
            "xsteps m;       list best/worst fits",
            "xsteps l;       list tree lengths",
            "xsteps u;       produce file of distinct trees",
            "xsteps w;       set character weights according to fits",
        ],
    },
    HelpTopic {
        name: "xx",
        summary: "display and modify diagnosed tree",
        details: &[
            "xx              display and modify diagnosed tree",
            "                number switchs to specific character",
            "                /[]+- follow ccode grammar",
            "                \\(branch1) (branch2) move branch",
            "                \\\\(branch) delete a branch",
            "                = save current settings and exit",
            "                ; exit without saving",
        ],
    },
    HelpTopic {
        name: "yama",
        summary: "return to dos",
        details: &["yama;           return to dos"],
    },
];

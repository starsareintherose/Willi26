#![forbid(unsafe_code)]

use std::fs;
use std::io;
use std::path::Path;

fn escape_roff(s: &str) -> String {
    s.replace('\\', "\\\\").replace('-', "\\-")
}

fn main() -> io::Result<()> {
    let manual = fs::read_to_string("docs/manual.md")?;

    let mut code_blocks = Vec::new();
    let mut in_code = false;
    let mut current = String::new();

    for line in manual.lines() {
        if line.trim() == "```" {
            if in_code {
                code_blocks.push(current.clone());
                current.clear();
                in_code = false;
            } else {
                in_code = true;
            }
            continue;
        }

        if in_code {
            current.push_str(line);
            current.push('\n');
        }
    }

    let commands = code_blocks.get(1).map(String::as_str).unwrap_or("");

    let mut out = String::new();

    out.push_str(".TH WILLI26 1 \"2026\" \"Willi26\" \"User Commands\"\n");
    out.push_str(".SH NAME\n");
    out.push_str("willi26 \\- A counterfeit of Hennig86\n");

    out.push_str(".SH SYNOPSIS\n");
    out.push_str(".B willi26\n");
    out.push_str("\n");

    out.push_str(".SH DESCRIPTION\n");
    out.push_str("Willi26 is an interactive parsimony phylogenetic analysis program written in safe Rust without third-party crates.\n");

    out.push_str(".SH COMMANDS\n");
    out.push_str(".nf\n");
    out.push_str(&escape_roff(commands));
    out.push_str(".fi\n");

    out.push_str(".SH AUTHOR\n");
    out.push_str("Guoyi Zhang\n");

    out.push_str(".SH COPYRIGHT\n");
    out.push_str("Copyright \\(co Guoyi Zhang 2026. All rights reserved.\n");

    fs::create_dir_all("man")?;
    fs::write(Path::new("man/willi26.1"), out)?;

    Ok(())
}

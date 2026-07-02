/*!
Module: binary entry point for the Willi26 command interpreter.
 */

#![forbid(unsafe_code)]

mod ast;
mod engines;
mod error;
mod help;
mod io;
mod lexer;
mod parser;
mod span;
mod state;
mod vm;

use std::env;
use std::path::Path;
/// selects interactive REPL mode when no script is provided, or runs a
/// batch/procedure file and reports any lenient command errors after
/// execution.

fn main() {
    let mut args = env::args().skip(1);

    /* Without a file argument the executable behaves like the original DOS tool
    and enters the interactive command loop. */
    let Some(first) = args.next() else {
        let mut vm = vm::Vm::new();
        if let Err(e) = vm.repl() {
            eprintln!("{e}");
            std::process::exit(1);
        }
        return;
    };

    /* With an argument, treat it as a batch/procedure script and keep reporting
    lenient command errors after the script finishes. */
    let mut vm = vm::Vm::new();
    match vm.run_file(Path::new(&first)) {
        Ok(errors) => {
            if !errors.is_empty() {
                eprintln!("completed with {} error(s):", errors.len());
                for e in errors {
                    eprintln!("{e}");
                }
            }
        }
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    }
}

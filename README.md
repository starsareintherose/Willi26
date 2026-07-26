# Willi26: An interactive parsimony phylogenetic analysis program 

An interactive parsimony phylogenetic analysis program, developed in safe Rust without third-party dependencies. Its command-driven interface follows theHennig86 tradition while providing an independent modern implementation.

It provides a command-driven workflow for reading discrete character matrices, searching for most-parsimonious trees, calculating consensus trees, and producing tree lists or simple tree plots.

## Compilation

Requirements:

* Rust toolchain
* Cargo

Build the release version:

```bash
cargo build --release
```

The executable will be generated at `./target/release/willi26` on Linux/macOS or `./target/release/willi26.exe` on Windows.

## Usage

For practical usage, workflows, and command examples, see [here](docs/practical.md).

For the command reference, see [here](docs/manual.md).


## Notes

Legacy DOS-related commands, e.g. `batch`, are accepted as placeholders but not implemented.

For reproducible analyses, keep the full command sequence in a procedure file and run it using the following syntax in willi26:

```text
procedure filename;
```

or run in terminal directly:

```bash
<willi26 binary> filename
```

## Development Documentation

Generate the Rust API documentation:

```bash
cargo doc
```

## License

Willi26 is licensed under the GNU Affero General Public License v3.0 or later.

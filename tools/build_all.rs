#![forbid(unsafe_code)]

use std::env;
use std::ffi::OsString;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const BIN_NAME: &str = "willi26";
const DIST_DIR: &str = "bin";
const OSXCROSS_BIN: &str = "/opt/osxcross/bin";

#[derive(Debug)]
struct BuildError {
    msg: String,
}

impl BuildError {
    fn new(msg: impl Into<String>) -> Self {
        Self { msg: msg.into() }
    }
}

impl From<io::Error> for BuildError {
    fn from(e: io::Error) -> Self {
        Self::new(e.to_string())
    }
}

type Result<T> = std::result::Result<T, BuildError>;

fn prepend_path(dir: &str) -> OsString {
    let old_path = env::var_os("PATH").unwrap_or_default();
    let mut new_path = OsString::from(dir);
    new_path.push(":");
    new_path.push(old_path);
    new_path
}

fn command_exists(program: &str, path_env: &OsString) -> bool {
    if program.contains('/') {
        return Path::new(program).exists();
    }

    let Some(paths) = env::split_paths(path_env).next() else {
        return false;
    };

    drop(paths);

    for dir in env::split_paths(path_env) {
        let candidate = dir.join(program);
        if candidate.exists() {
            return true;
        }
    }

    false
}

fn run_command(mut cmd: Command) -> Result<()> {

    let status = cmd
        .stdin(Stdio::inherit())
        .stdout(Stdio::inherit())
        .stderr(Stdio::inherit())
        .status()?;

    if !status.success() {
        return Err(BuildError::new(format!(
            "command failed with status: {}",
            status
        )));
    }

    Ok(())
}

fn run_command_allow_fail(mut cmd: Command) {

    match cmd
        .stdin(Stdio::inherit())
        .stdout(Stdio::inherit())
        .stderr(Stdio::inherit())
        .status()
    {
        Ok(status) if status.success() => {}
        Ok(status) => {
            eprintln!("warning: command failed but ignored: {}", status);
        }
        Err(e) => {
            eprintln!("warning: failed to run command but ignored: {}", e);
        }
    }
}

fn rustup_target_add(target: &str, path_env: &OsString) -> Result<()> {
    let mut cmd = Command::new("rustup");
    cmd.env("PATH", path_env);
    cmd.arg("target").arg("add").arg(target);
    run_command(cmd)
}

fn add_targets(path_env: &OsString) -> Result<()> {
    let targets = [
        "x86_64-unknown-linux-musl",
        "aarch64-unknown-linux-gnu",
        "x86_64-pc-windows-gnu",
        "x86_64-apple-darwin",
        "aarch64-apple-darwin",
    ];

    for target in targets {
        rustup_target_add(target, path_env)?;
    }

    Ok(())
}

fn target_env_name(target: &str) -> String {
    target
        .chars()
        .map(|c| match c {
            'a'..='z' => c.to_ascii_uppercase(),
            '-' => '_',
            other => other,
        })
        .collect()
}

fn strip_if_exists(
    strip_bin: &str,
    file_path: &Path,
    strip_args: &[&str],
    path_env: &OsString,
) {
    if command_exists(strip_bin, path_env) {
        println!("==> Stripping with {}", strip_bin);

        let mut cmd = Command::new(strip_bin);
        cmd.env("PATH", path_env);
        cmd.args(strip_args);
        cmd.arg(file_path);

        run_command_allow_fail(cmd);
    } else {
        println!("==> Strip tool not found: {}, skipping", strip_bin);
    }
}

fn human_size(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "K", "M", "G", "T"];

    let mut size = bytes as f64;
    let mut unit = 0usize;

    while size >= 1024.0 && unit < UNITS.len() - 1 {
        size /= 1024.0;
        unit += 1;
    }

    if unit == 0 {
        format!("{}{}", bytes, UNITS[unit])
    } else if size >= 10.0 {
        format!("{:.0}{}", size, UNITS[unit])
    } else {
        format!("{:.1}{}", size, UNITS[unit])
    }
}

fn show_binary(file_path: &Path, path_env: &OsString) -> Result<()> {
    if command_exists("file", path_env) {
        let mut cmd = Command::new("file");
        cmd.env("PATH", path_env);
        cmd.arg(file_path);
        run_command(cmd)?;
    } else {
        println!("==> command `file` not found, skipping file type display");
    }

    let meta = fs::metadata(file_path)?;
    println!(
        "{} {}",
        human_size(meta.len()),
        file_path.display()
    );

    Ok(())
}

fn copy_to_bin(out: &Path, bin: &Path) -> Result<()> {
    fs::copy(out, bin)?;
    Ok(())
}

fn cargo_build(
    target: &str,
    envs: &[(&str, &str)],
    path_env: &OsString,
) -> Result<()> {
    let mut cmd = Command::new("cargo");
    cmd.env("PATH", path_env);

    for (k, v) in envs {
        cmd.env(k, v);
    }

    cmd.arg("build")
        .arg("--release")
        .arg("--target")
        .arg(target);

    run_command(cmd)
}

fn build_linux(
    arch: &str,
    target: &str,
    linker: &str,
    strip_bin: &str,
    rustflags: &str,
    path_env: &OsString,
) -> Result<()> {
    println!("\n==> Building Linux {}", arch);

    let env_target = target_env_name(target);
    let linker_key = format!("CARGO_TARGET_{}_LINKER", env_target);

    let out = PathBuf::from(format!("target/{}/release/{}", target, BIN_NAME));
    let bin = PathBuf::from(format!("{}/{}-linux-{}", DIST_DIR, BIN_NAME, arch));

    cargo_build(
        target,
        &[
            (&linker_key, linker),
            ("RUSTFLAGS", rustflags),
        ],
        path_env,
    )?;

    strip_if_exists(strip_bin, &out, &[], path_env);

    copy_to_bin(&out, &bin)?;
    show_binary(&bin, path_env)?;

    Ok(())
}

fn build_windows_gnu(
    arch: &str,
    target: &str,
    linker: &str,
    strip_bin: &str,
    path_env: &OsString,
) -> Result<()> {
    println!("\n==> Building Windows {} static", arch);

    let env_target = target_env_name(target);
    let linker_key = format!("CARGO_TARGET_{}_LINKER", env_target);

    let out = PathBuf::from(format!("target/{}/release/{}.exe", target, BIN_NAME));
    let bin = PathBuf::from(format!("{}/{}-windows-{}.exe", DIST_DIR, BIN_NAME, arch));

    cargo_build(
        target,
        &[
            (&linker_key, linker),
            (
                "RUSTFLAGS",
                "-C target-feature=+crt-static -C link-args=-static -C link-arg=-s",
            ),
        ],
        path_env,
    )?;

    strip_if_exists(strip_bin, &out, &[], path_env);

    copy_to_bin(&out, &bin)?;
    show_binary(&bin, path_env)?;

    Ok(())
}

struct MacosBuild<'a> {
    arch: &'a str,
    target: &'a str,
    deployment_target: &'a str,
    cc_bin: &'a str,
    ar_bin: &'a str,
    strip_bin: &'a str,
}

fn build_macos_darwin(config: MacosBuild<'_>, path_env: &OsString) -> Result<()> {
    println!("\n==> Building macOS {}", config.arch);

    let out = PathBuf::from(format!("target/{}/release/{}", config.target, BIN_NAME));
    let bin = PathBuf::from(format!(
        "{}/{}-macos-{}",
        DIST_DIR, BIN_NAME, config.arch
    ));

    match config.target {
        "x86_64-apple-darwin" => {
            cargo_build(
                config.target,
                &[
                    ("MACOSX_DEPLOYMENT_TARGET", config.deployment_target),
                    ("CC_x86_64_apple_darwin", config.cc_bin),
                    ("AR_x86_64_apple_darwin", config.ar_bin),
                    ("CARGO_TARGET_X86_64_APPLE_DARWIN_LINKER", config.cc_bin),
                    ("CARGO_TARGET_X86_64_APPLE_DARWIN_AR", config.ar_bin),
                    ("RUSTFLAGS", "-C link-arg=-Wl,-dead_strip"),
                ],
                path_env,
            )?;
        }
        "aarch64-apple-darwin" => {
            cargo_build(
                config.target,
                &[
                    ("MACOSX_DEPLOYMENT_TARGET", config.deployment_target),
                    ("CC_aarch64_apple_darwin", config.cc_bin),
                    ("AR_aarch64_apple_darwin", config.ar_bin),
                    ("CARGO_TARGET_AARCH64_APPLE_DARWIN_LINKER", config.cc_bin),
                    ("CARGO_TARGET_AARCH64_APPLE_DARWIN_AR", config.ar_bin),
                    ("RUSTFLAGS", "-C link-arg=-Wl,-dead_strip"),
                ],
                path_env,
            )?;
        }
        other => {
            return Err(BuildError::new(format!(
                "unsupported macOS target: {}",
                other
            )));
        }
    }

    strip_if_exists(config.strip_bin, &out, &["-x"], path_env);

    copy_to_bin(&out, &bin)?;
    show_binary(&bin, path_env)?;

    Ok(())
}

fn list_bin() -> Result<()> {
    println!();
    println!("==> Finished");

    let mut entries = Vec::new();

    for entry in fs::read_dir(DIST_DIR)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_file() {
            let meta = fs::metadata(&path)?;
            entries.push((path, meta.len()));
        }
    }

    entries.sort_by(|a, b| a.0.cmp(&b.0));

    let total: u64 = entries.iter().map(|(_, size)| *size).sum();

    println!("total {}", human_size(total));

    for (path, size) in entries {
        println!("{:>6} {}", human_size(size), path.display());
    }

    Ok(())
}

fn main_result() -> Result<()> {
    let path_env = prepend_path(OSXCROSS_BIN);

    fs::create_dir_all(DIST_DIR)?;

    add_targets(&path_env)?;

    build_linux(
        "x86_64",
        "x86_64-unknown-linux-musl",
        "cc",
        "strip",
        "-C target-feature=+crt-static -C link-args=-static -C link-arg=-s",
        &path_env,
    )?;

   /* requirement: aarch64-linux-gnu-gcc */
    build_linux(
        "aarch64",
        "aarch64-unknown-linux-gnu",
        "aarch64-linux-gnu-gcc",
        "aarch64-linux-gnu-strip",
        "-C link-arg=-s",
        &path_env,
    )?;

    /* requirement: mingw-w64-gcc */
    build_windows_gnu(
        "x86_64",
        "x86_64-pc-windows-gnu",
        "x86_64-w64-mingw32-gcc",
        "x86_64-w64-mingw32-strip",
        &path_env,
    )?;

    /* requirement: osxcross */
    build_macos_darwin(
        MacosBuild {
            arch: "x86_64",
            target: "x86_64-apple-darwin",
            deployment_target: "10.13",
            cc_bin: "o64-clang",
            ar_bin: "x86_64-apple-darwin25.1-ar",
            strip_bin: "x86_64-apple-darwin25.1-strip",
        },
        &path_env,
    )?;

    /* requirement: osxcross */
    build_macos_darwin(
        MacosBuild {
            arch: "aarch64",
            target: "aarch64-apple-darwin",
            deployment_target: "11.0",
            cc_bin: "oa64-clang",
            ar_bin: "aarch64-apple-darwin25.1-ar",
            strip_bin: "aarch64-apple-darwin25.1-strip",
        },
        &path_env,
    )?;

    list_bin()?;

    Ok(())
}

fn main() {
    if let Err(e) = main_result() {
        eprintln!("error: {}", e.msg);
        std::process::exit(1);
    }
}

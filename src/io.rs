/*!Module: script-file loading and include-path resolution for batch and
procedure execution.

 */
use std::{
    fs,
    path::{Path, PathBuf},
};

use crate::error::{Error, Result};

/// Loaded script content with its display path for error messaging.
#[derive(Debug, Clone)]
pub struct Source {
    pub file: Option<String>,
    pub text: String,
    pub base_dir: PathBuf,
}
/// reads a script from disk and records its display path plus the base
/// directory used for relative includes.

pub fn load_source(path: &Path) -> Result<Source> {
    let text =
        fs::read_to_string(path).map_err(|e| Error::io(e, Some(path.display().to_string())))?;
    let base_dir = path.parent().unwrap_or_else(|| Path::new(".")).to_path_buf();
    Ok(Source { file: Some(path.display().to_string()), text, base_dir })
}
/// resolves a procedure/include path against the current source directory
/// unless it is already absolute.

pub fn resolve_included_path(base_dir: &Path, include: &Path) -> PathBuf {
    if include.is_absolute() { include.to_path_buf() } else { base_dir.join(include) }
}

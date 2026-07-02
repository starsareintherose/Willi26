/*!Module: character-coding configuration shared by search, diagnostics, and
the interactive tree editor.

 */
/// Per-character configuration: active/inactive, additive/non-additive,
/// weight.
#[derive(Debug, Clone)]
pub struct CharConfig {
    pub chars: Vec<CharSetting>, /* index = character number (0-based) */
}

/// One character's coding: active flag, additive flag, weight value.
#[derive(Debug, Clone, Copy)]
pub struct CharSetting {
    pub active: bool,
    pub additive: bool,
    pub weight: u32,
}

impl CharConfig {
    /// creates one default active, non-additive, unit-weight setting for each
    /// character in a dataset.
    pub fn new(nchar: usize) -> Self {
        Self { chars: vec![CharSetting { active: true, additive: false, weight: 1 }; nchar] }
    }
    /// returns the number of character settings stored in the configuration.

    pub fn len(&self) -> usize {
        self.chars.len()
    }
}

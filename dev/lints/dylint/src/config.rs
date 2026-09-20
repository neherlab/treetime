use serde::Deserialize;

/// Which logging framework to suggest in diagnostics.
///
/// Deserialized from `dylint.toml` as `"tracing"` or `"log"`.
/// Invalid values produce a serde error at config load time.
#[derive(Clone, Copy, Default, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum LogFramework {
    #[default]
    Tracing,
    Log,
}

impl LogFramework {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Tracing => "tracing",
            Self::Log => "log",
        }
    }
}

#[derive(Default, Deserialize)]
#[serde(default)]
pub struct DebugRemnantsConfig {
    /// Which logging framework to suggest: `"tracing"` (default) or `"log"`.
    pub suggested_framework: LogFramework,
}

#[derive(Deserialize)]
#[serde(default)]
pub struct SuggestBuilderConfig {
    pub threshold: usize,
    /// Derive names that exempt a struct from this lint.
    /// Matches the last path segment (e.g. `"Default"` matches both
    /// `#[derive(Default)]` and `#[derive(std::default::Default)]`).
    pub skip_derives: Vec<String>,
}

impl Default for SuggestBuilderConfig {
    fn default() -> Self {
        Self {
            threshold: 6,
            skip_derives: vec![
                "Default".into(),
                "Queryable".into(),
                "Insertable".into(),
                "Selectable".into(),
            ],
        }
    }
}

#[derive(Deserialize)]
#[serde(default)]
pub struct NeedlessBuilderConfig {
    pub threshold: usize,
}

impl Default for NeedlessBuilderConfig {
    fn default() -> Self {
        Self { threshold: 2 }
    }
}

/// Config for the `fallible_new` lint.
#[derive(Deserialize)]
#[serde(default)]
pub struct FallibleNewConfig {
    /// Also lint `fn new_*()` methods, not just `fn new()`.
    pub check_new_variants: bool,
}

impl Default for FallibleNewConfig {
    fn default() -> Self {
        Self {
            check_new_variants: true,
        }
    }
}

/// Config for the `file_too_long` lint.
#[derive(Deserialize)]
#[serde(default)]
pub struct FileLengthConfig {
    /// Maximum number of lines allowed in a source file before the lint fires.
    pub threshold: usize,
}

impl Default for FileLengthConfig {
    fn default() -> Self {
        Self { threshold: 1000 }
    }
}

/// Config for the `no_comments` and `doc_comment_limit` lints.
#[derive(Deserialize)]
#[serde(default)]
pub struct NoCommentsConfig {
    /// Item-level derives whose doc comments render as help or schema text and are kept.
    pub help_derives: Vec<String>,
    /// Enum derives whose per-variant doc comments render as help text and are kept.
    pub help_variant_derives: Vec<String>,
    /// Attribute paths whose doc comments render as help text and are kept.
    pub help_attrs: Vec<String>,
    pub max_chars: usize,
    pub max_lines: usize,
    pub max_paragraphs: usize,
}

impl Default for NoCommentsConfig {
    fn default() -> Self {
        Self {
            help_derives: vec![
                "Parser".into(),
                "Args".into(),
                "Subcommand".into(),
                "JsonSchema".into(),
                "ToSchema".into(),
                "IntoParams".into(),
                "ToResponse".into(),
            ],
            help_variant_derives: vec!["ValueEnum".into()],
            help_attrs: vec!["utoipa::path".into()],
            max_chars: 300,
            max_lines: 3,
            max_paragraphs: 2,
        }
    }
}

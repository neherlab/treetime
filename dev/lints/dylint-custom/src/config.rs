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

/// Config for the `proper_error_type` lint.
#[derive(Default, Deserialize)]
#[serde(default)]
pub struct ProperErrorTypeConfig {
    /// The error types the project expects in public signatures.
    pub style: ErrorStyle,
}

/// Error convention of a project, named in `proper_error_type` diagnostics.
///
/// Deserialized from `dylint.toml` as `"eyre"` or `"typed"`.
#[derive(Clone, Copy, Default, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ErrorStyle {
    /// `eyre::Report` in every crate.
    #[default]
    Eyre,
    /// Typed error enums (`thiserror`) in libraries, `eyre::Report` in binaries.
    Typed,
}

/// Config for the `prefer_error_macros` lint.
#[derive(Default, Deserialize)]
#[serde(default)]
pub struct PreferErrorMacrosConfig {
    /// Project macros that wrap `eyre!`; an `eyre!` expanded from one of them
    /// is generated code and is not reported.
    pub helper_macros: Vec<String>,
    /// Macro suggested in place of `eyre!`, with the same arguments as `eyre!`.
    pub report_macro: Option<String>,
    /// Macro suggested in place of `bail!`, as `return <error_macro>!(...)`.
    pub error_macro: Option<String>,
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
    /// Render sources: each names a tool surface and the doc comments it keeps.
    pub rendered: Vec<RenderSource>,
    /// Non-doc comment prefixes that are kept, matched against trimmed comment
    /// text, e.g. `"SPDX-License-Identifier"`. A line comment on the line right
    /// after a kept line comment is kept too, so a marker heads a multi-line note.
    pub allowed_comment_prefixes: Vec<String>,
    pub max_chars: usize,
    pub max_lines: usize,
    pub max_paragraphs: usize,
}

/// One tool surface that renders doc comments into user-facing output.
#[derive(Clone, Default, Deserialize)]
#[serde(default)]
pub struct RenderSource {
    /// Surface name, e.g. `"clap"`, `"schemars"`, `"tsify"`.
    pub name: String,
    /// Derive names that place an item on this surface (last path segment).
    pub derives: Vec<String>,
    /// Attribute paths that place an item on this surface (e.g. `"wasm_bindgen"`).
    pub attrs: Vec<String>,
    /// Which doc comments this surface renders: item, fields, variants, methods.
    pub renders: Vec<Granularity>,
    /// Attributes that remove a field or variant from this surface, so this
    /// surface does not keep its doc comment. Form `outer(word)`, e.g.
    /// `"arg(skip)"`; `word` matches a bare word inside the attribute list.
    pub skip_attrs: Vec<String>,
}

/// Doc-comment location that a render source can carry into output.
#[derive(Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Granularity {
    Item,
    Fields,
    Variants,
    /// The public associated functions of an `impl` block.
    Methods,
}

impl Default for NoCommentsConfig {
    fn default() -> Self {
        Self {
            rendered: vec![],
            allowed_comment_prefixes: vec![],
            max_chars: 300,
            max_lines: 3,
            max_paragraphs: 2,
        }
    }
}

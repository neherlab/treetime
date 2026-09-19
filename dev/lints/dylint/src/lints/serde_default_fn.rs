//! Flags function-based serde defaults (`#[serde(default = "path")]`), which
//! duplicate default values away from the type's own `Default` impl and drift
//! from it silently.

use rustc_hir::Attribute;
use rustc_lint::{LateContext, LateLintPass, LintContext as _};

use clippy_utils::diagnostics::span_lint_and_help;

rustc_session::declare_lint! {
    /// Flags `#[serde(default = "path")]`. A function-based default lives apart
    /// from the type's `Default` impl and diverges from it without warning; keep
    /// one source of defaults (`#[serde(default)]` plus `Default`/`SmartDefault`).
    pub SERDE_DEFAULT_FN,
    Warn,
    "`#[serde(default = \"...\")]` -- keep defaults in one place via Default/SmartDefault"
}

pub struct SerdeDefaultFn;

impl SerdeDefaultFn {
    pub const fn new() -> Self {
        Self
    }
}

rustc_session::impl_lint_pass!(SerdeDefaultFn => [SERDE_DEFAULT_FN]);

impl<'tcx> LateLintPass<'tcx> for SerdeDefaultFn {
    fn check_attribute(&mut self, cx: &LateContext<'tcx>, attr: &'tcx Attribute) {
        // Only source attributes carry a real span; synthesized attributes such
        // as `PreludeImport` panic on `.span()`.
        if !matches!(attr, Attribute::Unparsed(_)) {
            return;
        }
        let span = attr.span();
        if span.from_expansion() {
            return;
        }
        let Ok(snippet) = cx.sess().source_map().span_to_snippet(span) else {
            return;
        };
        let compact: String = snippet.chars().filter(|c| !c.is_whitespace()).collect();
        if compact.contains("serde(") && compact.contains("default=\"") {
            span_lint_and_help(
                cx,
                SERDE_DEFAULT_FN,
                span,
                "`#[serde(default = \"...\")]` points defaults at a free function",
                None,
                "express defaults once via `Default`/`SmartDefault` and use bare `#[serde(default)]`",
            );
        }
    }
}

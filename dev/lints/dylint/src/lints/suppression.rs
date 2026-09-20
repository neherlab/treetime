use clippy_utils::diagnostics::span_lint_and_help;
use clippy_utils::is_in_test;
use rustc_hir::{Attribute, Expr, HirId};
use rustc_lint::{LateContext, LateLintPass, LintContext as _};

/// Returns `true` if the expression is in a test context:
///
/// - **Test crate** -- compiled with `--test` (integration tests in `tests/`,
///   or `cargo test` on the main crate). Covers test helper functions that
///   don't carry `#[test]` themselves (e.g. `tests/common/mod.rs`).
/// - **Test function** -- `#[test]`, `#[tokio::test]`, `#[rstest]`, etc.
///   Detected via `clippy_utils::is_in_test` (checks `#[rustc_test_marker]`).
/// - **`#[cfg(test)]` module** -- also covered by `is_in_test`.
pub fn is_in_test_zone(cx: &LateContext<'_>, expr: &Expr<'_>) -> bool {
    is_hir_in_test_zone(cx, expr.hir_id)
}

pub fn is_hir_in_test_zone(cx: &LateContext<'_>, hir_id: HirId) -> bool {
    cx.sess().is_test_crate() || is_in_test(cx.tcx, hir_id)
}

rustc_session::declare_lint! {
    /// Flags an `#[allow(...)]` or `#[expect(...)]` lint suppression without a
    /// `reason = "..."`. Every suppression must justify itself, and `#[expect]`
    /// is preferred because rustc reports it when the lint no longer fires,
    /// catching suppressions that match no diagnostic.
    pub UNJUSTIFIED_SUPPRESSION,
    Warn,
    "lint suppression without a reason -- add `reason = \"...\"`, preferring `#[expect]`"
}

pub struct UnjustifiedSuppression;

impl UnjustifiedSuppression {
    pub const fn new() -> Self {
        Self
    }
}

rustc_session::impl_lint_pass!(UnjustifiedSuppression => [UNJUSTIFIED_SUPPRESSION]);

impl<'tcx> LateLintPass<'tcx> for UnjustifiedSuppression {
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
        let is_suppression = compact.contains("allow(") || compact.contains("expect(");
        if is_suppression && !compact.contains("reason=") {
            span_lint_and_help(
                cx,
                UNJUSTIFIED_SUPPRESSION,
                span,
                "lint suppression without a `reason`",
                None,
                "add `reason = \"...\"`, and prefer `#[expect(...)]` so rustc flags it once the lint stops firing",
            );
        }
    }
}

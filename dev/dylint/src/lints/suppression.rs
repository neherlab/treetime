use clippy_utils::is_in_test;
use rustc_hir::Expr;
use rustc_lint::{LateContext, LintContext as _};

/// Returns `true` if the expression is in a test context:
///
/// - **Test crate** -- compiled with `--test` (integration tests in `tests/`,
///   or `cargo test` on the main crate). Covers test helper functions that
///   don't carry `#[test]` themselves (e.g. `tests/common/mod.rs`).
/// - **Test function** -- `#[test]`, `#[tokio::test]`, `#[rstest]`, etc.
///   Detected via `clippy_utils::is_in_test` (checks `#[rustc_test_marker]`).
/// - **`#[cfg(test)]` module** -- also covered by `is_in_test`.
pub fn is_in_test_zone(cx: &LateContext<'_>, expr: &Expr<'_>) -> bool {
    cx.sess().is_test_crate() || is_in_test(cx.tcx, expr.hir_id)
}

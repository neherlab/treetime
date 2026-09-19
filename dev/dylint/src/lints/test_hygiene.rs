//! Test-only hygiene lints: printing to stdout, real sleeps, assertions inside
//! loops, and `#[rstest]` cases without `#[trace]`.

use rustc_data_structures::fx::FxHashSet;
use rustc_hir::intravisit::FnKind;
use rustc_hir::{Body, Expr, ExprKind, FnDecl, HirId, Node};
use rustc_lint::{LateContext, LateLintPass, LintContext as _};
use rustc_span::Span;
use rustc_span::def_id::LocalDefId;

use clippy_utils::diagnostics::span_lint_and_help;

use super::hir_refs::{PanicMacro, find_panic_macro, resolve_expr_def_id};
use super::suppression::is_in_test_zone;

rustc_session::declare_lint! {
    /// Flags `println!` / `print!` in test code. Test output belongs on stderr
    /// (`eprintln!`, `dbg!`) so test runners can capture and display it.
    pub TEST_STDOUT_PRINT,
    Warn,
    "printing to stdout in a test -- use eprintln!/dbg! or an assertion"
}

rustc_session::declare_lint! {
    /// Flags `thread::sleep` / `tokio::time::sleep` in test code. Real sleeps make
    /// tests slow and flaky; synchronize with channels, barriers, or joins.
    pub TEST_REAL_SLEEP,
    Warn,
    "real sleep in a test -- synchronize with a channel, barrier, or join instead"
}

rustc_session::declare_lint! {
    /// Flags an assertion inside a loop in test code. A failing case is reported
    /// once with no index; use a parameterized test (`#[rstest]` cases) instead.
    pub ASSERT_IN_LOOP,
    Warn,
    "assertion inside a loop -- use a parameterized test so each case is reported"
}

rustc_session::declare_lint! {
    /// Flags an `#[rstest]` function without `#[trace]`. Without it, a failing
    /// parameterized case does not print the input values that produced it.
    pub RSTEST_WITHOUT_TRACE,
    Warn,
    "`#[rstest]` without `#[trace]` -- add `#[trace]` so failing cases print inputs"
}

pub struct TestHygiene {
    seen_assert: FxHashSet<Span>,
}

impl TestHygiene {
    pub fn new() -> Self {
        Self {
            seen_assert: FxHashSet::default(),
        }
    }
}

rustc_session::impl_lint_pass!(TestHygiene => [TEST_STDOUT_PRINT, TEST_REAL_SLEEP, ASSERT_IN_LOOP, RSTEST_WITHOUT_TRACE]);

impl<'tcx> LateLintPass<'tcx> for TestHygiene {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        // stdout printing macros.
        if expr.span.from_expansion() {
            let expn = expr.span.ctxt().outer_expn_data();
            if let rustc_span::ExpnKind::Macro(_, name) = expn.kind
                && matches!(name.as_str(), "println" | "print")
                && is_in_test_zone(cx, expr)
            {
                span_lint_and_help(
                    cx,
                    TEST_STDOUT_PRINT,
                    expn.call_site,
                    format!("`{name}!` writes to stdout in a test"),
                    None,
                    "use `eprintln!`/`dbg!` (captured by the test runner) or assert on the value",
                );
            }

            // Assertions inside loops.
            if let Some((call_site, kind)) = find_panic_macro(expr.span)
                && matches!(kind, PanicMacro::Assert | PanicMacro::AssertEq | PanicMacro::AssertNe)
                && self.seen_assert.insert(call_site)
                && expr_in_loop(cx, expr.hir_id)
                && is_in_test_zone(cx, expr)
            {
                span_lint_and_help(
                    cx,
                    ASSERT_IN_LOOP,
                    call_site,
                    "assertion inside a loop",
                    None,
                    "use a parameterized test (`#[rstest]` with `#[case]`) so each case is reported separately",
                );
            }
            return;
        }

        // Real sleeps.
        if let ExprKind::Call(callee, _) = &expr.kind
            && let Some((def_id, _, _)) = resolve_expr_def_id(cx, callee)
        {
            let path = cx.tcx.def_path_str(def_id);
            if (path.ends_with("thread::sleep") || path.ends_with("time::sleep"))
                && is_in_test_zone(cx, expr)
            {
                span_lint_and_help(
                    cx,
                    TEST_REAL_SLEEP,
                    expr.span,
                    "real sleep in a test",
                    None,
                    "synchronize with a channel, barrier, or join handle instead of sleeping",
                );
            }
        }
    }

    fn check_fn(
        &mut self,
        cx: &LateContext<'tcx>,
        _kind: FnKind<'tcx>,
        _decl: &'tcx FnDecl<'tcx>,
        _body: &'tcx Body<'tcx>,
        span: Span,
        def_id: LocalDefId,
    ) {
        if span.from_expansion() {
            return;
        }
        let attrs = cx.tcx.hir_attrs(cx.tcx.local_def_id_to_hir_id(def_id));
        let snippets: Vec<String> = attrs
            .iter()
            .filter(|attr| matches!(attr, rustc_hir::Attribute::Unparsed(_)))
            .filter_map(|attr| cx.sess().source_map().span_to_snippet(attr.span()).ok())
            .collect();
        let has_rstest = snippets.iter().any(|s| s.contains("rstest"));
        let has_trace = snippets.iter().any(|s| s.contains("trace"));
        if has_rstest && !has_trace {
            span_lint_and_help(
                cx,
                RSTEST_WITHOUT_TRACE,
                span,
                "`#[rstest]` function without `#[trace]`",
                None,
                "add `#[trace]` so a failing parameterized case prints its input values",
            );
        }
    }
}

/// Returns `true` if a `loop`/`while`/`for` expression encloses `hir_id` within
/// the same function body.
fn expr_in_loop(cx: &LateContext<'_>, hir_id: HirId) -> bool {
    let mut id = hir_id;
    loop {
        match cx.tcx.parent_hir_node(id) {
            Node::Expr(e) => {
                if matches!(e.kind, ExprKind::Loop(..)) {
                    return true;
                }
                id = e.hir_id;
            }
            Node::Block(b) => id = b.hir_id,
            Node::Stmt(s) => id = s.hir_id,
            Node::LetStmt(l) => id = l.hir_id,
            _ => return false,
        }
    }
}

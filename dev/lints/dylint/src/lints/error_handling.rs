//! Error-handling defect lints: a `Result` bound to `_` and thrown away, and a
//! `Result` error collapsed into a default value.

use rustc_hir::{Expr, ExprKind, PatKind, Stmt, StmtKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::sym;

use clippy_utils::diagnostics::span_lint_and_help;

rustc_session::declare_lint! {
    /// Flags `let _ = <expr>` where the expression is a `Result`. Binding a
    /// fallible value to `_` drops the error without handling or propagating it.
    pub DISCARDED_RESULT,
    Warn,
    "`Result` bound to `_` -- handle the error or propagate it with `?`"
}

rustc_session::declare_lint! {
    /// Flags `.unwrap_or(...)`, `.unwrap_or_default()`, and `.unwrap_or_else(...)`
    /// applied to a `Result`. These discard the error and substitute a value,
    /// hiding the failure.
    pub DEFAULT_MASKS_ERROR,
    Warn,
    "`Result` error replaced by a default -- inspect or propagate the error instead"
}

pub struct ErrorHandling;

impl ErrorHandling {
    pub const fn new() -> Self {
        Self
    }
}

rustc_session::impl_lint_pass!(ErrorHandling => [DISCARDED_RESULT, DEFAULT_MASKS_ERROR]);

impl<'tcx> LateLintPass<'tcx> for ErrorHandling {
    fn check_stmt(&mut self, cx: &LateContext<'tcx>, stmt: &'tcx Stmt<'tcx>) {
        let StmtKind::Let(local) = &stmt.kind else {
            return;
        };
        if local.span.from_expansion() || !matches!(local.pat.kind, PatKind::Wild) {
            return;
        }
        let Some(init) = local.init else {
            return;
        };
        if is_result(cx, init) {
            span_lint_and_help(
                cx,
                DISCARDED_RESULT,
                stmt.span,
                "`Result` value is bound to `_` and its error discarded",
                None,
                "handle the error, propagate it with `?`, or assert success with `.expect(...)`",
            );
        }
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        if expr.span.from_expansion() {
            return;
        }
        let ExprKind::MethodCall(method, receiver, _, span) = &expr.kind else {
            return;
        };
        let name = method.ident.as_str();
        if !matches!(name, "unwrap_or" | "unwrap_or_default" | "unwrap_or_else") {
            return;
        }
        if is_result(cx, receiver) {
            span_lint_and_help(
                cx,
                DEFAULT_MASKS_ERROR,
                *span,
                format!("`.{name}()` replaces a `Result` error with a default value"),
                None,
                "match on the error, log it, or propagate it with `?` instead of discarding it",
            );
        }
    }
}

/// Returns `true` if the adjusted type of `expr` resolves to `Result<_, _>`.
fn is_result<'tcx>(cx: &LateContext<'tcx>, expr: &Expr<'tcx>) -> bool {
    let ty = cx.typeck_results().expr_ty(expr).peel_refs();
    if let ty::Adt(adt, _) = ty.kind() {
        return cx.tcx.is_diagnostic_item(sym::Result, adt.did());
    }
    false
}

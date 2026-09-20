use crate::adt_facts::result_err_ty;
use crate::baseline::emit;
use rustc_hir::{ExprKind, Stmt, StmtKind};
use rustc_lint::{LateContext, LateLintPass};

rustc_session::declare_lint! {
    /// Flags `.ok();` as a statement: it converts the `Result` to an
    /// `Option` and drops it, so the error disappears in a line that looks
    /// like handling. A deliberate discard is written `let _ = ...`, which
    /// states the intent and survives review.
    pub DISCARDED_ERROR,
    Warn,
    "statement-position .ok() drops the error in a line that looks like handling"
}

rustc_session::declare_lint_pass!(DiscardedError => [DISCARDED_ERROR]);

impl<'tcx> LateLintPass<'tcx> for DiscardedError {
    fn check_stmt(&mut self, cx: &LateContext<'tcx>, stmt: &'tcx Stmt<'tcx>) {
        let StmtKind::Semi(expr) = stmt.kind else {
            return;
        };
        let ExprKind::MethodCall(seg, recv, [], _) = expr.kind else {
            return;
        };
        if seg.ident.as_str() != "ok" {
            return;
        }
        let recv_ty = cx.typeck_results().expr_ty_adjusted(recv).peel_refs();
        let Some(err_ty) = result_err_ty(cx.tcx, recv_ty) else {
            return;
        };
        emit(
            cx,
            DISCARDED_ERROR,
            expr.span,
            format!(
                "`.ok();` converts this `Result` to an `Option` and drops it, so the `{err_ty}` error disappears in a line that looks like handling"
            ),
            "handle the `Err` or pass it on with `?`. If dropping it is intended, write `let _ = ...;` to say so",
        );
    }
}

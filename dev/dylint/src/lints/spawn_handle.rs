//! Flags a spawned task or thread whose `JoinHandle` is dropped immediately,
//! so the task is detached and its panic or result is lost.

use rustc_hir::{Expr, PatKind, Stmt, StmtKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;

use clippy_utils::diagnostics::span_lint_and_help;

rustc_session::declare_lint! {
    /// Flags a statement whose value is a `JoinHandle` that is then dropped:
    /// `tokio::spawn(...);` or `let _ = std::thread::spawn(...);`. The handle is
    /// the only way to observe the task's completion, panic, or return value.
    pub SPAWN_HANDLE_DROPPED,
    Warn,
    "spawned task handle dropped -- bind the JoinHandle and await or join it"
}

pub struct SpawnHandle;

impl SpawnHandle {
    pub const fn new() -> Self {
        Self
    }
}

rustc_session::impl_lint_pass!(SpawnHandle => [SPAWN_HANDLE_DROPPED]);

impl<'tcx> LateLintPass<'tcx> for SpawnHandle {
    fn check_stmt(&mut self, cx: &LateContext<'tcx>, stmt: &'tcx Stmt<'tcx>) {
        let dropped_expr = match &stmt.kind {
            // `spawn(...);` -- the handle is dropped at the end of the statement.
            StmtKind::Semi(expr) => *expr,
            // `let _ = spawn(...);` -- explicitly discarded.
            StmtKind::Let(local) if matches!(local.pat.kind, PatKind::Wild) => {
                let Some(init) = local.init else {
                    return;
                };
                init
            }
            _ => return,
        };

        if stmt.span.from_expansion() || !is_join_handle(cx, dropped_expr) {
            return;
        }

        span_lint_and_help(
            cx,
            SPAWN_HANDLE_DROPPED,
            stmt.span,
            "the `JoinHandle` of a spawned task is dropped immediately",
            None,
            "bind the handle and `.await` (tokio) or `.join()` (thread) it so a panic or result is not lost",
        );
    }
}

/// Returns `true` if the type of `expr` is a `JoinHandle` (tokio or std thread).
fn is_join_handle<'tcx>(cx: &LateContext<'tcx>, expr: &Expr<'tcx>) -> bool {
    let ty = cx.typeck_results().expr_ty(expr).peel_refs();
    if let ty::Adt(adt, _) = ty.kind() {
        return cx.tcx.item_name(adt.did()).as_str() == "JoinHandle";
    }
    false
}

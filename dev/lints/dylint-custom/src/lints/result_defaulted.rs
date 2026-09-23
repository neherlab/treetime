use clippy_utils::diagnostics::span_lint_and_help;
use rustc_hir::{Expr, ExprKind, PatKind};
use rustc_lint::{LateContext, LateLintPass};

use crate::lints::hir_refs::{receiver_is_result, result_error_carries_no_cause};
use crate::lints::suppression::is_in_test_zone;

rustc_session::declare_lint! {
    pub RESULT_DEFAULTED,
    Warn,
    "a `Result` error is replaced by a default value and its cause is lost"
}

pub struct ResultDefaulted;

impl ResultDefaulted {
  pub const fn new() -> Self {
    Self
  }
}

rustc_session::impl_lint_pass!(ResultDefaulted => [RESULT_DEFAULTED]);

impl<'tcx> LateLintPass<'tcx> for ResultDefaulted {
  fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
    if expr.span.from_expansion() || is_in_test_zone(cx, expr) {
      return;
    }
    let ExprKind::MethodCall(method, receiver, args, span) = &expr.kind else {
      return;
    };
    let name = method.ident.as_str();
    if !matches!(name, "unwrap_or" | "unwrap_or_else" | "unwrap_or_default" | "ok") {
      return;
    }
    if !receiver_is_result(cx, cx.typeck_results(), receiver) || result_error_carries_no_cause(cx, cx.typeck_results().expr_ty_adjusted(receiver).peel_refs()) {
      return;
    }
    if name == "unwrap_or_else" && closure_binds_error(cx, args) {
      return;
    }
    span_lint_and_help(
      cx,
      RESULT_DEFAULTED,
      *span,
      format!("`.{name}()` replaces the `Err` value with a default, so the cause of the failure is lost"),
      None,
      "handle the error (`match`, `?` with `wrap_err`), or record it before defaulting: \
             `unwrap_or_else(|error| { warn!(%error, ..); default })`",
    );
  }
}

fn closure_binds_error(cx: &LateContext<'_>, args: &[Expr<'_>]) -> bool {
  let [arg] = args else {
    return false;
  };
  let ExprKind::Closure(closure) = arg.kind else {
    return false;
  };
  let body = cx.tcx.hir_body(closure.body);
  body
    .params
    .first()
    .is_some_and(|param| matches!(param.pat.kind, PatKind::Binding(..)))
}

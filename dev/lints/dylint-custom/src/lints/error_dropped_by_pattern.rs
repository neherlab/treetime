use clippy_utils::diagnostics::span_lint_and_help;
use clippy_utils::is_lang_item_or_ctor;
use rustc_hir::def::Res;
use rustc_hir::{Expr, ExprKind, LangItem, LetStmt, Pat, PatKind, Stmt, StmtKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::Span;

use crate::lints::hir_refs::{expr_is_result, result_error_carries_no_cause};
use crate::lints::suppression::{is_hir_in_test_zone, is_in_test_zone};

rustc_session::declare_lint! {
    pub ERROR_DROPPED_BY_PATTERN,
    Warn,
    "a pattern on a `Result` discards the error"
}

const HELP: &str = "bind the error (`Err(error)`, or a `match` with both arms) and record why the operation failed";

pub struct ErrorDroppedByPattern;

impl ErrorDroppedByPattern {
  pub const fn new() -> Self {
    Self
  }
}

rustc_session::impl_lint_pass!(ErrorDroppedByPattern => [ERROR_DROPPED_BY_PATTERN]);

impl<'tcx> LateLintPass<'tcx> for ErrorDroppedByPattern {
  fn check_pat(&mut self, cx: &LateContext<'tcx>, pat: &'tcx Pat<'tcx>) {
    if pat.span.from_expansion() || is_hir_in_test_zone(cx, pat.hir_id) {
      return;
    }
    if result_ctor(cx, pat) == Some(ResultCtor::Err) && subpatterns_are_wild(pat) && error_has_cause(cx, pat) {
      emit(cx, pat.span, "`Err(_)` discards the error");
    }
  }

  fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
    if expr.span.from_expansion() || is_in_test_zone(cx, expr) {
      return;
    }
    let ExprKind::Let(let_expr) = expr.kind else {
      return;
    };
    if result_ctor(cx, let_expr.pat) == Some(ResultCtor::Ok)
      && expr_is_result(cx, let_expr.init)
      && error_has_cause(cx, let_expr.pat)
    {
      emit(
        cx,
        let_expr.pat.span,
        "`if let Ok(..)` discards the error when the pattern fails",
      );
    }
  }

  fn check_stmt(&mut self, cx: &LateContext<'tcx>, stmt: &'tcx Stmt<'tcx>) {
    if stmt.span.from_expansion() || is_hir_in_test_zone(cx, stmt.hir_id) {
      return;
    }
    let StmtKind::Let(LetStmt {
      pat,
      init: Some(init),
      els: Some(_),
      ..
    }) = stmt.kind
    else {
      return;
    };
    if result_ctor(cx, pat) == Some(ResultCtor::Ok) && expr_is_result(cx, init) && error_has_cause(cx, pat) {
      emit(
        cx,
        pat.span,
        "`let Ok(..) = .. else` discards the error in the `else` block",
      );
    }
  }
}

fn error_has_cause<'tcx>(cx: &LateContext<'tcx>, pat: &Pat<'tcx>) -> bool {
  !result_error_carries_no_cause(cx, cx.typeck_results().pat_ty(pat).peel_refs())
}

fn emit(cx: &LateContext<'_>, span: Span, msg: &'static str) {
  span_lint_and_help(cx, ERROR_DROPPED_BY_PATTERN, span, msg, None, HELP);
}

#[derive(PartialEq, Eq)]
enum ResultCtor {
  Ok,
  Err,
}

fn result_ctor(cx: &LateContext<'_>, pat: &Pat<'_>) -> Option<ResultCtor> {
  let PatKind::TupleStruct(qpath, ..) = &pat.kind else {
    return None;
  };
  let Res::Def(_, did) = cx.qpath_res(qpath, pat.hir_id) else {
    return None;
  };
  if is_lang_item_or_ctor(cx, did, LangItem::ResultOk) {
    Some(ResultCtor::Ok)
  } else if is_lang_item_or_ctor(cx, did, LangItem::ResultErr) {
    Some(ResultCtor::Err)
  } else {
    None
  }
}

fn subpatterns_are_wild(pat: &Pat<'_>) -> bool {
  let PatKind::TupleStruct(_, subpatterns, _) = &pat.kind else {
    return false;
  };
  subpatterns.iter().all(|sub| matches!(sub.kind, PatKind::Wild))
}

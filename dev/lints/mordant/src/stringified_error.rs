use crate::adt_facts::result_err_ty;
use crate::baseline::emit;
use crate::hir_shapes::callee_of;
use rustc_hir::{Expr, ExprKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;

rustc_session::declare_lint! {
    /// Flags a `.map_err(|e| e.to_string())` that turns a typed error into a
    /// `String`. From there on, callers cannot match on which failure it
    /// was. `stringly_error` flags the signature that demands this, and this
    /// lint flags the expression that does it.
    pub STRINGIFIED_ERROR,
    Warn,
    "typed error collapsed into a string"
}

rustc_session::declare_lint_pass!(StringifiedError => [STRINGIFIED_ERROR]);

impl<'tcx> LateLintPass<'tcx> for StringifiedError {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let ExprKind::MethodCall(seg, recv, [arg], _) = expr.kind else {
            return;
        };
        if seg.ident.as_str() != "map_err" {
            return;
        }
        let recv_ty = cx.typeck_results().expr_ty_adjusted(recv);
        let Some(err_ty) = result_err_ty(cx.tcx, recv_ty.peel_refs()) else {
            return;
        };
        // Already a string, or opaque to us: nothing is being destroyed.
        let ty::Adt(err_adt, _) = err_ty.peel_refs().kind() else {
            return;
        };
        if cx
            .tcx
            .is_lang_item(err_adt.did(), rustc_hir::LangItem::String)
        {
            return;
        }
        // The call must actually produce a string error, or nothing was
        // collapsed.
        let out_ty = cx.typeck_results().expr_ty(expr);
        let Some(out_err) = result_err_ty(cx.tcx, out_ty) else {
            return;
        };
        let is_string_out = match out_err.peel_refs().kind() {
            ty::Adt(a, _) => cx.tcx.is_lang_item(a.did(), rustc_hir::LangItem::String),
            _ => out_err.peel_refs().is_str(),
        };
        if !is_string_out {
            return;
        }
        let ExprKind::Closure(closure) = arg.kind else {
            return;
        };
        let body = cx.tcx.hir_body(closure.body);
        if closure_body_stringifies(cx, body.value) {
            emit(
                cx,
                STRINGIFIED_ERROR,
                expr.span,
                format!(
                    "this `map_err` turns a `{err_ty}` into a `String`. From here on, callers cannot match on which failure it was"
                ),
                format!(
                    "return `{err_ty}` itself, or an error enum that wraps it, and turn it into text only where it is shown"
                ),
            );
        }
    }
}

/// The closure body is exactly `param.to_string()`, `String::from(param)`, or a
/// `format!` invocation — a stringification and nothing else.
fn closure_body_stringifies<'tcx>(cx: &LateContext<'tcx>, body: &Expr<'tcx>) -> bool {
    // A `format!` body expands to a call, so the macro check runs first.
    let from_format = clippy_utils::macros::macro_backtrace(body.span)
        .next()
        .is_some_and(|mac| clippy_utils::macros::is_format_macro(cx, mac.def_id));
    if from_format {
        return true;
    }
    match body.kind {
        ExprKind::MethodCall(seg, recv, [], _) => {
            seg.ident.as_str() == "to_string" && is_param(recv)
        }
        ExprKind::Call(_, [inner]) => callee_of(cx, body).is_some_and(|c| {
            cx.tcx.def_path_str(c.def()) == "std::string::String::from" && is_param(inner)
        }),
        _ => false,
    }
}

fn is_param(expr: &Expr<'_>) -> bool {
    matches!(expr.kind, ExprKind::Path(_))
}

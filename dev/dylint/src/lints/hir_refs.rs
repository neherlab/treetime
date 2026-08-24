#![allow(
    clippy::wildcard_enum_match_arm,
    reason = "only specific variants are relevant"
)]

//! Shared helpers for resolving references and panic sites from HIR nodes.

use clippy_utils::macros::expn_backtrace;
use rustc_hir::def::Res;
use rustc_hir::{Body, Expr, ExprKind, HirId, Node};
use rustc_lint::LateContext;
use rustc_middle::ty::{self, TyCtxt};
use rustc_span::def_id::DefId;
use rustc_span::{ExpnKind, Span, sym};

/// Resolves the `DefId` of the item referenced by an expression, if any.
///
/// Handles path expressions, struct literals, and method calls.
pub fn resolve_expr_def_id(cx: &LateContext<'_>, expr: &Expr<'_>) -> Option<(DefId, HirId, Span)> {
    let qpath_def = |qpath| match cx.qpath_res(qpath, expr.hir_id) {
        Res::Def(_, def_id) => Some(def_id),
        _ => None,
    };
    let def_id = match &expr.kind {
        ExprKind::Path(qpath) => qpath_def(qpath),
        ExprKind::Struct(qpath, _, _) => qpath_def(qpath),
        ExprKind::MethodCall(..) => cx.typeck_results().type_dependent_def_id(expr.hir_id),
        _ => None,
    };
    def_id.map(|id| (id, expr.hir_id, expr.span))
}

/// Resolves the `DefId` of a type reference, if it refers to a named definition.
pub fn resolve_ty_def_id<'tcx>(
    cx: &LateContext<'tcx>,
    ty: &'tcx rustc_hir::Ty<'tcx, rustc_hir::AmbigArg>,
) -> Option<(DefId, HirId, Span)> {
    if let rustc_hir::TyKind::Path(ref qpath) = ty.kind
        && let Res::Def(_, def_id) = cx.qpath_res(qpath, ty.hir_id)
    {
        Some((def_id, ty.hir_id, ty.span))
    } else {
        None
    }
}

/// Returns `true` if the receiver type of a method call is `Option` or `Result`.
///
/// Accepts explicit `TypeckResults` because callers in `check_impl_item`
/// callbacks may not have body-level typeck results set on the `LateContext`.
pub fn receiver_is_option_or_result<'tcx>(
    cx: &LateContext<'tcx>,
    typeck: &rustc_middle::ty::TypeckResults<'tcx>,
    receiver: &Expr<'tcx>,
) -> bool {
    let recv_ty = typeck.expr_ty_adjusted(receiver).peel_refs();
    if let ty::Adt(adt, _) = recv_ty.kind() {
        let did = adt.did();
        return cx.tcx.is_diagnostic_item(sym::Option, did)
            || cx.tcx.is_diagnostic_item(sym::Result, did);
    }
    false
}

/// If `expr` is a panicking `.unwrap()` or `.expect()` on `Option`/`Result`,
/// returns the method-call span and a short description for diagnostics.
pub fn panicking_unwrap_or_expect<'tcx>(
    cx: &LateContext<'tcx>,
    typeck: &rustc_middle::ty::TypeckResults<'tcx>,
    expr: &Expr<'tcx>,
) -> Option<(Span, &'static str)> {
    let ExprKind::MethodCall(method, receiver, _, span) = &expr.kind else {
        return None;
    };
    let desc = match method.ident.as_str() {
        "unwrap" => ".unwrap()",
        "expect" => ".expect()",
        _ => return None,
    };
    receiver_is_option_or_result(cx, typeck, receiver).then_some((*span, desc))
}

/// If `expr` is a closure that is immediately invoked (e.g. `(|| panic!())()`),
/// returns its body for the caller to walk. Returns `None` for closures stored
/// in a field, passed as a callback, returned, etc. -- those don't run eagerly.
pub fn iife_closure_body<'tcx>(tcx: TyCtxt<'tcx>, expr: &Expr<'_>) -> Option<&'tcx Body<'tcx>> {
    let ExprKind::Closure(closure) = expr.kind else {
        return None;
    };
    let is_iife = matches!(
        tcx.parent_hir_node(expr.hir_id),
        Node::Expr(Expr {
            kind: ExprKind::Call(callee, _),
            ..
        }) if callee.hir_id == expr.hir_id
    );
    is_iife.then(|| tcx.hir_body(closure.body))
}

/// Which panic-family macro was detected.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum PanicMacro {
    Panic,
    Unreachable,
    Assert,
    AssertEq,
    AssertNe,
}

impl PanicMacro {
    /// Human-readable label for diagnostics.
    pub fn desc(self) -> &'static str {
        match self {
            Self::Panic => "panic!()",
            Self::Unreachable => "unreachable!()",
            Self::Assert => "assert!()",
            Self::AssertEq => "assert_eq!()",
            Self::AssertNe => "assert_ne!()",
        }
    }
}

/// Checks if a span originates from a panic-related macro, walking up the
/// expansion chain to handle cases like `panic!` expanding through internal
/// macros (`panic_fmt`, `panic_2021`, etc.).
pub fn find_panic_macro(span: Span) -> Option<(Span, PanicMacro)> {
    expn_backtrace(span).find_map(|(_, data)| {
        let ExpnKind::Macro(_, name) = data.kind else {
            return None;
        };
        let kind = match name.as_str() {
            "panic" => PanicMacro::Panic,
            "unreachable" => PanicMacro::Unreachable,
            "assert" => PanicMacro::Assert,
            "assert_eq" => PanicMacro::AssertEq,
            "assert_ne" => PanicMacro::AssertNe,
            _ => return None,
        };
        Some((data.call_site, kind))
    })
}

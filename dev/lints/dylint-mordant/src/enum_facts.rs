//! Shared analysis for the enum lints: which variant a resolution, a
//! constructor literal or a pattern head names, which of those variants
//! belong to a crate-private enum, and whether a match arm is a panic rather
//! than any other divergence. Every lint that asks "which variant is this"
//! goes through here, so the answer is the same in all of them.

use rustc_hir::def::{CtorKind, CtorOf, DefKind, Res};
use rustc_hir::def_id::DefId;
use rustc_hir::{Expr, ExprKind, Pat, PatExpr, PatExprKind, PatKind, QPath};
use rustc_lint::LateContext;

use crate::hir_shapes::peel_blocks_unsafe;

/// The variant `res` names, whether spelled as the variant itself (unit and
/// struct variants) or as its constructor (tuple variants).
pub(crate) fn variant_of_res(cx: &LateContext<'_>, res: Res) -> Option<DefId> {
    match res {
        Res::Def(DefKind::Variant, id) => Some(id),
        Res::Def(DefKind::Ctor(CtorOf::Variant, _), id) => Some(cx.tcx.parent(id)),
        _ => None,
    }
}

/// The path at the head of a pattern that names something: `E::V`,
/// `E::V(..)` or `E::V { .. }`. Bindings, wildcards, literals, ranges, ors,
/// tuples and the rest have no head path.
pub(crate) fn pat_head_qpath<'h>(pat: &'h Pat<'h>) -> Option<&'h QPath<'h>> {
    match &pat.kind {
        PatKind::TupleStruct(qpath, ..) | PatKind::Struct(qpath, ..) => Some(qpath),
        PatKind::Expr(PatExpr {
            kind: PatExprKind::Path(qpath),
            ..
        }) => Some(qpath),
        _ => None,
    }
}

/// The enum owning `variant`, when that enum is defined in this crate and not
/// reachable from outside it; a variant an outside crate could name is not
/// this crate's to account for.
pub(crate) fn private_enum_of(cx: &LateContext<'_>, variant: DefId) -> Option<DefId> {
    let enum_did = cx.tcx.parent(variant);
    let local = enum_did.as_local()?;
    if cx.effective_visibilities.is_exported(local) {
        return None;
    }
    Some(enum_did)
}

/// The variant whose VALUE `e` is: `E::Unit`, `E::Tuple(..)` or
/// `E::Struct { .. }` written out. None for anything else, including a bare
/// `E::Tuple` passed along as a function value — that expression is a
/// constructor, not a value of the enum, so it passes no variant anywhere.
pub(crate) fn ctor_literal_variant(cx: &LateContext<'_>, e: &Expr<'_>) -> Option<DefId> {
    let res = match &e.kind {
        ExprKind::Call(callee, _) => {
            let ExprKind::Path(qpath) = &callee.kind else {
                return None;
            };
            cx.qpath_res(qpath, callee.hir_id)
        }
        ExprKind::Path(qpath) => match cx.qpath_res(qpath, e.hir_id) {
            Res::Def(DefKind::Ctor(CtorOf::Variant, CtorKind::Fn), _) => return None,
            res => res,
        },
        ExprKind::Struct(qpath, ..) => cx.qpath_res(qpath, e.hir_id),
        _ => return None,
    };
    variant_of_res(cx, res)
}

/// The variant `e` constructs, now or later: everything
/// [`ctor_literal_variant`] accepts, plus a bare `E::Tuple` handed to
/// `map`/`map_err`/`ok_or_else`, which builds the variant when it is called.
/// This is the question "does the crate ever make one of these"; the value
/// question above is the one for "what does this argument pass".
pub(crate) fn constructed_variant(cx: &LateContext<'_>, e: &Expr<'_>) -> Option<DefId> {
    if let ExprKind::Path(qpath) = &e.kind {
        return variant_of_res(cx, cx.qpath_res(qpath, e.hir_id));
    }
    ctor_literal_variant(cx, e)
}

/// The variant a match-arm pattern names at its head.
pub(crate) fn arm_variant(cx: &LateContext<'_>, pat: &Pat<'_>) -> Option<DefId> {
    let qpath = pat_head_qpath(pat)?;
    variant_of_res(cx, cx.qpath_res(qpath, pat.hir_id))
}

/// A diverging arm that is a panic, not a `return`/`continue`: never-typed
/// AND rooted in a panic-family macro.
pub(crate) fn is_panic_arm(cx: &LateContext<'_>, body: &Expr<'_>) -> bool {
    if !cx.typeck_results().expr_ty(body).is_never() {
        return false;
    }
    let inner = peel_blocks_unsafe(body);
    clippy_utils::macros::macro_backtrace(inner.span).any(|mac| {
        matches!(
            cx.tcx.item_name(mac.def_id).as_str(),
            "panic" | "unreachable" | "todo" | "unimplemented"
        )
    })
}

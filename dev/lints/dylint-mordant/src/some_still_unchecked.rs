use clippy_utils::source::snippet;
use clippy_utils::visitors::is_local_used;
use clippy_utils::{SpanlessEq, is_lang_item_or_ctor, is_refutable};
use rustc_hir::def::Res;
use rustc_hir::{
    Arm, BinOpKind, Expr, ExprKind, HirId, LangItem, LetExpr, MatchSource, Pat, PatKind,
};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::Ty;
use rustc_span::{Span, SyntaxContext, sym};

use crate::baseline::emit_with_note;
use crate::enum_facts::pat_head_qpath;
use crate::hir_shapes::sole_expr;

rustc_session::declare_lint! {
    /// Flags a `match` or `if let` over an `Option` or `Result` whose
    /// present case is only taken under a further condition on the value,
    /// with the value that fails it treated exactly like the absent case:
    /// `Some(x) if x.ready() => go(x), _ => wait()`, `if let Some(x) = next
    /// && x.ready() { .. } else { .. }`, or the same with the condition as an
    /// inner `if` whose `else` repeats the outer one. A `Some` that fails
    /// the condition is treated exactly like `None`, so at this point `Some`
    /// alone does not mean usable, and the condition saying which is written
    /// here, at a consumer, rather than where the value is produced. Every
    /// other consumer has to repeat it or takes the unwanted `Some` as a
    /// valid one. Filtering at the source, or a payload type that cannot
    /// hold the unwanted value, makes `Some` mean present.
    ///
    /// Stays quiet when the condition never reads what the pattern bound (a
    /// guard on unrelated state is a different shape), when the failing
    /// value and the absent case are handled by different code, when a
    /// second `Some` arm does its own work, when the fallback does nothing
    /// (`_ => {}`, or an `if let` chain with no `else`), on any enum but
    /// `Option` and `Result`, and inside macros; `matches!(v, Some(x) if ..)`
    /// and `.is_some_and(..)` are conditions, not handling, and are left
    /// alone.
    pub SOME_STILL_UNCHECKED,
    Warn,
    "a guarded `Some` arm whose failures fall through to the `None` handling"
}

rustc_session::declare_lint_pass!(SomeStillUnchecked => [SOME_STILL_UNCHECKED]);

/// The two types this lint reads, each with a present and an absent variant.
#[derive(Clone, Copy)]
enum Wrapper {
    Option,
    Result,
}

impl Wrapper {
    fn of(cx: &LateContext<'_>, ty: Ty<'_>) -> Option<Self> {
        let did = ty.peel_refs().ty_adt_def()?.did();
        if cx.tcx.is_diagnostic_item(sym::Option, did) {
            Some(Wrapper::Option)
        } else if cx.tcx.is_diagnostic_item(sym::Result, did) {
            Some(Wrapper::Result)
        } else {
            None
        }
    }

    fn present(self) -> LangItem {
        match self {
            Wrapper::Option => LangItem::OptionSome,
            Wrapper::Result => LangItem::ResultOk,
        }
    }

    fn absent(self) -> LangItem {
        match self {
            Wrapper::Option => LangItem::OptionNone,
            Wrapper::Result => LangItem::ResultErr,
        }
    }

    /// How the present variant is spelled in a message.
    fn present_name(self) -> &'static str {
        match self {
            Wrapper::Option => "Some",
            Wrapper::Result => "Ok",
        }
    }

    fn absent_name(self) -> &'static str {
        match self {
            Wrapper::Option => "None",
            Wrapper::Result => "Err",
        }
    }

    fn message(self, cond: &str) -> String {
        let (present, absent) = (self.present_name(), self.absent_name());
        let article = match self {
            Wrapper::Option => "a",
            Wrapper::Result => "an",
        };
        format!(
            "{article} `{present}` that fails `{cond}` is treated exactly like `{absent}` here. `{present}` alone does not mean usable, and every consumer has to repeat the check"
        )
    }

    fn note(self) -> String {
        format!(
            "the failing `{}` ends up here, with `{}`",
            self.present_name(),
            self.absent_name()
        )
    }

    fn help(self, scrut: &str, cond: &str) -> String {
        match self {
            Wrapper::Option => format!(
                "filter where `{scrut}` is produced, with `.filter(..)` or a payload type that cannot fail `{cond}`"
            ),
            Wrapper::Result => format!(
                "reject that case where `{scrut}` is produced, with an error for it or a payload type that cannot fail `{cond}`"
            ),
        }
    }
}

/// Whether the pattern's head names the lang variant `item` (`Some(..)`,
/// `None`, `Ok(..)`, `Err(..)`), however the path to it is spelled.
fn head_is(cx: &LateContext<'_>, pat: &Pat<'_>, item: LangItem) -> bool {
    let Some(qpath) = pat_head_qpath(pat) else {
        return false;
    };
    match cx.qpath_res(qpath, pat.hir_id) {
        Res::Def(_, did) => is_lang_item_or_ctor(cx, did, item),
        _ => false,
    }
}

/// The locals bound by `Some(p)` / `Ok(p)` when `p` itself cannot fail, so
/// the guard is the only thing standing between this arm and a present
/// value. None for any other pattern.
fn present_bindings(cx: &LateContext<'_>, pat: &Pat<'_>, w: Wrapper) -> Option<Vec<HirId>> {
    let PatKind::TupleStruct(_, [inner], _) = pat.kind else {
        return None;
    };
    if !head_is(cx, pat, w.present()) || is_refutable(cx, inner) {
        return None;
    }
    let mut ids = Vec::new();
    inner.each_binding(|_, id, _, _| ids.push(id));
    (!ids.is_empty()).then_some(ids)
}

fn reads_any<'tcx>(cx: &LateContext<'tcx>, cond: &'tcx Expr<'tcx>, ids: &[HirId]) -> bool {
    ids.iter().any(|id| is_local_used(cx, cond, *id))
}

/// A pattern that takes whatever the guarded arm let through without
/// looking inside it: `_`, an unread binding, `None`, `Err(..)`, `Some(_)`,
/// or an `|` of those.
fn covers_rest<'tcx>(
    cx: &LateContext<'tcx>,
    pat: &Pat<'_>,
    body: &'tcx Expr<'tcx>,
    w: Wrapper,
) -> bool {
    match pat.kind {
        PatKind::Wild => true,
        PatKind::Binding(_, id, _, None) => !is_local_used(cx, body, id),
        PatKind::Or(pats) => pats.iter().all(|p| covers_rest(cx, p, body, w)),
        PatKind::TupleStruct(_, [inner], _) if head_is(cx, pat, w.present()) => {
            matches!(inner.kind, PatKind::Wild | PatKind::Binding(.., None))
        }
        _ => head_is(cx, pat, w.absent()),
    }
}

/// Whether a rest pattern receives present values (as opposed to naming
/// only the absent variant), so the note can point at the arm the failing
/// `Some` actually reaches.
fn takes_present(cx: &LateContext<'_>, pat: &Pat<'_>, w: Wrapper) -> bool {
    match pat.kind {
        PatKind::Wild | PatKind::Binding(.., None) => true,
        PatKind::Or(pats) => pats.iter().any(|p| takes_present(cx, p, w)),
        _ => head_is(cx, pat, w.present()),
    }
}

/// `{}`, `()` and nests of them: a fallback that handles nothing.
fn does_nothing(e: &Expr<'_>) -> bool {
    match e.kind {
        ExprKind::Block(b, _) => b.stmts.is_empty() && b.expr.is_none_or(does_nothing),
        ExprKind::Tup([]) => true,
        ExprKind::DropTemps(inner) => does_nothing(inner),
        _ => false,
    }
}

/// The operands of an `&&` chain, left to right.
fn and_operands<'h>(e: &'h Expr<'h>, out: &mut Vec<&'h Expr<'h>>) {
    match e.kind {
        ExprKind::Binary(op, l, r) if op.node == BinOpKind::And => {
            and_operands(l, out);
            and_operands(r, out);
        }
        ExprKind::DropTemps(inner) => and_operands(inner, out),
        _ => out.push(e),
    }
}

fn has_let(operands: &[&Expr<'_>]) -> bool {
    operands.iter().any(|e| matches!(e.kind, ExprKind::Let(..)))
}

struct Finding {
    w: Wrapper,
    span: Span,
    /// The `Option`/`Result` being matched, for the help.
    scrut: Span,
    cond: Span,
    lands: Span,
}

fn report(cx: &LateContext<'_>, f: Finding) {
    let cond = snippet(cx, f.cond, "..");
    let scrut = snippet(cx, f.scrut, "..");
    emit_with_note(
        cx,
        SOME_STILL_UNCHECKED,
        f.span,
        f.w.message(&cond),
        f.lands,
        f.w.note(),
        f.w.help(&scrut, &cond),
    );
}

fn check_match<'tcx>(
    cx: &LateContext<'tcx>,
    scrut: &'tcx Expr<'tcx>,
    arms: &'tcx [Arm<'tcx>],
) -> Option<Finding> {
    let w = Wrapper::of(cx, cx.typeck_results().expr_ty(scrut))?;
    let mut guarded = None;
    let mut rest = Vec::new();
    for arm in arms {
        match arm.guard {
            Some(guard) => {
                if guarded.replace((arm, guard)).is_some() {
                    return None;
                }
            }
            None => rest.push(arm),
        }
    }
    let (arm, cond) = guarded?;
    let ids = present_bindings(cx, arm.pat, w)?;
    let mut cond_parts = Vec::new();
    and_operands(cond, &mut cond_parts);
    if has_let(&cond_parts) || !reads_any(cx, cond, &ids) {
        return None;
    }
    let (first, others) = rest.split_first()?;
    if does_nothing(first.body) || !rest.iter().all(|a| covers_rest(cx, a.pat, a.body, w)) {
        return None;
    }
    let mut eq = SpanlessEq::new(cx);
    if !others
        .iter()
        .all(|a| eq.eq_expr(SyntaxContext::root(), first.body, a.body))
    {
        return None;
    }
    let lands = rest
        .iter()
        .find(|a| takes_present(cx, a.pat, w))
        .unwrap_or(first);
    Some(Finding {
        w,
        span: arm.pat.span.to(cond.span),
        scrut: scrut.span,
        cond: cond.span,
        lands: lands.pat.span,
    })
}

/// `let Some(x) = init` at the head of a condition: the wrapper it opens,
/// the locals it binds, and where `init` is.
fn present_let<'tcx>(
    cx: &LateContext<'tcx>,
    e: &Expr<'tcx>,
) -> Option<(Wrapper, Vec<HirId>, Span)> {
    let ExprKind::Let(LetExpr { pat, init, .. }) = e.kind else {
        return None;
    };
    let w = Wrapper::of(cx, cx.typeck_results().expr_ty(init))?;
    Some((w, present_bindings(cx, pat, w)?, init.span))
}

fn check_if<'tcx>(
    cx: &LateContext<'tcx>,
    cond: &'tcx Expr<'tcx>,
    then: &'tcx Expr<'tcx>,
    els: &'tcx Expr<'tcx>,
) -> Option<Finding> {
    if does_nothing(els) {
        return None;
    }
    let mut parts = Vec::new();
    and_operands(cond, &mut parts);
    let (head, tail) = parts.split_first()?;
    let (w, ids, scrut) = present_let(cx, head)?;
    if let (Some(first), Some(last)) = (tail.first(), tail.last()) {
        // `if let Some(x) = v && COND { A } else { B }`
        if has_let(tail) || !tail.iter().any(|e| reads_any(cx, e, &ids)) {
            return None;
        }
        return Some(Finding {
            w,
            span: cond.span,
            scrut,
            cond: first.span.to(last.span),
            lands: els.span,
        });
    }
    // `if let Some(x) = v { if COND { A } else { B } } else { B }`
    let ExprKind::Block(block, None) = then.kind else {
        return None;
    };
    let ExprKind::If(inner_cond, _, Some(inner_els)) = sole_expr(block)?.kind else {
        return None;
    };
    let mut inner_parts = Vec::new();
    and_operands(inner_cond, &mut inner_parts);
    if has_let(&inner_parts)
        || !reads_any(cx, inner_cond, &ids)
        || !SpanlessEq::new(cx).eq_expr(SyntaxContext::root(), inner_els, els)
    {
        return None;
    }
    Some(Finding {
        w,
        span: head.span,
        scrut,
        cond: inner_cond.span,
        lands: inner_els.span,
    })
}

impl<'tcx> LateLintPass<'tcx> for SomeStillUnchecked {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        if expr.span.from_expansion() {
            return;
        }
        let finding = match expr.kind {
            ExprKind::Match(scrut, arms, MatchSource::Normal) => check_match(cx, scrut, arms),
            ExprKind::If(cond, then, Some(els)) => check_if(cx, cond, then, els),
            _ => None,
        };
        if let Some(f) = finding {
            report(cx, f);
        }
    }
}

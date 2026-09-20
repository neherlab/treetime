use std::collections::{HashMap, HashSet};

use clippy_utils::consts::ConstEvalCtxt;
use clippy_utils::res::MaybeResPath;
use clippy_utils::source::snippet;
use rustc_hir::def::{DefKind, Res};
use rustc_hir::def_id::DefId;
use rustc_hir::{BinOpKind, Expr, ExprKind, HirId, Pat, PatKind, UnOp};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::{self, Ty};
use rustc_span::{Span, Symbol, sym};

use crate::adt_facts::has_fixed_repr;
use crate::baseline::emit_with_note;
use crate::hir_shapes::{Callee, callee_of};

rustc_session::declare_lint! {
    /// Flags an integer place (a struct field, or a local or parameter
    /// within one function) that is cut down with a bare `as` at one site,
    /// which wraps silently, while another site converts the same place into
    /// a type no wider with `u32::try_from(x)` or `x.try_into()` because it
    /// might not fit. The place's declared type
    /// is wider than what its readers need, so the range check lives in
    /// whichever reader remembered it.
    ///
    /// A field is one place whatever it is read through, keyed by the struct
    /// that declares it; a local is a place only inside its own body. Only
    /// bare places count: `(x & 0xff) as u8`, `x.min(MAX) as u32` and casts
    /// of call results are computations, not the place, and `x as u32 &
    /// mask` keeps low bits on purpose. An `as` in a function that also
    /// compares the place against a constant (`x <= u32::MAX as usize`, a
    /// range pattern) is that function's checked form and stays quiet, but
    /// it is not evidence against other sites: what a comparison guards is
    /// not in its syntax. `usize` and `isize` have the target's pointer
    /// width, so `u64 as usize` on a 64-bit target is not a narrowing. A
    /// checked site only condemns `as` casts to a type at most as wide as
    /// its own target: `u8::try_from(c)` inside an arm that already matched
    /// a letter says nothing about `c as u32` elsewhere. A field of a
    /// `repr(C)`, packed or transparent struct has its width fixed by that
    /// layout and is not read; neither are sites inside macro expansions.
    pub NARROWED_TWO_WAYS,
    Warn,
    "an integer place range-checked at one narrowing and truncated with `as` at another"
}

/// The place a conversion reads: a field of whichever struct declares it, or
/// a local binding.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Place {
    Field(DefId, Symbol),
    Local(HirId),
}

struct Site {
    span: Span,
    /// The fn (closures folded into their parent) the site sits in, for the
    /// compare-then-cast form.
    body: DefId,
    dst_bits: u64,
    /// `try_from` / `try_into`: the evidence that the place may not fit.
    checked: bool,
    shown: String,
    src: String,
    dst: String,
}

#[derive(Default)]
pub struct NarrowedTwoWays {
    sites: HashMap<Place, Vec<Site>>,
    /// (fn, place): the fn compares the place against a constant, so its
    /// `as` casts of that place are excused (not reported, not evidence).
    range_checked: HashSet<(DefId, Place)>,
}

rustc_session::impl_lint_pass!(NarrowedTwoWays => [NARROWED_TWO_WAYS]);

#[derive(Clone, Copy, PartialEq, Eq)]
struct IntLayout {
    bits: u64,
    signed: bool,
}

/// `usize`/`isize` take the target's pointer width: whether `u64 as usize`
/// narrows is a question about the target being linted, not about all of
/// them.
fn int_layout(cx: &LateContext<'_>, ty: Ty<'_>) -> Option<IntLayout> {
    let ptr_bits = || cx.tcx.data_layout.pointer_size().bits();
    match ty.kind() {
        ty::Int(i) => Some(IntLayout {
            bits: i.bit_width().unwrap_or_else(ptr_bits),
            signed: true,
        }),
        ty::Uint(u) => Some(IntLayout {
            bits: u.bit_width().unwrap_or_else(ptr_bits),
            signed: false,
        }),
        _ => None,
    }
}

/// Some value of `src` has no representation in `dst`.
fn lossy(src: IntLayout, dst: IntLayout) -> bool {
    let preserving = (src.signed == dst.signed && dst.bits >= src.bits)
        || (!src.signed && dst.signed && dst.bits > src.bits);
    !preserving
}

/// The integer place `e` reads, through `&`, `*`, HIR temporaries and inner
/// integer-to-integer casts (which change representation, not which place is
/// read), with the place's own layout and its source text to show for it.
fn place_of<'tcx>(
    cx: &LateContext<'tcx>,
    mut e: &'tcx Expr<'tcx>,
) -> Option<(Place, IntLayout, Ty<'tcx>, String)> {
    let typeck = cx.typeck_results();
    loop {
        match e.kind {
            ExprKind::Cast(inner, _)
                if int_layout(cx, typeck.expr_ty(inner)).is_some()
                    && int_layout(cx, typeck.expr_ty(e)).is_some() =>
            {
                e = inner;
            }
            ExprKind::DropTemps(inner)
            | ExprKind::AddrOf(_, _, inner)
            | ExprKind::Unary(UnOp::Deref, inner) => e = inner,
            _ => break,
        }
    }
    let ty = typeck.expr_ty(e).peel_refs();
    let layout = int_layout(cx, ty)?;
    let shown = snippet(cx, e.span, "..").into_owned();
    match e.kind {
        ExprKind::Field(base, ident) => {
            let adt = typeck.expr_ty_adjusted(base).peel_refs().ty_adt_def()?;
            // A layout-fixed struct cannot redeclare the field narrow.
            if has_fixed_repr(adt) {
                return None;
            }
            Some((Place::Field(adt.did(), ident.name), layout, ty, shown))
        }
        _ => {
            let local = e.res_local_id()?;
            Some((Place::Local(local), layout, ty, shown))
        }
    }
}

fn has_range_pat(pat: &Pat<'_>) -> bool {
    let mut found = false;
    pat.walk_always(|p| found |= matches!(p.kind, PatKind::Range(..)));
    found
}

/// `T::try_from(x)` / `x.try_into()`: the operand and the `T` it is checked
/// into.
fn checked_conversion<'tcx>(
    cx: &LateContext<'tcx>,
    e: &'tcx Expr<'tcx>,
) -> Option<(&'tcx Expr<'tcx>, Ty<'tcx>)> {
    // By name, not `clippy_utils::sym::try_from_fn`: that is one of clippy's
    // extra symbols, interned at an index the dylint driver need not share.
    let operand = match callee_of(cx, e)? {
        Callee::Path { def, args: [arg] }
            if cx
                .tcx
                .get_diagnostic_name(def)
                .is_some_and(|n| n.as_str() == "try_from_fn") =>
        {
            arg
        }
        Callee::Method {
            def,
            recv,
            args: [],
        } if cx.tcx.item_name(def) == sym::try_into
            && cx
                .tcx
                .opt_parent(def)
                .is_some_and(|t| cx.tcx.is_diagnostic_item(sym::TryInto, t)) =>
        {
            recv
        }
        _ => return None,
    };
    let ty::Adt(adt, args) = cx.typeck_results().expr_ty(e).kind() else {
        return None;
    };
    if !cx.tcx.is_diagnostic_item(sym::Result, adt.did()) {
        return None;
    }
    Some((operand, args.types().next()?))
}

/// A comparand fixed at compile time: a literal, a named constant
/// (`u32::MAX`, `LIMIT`), a `const fn` of constants (`size_of::<T>()`), or
/// arithmetic and casts over those. `ConstEvalCtxt` alone stops at `as`.
fn is_constant<'tcx>(cx: &LateContext<'tcx>, e: &'tcx Expr<'tcx>) -> bool {
    match e.kind {
        ExprKind::Lit(_) => true,
        ExprKind::Cast(inner, _) | ExprKind::DropTemps(inner) | ExprKind::Unary(_, inner) => {
            is_constant(cx, inner)
        }
        ExprKind::Binary(_, l, r) => is_constant(cx, l) && is_constant(cx, r),
        ExprKind::Path(ref qpath) => matches!(
            cx.qpath_res(qpath, e.hir_id),
            Res::Def(
                DefKind::Const { .. } | DefKind::AssocConst { .. } | DefKind::ConstParam,
                _
            )
        ),
        ExprKind::Call(..) => {
            matches!(callee_of(cx, e), Some(Callee::Path { def, args })
                if cx.tcx.is_const_fn(def) && args.iter().all(|a| is_constant(cx, a)))
        }
        _ => ConstEvalCtxt::new(cx).eval(e).is_some(),
    }
}

fn enclosing_fn(cx: &LateContext<'_>, hir_id: HirId) -> DefId {
    let owner = cx.tcx.hir_enclosing_body_owner(hir_id).to_def_id();
    cx.tcx.typeck_root_def_id(owner)
}

impl NarrowedTwoWays {
    fn record<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        site: &'tcx Expr<'tcx>,
        operand: &'tcx Expr<'tcx>,
        dst_ty: Ty<'tcx>,
        checked: bool,
    ) {
        let operand_ty = cx.typeck_results().expr_ty(operand).peel_refs();
        let (Some(src), Some(dst)) = (int_layout(cx, operand_ty), int_layout(cx, dst_ty)) else {
            return;
        };
        // A conversion that cannot fail is neither a truncation nor a check.
        if !lossy(src, dst) {
            return;
        }
        let Some((place, place_layout, place_ty, shown)) = place_of(cx, operand) else {
            return;
        };
        // `x as usize as u32` on a `u32` is `x` again: what the site does to
        // the place is a question about the place's type, not the operand's.
        if !lossy(place_layout, dst) {
            return;
        }
        self.sites.entry(place).or_default().push(Site {
            span: site.span,
            body: enclosing_fn(cx, site.hir_id),
            dst_bits: dst.bits,
            checked,
            shown,
            src: place_ty.to_string(),
            dst: dst_ty.to_string(),
        });
    }

    fn mark_range_checked<'tcx>(&mut self, cx: &LateContext<'tcx>, at: HirId, e: &'tcx Expr<'tcx>) {
        if let Some((place, ..)) = place_of(cx, e) {
            self.range_checked.insert((enclosing_fn(cx, at), place));
        }
    }
}

impl<'tcx> LateLintPass<'tcx> for NarrowedTwoWays {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        match expr.kind {
            ExprKind::Cast(operand, _) if !expr.span.from_expansion() => {
                // `x as u32 & (BITS - 1)`: the low bits are what is wanted.
                if let Some(parent) = clippy_utils::get_parent_expr(cx, expr)
                    && let ExprKind::Binary(op, ..) = parent.kind
                    && op.node == BinOpKind::BitAnd
                {
                    return;
                }
                let dst_ty = cx.typeck_results().expr_ty(expr);
                self.record(cx, expr, operand, dst_ty, false);
            }
            // `x < LIMIT`, `u32::MAX as usize >= x`: an explicit range test
            // of the place against something constant. `i < self.len` bounds
            // `i`, not `self.len`, which is why the other side must be one.
            ExprKind::Binary(op, l, r)
                if matches!(
                    op.node,
                    BinOpKind::Lt | BinOpKind::Le | BinOpKind::Gt | BinOpKind::Ge
                ) =>
            {
                if is_constant(cx, r) {
                    self.mark_range_checked(cx, expr.hir_id, l);
                }
                if is_constant(cx, l) {
                    self.mark_range_checked(cx, expr.hir_id, r);
                }
            }
            // `match x { 0..=9 => .., _ => .. }` and `if let 0..=9 = x` test
            // the range too.
            ExprKind::Match(scrut, arms, _) if arms.iter().any(|arm| has_range_pat(arm.pat)) => {
                self.mark_range_checked(cx, expr.hir_id, scrut);
            }
            ExprKind::Let(l) if has_range_pat(l.pat) => {
                self.mark_range_checked(cx, expr.hir_id, l.init);
            }
            _ if !expr.span.from_expansion() => {
                if let Some((operand, dst_ty)) = checked_conversion(cx, expr) {
                    self.record(cx, expr, operand, dst_ty, true);
                }
            }
            _ => {}
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let mut findings: Vec<(Span, Span, String, String)> = Vec::new();
        for (place, sites) in &self.sites {
            // The widest checked target: a check into `i64` says the value
            // may exceed even that, so every narrower `as` is condemned; a
            // check into `u8` says nothing about an `as u32`.
            let Some(check) = sites
                .iter()
                .filter(|s| s.checked)
                .max_by_key(|s| (s.dst_bits, std::cmp::Reverse(s.span.lo())))
            else {
                continue;
            };
            for site in sites {
                if site.checked
                    || site.dst_bits > check.dst_bits
                    || self.range_checked.contains(&(site.body, *place))
                {
                    continue;
                }
                findings.push((
                    site.span,
                    check.span,
                    format!(
                        "`{shown}` is a `{}` cut down to `{dst}` with `as` here, which wraps silently. Another place converts the same `{shown}` to `{}` with `try_from` because it might not fit",
                        site.src,
                        check.dst,
                        shown = site.shown,
                        dst = site.dst,
                    ),
                    format!(
                        "store `{}` as `{dst}` (or a newtype checked once at construction) so nobody has to narrow it. At least use `{dst}::try_from` here as well",
                        site.shown,
                        dst = site.dst,
                    ),
                ));
            }
        }
        // `sites` is a HashMap; report in source order.
        findings.sort_by_key(|(span, ..)| span.lo());
        for (span, check, msg, help) in findings {
            emit_with_note(
                cx,
                NARROWED_TWO_WAYS,
                span,
                msg,
                check,
                "the checked conversion",
                help,
            );
        }
    }
}

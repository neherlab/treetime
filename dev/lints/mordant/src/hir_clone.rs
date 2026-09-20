//! Structural identity of HIR: a hash to bucket candidates and an equality
//! to confirm them, both blind to spans, `HirId`s and binding names, both
//! reading paths by what they resolve to. Thin over
//! `clippy_utils::hir_utils`; what lives here is the pairing of local
//! bindings across two bodies or two expressions, which `SpanlessEq` leaves
//! to its caller (a `Res::Local` on the left equals one on the right when
//! the two are pre-mapped or identical), the comparison of parameter
//! patterns, which `SpanlessEq::eq_body` skips, and never
//! `deny_side_effects`, since that makes every method call unequal to
//! itself.
//!
//! Hash where typeck results exist (`check_fn`, `check_expr`); the equality
//! functions work from `check_crate_post`, where the context has no
//! enclosing body, by fetching each side's typeck results themselves.

use clippy_utils::{SpanlessEq, SpanlessHash};
use rustc_hir::def::Res;
use rustc_hir::def_id::LocalDefId;
use rustc_hir::intravisit::{Visitor, walk_expr, walk_pat};
use rustc_hir::{BodyId, Expr, ExprKind, HirId, HirIdMap, HirIdSet, Pat, PatKind, QPath};
use rustc_lint::LateContext;
use rustc_span::SyntaxContext;

/// Structural hash of a fn body: spans, `HirId`s and binding names do not
/// contribute; resolved paths, literals, field and method names, operators
/// and shape do. Parameter patterns do not contribute either, so two bodies
/// `bodies_equal` tells apart by their parameters can share a hash.
pub(crate) fn body_hash(cx: &LateContext<'_>, body: BodyId) -> u64 {
    let mut h = SpanlessHash::new(cx).paths_by_resolution();
    h.hash_body(body);
    h.finish()
}

/// Two parameter patterns destructure the same way: same shape at every
/// level, same constructor, same fields in the same order, same binding
/// modes. Each binding on the left is paired with the binding at its place
/// on the right in `locals`. `SpanlessEq` compares a body's value but not
/// its parameters, so without this `(a, _)` and `(_, b)` of one tuple type
/// would pair `a` with `b`. Or-patterns, ranges and literals are refutable
/// and cannot be a whole parameter; anything else unlisted is unequal.
fn eq_param_pat(
    cx: &LateContext<'_>,
    locals: &mut HirIdMap<HirId>,
    l: &Pat<'_>,
    r: &Pat<'_>,
) -> bool {
    let same_ctor = |lq: &QPath<'_>, rq: &QPath<'_>| {
        let res = cx.qpath_res(lq, l.hir_id);
        res != Res::Err && res == cx.qpath_res(rq, r.hir_id)
    };
    match (&l.kind, &r.kind) {
        (PatKind::Wild, PatKind::Wild) => true,
        (PatKind::Binding(lm, lid, _, lsub), PatKind::Binding(rm, rid, _, rsub)) => {
            let eq = lm == rm && eq_opt_pat(cx, locals, *lsub, *rsub);
            if eq {
                locals.insert(*lid, *rid);
            }
            eq
        }
        (PatKind::Tuple(ls, ld), PatKind::Tuple(rs, rd)) => {
            ld == rd && eq_param_pats(cx, locals, ls, rs)
        }
        (PatKind::TupleStruct(lq, ls, ld), PatKind::TupleStruct(rq, rs, rd)) => {
            ld == rd && same_ctor(lq, rq) && eq_param_pats(cx, locals, ls, rs)
        }
        (PatKind::Struct(lq, lfs, lrest), PatKind::Struct(rq, rfs, rrest)) => {
            lrest.is_some() == rrest.is_some()
                && same_ctor(lq, rq)
                && lfs.len() == rfs.len()
                && lfs.iter().zip(*rfs).all(|(lf, rf)| {
                    lf.ident.name == rf.ident.name && eq_param_pat(cx, locals, lf.pat, rf.pat)
                })
        }
        (PatKind::Ref(lp, lpin, lm), PatKind::Ref(rp, rpin, rm)) => {
            lpin == rpin && lm == rm && eq_param_pat(cx, locals, lp, rp)
        }
        (PatKind::Box(lp), PatKind::Box(rp)) | (PatKind::Deref(lp), PatKind::Deref(rp)) => {
            eq_param_pat(cx, locals, lp, rp)
        }
        (PatKind::Slice(lb, lm, la), PatKind::Slice(rb, rm, ra)) => {
            eq_param_pats(cx, locals, lb, rb)
                && eq_opt_pat(cx, locals, *lm, *rm)
                && eq_param_pats(cx, locals, la, ra)
        }
        _ => false,
    }
}

fn eq_param_pats(
    cx: &LateContext<'_>,
    locals: &mut HirIdMap<HirId>,
    ls: &[Pat<'_>],
    rs: &[Pat<'_>],
) -> bool {
    ls.len() == rs.len()
        && ls
            .iter()
            .zip(rs)
            .all(|(l, r)| eq_param_pat(cx, locals, l, r))
}

fn eq_opt_pat(
    cx: &LateContext<'_>,
    locals: &mut HirIdMap<HirId>,
    l: Option<&Pat<'_>>,
    r: Option<&Pat<'_>>,
) -> bool {
    match (l, r) {
        (None, None) => true,
        (Some(l), Some(r)) => eq_param_pat(cx, locals, l, r),
        _ => false,
    }
}

/// Two fn bodies are the same computation up to renaming of parameters and
/// locals. Callers must have checked the two signatures are equal first
/// (`fn_sigs_equal`): method calls compare by name, so `.len()` on two
/// different receiver types is the same call to this function. Two bodies
/// containing a closure are never equal (`SpanlessEq` refuses closures).
pub(crate) fn bodies_equal(cx: &LateContext<'_>, l: BodyId, r: BodyId) -> bool {
    let (lp, rp) = (cx.tcx.hir_body(l).params, cx.tcx.hir_body(r).params);
    let mut locals = HirIdMap::default();
    if lp.len() != rp.len()
        || !lp
            .iter()
            .zip(rp)
            .all(|(l, r)| eq_param_pat(cx, &mut locals, l.pat, r.pat))
    {
        return false;
    }
    let mut eq = SpanlessEq::new(cx).paths_by_resolution();
    let mut ie = eq.inter_expr(SyntaxContext::root());
    ie.locals = locals;
    ie.eq_body(l, r)
}

/// Erased signatures equal: same arity, each input and the output the same
/// `Ty` once late-bound regions are erased, and the same where-clauses. For
/// methods input 0 is `Self`, which is what keeps `Foo::is_empty` and
/// `Bar::is_empty` apart. The where-clauses matter because a method call in
/// the body compares by name: `t.go()` under `T: X` and under `T: Y` are two
/// different functions spelled alike.
pub(crate) fn fn_sigs_equal(cx: &LateContext<'_>, l: LocalDefId, r: LocalDefId) -> bool {
    let sig = |d: LocalDefId| {
        cx.tcx.instantiate_bound_regions_with_erased(
            cx.tcx
                .fn_sig(d.to_def_id())
                .instantiate_identity()
                .skip_normalization(),
        )
    };
    let bounds = |d: LocalDefId| {
        cx.tcx
            .predicates_of(d.to_def_id())
            .instantiate_identity(cx.tcx)
            .predicates
    };
    sig(l).inputs_and_output == sig(r).inputs_and_output && bounds(l) == bounds(r)
}

/// Structural hash of one expression, for bucketing across bodies. Every
/// local hashes alike, so which binding an arm reads never separates two
/// buckets; `exprs_equal` decides that.
pub(crate) fn expr_hash(cx: &LateContext<'_>, e: &Expr<'_>) -> u64 {
    let mut h = SpanlessHash::new(cx).paths_by_resolution();
    h.hash_expr(e);
    h.finish()
}

/// The locals `e` reads but does not bind, in order of first use.
fn free_locals(e: &Expr<'_>) -> Vec<HirId> {
    struct V {
        bound: HirIdSet,
        free: Vec<HirId>,
    }
    impl<'tcx> Visitor<'tcx> for V {
        fn visit_pat(&mut self, p: &'tcx Pat<'tcx>) {
            if let PatKind::Binding(_, id, ..) = p.kind {
                self.bound.insert(id);
            }
            walk_pat(self, p);
        }
        fn visit_expr(&mut self, e: &'tcx Expr<'tcx>) {
            if let ExprKind::Path(QPath::Resolved(None, path)) = e.kind
                && let Res::Local(id) = path.res
                && !self.bound.contains(&id)
                && !self.free.contains(&id)
            {
                self.free.push(id);
            }
            walk_expr(self, e);
        }
    }
    let mut v = V {
        bound: HirIdSet::default(),
        free: Vec::new(),
    };
    v.visit_expr(e);
    v.free
}

/// Two expressions (each given with its body owner, which may be the same)
/// are the same computation up to renaming: the locals each reads from
/// outside itself are paired in order of first use and must have the same
/// type in their own body, and under that pairing the two are structurally
/// equal. When both come from one body a local read by both sides can only
/// pair with itself: `SpanlessEq` accepts an identical local regardless of
/// the pairing, so any other pairing of it would not be a renaming.
pub(crate) fn exprs_equal(
    cx: &LateContext<'_>,
    (l_owner, l): (LocalDefId, &Expr<'_>),
    (r_owner, r): (LocalDefId, &Expr<'_>),
) -> bool {
    let (lf, rf) = (free_locals(l), free_locals(r));
    if lf.len() != rf.len()
        || lf
            .iter()
            .zip(&rf)
            .any(|(a, b)| a != b && (rf.contains(a) || lf.contains(b)))
    {
        return false;
    }
    let (lt, rt) = (cx.tcx.typeck(l_owner), cx.tcx.typeck(r_owner));
    let erased = |ty| cx.tcx.erase_and_anonymize_regions(ty);
    for (&a, &b) in lf.iter().zip(&rf) {
        match (lt.node_type_opt(a), rt.node_type_opt(b)) {
            (Some(ta), Some(tb)) if erased(ta) == erased(tb) => {}
            _ => return false,
        }
    }
    let mut eq = SpanlessEq::new(cx).paths_by_resolution();
    let mut ie = eq.inter_expr(SyntaxContext::root());
    ie.locals.extend(lf.into_iter().zip(rf));
    ie.eq_expr(l, r)
}

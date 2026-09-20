use std::collections::HashMap;

use rustc_hir::def_id::{DefId, LocalDefId};
use rustc_hir::{Expr, ExprKind, HirId, MatchSource, PatKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::{Span, sym};

use crate::adt_facts::inside_own_trait_impl;
use crate::baseline::emit_with_note;
use crate::hir_clone::{expr_hash, exprs_equal};

rustc_session::declare_lint! {
    /// Flags a `match` over an enum that repeats another one somewhere else
    /// in the crate arm for arm: same scrutinee type, same patterns, same arm
    /// bodies up to the names of the locals they read. The mapping exists
    /// twice, which is a method the enum does not have. A change to one copy
    /// will miss the other, and a variant added later is handled in
    /// whichever copy the author remembers.
    ///
    /// Only matches with at least two arms that name a pattern count (a
    /// `_`/binding catch-all plus one arm is a test, not a table), and only
    /// on enums outside the standard library, since a repeated `match` on
    /// `Option` is an idiom rather than a table. Matches inside the enum's
    /// own trait impls (`Display` beside `Debug`) stay quiet, as do matches
    /// produced by macros. Two copies whose free locals differ in type, or
    /// that read a different number of them, are different code and stay
    /// quiet; so does a copy nested inside a larger copy already reported.
    /// Two blind spots keep it quiet on some true copies: arms using `A | B`
    /// patterns are never confirmed equal, and a macro call in an arm
    /// (`format!(..)`) compares by its tokens, so a renamed local inside its
    /// arguments makes the arms differ.
    pub SAME_MATCH_TWICE,
    Warn,
    "the same match over one enum written out in two places"
}

struct Site {
    hash: u64,
    owner: LocalDefId,
    expr: HirId,
    span: Span,
}

#[derive(Default)]
pub struct SameMatchTwice {
    /// Scrutinee enum -> every counted match on it, in visit order.
    sites: HashMap<DefId, Vec<Site>>,
}

rustc_session::impl_lint_pass!(SameMatchTwice => [SAME_MATCH_TWICE]);

impl<'tcx> LateLintPass<'tcx> for SameMatchTwice {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let ExprKind::Match(scrut, arms, MatchSource::Normal) = expr.kind else {
            return;
        };
        if expr.span.from_expansion() {
            return;
        }
        let counted = arms
            .iter()
            .filter(|a| !matches!(a.pat.kind, PatKind::Wild | PatKind::Binding(.., None)))
            .count();
        if counted < 2 {
            return;
        }
        let Some(adt) = cx.typeck_results().expr_ty(scrut).peel_refs().ty_adt_def() else {
            return;
        };
        if !adt.is_enum()
            || matches!(
                cx.tcx.crate_name(adt.did().krate),
                sym::core | sym::alloc | sym::std
            )
            || inside_own_trait_impl(cx, expr.hir_id, adt.did())
        {
            return;
        }
        self.sites.entry(adt.did()).or_default().push(Site {
            hash: expr_hash(cx, expr),
            owner: cx.tcx.hir_enclosing_body_owner(expr.hir_id),
            expr: expr.hir_id,
            span: expr.span,
        });
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        // (later copy, earlier copy, enum)
        let mut findings: Vec<(Span, Span, DefId)> = Vec::new();
        for (adt, mut sites) in self.sites.drain() {
            sites.sort_by_key(|s| s.span.lo());
            let mut buckets: HashMap<u64, Vec<&Site>> = HashMap::new();
            for (i, site) in sites.iter().enumerate() {
                // A macro that expands its argument twice (`log!`) yields two
                // matches from one piece of source; that is one copy.
                if sites[..i].iter().any(|e| e.span.source_equal(site.span)) {
                    continue;
                }
                let bucket = buckets.entry(site.hash).or_default();
                let this = (site.owner, cx.tcx.hir_expect_expr(site.expr));
                let earlier = bucket
                    .iter()
                    .find(|e| exprs_equal(cx, (e.owner, cx.tcx.hir_expect_expr(e.expr)), this));
                match earlier {
                    Some(e) => findings.push((site.span, e.span, adt)),
                    // Only distinct shapes stand as candidates, so a third
                    // copy is reported against the first, once.
                    None => bucket.push(site),
                }
            }
        }
        findings.sort_by_key(|(span, ..)| span.lo());
        let mut reported: Vec<Span> = Vec::new();
        for (span, earlier, adt) in findings {
            // A match inside a repeated match is repeated too; the outer
            // report already covers it.
            if reported.iter().any(|r| r.contains(span)) {
                continue;
            }
            reported.push(span);
            let name = cx.tcx.def_path_str(adt);
            emit_with_note(
                cx,
                SAME_MATCH_TWICE,
                span,
                format!(
                    "this `match` on `{name}` repeats an earlier one arm for arm. The mapping exists twice, and a change to one copy will miss the other"
                ),
                earlier,
                "the other copy",
                format!("move the mapping into a method on `{name}` and call it from both places"),
            );
        }
    }
}

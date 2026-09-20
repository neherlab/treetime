use clippy_utils::res::MaybeResPath;
use clippy_utils::source::snippet_opt;
use rustc_ast::LitKind;
use rustc_errors::Applicability;
use rustc_hir::def_id::DefId;
use rustc_hir::{Arm, Expr, ExprKind, HirId, MatchSource, Pat, PatKind, StmtKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::Span;

use crate::MordantConfig;
use crate::baseline::emit_hir_then;
use crate::enum_facts::{pat_head_qpath, variant_of_res};
use crate::hir_shapes::{callee_of, peel_blocks_unsafe};

rustc_session::declare_lint! {
    /// Flags `_` (or a catch-all binding) matching over a small crate-local
    /// enum. The catch-all will also match any variant added to the enum
    /// later, so a new variant lands in it silently instead of failing to
    /// compile. Listing the remaining variants keeps exhaustiveness checking
    /// alive.
    ///
    /// Silent on an arm that answers "not this shape" (`None`, `false`, an
    /// empty slice, string or collection, or a `return` of one of those, with
    /// or without braces), and on a binding arm whose whole body is another
    /// `match` on that binding: the inner match is the dispatch, and is
    /// judged on its own.
    pub WILDCARD_OVER_OWN_ENUM,
    Warn,
    "wildcard arm over a small crate-local enum"
}

pub struct WildcardOverOwnEnum {
    pub config: &'static MordantConfig,
}

rustc_session::impl_lint_pass!(WildcardOverOwnEnum => [WILDCARD_OVER_OWN_ENUM]);

fn is_negative_extractor<'tcx>(cx: &LateContext<'tcx>, body: &Expr<'tcx>) -> bool {
    match body.kind {
        // `""` and `b""` are the empty-slice answer spelled as a literal.
        ExprKind::Lit(lit) => match lit.node {
            LitKind::Bool(b) => !b,
            LitKind::Str(s, _) => s.as_str().is_empty(),
            LitKind::ByteStr(s, _) => s.as_byte_str().is_empty(),
            _ => false,
        },
        ExprKind::Path(_) => clippy_utils::is_none_expr(cx, body),
        // `&[]` and `[]`: the empty-slice answer.
        ExprKind::AddrOf(_, _, inner) => is_negative_extractor(cx, inner),
        ExprKind::Array(elems) => elems.is_empty(),
        // `Vec::new()` / `String::new()`: an empty collection, constructed
        // fresh, is an answer with no content, not behavior.
        ExprKind::Call(_, []) => callee_of(cx, body).is_some_and(|c| {
            let path = cx.tcx.def_path_str(c.def());
            path.ends_with("Vec::<T>::new")
                || path.ends_with("Vec::new")
                || path.ends_with("String::new")
        }),
        // `return None` / `return false`: the early-exit spelling of the same
        // empty answers.
        ExprKind::Ret(Some(inner)) => is_negative_extractor(cx, inner),
        // `{ None }` and `{ return None; }`: the same answers inside the
        // braces rustfmt or a comment puts around them. Any other statement
        // is behavior the arm gives future variants.
        ExprKind::Block(block, None) => match (block.stmts, block.expr) {
            ([], Some(tail)) => is_negative_extractor(cx, tail),
            ([stmt], None) => match stmt.kind {
                StmtKind::Semi(e) | StmtKind::Expr(e) => {
                    matches!(e.kind, ExprKind::Ret(Some(_))) && is_negative_extractor(cx, e)
                }
                _ => false,
            },
            _ => false,
        },
        _ => false,
    }
}

/// A binding arm whose whole body is another `match` on the binding passes
/// the dispatch on; the inner match is where the variants are or are not
/// listed, and this lint judges it on its own.
fn redispatches(cx: &LateContext<'_>, body: &Expr<'_>, binding: HirId) -> bool {
    let ExprKind::Match(scrut, _, MatchSource::Normal) = peel_blocks_unsafe(body).kind else {
        return false;
    };
    clippy_utils::peel_ref_operators(cx, scrut).res_local_id() == Some(binding)
}

/// Every variant the non-catch-all arms cover, or None when an arm's shape is
/// beyond this analysis (then no fix is offered). `qspan` is a variant path
/// span to copy the file's path style from.
fn covered_variants(
    cx: &LateContext<'_>,
    arms: &[Arm<'_>],
    skip: &Pat<'_>,
) -> Option<(Vec<DefId>, Option<Span>)> {
    let mut covered = Vec::new();
    let mut qspan = None;
    for arm in arms {
        if arm.pat.hir_id == skip.hir_id {
            continue;
        }
        // A guarded arm covers nothing: its variant still reaches the
        // catch-all when the guard fails.
        if arm.guard.is_some() {
            continue;
        }
        collect(cx, arm.pat, &mut covered, &mut qspan)?;
    }
    return Some((covered, qspan));

    fn collect(
        cx: &LateContext<'_>,
        pat: &Pat<'_>,
        covered: &mut Vec<DefId>,
        qspan: &mut Option<Span>,
    ) -> Option<()> {
        if let Some(qpath) = pat_head_qpath(pat) {
            let variant = variant_of_res(cx, cx.qpath_res(qpath, pat.hir_id))?;
            covered.push(variant);
            qspan.get_or_insert(qpath.span());
            return Some(());
        }
        let PatKind::Or(pats) = pat.kind else {
            return None;
        };
        for p in pats {
            collect(cx, p, covered, qspan)?;
        }
        Some(())
    }
}

/// `Enum::Variant`, `Enum::Variant(..)`, or `Enum::Variant { .. }`, spelled
/// with the same path prefix the match's other arms use.
fn render_variant(_cx: &LateContext<'_>, prefix: &str, variant: &ty::VariantDef) -> String {
    let name = variant.name;
    match variant.ctor_kind() {
        Some(rustc_hir::def::CtorKind::Const) => format!("{prefix}{name}"),
        Some(rustc_hir::def::CtorKind::Fn) => format!("{prefix}{name}(..)"),
        None => format!("{prefix}{name} {{ .. }}"),
    }
}

impl<'tcx> LateLintPass<'tcx> for WildcardOverOwnEnum {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let ExprKind::Match(scrut, arms, MatchSource::Normal) = expr.kind else {
            return;
        };
        // A single-arm match is a destructuring, not a dispatch.
        if arms.len() < 2 {
            return;
        }
        let ty::Adt(adt, _) = cx
            .typeck_results()
            .expr_ty_adjusted(scrut)
            .peel_refs()
            .kind()
        else {
            return;
        };
        if !adt.is_enum() || !adt.did().is_local() {
            return;
        }
        // No `#[non_exhaustive]` exemption, deliberately. `is_local()` above
        // means the enum is defined in the crate being compiled, and a
        // `LateLintPass` only ever walks that same crate's bodies — so every
        // match this lint can reach is a match in the defining crate, where
        // `#[non_exhaustive]` has no effect and rustc checks exhaustiveness
        // normally. It constrains downstream crates only. An earlier exemption
        // here read "already opted out of exhaustiveness", which was false for
        // every match the lint sees, and it silenced the lint on exactly the
        // enums annotated *because* the author expects new variants.
        let n = adt.variants().len();
        if n > self.config.wildcard_over_own_enum_max_variants {
            return;
        }
        for arm in arms {
            if arm.guard.is_some() {
                continue;
            }
            let binding = match arm.pat.kind {
                PatKind::Wild => None,
                PatKind::Binding(_, id, _, None) => Some(id),
                _ => continue,
            };
            // `_ => None` and `_ => false` are extractors asking "is it this
            // one shape?" — future variants are correctly not that shape, so
            // absorption is the right behavior, not a hazard.
            if is_negative_extractor(cx, arm.body) {
                continue;
            }
            if binding.is_some_and(|b| redispatches(cx, arm.body, b)) {
                continue;
            }
            // The fix replaces the catch-all with the uncovered variants,
            // spelled with the same path prefix the sibling arms use. Offered
            // only when every sibling arm's coverage is provable.
            let fix = covered_variants(cx, arms, arm.pat).and_then(|(covered, qspan)| {
                let missing: Vec<&ty::VariantDef> = adt
                    .variants()
                    .iter()
                    .filter(|v| !covered.contains(&v.def_id))
                    .collect();
                if missing.is_empty() {
                    return None;
                }
                let prefix =
                    qspan
                        .and_then(|s| snippet_opt(cx, s))
                        .map(|snip| match snip.rfind("::") {
                            Some(i) => snip[..i + 2].to_string(),
                            None => String::new(),
                        })?;
                let list = missing
                    .iter()
                    .map(|v| render_variant(cx, &prefix, v))
                    .collect::<Vec<_>>()
                    .join(" | ");
                Some(match arm.pat.kind {
                    PatKind::Binding(_, _, ident, None) => format!("{ident} @ ({list})"),
                    _ => list,
                })
            });
            // Emitted against the ARM's hir id, so an #[allow] placed on the
            // arm itself is honored; a plain span emission would only see the
            // enclosing match's lint level.
            let name = cx.tcx.def_path_str(adt.did());
            emit_hir_then(
                cx,
                WILDCARD_OVER_OWN_ENUM,
                arm.hir_id,
                arm.pat.span,
                format!(
                    "this catch-all will also match any variant added to `{name}` later ({n} today). A new variant lands here silently instead of failing to compile"
                ),
                |diag| {
                    diag.span_note(
                        cx.tcx.def_span(adt.did()),
                        format!("`{name}` is defined here"),
                    );
                    let wildcard = snippet_opt(cx, arm.pat.span).unwrap_or_else(|| "_".into());
                    let msg = format!("list the remaining variants instead of `{wildcard}`");
                    match &fix {
                        Some(sugg) => {
                            diag.span_suggestion_verbose(
                                arm.pat.span,
                                msg,
                                sugg,
                                Applicability::MachineApplicable,
                            );
                        }
                        None => {
                            diag.help(msg);
                        }
                    }
                },
            );
        }
    }
}

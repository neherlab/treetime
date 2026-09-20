use clippy_utils::visitors::for_each_expr;
use rustc_hir::{Block, Expr, ExprKind, Pat, QPath, StmtKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::sym;

use crate::baseline::emit_with_note;
use crate::hir_shapes::{FieldChain, dotted, field_chain, stmt_expr};

rustc_session::declare_lint! {
    /// Flags `map.get(&k).unwrap()` right after a `map.insert(k, ..)` put
    /// the key there, with nothing in between that could touch the map or the
    /// key: no calls, no assignments to either. The lookup and its panic are
    /// redundant. Keeping the inserted value (or the entry API) removes them.
    ///
    /// Silent once a `let` in between rebinds the name the map or key was
    /// spelled with, and treats `-k` / `!k` as a different key from `k`.
    pub INSERT_THEN_UNWRAP,
    Warn,
    "unwrap of a lookup proven by an insert just above"
}

rustc_session::declare_lint_pass!(InsertThenUnwrap => [INSERT_THEN_UNWRAP]);

/// A stable textual identity for the small expressions worth tracking:
/// `self.a.b` chains, plain locals, and literals, seen through `&` and `*`
/// (which name the same value); `-k` and `!k` are different values from `k`.
/// Anything else is `None`, and untrackable means unprovable means silent.
fn identity(e: &Expr<'_>) -> Option<String> {
    let FieldChain { root, fields } = field_chain(e);
    let head = match &root.kind {
        ExprKind::Path(QPath::Resolved(None, p)) if p.segments.len() == 1 => {
            p.segments[0].ident.name.to_string()
        }
        ExprKind::Lit(lit) => format!("lit:{:?}", lit.node),
        _ => return None,
    };
    Some(dotted(head, &fields))
}

/// The local an identity is rooted at: `k` for `k` and `k.id`, `None` for
/// `self` chains and literals, which no `let` can rebind.
fn root_local(id: &str) -> Option<&str> {
    let root = id.split('.').next()?;
    (root != "self" && !root.starts_with("lit:")).then_some(root)
}

/// A `let` between the insert and the lookup that rebinds the map's or the
/// key's root name makes the later spelling name a different value.
fn rebinds(pat: &Pat<'_>, roots: &[&str]) -> bool {
    let mut hit = false;
    pat.each_binding(|_, _, _, ident| hit |= roots.contains(&ident.as_str()));
    hit
}

fn is_map(cx: &LateContext<'_>, recv: &Expr<'_>) -> bool {
    let ty::Adt(adt, _) = cx
        .typeck_results()
        .expr_ty_adjusted(recv)
        .peel_refs()
        .kind()
    else {
        return false;
    };
    cx.tcx.is_diagnostic_item(sym::HashMap, adt.did())
        || cx.tcx.is_diagnostic_item(sym::BTreeMap, adt.did())
}

/// `map.insert(k, _)` as (map identity, key identity).
fn insert_of(cx: &LateContext<'_>, e: &Expr<'_>) -> Option<(String, String)> {
    let ExprKind::MethodCall(seg, recv, [k, _], _) = &e.kind else {
        return None;
    };
    if seg.ident.as_str() != "insert" || !is_map(cx, recv) {
        return None;
    }
    Some((identity(recv)?, identity(k)?))
}

/// `map.get(&k).unwrap()` (or `.expect(..)`) as (map, key, whole span).
fn get_unwrap_of<'tcx>(
    cx: &LateContext<'tcx>,
    e: &'tcx Expr<'tcx>,
) -> Option<(String, String, String, rustc_span::Span)> {
    let ExprKind::MethodCall(useg, urecv, _, _) = &e.kind else {
        return None;
    };
    if !matches!(useg.ident.as_str(), "unwrap" | "expect") {
        return None;
    }
    let ExprKind::MethodCall(gseg, grecv, [k], _) = &urecv.kind else {
        return None;
    };
    if !matches!(gseg.ident.as_str(), "get" | "get_mut") || !is_map(cx, grecv) {
        return None;
    }
    let shown = clippy_utils::source::snippet_opt(cx, k.span)
        .unwrap_or_else(|| "..".to_string())
        .trim_start_matches('&')
        .to_string();
    Some((identity(grecv)?, identity(k)?, shown, e.span))
}

/// True when this statement could invalidate a tracked (map, key) fact:
/// any call at all, or any assignment. Purity is not analyzed; anything
/// callable is assumed able to touch the map, which can only cause silence.
fn may_disturb<'tcx>(cx: &LateContext<'tcx>, e: &'tcx Expr<'tcx>) -> bool {
    let mut disturbs = false;
    for_each_expr(cx, e, |inner: &Expr<'tcx>| {
        match &inner.kind {
            ExprKind::MethodCall(..) | ExprKind::Call(..) => {
                // The proven lookup itself is a call; the caller filters it
                // out by checking lookups before disturbance.
                disturbs = true;
            }
            ExprKind::Assign(..) | ExprKind::AssignOp(..) => disturbs = true,
            _ => {}
        }
        std::ops::ControlFlow::<()>::Continue(())
    });
    disturbs
}

impl<'tcx> LateLintPass<'tcx> for InsertThenUnwrap {
    fn check_block(&mut self, cx: &LateContext<'tcx>, block: &'tcx Block<'tcx>) {
        for (i, stmt) in block.stmts.iter().enumerate() {
            let (StmtKind::Expr(e) | StmtKind::Semi(e)) = stmt.kind else {
                continue;
            };
            let Some((map, key)) = insert_of(cx, e) else {
                continue;
            };
            let roots: Vec<&str> = [&map, &key]
                .into_iter()
                .filter_map(|id| root_local(id))
                .collect();
            for later in &block.stmts[i + 1..] {
                let le = stmt_expr(later);
                // The lookup must be the statement's own top-level shape (or
                // the let initializer); a lookup buried under other calls is
                // checked as a disturbance instead.
                if let Some(le) = le
                    && let Some((gmap, gkey, shown, at)) = get_unwrap_of(cx, le)
                    && gmap == map
                    && gkey == key
                {
                    let map_shown = map.strip_prefix("self.").unwrap_or(&map);
                    emit_with_note(
                        cx,
                        INSERT_THEN_UNWRAP,
                        at,
                        format!(
                            "this looks up `{shown}` in `{map_shown}` and unwraps it right after the `insert` above put it there. Nothing in between could remove it, so the lookup and its panic are redundant"
                        ),
                        e.span,
                        "the insert",
                        format!(
                            "keep the value you inserted, or use `{map_shown}.entry({shown})` which hands it back"
                        ),
                    );
                    return;
                }
                // The initializer still saw the names as the insert did; the
                // pattern takes effect for the statements after it.
                if let StmtKind::Let(l) = later.kind
                    && rebinds(l.pat, &roots)
                {
                    break;
                }
                if le.is_some_and(|le| may_disturb(cx, le)) {
                    break;
                }
            }
        }
    }
}

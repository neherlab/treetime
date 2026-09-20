use std::collections::{HashMap, HashSet};

use clippy_utils::visitors::{Descend, for_each_expr};
use rustc_hir::def_id::DefId;
use rustc_hir::intravisit::FnKind;
use rustc_hir::{Block, Body, Expr, ExprKind, FnDecl, HirId, StmtKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::def_id::LocalDefId;
use rustc_span::{Span, Symbol};

use crate::adt_facts::impl_self_adt;
use crate::baseline::emit_hir_then;
use crate::hir_shapes::{
    Callee, SelfField, callee_of, ends_in_return, is_self_path, peel_not, self_field, stmt_expr,
};

rustc_session::declare_lint! {
    /// Flags a call that only runs when `self.can_x()` returns true, but
    /// changes a field of `self` that `can_x` never reads, directly or
    /// through anything it calls. The check cannot tell whether the call is
    /// safe. The guard/mutator pair drifting apart is how a live connection
    /// gets evicted by a donation its guard approved.
    ///
    /// Only calls at the guard's own level count as its actions. A call
    /// nested under a further `if`, `match` or loop is gated by that
    /// condition; typically it is a state transition of its own that happens
    /// to sit in a method the guard opened. The exception is a further guard:
    /// a call approved by `can_x()` and then `can_y()` is judged against what
    /// the two read together, and is reported once, naming both.
    ///
    /// A guard is followed through the fields it reads and the same-type
    /// methods it calls. One that hands `self` to anything else (a free
    /// function, a trait method) reads who knows what, and gates nothing as
    /// far as this lint is concerned: silence, not a guess.
    pub GUARD_BLIND_TO_ACTION,
    Warn,
    "guarded call touches state its guard never reads"
}

#[derive(Default)]
struct MethodFacts {
    /// `self.field` occurrences, read or written.
    touched: HashSet<Symbol>,
    /// Same-type inherent methods called on `self`.
    calls: HashSet<DefId>,
    /// `self` handed whole to something the two sets above do not follow (a
    /// free function, a trait method, a closure): what gets read behind it
    /// is unknown, so this method's coverage is not computable.
    escapes: bool,
}

/// One guarded call site and every guard that approved it.
struct Gate {
    guards: Vec<DefId>,
    action: DefId,
    span: Span,
}

#[derive(Default)]
pub struct GuardBlindToAction {
    facts: HashMap<DefId, MethodFacts>,
    /// In the order the sites were first seen, so findings come out in
    /// source order; `site_index` folds a second guard into the same entry.
    gates: Vec<Gate>,
    site_index: HashMap<Span, usize>,
}

rustc_session::impl_lint_pass!(GuardBlindToAction => [GUARD_BLIND_TO_ACTION]);

/// The inherent method a `self.m(..)` call resolves to, with its self type,
/// when that type is a crate-local ADT.
fn self_method_call<'tcx>(cx: &LateContext<'tcx>, e: &Expr<'tcx>) -> Option<(DefId, DefId)> {
    let Some(Callee::Method {
        def: method, recv, ..
    }) = callee_of(cx, e)
    else {
        return None;
    };
    if !is_self_path(recv) {
        return None;
    }
    let impl_did = cx.tcx.impl_of_assoc(method)?;
    if !matches!(
        cx.tcx.def_kind(impl_did),
        rustc_hir::def::DefKind::Impl { of_trait: false }
    ) {
        return None;
    }
    let adt = impl_self_adt(cx, impl_did)?;
    adt.did().is_local().then(|| (method, adt.did()))
}

/// A guard is a permission-flavored boolean method: the names that promise
/// "this action is safe to take".
fn is_guard_name(cx: &LateContext<'_>, method: DefId) -> bool {
    let name = cx.tcx.item_name(method);
    let n = name.as_str();
    n.starts_with("can_") || n.starts_with("may_") || n.starts_with("check_")
}

/// `if self.can_y() { .. }` with no else: the then-branch is still the outer
/// guard's territory, only more narrowly approved. Any other condition, or a
/// denial branch, is a decision of its own.
fn positive_guard_if<'tcx>(
    cx: &LateContext<'tcx>,
    e: &'tcx Expr<'tcx>,
) -> Option<(DefId, DefId, &'tcx Expr<'tcx>)> {
    let ExprKind::If(cond, then, None) = e.kind else {
        return None;
    };
    let (gexpr, negated) = peel_not(cond);
    let (guard, adt) = self_method_call(cx, gexpr)?;
    (!negated && is_guard_name(cx, guard)).then_some((guard, adt, then))
}

impl GuardBlindToAction {
    fn note_gate(&mut self, guards: &[DefId], action: DefId, span: Span) {
        match self.site_index.get(&span) {
            Some(&i) => {
                let known = &mut self.gates[i].guards;
                for &g in guards {
                    if !known.contains(&g) {
                        known.push(g);
                    }
                }
            }
            None => {
                self.site_index.insert(span, self.gates.len());
                self.gates.push(Gate {
                    guards: guards.to_vec(),
                    action,
                    span,
                });
            }
        }
    }

    /// Every same-type call in `scope` is approved by `guard`; below a further
    /// positive guard it is approved by both.
    fn collect_actions<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        scope: &'tcx Expr<'tcx>,
        guard: DefId,
        adt: DefId,
    ) {
        let mut sites: Vec<(DefId, Span, Vec<DefId>)> = Vec::new();
        let mut pending = vec![(scope, vec![guard])];
        while let Some((scope, guards)) = pending.pop() {
            for_each_expr(cx, scope, |e: &'tcx Expr<'tcx>| {
                if let Some((inner, _, then)) = positive_guard_if(cx, e) {
                    let mut below = guards.clone();
                    if !below.contains(&inner) {
                        below.push(inner);
                    }
                    pending.push((then, below));
                    return std::ops::ControlFlow::<(), Descend>::Continue(Descend::No);
                }
                // Anything below its own condition answers to that condition. (The
                // then-branch form hands in the guard's own block, which is not one.)
                if matches!(
                    e.kind,
                    ExprKind::If(..)
                        | ExprKind::Match(..)
                        | ExprKind::Loop(..)
                        | ExprKind::Closure(..)
                ) {
                    return std::ops::ControlFlow::<(), Descend>::Continue(Descend::No);
                }
                if let Some((m, m_adt)) = self_method_call(cx, e)
                    && m_adt == adt
                    && !guards.contains(&m)
                {
                    sites.push((m, e.span, guards.clone()));
                }
                std::ops::ControlFlow::<(), Descend>::Continue(Descend::Yes)
            });
        }
        for (action, span, guards) in sites {
            self.note_gate(&guards, action, span);
        }
    }
}

impl GuardBlindToAction {
    /// The fields `guard` reads, directly or through the same-type methods
    /// it calls, sorted by name; None when it hands `self` to something the
    /// walk cannot follow, so its coverage is not computable.
    fn coverage(&self, guard: DefId) -> Option<Vec<Symbol>> {
        let mut read: HashSet<Symbol> = HashSet::new();
        let mut queue = vec![guard];
        let mut seen: HashSet<DefId> = HashSet::new();
        while let Some(m) = queue.pop() {
            if !seen.insert(m) {
                continue;
            }
            if let Some(f) = self.facts.get(&m) {
                if f.escapes {
                    return None;
                }
                read.extend(f.touched.iter().copied());
                queue.extend(f.calls.iter().copied());
            }
        }
        let mut read: Vec<Symbol> = read.into_iter().collect();
        read.sort_by_key(|s| s.as_str().to_owned());
        Some(read)
    }
}

impl<'tcx> LateLintPass<'tcx> for GuardBlindToAction {
    fn check_fn(
        &mut self,
        cx: &LateContext<'tcx>,
        kind: FnKind<'tcx>,
        _decl: &'tcx FnDecl<'tcx>,
        body: &'tcx Body<'tcx>,
        _span: Span,
        def_id: LocalDefId,
    ) {
        if matches!(kind, FnKind::Closure) {
            return;
        }
        // Facts: which self-fields this method touches, and what it calls.
        let method = def_id.to_def_id();
        let mut facts = MethodFacts::default();
        // The walk is pre-order, so a `self` that is the base of a field or
        // the receiver of a followed call is marked before it is visited; any
        // other `self` reaching the walk is one that got away.
        let mut followed: HashSet<HirId> = HashSet::new();
        for_each_expr(cx, body.value, |e: &Expr<'_>| {
            if let Some(SelfField { base, ident }) = self_field(e) {
                facts.touched.insert(ident.name);
                followed.insert(base.hir_id);
            } else if let Some((m, _)) = self_method_call(cx, e)
                && let ExprKind::MethodCall(_, recv, ..) = e.kind
            {
                facts.calls.insert(m);
                followed.insert(recv.hir_id);
            } else if is_self_path(e) && !followed.contains(&e.hir_id) {
                facts.escapes = true;
            }
            std::ops::ControlFlow::<()>::Continue(())
        });
        self.facts.insert(method, facts);
    }

    fn check_block(&mut self, cx: &LateContext<'tcx>, block: &'tcx Block<'tcx>) {
        // `if self.can_x() { self.y() }` as the block's tail gates the same
        // way it does as a statement; with nothing after it, only that form.
        if let Some(tail) = block.expr
            && let Some((guard, adt, then)) = positive_guard_if(cx, tail)
        {
            self.collect_actions(cx, then, guard, adt);
        }
        // `if !self.can_x() { return } ... self.y()` — the guard covers the
        // rest of the block.
        for (i, stmt) in block.stmts.iter().enumerate() {
            let (StmtKind::Expr(e) | StmtKind::Semi(e)) = stmt.kind else {
                continue;
            };
            let ExprKind::If(cond, then, None) = e.kind else {
                continue;
            };
            let (gexpr, negated) = peel_not(cond);
            let Some((guard, adt)) = self_method_call(cx, gexpr) else {
                continue;
            };
            if !is_guard_name(cx, guard) {
                continue;
            }
            if negated && ends_in_return(then) {
                // `let` initializers count: the real-world shape binds the
                // result, `let connections = self.detach_fds(&victim);`.
                let later = block.stmts[i + 1..].iter().filter_map(stmt_expr);
                for e in later.chain(block.expr) {
                    self.collect_actions(cx, e, guard, adt);
                }
            } else if !negated {
                // `if self.can_x() { self.y() }` — the then-branch is gated.
                self.collect_actions(cx, then, guard, adt);
            }
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        for gate in &self.gates {
            let Some(action) = self.facts.get(&gate.action) else {
                continue;
            };
            // Only mutations: a guarded read cannot invalidate anything. The
            // honest mutator signal is the receiver: `&mut self`.
            let sig = cx
                .tcx
                .fn_sig(gate.action)
                .instantiate_identity()
                .skip_binder();
            let mutates = sig
                .inputs()
                .first()
                .is_some_and(|t| matches!(t.kind(), ty::Ref(_, _, m) if m.is_mut()));
            if !mutates {
                continue;
            }
            // A guard's coverage is everything it touches transitively
            // through same-type calls, so a guard delegating to helpers is
            // never accused falsely; one delegating to something the walk
            // cannot follow has no computable coverage, and is not judged.
            let mut guards: Vec<(String, DefId, Vec<Symbol>)> = Vec::new();
            for &g in &gate.guards {
                let Some(reads) = self.coverage(g) else {
                    break;
                };
                guards.push((format!("`{}`", cx.tcx.item_name(g)), g, reads));
            }
            if guards.len() != gate.guards.len() {
                continue;
            }
            let covered: HashSet<Symbol> = guards
                .iter()
                .flat_map(|(.., r)| r.iter().copied())
                .collect();
            let mut missed: Vec<&Symbol> = action
                .touched
                .iter()
                .filter(|f| !covered.contains(f))
                .collect();
            if missed.is_empty() {
                continue;
            }
            missed.sort_by_key(|s| s.as_str().to_owned());
            let fields: Vec<String> = missed.iter().map(|s| format!("`{s}`")).collect();
            let fields = fields.join(", ");
            guards.sort_by(|a, b| a.0.cmp(&b.0));
            let names: Vec<&str> = guards.iter().map(|(n, ..)| n.as_str()).collect();
            let (named, returns, blind, checks, either) = match names.as_slice() {
                [one] => (
                    (*one).to_owned(),
                    "returns",
                    format!("{one} never reads"),
                    "The check",
                    (*one).to_owned(),
                ),
                [a, b] => (
                    format!("{a} and {b}"),
                    "return",
                    format!("neither {a} nor {b} reads"),
                    "The checks",
                    format!("{a} or {b}"),
                ),
                [head @ .., last] => (
                    format!("{} and {last}", head.join(", ")),
                    "return",
                    format!("none of {} and {last} reads", head.join(", ")),
                    "The checks",
                    format!("{} or {last}", head.join(", ")),
                ),
                [] => continue,
            };
            let action = cx.tcx.item_name(gate.action);
            emit_hir_then(
                cx,
                GUARD_BLIND_TO_ACTION,
                cx.last_node_with_lint_attrs,
                gate.span,
                format!(
                    "`{action}` only runs when {named} {returns} true, but it changes {fields}, which {blind}. {checks} cannot tell whether `{action}` is safe"
                ),
                |diag| {
                    // What each guard does read: the other half of the
                    // mismatch, shown at the guard.
                    for (name, g, reads) in &guards {
                        let read: Vec<String> = reads.iter().map(|s| format!("`{s}`")).collect();
                        let what = if read.is_empty() {
                            "no field of `self`".to_owned()
                        } else {
                            format!("only {}", read.join(", "))
                        };
                        diag.span_note(cx.tcx.def_span(*g), format!("{name} reads {what}"));
                    }
                    diag.help(format!(
                        "make {either} look at {fields} too, or move the {fields} work out from under this check"
                    ));
                },
            );
        }
    }
}

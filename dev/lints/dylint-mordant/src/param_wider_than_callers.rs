use std::collections::{HashMap, HashSet};

use clippy_utils::res::MaybeResPath;
use clippy_utils::visitors::for_each_expr;
use rustc_hir::def::{DefKind, Res};
use rustc_hir::def_id::DefId;
use rustc_hir::intravisit::FnKind;
use rustc_hir::{Body, Expr, ExprKind, FnDecl, HirId, PatKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::Span;
use rustc_span::def_id::LocalDefId;

use crate::baseline::emit_with_note;
use crate::hir_shapes::{Callee, callee_of};

rustc_session::declare_lint! {
    /// Flags a match arm that panics on an enum variant no existing call
    /// passes: the function takes the full enum, panics on `E::C`, and every
    /// call in the crate provably passes a different variant (constructor
    /// literals only, anything else makes the set unknowable and the lint
    /// silent. A method's receiver is the value sent to `self`, so `m.run()`
    /// makes `run`'s domain unknowable while `Mode::Off.run()` passes `Off`).
    /// The parameter type allows a case no caller uses. Narrowing the type
    /// turns it into a compile error for future callers.
    pub PARAM_WIDER_THAN_CALLERS,
    Warn,
    "panicking arm for a variant no existing call site passes"
}

#[derive(Default)]
struct CallFacts {
    passed: HashSet<DefId>,
    sites: usize,
    /// The first argument recorded, for the note.
    first: Option<Span>,
    unknown: bool,
}

#[derive(Default)]
pub struct ParamWiderThanCallers {
    /// (fn, param idx) -> variants with a diverging panic arm in the fn body.
    panicking: HashMap<(DefId, usize), Vec<(DefId, Span)>>,
    calls: HashMap<(DefId, usize), CallFacts>,
    /// Functions referenced other than by direct call: indirect callers are
    /// invisible, so their call-site sets are unknowable.
    poisoned: HashSet<DefId>,
}

rustc_session::impl_lint_pass!(ParamWiderThanCallers => [PARAM_WIDER_THAN_CALLERS]);

impl ParamWiderThanCallers {
    /// One call of `def` sends `value` to its body param `i`.
    fn record_arg<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        def: DefId,
        i: usize,
        value: &Expr<'tcx>,
    ) {
        let facts = self.calls.entry((def, i)).or_default();
        facts.sites += 1;
        facts.first.get_or_insert(value.span);
        match crate::enum_facts::ctor_literal_variant(cx, value) {
            Some(v) => {
                facts.passed.insert(v);
            }
            None => facts.unknown = true,
        }
    }
}

impl<'tcx> LateLintPass<'tcx> for ParamWiderThanCallers {
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
        // External callers are invisible; only a crate-private fn has a
        // knowable call-site set.
        if cx.effective_visibilities.is_exported(def_id) {
            return;
        }
        // Param locals whose type is a crate-local enum.
        let params: Vec<Option<HirId>> = body
            .params
            .iter()
            .map(|param| {
                let PatKind::Binding(_, hir_id, _, None) = param.pat.kind else {
                    return None;
                };
                let pty = cx.typeck_results().pat_ty(param.pat).peel_refs();
                matches!(pty.kind(), ty::Adt(adt, _) if adt.is_enum() && adt.did().is_local())
                    .then_some(hir_id)
            })
            .collect();
        if params.iter().all(Option::is_none) {
            return;
        }
        // Panicking arms in matches whose scrutinee is exactly the param.
        let fn_def = def_id.to_def_id();
        for_each_expr(cx, body.value, |e: &Expr<'tcx>| {
            if let ExprKind::Match(scrut, arms, _) = e.kind
                && let Some(scrut_local) = scrut.res_local_id()
                && let Some(idx) = params.iter().position(|p| *p == Some(scrut_local))
            {
                for arm in arms {
                    if arm.guard.is_none()
                        && let Some(variant) = crate::enum_facts::arm_variant(cx, arm.pat)
                        && crate::enum_facts::is_panic_arm(cx, arm.body)
                    {
                        self.panicking
                            .entry((fn_def, idx))
                            .or_default()
                            .push((variant, arm.span));
                    }
                }
            }
            std::ops::ControlFlow::<()>::Continue(())
        });
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        match callee_of(cx, expr) {
            Some(Callee::Path { def, args }) => {
                if !def.is_local()
                    || !matches!(cx.tcx.def_kind(def), DefKind::Fn | DefKind::AssocFn)
                {
                    return;
                }
                for (i, a) in args.iter().enumerate() {
                    self.record_arg(cx, def, i, a);
                }
            }
            Some(Callee::Method { def, recv, args }) => {
                if !def.is_local() {
                    return;
                }
                // The receiver is the body's param 0 (`self`), so it is a
                // call-site value for that param like any other: `m.run()`
                // sends whatever `m` holds, `Mode::Off.run()` sends `Off`.
                self.record_arg(cx, def, 0, recv);
                for (i, a) in args.iter().enumerate() {
                    self.record_arg(cx, def, i + 1, a);
                }
            }
            // A bare reference to a local fn (fn pointer, higher-order use):
            // its future call sites are invisible. The callee position of a
            // direct call is not a bare reference; that call is counted above.
            None => {
                let ExprKind::Path(qpath) = &expr.kind else {
                    return;
                };
                if matches!(
                    clippy_utils::get_parent_expr(cx, expr),
                    Some(Expr { kind: ExprKind::Call(callee, _), .. }) if callee.hir_id == expr.hir_id
                ) {
                    return;
                }
                if let Res::Def(DefKind::Fn | DefKind::AssocFn, def) =
                    cx.qpath_res(qpath, expr.hir_id)
                    && def.is_local()
                {
                    self.poisoned.insert(def);
                }
            }
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let mut findings: Vec<(Span, String, Span, String)> = Vec::new();
        for (key @ (fn_def, _), arms) in &self.panicking {
            if self.poisoned.contains(fn_def) {
                continue;
            }
            let Some(calls) = self.calls.get(key) else {
                // Never called: nothing provable about its domain.
                continue;
            };
            let Some(first) = calls.first else {
                continue;
            };
            if calls.unknown {
                continue;
            }
            for (variant, span) in arms {
                if calls.passed.contains(variant) {
                    continue;
                }
                let mut passed: Vec<String> = calls
                    .passed
                    .iter()
                    .map(|v| format!("`{}`", cx.tcx.item_name(*v)))
                    .collect();
                passed.sort();
                let fn_name = cx.tcx.item_name(*fn_def);
                let variant = cx.tcx.item_name(*variant);
                findings.push((
                    *span,
                    format!(
                        "this arm panics on `{variant}`, but all {} calls to `{fn_name}` pass only {}. The parameter type allows a case no caller uses",
                        calls.sites,
                        passed.join(", "),
                    ),
                    first,
                    format!(
                        "narrow `{fn_name}`'s parameter to a type that cannot be `{variant}`, such as a smaller enum or the payload itself. A caller passing `{variant}` then fails to compile"
                    ),
                ));
            }
        }
        // `panicking` is a HashMap; report in source order.
        findings.sort_by_key(|(span, ..)| span.lo());
        for (span, msg, first, help) in findings {
            emit_with_note(
                cx,
                PARAM_WIDER_THAN_CALLERS,
                span,
                msg,
                first,
                "one of the calls",
                help,
            );
        }
    }
}

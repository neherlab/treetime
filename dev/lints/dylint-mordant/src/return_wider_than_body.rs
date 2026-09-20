use std::collections::{HashMap, HashSet};

use rustc_hir::def::DefKind;
use rustc_hir::def_id::DefId;
use rustc_hir::intravisit::FnKind;
use rustc_hir::{Body, Expr, ExprKind, FnDecl};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::Span;
use rustc_span::def_id::LocalDefId;

use crate::baseline::emit_with_note;
use crate::enum_facts::{arm_variant, is_panic_arm};
use crate::hir_shapes::{Callee, callee_of};
use crate::variant_flow::returned_variants;

rustc_session::declare_lint! {
    /// Flags a match arm that panics on a variant the matched call never
    /// returns: the callee's return type is a crate-local enum, every value
    /// its body returns is a constructor literal, and the panicked-on variant
    /// is not among them. The return type is wider than what the function
    /// produces. Narrowing it deletes the arm at compile time.
    ///
    /// The return set comes from MIR dataflow: constructor aggregates traced
    /// through plain copies between locals, so `let t = ...; t` and branches
    /// resolve. Anything untraceable — a parameter, a call result, a
    /// projection — makes the set unknowable and the function is skipped.
    pub RETURN_WIDER_THAN_BODY,
    Warn,
    "panicking arm for a variant the callee never constructs"
}

#[derive(Default)]
pub struct ReturnWiderThanBody {
    /// fn -> its enum and the provably-complete set of returned variants.
    returns: HashMap<DefId, (DefId, HashSet<rustc_abi::VariantIdx>)>,
    /// (callee, panicked variant, arm span, caller-visible name) discovered
    /// at match sites, resolved against `returns` at the end.
    suspects: Vec<(DefId, DefId, Span)>,
}

rustc_session::impl_lint_pass!(ReturnWiderThanBody => [RETURN_WIDER_THAN_BODY]);

impl<'tcx> LateLintPass<'tcx> for ReturnWiderThanBody {
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
        // Only fns returning a crate-local enum directly.
        let ret_ty = cx.typeck_results().expr_ty(body.value);
        let ty::Adt(adt, _) = ret_ty.kind() else {
            return;
        };
        if !adt.is_enum() || !adt.did().is_local() {
            return;
        }
        if let Some(set) = returned_variants(cx, def_id, adt.did()) {
            self.returns.insert(def_id.to_def_id(), (adt.did(), set));
        }
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        // `match f(..) { ... E::V => panic!(), ... }` with the call as the
        // direct scrutinee: no aliasing to reason about.
        let ExprKind::Match(scrut, arms, _) = expr.kind else {
            return;
        };
        let Some(Callee::Path { def, .. }) = callee_of(cx, scrut) else {
            return;
        };
        if !def.is_local() || !matches!(cx.tcx.def_kind(def), DefKind::Fn | DefKind::AssocFn) {
            return;
        }
        for arm in arms {
            if arm.guard.is_some() {
                continue;
            }
            if let Some(variant) = arm_variant(cx, arm.pat)
                && is_panic_arm(cx, arm.body)
            {
                self.suspects.push((def, variant, arm.span));
            }
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        for (callee, variant, span) in &self.suspects {
            let Some((enum_did, returned)) = self.returns.get(callee) else {
                continue;
            };
            let adt = cx.tcx.adt_def(*enum_did);
            let Some((vidx, _)) = adt
                .variants()
                .iter_enumerated()
                .find(|(_, v)| v.def_id == *variant)
            else {
                continue;
            };
            if returned.contains(&vidx) {
                continue;
            }
            let mut names: Vec<String> = returned
                .iter()
                .map(|v| format!("`{}`", adt.variant(*v).name))
                .collect();
            names.sort();
            let fn_name = cx.tcx.item_name(*callee);
            let variant = cx.tcx.item_name(*variant);
            let enum_name = cx.tcx.item_name(*enum_did);
            emit_with_note(
                cx,
                RETURN_WIDER_THAN_BODY,
                *span,
                format!(
                    "this arm panics on `{variant}`, but `{fn_name}` only ever returns {}. The return type `{enum_name}` is wider than what the function produces",
                    names.join(", "),
                ),
                cx.tcx.def_span(*callee),
                format!("`{fn_name}`'s signature"),
                format!(
                    "give `{fn_name}` a return type without `{variant}`. This arm then cannot be written"
                ),
            );
        }
    }
}

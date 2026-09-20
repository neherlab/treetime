use std::collections::{HashMap, HashSet};

use crate::adt_facts::{field_ty, struct_field};
use crate::baseline::emit_with_note;
use crate::hir_shapes::{SelfField, assigned_adt_field, ends_in_return, peel_not, self_field};
use rustc_hir::def_id::DefId;
use rustc_hir::intravisit::FnKind;
use rustc_hir::{Body, Expr, ExprKind, FnDecl, Mutability, Stmt, StmtKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::def_id::LocalDefId;
use rustc_span::{Span, Symbol};

rustc_session::declare_lint! {
    /// Flags a bool field that two or more methods begin by testing and
    /// returning early on: `if self.flag { return ... }` as the first
    /// statement. Whether those methods may be called yet is checked at
    /// runtime, one method at a time, and only where someone remembered. A
    /// type per state checks it at compile time everywhere.
    ///
    /// The field must also be written somewhere after construction: a flag
    /// that is only ever set in a literal (`is_server`, `minify`) is a role
    /// or an option, and the methods bailing on it are not enforcing an
    /// order.
    pub RUNTIME_TYPESTATE,
    Warn,
    "bool field tested at the top of several methods to order calls at runtime"
}

#[derive(Default)]
pub struct RuntimeTypestate {
    /// (struct, field) -> how many methods bail on it at entry, and the
    /// first such test.
    guards: HashMap<(DefId, Symbol), (usize, Span)>,
    /// Fields assigned, compound-assigned, or mutably borrowed anywhere in
    /// the crate: the ones whose value is a state rather than a setting.
    written: HashSet<(DefId, Symbol)>,
}

rustc_session::impl_lint_pass!(RuntimeTypestate => [RUNTIME_TYPESTATE]);

/// `self.field` where `self` is the literal receiver.
fn self_bool_field<'tcx>(
    cx: &LateContext<'tcx>,
    e: &'tcx Expr<'tcx>,
) -> Option<(ty::AdtDef<'tcx>, Symbol)> {
    let SelfField { base, ident } = self_field(e)?;
    let ty::Adt(adt, _) = cx.typeck_results().expr_ty(base).peel_refs().kind() else {
        return None;
    };
    if !adt.did().is_local() || !adt.is_struct() {
        return None;
    }
    let is_bool = struct_field(*adt, ident.name).is_some_and(|f| field_ty(cx, f).is_bool());
    is_bool.then_some((*adt, ident.name))
}

/// The struct field a written place denotes, for any receiver: `x.flag`,
/// `(*this).flag`, `self.inner.flag`.
fn written_field<'tcx>(cx: &LateContext<'tcx>, place: &'tcx Expr<'tcx>) -> Option<(DefId, Symbol)> {
    let (adt, ident, _) = assigned_adt_field(cx, place)?;
    (adt.is_struct() && adt.did().is_local()).then(|| (adt.did(), ident.name))
}

impl<'tcx> LateLintPass<'tcx> for RuntimeTypestate {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let place = match expr.kind {
            ExprKind::Assign(place, ..) | ExprKind::AssignOp(_, place, _) => place,
            ExprKind::AddrOf(_, Mutability::Mut, place) => place,
            _ => return,
        };
        if let Some(key) = written_field(cx, place) {
            self.written.insert(key);
        }
    }

    fn check_fn(
        &mut self,
        cx: &LateContext<'tcx>,
        kind: FnKind<'tcx>,
        _decl: &'tcx FnDecl<'tcx>,
        body: &'tcx Body<'tcx>,
        _span: Span,
        _def_id: LocalDefId,
    ) {
        if matches!(kind, FnKind::Closure) {
            return;
        }
        let ExprKind::Block(block, _) = body.value.kind else {
            return;
        };
        let Some(Stmt {
            kind: StmtKind::Expr(first) | StmtKind::Semi(first),
            ..
        }) = block.stmts.first()
        else {
            return;
        };
        let ExprKind::If(cond, then, _) = first.kind else {
            return;
        };
        let Some((adt, field)) = self_bool_field(cx, peel_not(cond).0) else {
            return;
        };
        if ends_in_return(then) {
            self.guards
                .entry((adt.did(), field))
                .or_insert((0, cond.span))
                .0 += 1;
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        for ((did, field), (count, first)) in &self.guards {
            if *count < 2 || !self.written.contains(&(*did, *field)) {
                continue;
            }
            let Some(fdef) = struct_field(cx.tcx.adt_def(*did), *field) else {
                continue;
            };
            let owner = cx.tcx.def_path_str(*did);
            emit_with_note(
                cx,
                RUNTIME_TYPESTATE,
                cx.tcx.def_span(fdef.did),
                format!(
                    "{count} methods of `{owner}` begin by testing `{field}` and returning early. Whether they may be called yet is checked at runtime, one method at a time"
                ),
                *first,
                "the test at the top of one of them",
                format!(
                    "split `{owner}` into a type per state and put those {count} methods only on the one `{field}` allows. Calling them too early becomes a compile error"
                ),
            );
        }
    }
}

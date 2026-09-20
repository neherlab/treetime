use std::collections::HashMap;

use rustc_hir::def_id::DefId;
use rustc_hir::{Expr, ExprKind, Node};
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::Symbol;

use crate::adt_facts::{field_ty, is_option_ty, private_local_struct, struct_field};
use crate::baseline::emit_with_note;

rustc_session::declare_lint! {
    /// Flags an `Option` field that is unwrapped at all of its reads (two or
    /// more `unwrap`/`expect`s) and `None` is never handled, so `None` can
    /// only crash. That usually means a two-phase object:
    /// `None` between construction and setup. Storing the value directly and
    /// constructing the struct once it is available, or splitting it into a
    /// before type and an after type, deletes every one of those panics.
    pub ALWAYS_UNWRAPPED_OPTION,
    Warn,
    "Option field whose None is never handled by any reader"
}

#[derive(Default)]
struct FieldFacts {
    panicking: Vec<rustc_span::Span>,
    handled: usize,
}

#[derive(Default)]
pub struct AlwaysUnwrappedOption {
    fields: HashMap<(DefId, Symbol), FieldFacts>,
}

rustc_session::impl_lint_pass!(AlwaysUnwrappedOption => [ALWAYS_UNWRAPPED_OPTION]);

/// The crate-private local struct owning this field access, when the field is
/// an `Option`.
fn option_field_of<'tcx>(
    cx: &LateContext<'tcx>,
    base: &Expr<'tcx>,
    field: Symbol,
) -> Option<DefId> {
    let adt = private_local_struct(cx, cx.typeck_results().expr_ty_adjusted(base))?;
    let is_option = struct_field(adt, field).is_some_and(|f| is_option_ty(cx, field_ty(cx, f)));
    is_option.then(|| adt.did())
}

enum ReadKind {
    Panicking,
    Handled,
    Write,
}

/// Walk up from the field access through transparent adapters to see what
/// finally consumes it. Anything that is not literally `unwrap`/`expect` is
/// treated as handling `None` — the conservative direction: uncounted
/// handling produces silence, never a finding.
fn classify<'tcx>(cx: &LateContext<'tcx>, field_expr: &'tcx Expr<'tcx>) -> ReadKind {
    let mut current = field_expr.hir_id;
    for (parent_id, node) in cx.tcx.hir_parent_iter(field_expr.hir_id).take(6) {
        let Node::Expr(parent) = node else {
            return ReadKind::Handled;
        };
        match &parent.kind {
            ExprKind::MethodCall(seg, recv, _, _) if recv.hir_id == current => {
                match seg.ident.as_str() {
                    "unwrap" | "expect" => return ReadKind::Panicking,
                    // Transparent adapters: keep looking at what consumes the
                    // adapted value.
                    "as_ref" | "as_mut" | "clone" => {
                        current = parent.hir_id;
                    }
                    _ => return ReadKind::Handled,
                }
            }
            ExprKind::AddrOf(_, _, inner) if inner.hir_id == current => {
                current = parent.hir_id;
            }
            ExprKind::Assign(lhs, ..) | ExprKind::AssignOp(_, lhs, _) if lhs.hir_id == current => {
                return ReadKind::Write;
            }
            _ => return ReadKind::Handled,
        }
        let _ = parent_id;
    }
    ReadKind::Handled
}

impl<'tcx> LateLintPass<'tcx> for AlwaysUnwrappedOption {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let ExprKind::Field(base, ident) = expr.kind else {
            return;
        };
        let Some(adt_did) = option_field_of(cx, base, ident.name) else {
            return;
        };
        let facts = self.fields.entry((adt_did, ident.name)).or_default();
        match classify(cx, expr) {
            ReadKind::Panicking => facts.panicking.push(expr.span),
            ReadKind::Handled => facts.handled += 1,
            ReadKind::Write => {}
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        for ((adt, field), facts) in &self.fields {
            if facts.handled > 0 || facts.panicking.len() < 2 {
                continue;
            }
            let Some(fdef) = struct_field(cx.tcx.adt_def(*adt), *field) else {
                continue;
            };
            let owner = cx.tcx.item_name(*adt);
            emit_with_note(
                cx,
                ALWAYS_UNWRAPPED_OPTION,
                cx.tcx.def_span(fdef.did),
                format!(
                    "`{owner}.{field}` is unwrapped at all {} of its reads and `None` is never handled, so `None` can only crash",
                    facts.panicking.len(),
                ),
                facts.panicking[0],
                "one of the reads",
                format!(
                    "store the `{field}` value itself and build `{owner}` once you have it. Or split `{owner}` into a type without `{field}` and one with it"
                ),
            );
        }
    }
}

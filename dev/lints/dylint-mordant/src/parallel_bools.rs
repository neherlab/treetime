use std::collections::HashMap;

use crate::adt_facts::{field_ty, private_local_struct, struct_field};
use crate::baseline::{emit_with_note, join};
use crate::hir_shapes::assigned_field;
use clippy_utils::get_enclosing_block;
use rustc_hir::def_id::DefId;
use rustc_hir::{Expr, ExprKind, HirId};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::{Span, Symbol};

rustc_session::declare_lint! {
    /// Flags bool fields that are always assigned together: every write to
    /// one sits in the same block as a write to the other, across at least two
    /// functions. Together they hold one state, and the struct lets them
    /// drift apart. One enum field naming those states cannot.
    ///
    /// Only fires on structs private to the crate. One lone write to either
    /// field anywhere — `s.f = ..`, `s.f |= ..`, or either through a `Box`,
    /// guard or other `Deref` container — disproves the pairing and silences
    /// the lint.
    pub PARALLEL_BOOLS,
    Warn,
    "bool fields only ever assigned together"
}

/// One assignment to a bool field: the body and block it sits in, and the
/// assignment itself.
struct Write {
    body: DefId,
    block: HirId,
    at: Span,
}

#[derive(Default)]
pub struct ParallelBools {
    /// (struct, field) -> the sites that assign it.
    writes: HashMap<(DefId, Symbol), Vec<Write>>,
}

rustc_session::impl_lint_pass!(ParallelBools => [PARALLEL_BOOLS]);

/// The crate-private local struct behind `ty`, if it has 2+ bool fields.
fn relevant_struct<'tcx>(cx: &LateContext<'tcx>, ty: ty::Ty<'tcx>) -> Option<ty::AdtDef<'tcx>> {
    let adt = private_local_struct(cx, ty)?;
    let bools = adt
        .non_enum_variant()
        .fields
        .iter()
        .filter(|f| field_ty(cx, f).is_bool())
        .count();
    (bools >= 2).then_some(adt)
}

impl<'tcx> LateLintPass<'tcx> for ParallelBools {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let (ExprKind::Assign(place, _, _) | ExprKind::AssignOp(_, place, _)) = expr.kind else {
            return;
        };
        let Some((base, ident, _)) = assigned_field(place) else {
            return;
        };
        // The adjusted type: a write through a `Box`, a guard or any other
        // `Deref` container is a write to the struct behind it, and one such
        // lone write disproves the pairing like any other.
        let Some(adt) = relevant_struct(cx, cx.typeck_results().expr_ty_adjusted(base)) else {
            return;
        };
        let field_is_bool =
            struct_field(adt, ident.name).is_some_and(|f| field_ty(cx, f).is_bool());
        if !field_is_bool {
            return;
        }
        let Some(block) = get_enclosing_block(cx, expr.hir_id) else {
            return;
        };
        let body = cx.tcx.hir_enclosing_body_owner(expr.hir_id).to_def_id();
        self.writes
            .entry((adt.did(), ident.name))
            .or_default()
            .push(Write {
                body,
                block: block.hir_id,
                at: expr.span,
            });
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        // Group fields per struct, then find groups whose write-site sets are
        // identical. Sites are compared as sets of blocks.
        type Sites = std::collections::HashSet<(DefId, HirId)>;
        let mut per_struct: HashMap<DefId, Vec<(Symbol, Sites, Span)>> = HashMap::new();
        for ((did, field), writes) in self.writes.drain() {
            let Some(first) = writes.iter().map(|w| w.at).min() else {
                continue;
            };
            let sites: Sites = writes.into_iter().map(|w| (w.body, w.block)).collect();
            per_struct
                .entry(did)
                .or_default()
                .push((field, sites, first));
        }
        for (did, mut fields) in per_struct {
            if fields.len() < 2 {
                continue;
            }
            fields.sort_by_key(|(f, ..)| f.as_str().to_owned());
            // Partition by identical site-sets.
            let mut groups: Vec<(&Sites, Vec<Symbol>, Span)> = Vec::new();
            for (field, sites, first) in &fields {
                match groups.iter_mut().find(|(s, ..)| *s == sites) {
                    Some((_, members, _)) => members.push(*field),
                    None => groups.push((sites, vec![*field], *first)),
                }
            }
            for (sites, members, first) in groups {
                if members.len() < 2 || sites.len() < 2 {
                    continue;
                }
                let bodies: std::collections::HashSet<DefId> =
                    sites.iter().map(|(b, _)| *b).collect();
                if bodies.len() < 2 {
                    continue;
                }
                let names: Vec<String> = members.iter().map(|s| format!("`{s}`")).collect();
                let anded = join(&names, "and");
                let names = names.join(", ");
                emit_with_note(
                    cx,
                    PARALLEL_BOOLS,
                    cx.tcx.def_span(did),
                    format!(
                        "{anded} of `{}` are always assigned together ({} blocks in {} functions). Together they hold one state, and the struct lets them drift apart",
                        cx.tcx.def_path_str(did),
                        sites.len(),
                        bodies.len(),
                    ),
                    first,
                    "one of the blocks that sets both",
                    format!(
                        "replace {names} with one enum field whose variants are the states those blocks set"
                    ),
                );
            }
        }
    }
}

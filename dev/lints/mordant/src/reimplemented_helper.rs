use std::collections::HashMap;
use std::ops::ControlFlow;

use clippy_utils::visitors::for_each_expr_without_closures;
use rustc_hir::def::DefKind;
use rustc_hir::def_id::LocalDefId;
use rustc_hir::intravisit::FnKind;
use rustc_hir::{Body, FnDecl};
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::Span;

use crate::MordantConfig;
use crate::baseline::emit_with_note;
use crate::hir_clone::{bodies_equal, body_hash, fn_sigs_equal};

rustc_session::declare_lint! {
    /// Flags a function with the same signature and the same body as
    /// another function in the crate, apart from local names: the same
    /// parameter and return types and bounds, parameters destructured the
    /// same way, and the same computation. One helper exists twice under two
    /// names, nothing ties the copies together, and a fix to one will miss
    /// the other.
    ///
    /// Bodies smaller than `reimplemented-helper-min-nodes` expression nodes
    /// (default 12) are not compared, so one-line accessors and constructors
    /// never pair up. Signatures are compared as types, so two methods that
    /// read the same field names off different `Self` types are different
    /// functions. A body containing a closure is never matched (closures are
    /// not compared structurally), and neither is a function a macro wrote.
    /// Two methods of one trait impl (`grow` and `shrink` both forwarding to
    /// one `remap`) stay quiet: the trait requires both to exist, so there is
    /// no copy to delete.
    pub REIMPLEMENTED_HELPER,
    Warn,
    "a function whose signature and body repeat another function in the crate"
}

struct FnFact {
    def: LocalDefId,
    span: Span,
}

pub struct ReimplementedHelper {
    min_nodes: usize,
    /// Body hash -> the functions with that hash, in visit order.
    fns: HashMap<u64, Vec<FnFact>>,
}

rustc_session::impl_lint_pass!(ReimplementedHelper => [REIMPLEMENTED_HELPER]);

impl ReimplementedHelper {
    pub fn new(config: &MordantConfig) -> Self {
        Self {
            min_nodes: config.reimplemented_helper_min_nodes,
            fns: HashMap::new(),
        }
    }
}

fn expr_nodes(body: &Body<'_>) -> usize {
    let mut n = 0usize;
    for_each_expr_without_closures(body.value, |_| {
        n += 1;
        ControlFlow::<()>::Continue(())
    });
    n
}

/// Both are items of one trait impl block, which the trait obliges to
/// define each of them.
fn same_trait_impl(cx: &LateContext<'_>, l: LocalDefId, r: LocalDefId) -> bool {
    let parent = cx.tcx.local_parent(l);
    parent == cx.tcx.local_parent(r)
        && matches!(cx.tcx.def_kind(parent), DefKind::Impl { of_trait: true })
}

impl<'tcx> LateLintPass<'tcx> for ReimplementedHelper {
    fn check_fn(
        &mut self,
        cx: &LateContext<'tcx>,
        kind: FnKind<'tcx>,
        _decl: &'tcx FnDecl<'tcx>,
        body: &'tcx Body<'tcx>,
        span: Span,
        def_id: LocalDefId,
    ) {
        if matches!(kind, FnKind::Closure)
            || span.from_expansion()
            || body.value.span.from_expansion()
            || expr_nodes(body) < self.min_nodes
        {
            return;
        }
        self.fns
            .entry(body_hash(cx, body.id()))
            .or_default()
            .push(FnFact { def: def_id, span });
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        // (copy, original)
        let mut findings: Vec<(&FnFact, &FnFact)> = Vec::new();
        for bucket in self.fns.values_mut() {
            if bucket.len() < 2 {
                continue;
            }
            bucket.sort_by_key(|f| f.span.lo());
            let mut distinct: Vec<&FnFact> = Vec::new();
            for f in bucket.iter() {
                let body = cx.tcx.hir_body_owned_by(f.def).id();
                let original = distinct.iter().find(|o| {
                    !same_trait_impl(cx, o.def, f.def)
                        && fn_sigs_equal(cx, o.def, f.def)
                        && bodies_equal(cx, cx.tcx.hir_body_owned_by(o.def).id(), body)
                });
                match original {
                    Some(o) => findings.push((f, o)),
                    None => distinct.push(f),
                }
            }
        }
        findings.sort_by_key(|(f, _)| f.span.lo());
        for (copy, original) in findings {
            let (copy_name, original_name) = (
                cx.tcx.def_path_str(copy.def),
                cx.tcx.def_path_str(original.def),
            );
            emit_with_note(
                cx,
                REIMPLEMENTED_HELPER,
                cx.tcx.def_span(copy.def),
                format!(
                    "`{copy_name}` has the same signature and the same body as `{original_name}`, apart from local names. A fix to one will miss the other"
                ),
                cx.tcx.def_span(original.def),
                format!("`{original_name}`, the other copy"),
                format!(
                    "delete `{copy_name}` and call `{original_name}` instead, or the other way round"
                ),
            );
        }
    }
}

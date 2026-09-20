//! Finds generic functions where a large part of the body does not use the parameters.

mod borrows;
mod classify;
mod describe;
mod instances;
mod region;
mod region_io;
mod source_span;

use rustc_hir::def_id::LocalDefId;
use rustc_hir::{ClosureKind, ExprKind};
use rustc_index::bit_set::DenseBitSet;
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::middle::codegen_fn_attrs::CodegenFnAttrFlags;
use rustc_middle::ty::TyCtxt;

use crate::MordantConfig;
use crate::hir_shapes::owns_signature;
use crate::mir_flow::mir_for;

use classify::classify;
use describe::{Finding, conversion_help, only_conversions, render_locals, report};
use instances::count_instantiations;
use region::best_shared_part;
use source_span::source_span;

rustc_session::declare_lint! {
    /// Finds a generic function where a large part of the body is the same
    /// every time it is compiled. rustc compiles the whole body once for
    /// each distinct set of concrete arguments this crate uses the function
    /// with. That repeats the statements that do not use the parameters.
    /// Moving them into a separate non-generic function removes the copies
    /// only if the compiler keeps it out of line. It usually inlines a small
    /// one back, and release links that fold identical functions already
    /// merge some copies. So measure each finding, largest first.
    ///
    /// Size is counted in statements of MIR, the compiler's form of the
    /// body before optimization, where one source line is usually several.
    /// The part is made of whole basic blocks: runs of statements that end
    /// where control branches, calls or returns. It meets three conditions.
    /// Nothing in it uses a parameter or a value whose type contains one.
    /// Control enters it at one point and leaves it at one point, or it
    /// ends at the return. Every value it reads from earlier code or
    /// produces for later code has a type without parameters. Those values
    /// become the new function's arguments and results. The finding lists
    /// them and points at one call site. A part in a loop is not reported.
    ///
    /// Runs only with `generic-body-not-generic-enabled = true` in
    /// `dylint.toml`. `generic-body-not-generic-min-statements` (default 24) is
    /// the smallest part reported, in hand-written statements. Code from a
    /// macro adds to the printed size only. The crate must use the function
    /// with `generic-body-not-generic-min-instantiations` (default 2) argument
    /// sets. Calls, uses as a value, and inserted `Deref` and `Drop` calls
    /// count, followed through generic callers. Calls through `dyn`, function
    /// pointers, other crates and compile-time evaluation are not seen.
    ///
    /// Not reported: a part entered under a branch that a const parameter
    /// decides, or one that would take a value computed from a const
    /// parameter. Nor are `async`, `gen`, `#[inline(always)]`,
    /// `#[track_caller]`, `#[target_feature]` and `#[naked]` functions.
    pub(crate) GENERIC_BODY_NOT_GENERIC,
    Warn,
    "a generic fn compiled several times where a large part of the body does not use its parameters"
}

pub(crate) struct GenericBodyNotGeneric {
    min_statements: usize,
    min_instantiations: usize,
}

rustc_session::impl_lint_pass!(GenericBodyNotGeneric => [GENERIC_BODY_NOT_GENERIC]);

impl GenericBodyNotGeneric {
    pub(crate) fn new(config: &MordantConfig) -> Self {
        Self {
            min_statements: config.generic_body_not_generic_min_statements,
            min_instantiations: config.generic_body_not_generic_min_instantiations,
        }
    }
}

/// A fn or method with a type or const parameter in scope. Not a closure.
fn generic_fn(tcx: TyCtxt<'_>, def: LocalDefId) -> bool {
    tcx.def_kind(def).is_fn_like()
        && !tcx.is_closure_like(def.to_def_id())
        && tcx.generics_of(def).requires_monomorphization(tcx)
}

/// False for `#[inline(always)]`, `#[track_caller]`, `#[target_feature]`
/// and `#[naked]`: a call added inside would change cost or behaviour.
fn splittable(tcx: TyCtxt<'_>, def: LocalDefId) -> bool {
    let attrs = tcx.codegen_fn_attrs(def);
    !attrs.inline.always()
        && !attrs
            .flags
            .intersects(CodegenFnAttrFlags::TRACK_CALLER | CodegenFnAttrFlags::NAKED)
        && attrs.target_features.is_empty()
}

/// Generic, `splittable`, not written by a macro, and not `async` or `gen`
/// (their written body is a coroutine, which `generic_fn` excludes).
fn candidate_fn(tcx: TyCtxt<'_>, def: LocalDefId) -> bool {
    if !generic_fn(tcx, def) || !splittable(tcx, def) || tcx.def_span(def).from_expansion() {
        return false;
    }
    let body = tcx.hir_body_owned_by(def).value;
    !body.span.from_expansion()
        && !matches!(body.kind, ExprKind::Closure(c) if matches!(c.kind, ClosureKind::Coroutine(_)))
}

impl<'tcx> LateLintPass<'tcx> for GenericBodyNotGeneric {
    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let tcx = cx.tcx;
        let candidates: Vec<LocalDefId> = tcx
            .hir_body_owners()
            .filter(|&def| candidate_fn(tcx, def))
            .collect();
        if candidates.is_empty() {
            return;
        }
        // Count uses first, then classify only functions used with enough argument sets.
        let counts = count_instantiations(tcx);
        let sets = |def: LocalDefId| counts.concrete.get(&def).map_or(0, |s| s.len());
        let mut findings: Vec<Finding> = candidates
            .into_iter()
            .filter(|&def| sets(def) >= self.min_instantiations)
            .filter_map(|def| {
                let body = mir_for(tcx, def)?;
                let (facts, total) = classify(tcx, &body);
                let part = best_shared_part(tcx, &body, &facts, self.min_statements)?;
                let mut after = DenseBitSet::new_filled(body.basic_blocks.len());
                after.subtract(&part.blocks);
                let site = source_span(tcx, def, &body, &part.blocks);
                Some(Finding {
                    def,
                    site,
                    size: part.size,
                    reads: render_locals(cx, def, &body, &part.params, &part.blocks),
                    produces: render_locals(cx, def, &body, &part.returns, &after),
                    other_parts: part.other_parts,
                    total,
                    // Suggest a new signature only where this crate may change it.
                    signature_help: owns_signature(cx, def)
                        .then(|| conversion_help(&body, &only_conversions(tcx, &body)))
                        .flatten(),
                })
            })
            .collect();
        findings.sort_by_key(|s| tcx.def_span(s.def).lo());
        for finding in &findings {
            let site = counts.first_site.get(&finding.def).copied();
            report(cx, finding, sets(finding.def), site, self.min_statements);
        }
    }
}

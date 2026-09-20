//! Counts, for each generic function, the distinct sets of concrete arguments
//! this crate uses it with.

use std::collections::VecDeque;
use std::ops::ControlFlow;

use clippy_utils::visitors::for_each_expr_without_closures;
use rustc_data_structures::fx::{FxHashMap, FxHashSet, FxIndexMap, FxIndexSet};
use rustc_data_structures::stack::ensure_sufficient_stack;
use rustc_hir::def::DefKind;
use rustc_hir::def_id::{DefId, LocalDefId};
use rustc_hir::{ConstContext, Expr, ExprKind, LangItem, Node};
use rustc_middle::mir::TerminatorKind;
use rustc_middle::ty::adjustment::{Adjust, DerefAdjustKind};
use rustc_middle::ty::{
    self, GenericArgsRef, Ty, TyCtxt, TypeSuperVisitable, TypeVisitable, TypeVisitableExt,
    TypeVisitor, TypeckResults,
};
use rustc_span::Span;

use super::generic_fn;
use crate::mir_flow::mir_for;

/// A use of a generic fn, found in `caller` (a typeck root), whose arguments
/// still contain `caller`'s parameters. Dropping a `T` is a use of `drop_glue::<T>`.
struct DependentUse<'tcx> {
    caller: LocalDefId,
    callee: DefId,
    args: GenericArgsRef<'tcx>,
    site: Span,
}

type Instantiation<'tcx> = (LocalDefId, GenericArgsRef<'tcx>);

#[derive(Default)]
pub(super) struct InstantiationCounts<'tcx> {
    pub(super) concrete: FxIndexMap<LocalDefId, FxIndexSet<GenericArgsRef<'tcx>>>,
    pub(super) first_site: FxHashMap<LocalDefId, (Span, GenericArgsRef<'tcx>)>,
    dependent_uses: Vec<DependentUse<'tcx>>,
    queue: VecDeque<Instantiation<'tcx>>,
    drop_fn: Option<DefId>,
}

impl<'tcx> InstantiationCounts<'tcx> {
    fn add_use(
        &mut self,
        tcx: TyCtxt<'tcx>,
        caller: LocalDefId,
        callee: DefId,
        args: GenericArgsRef<'tcx>,
        site: Span,
    ) {
        if args.has_non_region_param() {
            self.dependent_uses.push(DependentUse {
                caller,
                callee,
                args,
                site,
            });
        } else {
            self.record(tcx, callee, args, site);
        }
    }

    /// `drop_glue::<T>` is one `Drop::drop` use per local `Drop` impl inside `T`.
    fn record(&mut self, tcx: TyCtxt<'tcx>, callee: DefId, args: GenericArgsRef<'tcx>, site: Span) {
        if !tcx.is_lang_item(callee, LangItem::DropGlue) {
            self.count(tcx, callee, args, site);
        } else if let Some(drop_fn) = self.drop_fn {
            for part in types_dropped_by_local_impl(tcx, args.type_at(0)) {
                self.count(tcx, drop_fn, tcx.mk_args(&[part.into()]), site);
            }
        }
    }

    fn count(&mut self, tcx: TyCtxt<'tcx>, callee: DefId, args: GenericArgsRef<'tcx>, site: Span) {
        let Some((item, args)) = resolve(tcx, callee, args) else {
            return;
        };
        let lists = self.concrete.entry(item).or_default();
        if lists.insert(args) && lists.len() <= MAX_ARGUMENT_LISTS_PER_FN {
            self.queue.push_back((item, args));
        }
        let first = self.first_site.entry(item).or_insert((site, args));
        if site.lo() < first.0.lo() {
            *first = (site, args);
        }
    }
}

/// The local generic fn whose body is compiled for this use, and the
/// arguments as that body writes them. `dyn`, closure shims and foreign fns give `None`.
fn resolve<'tcx>(
    tcx: TyCtxt<'tcx>,
    callee: DefId,
    args: GenericArgsRef<'tcx>,
) -> Option<Instantiation<'tcx>> {
    // A foreign trait's item can still resolve to an impl in this crate.
    if !callee.is_local() && tcx.trait_of_assoc(callee).is_none() {
        return None;
    }
    let env = ty::TypingEnv::fully_monomorphized();
    let args = tcx
        .try_normalize_erasing_regions(env, ty::Unnormalized::new_wip(args))
        .ok()?;
    if args.has_non_region_param() {
        return None;
    }
    let instance = ty::Instance::try_resolve(tcx, env, callee, args).ok()??;
    let ty::InstanceKind::Item(item) = instance.def else {
        return None;
    };
    let item = item.as_local()?;
    if !generic_fn(tcx, item) {
        return None;
    }
    let args = tcx.erase_and_anonymize_regions(instance.args);
    (args.len() == tcx.generics_of(item).count() && !args.has_non_region_param())
        .then_some((item, args))
}

/// `Drop::drop`, or `None` when no `Drop` impl in this crate is generic: no drop counts then.
fn drop_fn(tcx: TyCtxt<'_>) -> Option<DefId> {
    let drop_trait = tcx.lang_items().drop_trait()?;
    let generic = |&imp: &LocalDefId| tcx.generics_of(imp).requires_monomorphization(tcx);
    if !tcx.local_trait_impls(drop_trait).iter().any(generic) {
        return None;
    }
    tcx.associated_items(drop_trait)
        .in_definition_order()
        .find(|item| item.is_fn())
        .map(|item| item.def_id)
}

/// The types inside a `ty` value whose local `Drop` impl runs with it. Not
/// behind a pointer, `ManuallyDrop`, `union`, `dyn`, or a foreign `Drop`.
fn types_dropped_by_local_impl<'tcx>(tcx: TyCtxt<'tcx>, ty: Ty<'tcx>) -> Vec<Ty<'tcx>> {
    let env = ty::TypingEnv::fully_monomorphized();
    let mut local = Vec::new();
    let Ok(ty) = tcx.try_normalize_erasing_regions(env, ty::Unnormalized::new_wip(ty)) else {
        return local;
    };
    let mut seen: FxHashSet<Ty<'tcx>> = FxHashSet::default();
    let mut stack = vec![ty];
    while let Some(ty) = stack.pop() {
        if !seen.insert(ty) || !ty.needs_drop(tcx, env) {
            continue;
        }
        match *ty.kind() {
            ty::Adt(def, args) => {
                if def.is_manually_drop() {
                    continue;
                }
                if def.destructor(tcx).is_some_and(|d| d.did.is_local()) {
                    local.push(ty);
                }
                if def.is_union() {
                    continue;
                }
                // A `Box`'s fields are a pointer. What it drops is its type argument.
                if def.is_box() {
                    stack.extend(args.types());
                }
                stack.extend(def.all_fields().filter_map(|field| {
                    tcx.try_normalize_erasing_regions(env, field.ty(tcx, args))
                        .ok()
                }));
            }
            ty::Array(elem, _) | ty::Slice(elem) | ty::Pat(elem, _) => stack.push(elem),
            ty::Tuple(tys) => stack.extend(tys),
            ty::Closure(_, args) => stack.extend(args.as_closure().upvar_tys()),
            ty::CoroutineClosure(_, args) => {
                stack.extend(args.as_coroutine_closure().upvar_tys());
            }
            ty::Coroutine(def, args) => {
                stack.extend(args.as_coroutine().upvar_tys());
                if let Some(layout) = tcx.mir_coroutine_witnesses(def) {
                    stack.extend(layout.field_tys.iter().filter_map(|saved| {
                        tcx.try_instantiate_and_normalize_erasing_regions(
                            args,
                            env,
                            ty::EarlyBinder::bind(saved.ty),
                        )
                        .ok()
                    }));
                }
            }
            _ => {}
        }
    }
    local
}

/// Every fn `e` uses and each overloaded auto-deref on it. A path's `FnDef`
/// type has its arguments even where the node has none (`for` desugaring).
fn fn_uses<'tcx>(
    tcx: TyCtxt<'tcx>,
    typeck: &TypeckResults<'tcx>,
    e: &Expr<'_>,
    mut each: impl FnMut(DefId, GenericArgsRef<'tcx>, Span),
) {
    let Some(ty) = typeck.expr_ty_opt(e) else {
        return;
    };
    let named = match e.kind {
        ExprKind::Path(..) => match *ty.kind() {
            ty::FnDef(def, args) => Some((def, args)),
            _ => None,
        },
        _ => typeck
            .type_dependent_def_id(e.hir_id)
            .map(|def| (def, typeck.node_args(e.hir_id))),
    };
    // Constructors are `FnDef`s too. `resolve` needs one argument per parameter.
    if let Some((def, args)) = named
        && matches!(tcx.def_kind(def), DefKind::Fn | DefKind::AssocFn)
        && args.len() == tcx.generics_of(def).count()
    {
        each(def, args, use_site(tcx, e));
    }
    let mut source = ty;
    for adjustment in typeck.expr_adjustments(e) {
        if let Adjust::Deref(DerefAdjustKind::Overloaded(deref)) = adjustment.kind {
            each(
                deref.method_call(tcx),
                tcx.mk_args(&[source.into()]),
                e.span,
            );
        }
        source = adjustment.target;
    }
}

fn use_site(tcx: TyCtxt<'_>, e: &Expr<'_>) -> Span {
    match tcx.parent_hir_node(e.hir_id) {
        Node::Expr(
            call @ Expr {
                kind: ExprKind::Call(callee, _),
                ..
            },
        ) if callee.hir_id == e.hir_id => call.span,
        _ => e.span,
    }
}

/// Not rustc's `collect_and_partition_mono_items`: that builds optimized
/// MIR, which steals the body every MIR lint here reads.
pub(super) fn count_instantiations<'tcx>(tcx: TyCtxt<'tcx>) -> InstantiationCounts<'tcx> {
    let mut counts = InstantiationCounts {
        drop_fn: drop_fn(tcx),
        ..InstantiationCounts::default()
    };
    for owner in tcx.hir_body_owners() {
        if !tcx.has_typeck_results(owner) {
            continue;
        }
        // Evaluated at compile time: what they call is not in the binary.
        if let Some(ConstContext::Const { .. } | ConstContext::Static(_)) =
            tcx.hir_body_const_context(owner)
        {
            continue;
        }
        let caller = tcx.typeck_root_def_id_local(owner);
        let typeck = tcx.typeck(owner);
        let body = tcx.hir_body_owned_by(owner);
        for_each_expr_without_closures(body.value, |e| {
            fn_uses(tcx, typeck, e, |callee, args, site| {
                counts.add_use(tcx, caller, callee, args, site);
            });
            ControlFlow::<()>::Continue(())
        });
        // Drops have no HIR node, so read MIR `Drop` terminators. Slow, so only if one counts.
        if counts.drop_fn.is_some()
            && let Some(mir) = mir_for(tcx, owner)
            && let Some(drop_glue) = tcx.lang_items().drop_glue_fn()
        {
            for data in mir.basic_blocks.iter() {
                if let Some(terminator) = &data.terminator
                    && let TerminatorKind::Drop { place, .. } = terminator.kind
                {
                    let dropped = place.ty(&mir.local_decls, tcx).ty;
                    let site = mir.local_decls[place.local].source_info.span;
                    counts.add_use(tcx, caller, drop_glue, tcx.mk_args(&[dropped.into()]), site);
                }
            }
        }
    }
    propagate(tcx, &mut counts);
    counts
}

/// Reached only by `f::<T>` calling `f::<A<T>>` and `f::<B<T>>`, which does not build.
const MAX_ARGUMENT_LISTS_PER_FN: usize = 1 << 16;

/// Substitutes each queued argument list into that fn's `DependentUse`s and
/// records the results. Types past the recursion limit stop `f::<W<T>>`.
fn propagate<'tcx>(tcx: TyCtxt<'tcx>, counts: &mut InstantiationCounts<'tcx>) {
    let limit = tcx.recursion_limit().0;
    let dependent_uses = std::mem::take(&mut counts.dependent_uses);
    let mut uses: FxHashMap<LocalDefId, Vec<&DependentUse<'tcx>>> = FxHashMap::default();
    for edge in &dependent_uses {
        uses.entry(edge.caller).or_default().push(edge);
    }
    let env = ty::TypingEnv::fully_monomorphized();
    let mut type_depth = TypeDepth::default();
    while let Some((caller, caller_args)) = counts.queue.pop_front() {
        let Some(edges) = uses.get(&caller) else {
            continue;
        };
        for edge in edges {
            let Ok(args) = tcx.try_instantiate_and_normalize_erasing_regions(
                caller_args,
                env,
                ty::EarlyBinder::bind(edge.args),
            ) else {
                continue;
            };
            if args.has_non_region_param() || type_depth.of(args) > limit {
                continue;
            }
            counts.record(tcx, edge.callee, args, edge.site);
        }
    }
}

/// Measures how deeply types nest: `u8` is 1, `&[Vec<u8>]` 4. Each distinct
/// type is measured once, since `f::<(T, T)>` doubles the type per step.
#[derive(Default)]
struct TypeDepth<'tcx> {
    depth: usize,
    measured: FxHashMap<Ty<'tcx>, usize>,
}

impl<'tcx> TypeDepth<'tcx> {
    fn of(&mut self, args: GenericArgsRef<'tcx>) -> usize {
        self.depth = 0;
        args.visit_with(self);
        self.depth
    }
}

impl<'tcx> TypeVisitor<TyCtxt<'tcx>> for TypeDepth<'tcx> {
    fn visit_ty(&mut self, ty: Ty<'tcx>) {
        let depth = match self.measured.get(&ty) {
            Some(&depth) => depth,
            None => {
                let outer = std::mem::take(&mut self.depth);
                ensure_sufficient_stack(|| ty.super_visit_with(self));
                let depth = std::mem::replace(&mut self.depth, outer) + 1;
                self.measured.insert(ty, depth);
                depth
            }
        };
        self.depth = self.depth.max(depth);
    }
}

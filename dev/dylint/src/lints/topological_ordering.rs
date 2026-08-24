#![allow(
    clippy::indexing_slicing,
    reason = "graph algorithm indices are always in-bounds"
)]

//! Lint enforcing the project's item ordering convention within a module.
//!
//! The TreeTime style reads top-to-bottom: a referencing item (caller, user of
//! a type) appears *before* the item it references (callee, dependency).  Public
//! API sits near the top and its private helpers follow below.  This is the
//! reverse of a classic callee-first topological sort.
//!
//! `const` / `static` items are exempt as ordering *targets*: they are placed at
//! the top of the module by convention, so a function referencing a `const`
//! declared above it is never flagged.
//!
//! Cycles (mutual references) are handled by collapsing strongly connected
//! components into single nodes: items within an SCC are unordered relative to
//! each other, but the SCC as a whole is ordered relative to outside items.

use clippy_utils::diagnostics::span_lint_hir_and_then;
use clippy_utils::is_in_test;
use rustc_data_structures::fx::{FxHashMap, FxHashSet};
use rustc_hir as hir;
use rustc_hir::def::DefKind;
use rustc_lint::{LateContext, LateLintPass, LintContext as _};
use rustc_middle::ty::TyCtxt;
use rustc_span::Span;
use rustc_span::def_id::LocalDefId;

use super::hir_refs;

rustc_session::declare_lint! {
    /// Flags items that appear out of the project's ordering convention within a
    /// module.
    ///
    /// Caller-first: an item should appear before the items it references, so
    /// the module reads top-to-bottom from public API down to private helpers.
    ///
    /// Cycles are tolerated: items in a strongly connected component are
    /// unordered relative to each other.
    pub TOPOLOGICAL_ORDERING,
    Warn,
    "items are not ordered callers-before-callees in this module"
}

fn is_relevant_item_kind(kind: &hir::ItemKind<'_>) -> bool {
    matches!(
        kind,
        hir::ItemKind::Fn { .. }
            | hir::ItemKind::Struct(..)
            | hir::ItemKind::Enum(..)
            | hir::ItemKind::Trait { .. }
            | hir::ItemKind::TyAlias(..)
            | hir::ItemKind::Const(..)
            | hir::ItemKind::Static(..)
            | hir::ItemKind::Impl(..)
    )
}

fn item_display_name(cx: &LateContext<'_>, item: &hir::Item<'_>) -> String {
    let def_id = item.owner_id.to_def_id();
    let kind = cx.tcx.def_kind(def_id);

    if let DefKind::Impl { .. } = kind {
        if let hir::ItemKind::Impl(impl_block) = &item.kind {
            let ty_name = cx
                .sess()
                .source_map()
                .span_to_snippet(impl_block.self_ty.span)
                .unwrap_or_else(|_| "?".into());
            return if let Some(trait_ref) = &impl_block.of_trait {
                let trait_name = cx
                    .sess()
                    .source_map()
                    .span_to_snippet(trait_ref.trait_ref.path.span)
                    .unwrap_or_else(|_| "..".into());
                format!("impl {trait_name} for {ty_name}")
            } else {
                format!("impl {ty_name}")
            };
        }
    }

    let name = cx
        .tcx
        .opt_item_name(def_id)
        .map_or_else(|| "?".into(), |s| s.to_string());

    #[expect(
        clippy::wildcard_enum_match_arm,
        reason = "only a few kinds get prefixes"
    )]
    let prefix = match kind {
        DefKind::Fn => "fn",
        DefKind::Struct => "struct",
        DefKind::Enum => "enum",
        DefKind::Trait => "trait",
        DefKind::TyAlias => "type",
        DefKind::Const { .. } => "const",
        DefKind::Static { .. } => "static",
        _ => "",
    };

    if prefix.is_empty() {
        name
    } else {
        format!("{prefix} {name}")
    }
}

fn resolve_self_ty_def_id(cx: &LateContext<'_>, impl_block: &hir::Impl<'_>) -> Option<LocalDefId> {
    // Cannot reuse `hir_refs::resolve_ty_def_id` due to `AmbigArg` type mismatch
    // on `Impl::self_ty`.
    if let hir::TyKind::Path(ref qpath) = impl_block.self_ty.kind
        && let hir::def::Res::Def(_, def_id) = cx.qpath_res(qpath, impl_block.self_ty.hir_id)
    {
        def_id.as_local()
    } else {
        None
    }
}

fn find_module_item(
    tcx: TyCtxt<'_>,
    def_id: LocalDefId,
    map: &FxHashMap<LocalDefId, (LocalDefId, usize)>,
) -> Option<(LocalDefId, usize)> {
    if let Some(&result) = map.get(&def_id) {
        return Some(result);
    }
    let mut current = def_id.to_def_id();
    loop {
        let parent = tcx.opt_parent(current);
        let Some(parent) = parent else {
            return None;
        };
        if let Some(local) = parent.as_local() {
            if let Some(&result) = map.get(&local) {
                return Some(result);
            }
        }
        current = parent;
    }
}

struct ModuleItem {
    def_id: LocalDefId,
    span: Span,
    /// Tight span covering just the declaration head (e.g. `fn callee` or
    /// `struct Config`), used for "defined here" labels.
    ident_span: Span,
    display_name: String,
    /// `true` for `const` / `static` items, which are exempt as ordering targets
    /// (they live at the top of the module by convention).
    is_const_or_static: bool,
    /// When the impl's self type is defined in this same module, the `LocalDefId`
    /// of the self type. Used to group the impl with its type for ordering.
    impl_self_ty: Option<LocalDefId>,
}

/// Remap refs so impl items share their self-type's index (merging edges from
/// `impl T` into edges from `T`), then drop self-loops.
fn remap_impl_refs(
    items: &[ModuleItem],
    def_id_to_idx: &FxHashMap<LocalDefId, usize>,
    refs: &[(usize, usize, Span)],
) -> Vec<(usize, usize, Span)> {
    let remap = |idx: usize| -> usize {
        items[idx]
            .impl_self_ty
            .and_then(|self_ty| def_id_to_idx.get(&self_ty).copied())
            .unwrap_or(idx)
    };

    refs.iter()
        .map(|&(from, to, span)| (remap(from), remap(to), span))
        .filter(|(from, to, _)| from != to)
        .collect()
}

fn build_adj_list(refs: &[(usize, usize, Span)], n: usize) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); n];
    for &(from, to, _) in refs {
        adj[from].push(to);
    }
    for list in &mut adj {
        list.sort_unstable();
        list.dedup();
    }
    adj
}

struct TarjanState {
    index_counter: usize,
    stack: Vec<usize>,
    on_stack: Vec<bool>,
    index: Vec<usize>,
    lowlink: Vec<usize>,
    sccs: Vec<Vec<usize>>,
}

fn strongconnect(v: usize, adj: &[Vec<usize>], state: &mut TarjanState) {
    state.index[v] = state.index_counter;
    state.lowlink[v] = state.index_counter;
    state.index_counter += 1;
    state.stack.push(v);
    state.on_stack[v] = true;

    for &w in &adj[v] {
        if state.index[w] == usize::MAX {
            strongconnect(w, adj, state);
            state.lowlink[v] = state.lowlink[v].min(state.lowlink[w]);
        } else if state.on_stack[w] {
            state.lowlink[v] = state.lowlink[v].min(state.index[w]);
        }
    }

    if state.lowlink[v] == state.index[v] {
        let mut scc = Vec::new();
        while let Some(w) = state.stack.pop() {
            state.on_stack[w] = false;
            scc.push(w);
            if w == v {
                break;
            }
        }
        state.sccs.push(scc);
    }
}

/// Returns a vector mapping each item index to its SCC id.
fn compute_sccs(adj: &[Vec<usize>], n: usize) -> Vec<usize> {
    let mut state = TarjanState {
        index_counter: 0,
        stack: Vec::new(),
        on_stack: vec![false; n],
        index: vec![usize::MAX; n],
        lowlink: vec![0; n],
        sccs: Vec::new(),
    };
    for v in 0..n {
        if state.index[v] == usize::MAX {
            strongconnect(v, adj, &mut state);
        }
    }
    let mut item_to_scc = vec![0usize; n];
    for (scc_idx, scc) in state.sccs.iter().enumerate() {
        for &item_idx in scc {
            item_to_scc[item_idx] = scc_idx;
        }
    }
    item_to_scc
}

struct OrderingViolation {
    item_def_id: LocalDefId,
    ref_span: Span,
    item_name: String,
    /// `(referenced item name, ref_span, def_span)`
    witnesses: Vec<(String, Span, Span)>,
}

struct GroupingViolation {
    impl_def_id: LocalDefId,
    impl_span: Span,
    impl_name: String,
    type_name: String,
    type_span: Span,
}

fn find_ordering_violations(
    items: &[ModuleItem],
    refs: &[(usize, usize, Span)],
    item_to_scc: &[usize],
) -> Vec<OrderingViolation> {
    let mut violation_map: FxHashMap<usize, Vec<(String, Span, Span)>> = FxHashMap::default();

    for &(from, to, ref_span) in refs {
        if item_to_scc[from] == item_to_scc[to] {
            continue;
        }

        // `const` / `static` targets are placed at the top of the module by
        // convention; referencing them from below is expected, not a violation.
        if items[to].is_const_or_static {
            continue;
        }

        // Project convention is caller-first: the referencing item (`from`) must
        // appear before the item it references (`to`). A violation is a
        // referencing item that appears *after* its referenced item.
        if from > to {
            violation_map
                .entry(from)
                .or_default()
                .push((items[to].display_name.clone(), ref_span, items[to].ident_span));
        }
    }

    let mut violations: Vec<OrderingViolation> = violation_map
        .into_iter()
        .map(|(from_idx, witnesses)| OrderingViolation {
            item_def_id: items[from_idx].def_id,
            ref_span: witnesses[0].1,
            item_name: items[from_idx].display_name.clone(),
            witnesses,
        })
        .collect();
    violations.sort_by_key(|v| v.ref_span.lo());
    violations
}

fn check_impl_grouping(
    items: &[ModuleItem],
    def_id_to_idx: &FxHashMap<LocalDefId, usize>,
) -> Vec<GroupingViolation> {
    let mut violations = Vec::new();

    for (impl_pos, item) in items.iter().enumerate() {
        let Some(self_ty_def_id) = item.impl_self_ty else {
            continue;
        };
        let Some(&type_idx) = def_id_to_idx.get(&self_ty_def_id) else {
            continue;
        };
        let type_item = &items[type_idx];

        let (lo, hi) = if impl_pos > type_idx {
            (type_idx, impl_pos)
        } else {
            (impl_pos, type_idx)
        };
        let has_unrelated = items[lo + 1..hi]
            .iter()
            .any(|i| i.def_id != self_ty_def_id && i.impl_self_ty != Some(self_ty_def_id));

        if has_unrelated {
            violations.push(GroupingViolation {
                impl_def_id: item.def_id,
                impl_span: item.span,
                impl_name: item.display_name.clone(),
                type_name: type_item
                    .display_name
                    .split_once(' ')
                    .map_or(type_item.display_name.as_str(), |(_, name)| name)
                    .to_string(),
                type_span: type_item.span,
            });
        }
    }

    violations
}

fn emit_module_diagnostic(
    cx: &LateContext<'_>,
    ordering_violations: &[OrderingViolation],
    grouping_violations: &[GroupingViolation],
) {
    if !ordering_violations.is_empty() {
        let first = &ordering_violations[0];
        let hir_id = cx.tcx.local_def_id_to_hir_id(first.item_def_id);

        span_lint_hir_and_then(
            cx,
            TOPOLOGICAL_ORDERING,
            hir_id,
            first.ref_span,
            "items are not ordered callers-before-callees in this module",
            |diag| {
                let mut seen_defs = FxHashSet::default();
                for violation in ordering_violations {
                    for (name, ref_span, def_span) in &violation.witnesses {
                        diag.span_label(
                            *ref_span,
                            format!(
                                "`{}` references `{name}` but appears after it",
                                violation.item_name,
                            ),
                        );
                        if seen_defs.insert(*def_span) {
                            diag.span_label(
                                *def_span,
                                format!("`{name}` defined here"),
                            );
                        }
                    }
                }
                diag.help(
                    "reorder items so referencing items appear before the items they reference (callers before callees)",
                );
            },
        );
    }

    for violation in grouping_violations {
        let hir_id = cx.tcx.local_def_id_to_hir_id(violation.impl_def_id);

        span_lint_hir_and_then(
            cx,
            TOPOLOGICAL_ORDERING,
            hir_id,
            violation.impl_span,
            format!(
                "`{}` is separated from its type definition",
                violation.impl_name
            ),
            |diag| {
                diag.span_label(
                    violation.type_span,
                    format!("`{}` defined here", violation.type_name),
                );
                // Suppress when ordering violations exist: the reorder
                // fix will implicitly place impls next to their types.
                if ordering_violations.is_empty() {
                    diag.help("move the impl block adjacent to its type definition");
                }
            },
        );
    }
}

/// A reference edge from the containing owner to the referenced target,
/// collected eagerly in `check_expr` / `check_ty`. Owner -> module-item
/// resolution is deferred to `check_crate_post` (items are not yet collected
/// when the ref is seen).
struct RawRef {
    source_owner: LocalDefId,
    target: LocalDefId,
    ref_span: Span,
}

pub struct TopologicalOrdering {
    modules: FxHashMap<LocalDefId, Vec<ModuleItem>>,
    raw_refs: Vec<RawRef>,
    /// Short-circuits all callbacks when the lint is `#[allow]`'d at crate level.
    enabled: bool,
}

impl TopologicalOrdering {
    pub fn new() -> Self {
        Self {
            modules: FxHashMap::default(),
            raw_refs: Vec::new(),
            enabled: false,
        }
    }

    fn record_ref(&mut self, owner: LocalDefId, def_id: rustc_span::def_id::DefId, span: Span) {
        if let Some(target) = def_id.as_local() {
            self.raw_refs.push(RawRef {
                source_owner: owner,
                target,
                ref_span: span,
            });
        }
    }
}

rustc_session::impl_lint_pass!(TopologicalOrdering => [TOPOLOGICAL_ORDERING]);

impl<'tcx> LateLintPass<'tcx> for TopologicalOrdering {
    fn check_crate(&mut self, cx: &LateContext<'tcx>) {
        self.enabled =
            !clippy_utils::is_lint_allowed(cx, TOPOLOGICAL_ORDERING, rustc_hir::CRATE_HIR_ID);
    }

    fn check_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx hir::Item<'tcx>) {
        if !self.enabled {
            return;
        }
        if item.span.from_expansion() {
            return;
        }
        if !is_relevant_item_kind(&item.kind) {
            return;
        }

        let item_def_id = item.owner_id.def_id;

        let parent_def_id = cx.tcx.parent(item_def_id.to_def_id());
        let Some(parent_local) = parent_def_id.as_local() else {
            return;
        };
        if cx.tcx.def_kind(parent_local) != DefKind::Mod {
            return;
        }

        if is_in_test(cx.tcx, item.hir_id()) {
            return;
        }

        let display_name = item_display_name(cx, item);

        let is_const_or_static =
            matches!(item.kind, hir::ItemKind::Const(..) | hir::ItemKind::Static(..));

        let impl_self_ty = if let hir::ItemKind::Impl(impl_block) = &item.kind {
            resolve_self_ty_def_id(cx, impl_block)
        } else {
            None
        };

        let ident_span = cx
            .tcx
            .def_ident_span(item_def_id.to_def_id())
            .map_or(item.span, |id_sp| item.span.with_hi(id_sp.hi()));

        self.modules.entry(parent_local).or_default().push(ModuleItem {
            def_id: item_def_id,
            span: item.span,
            ident_span,
            display_name,
            is_const_or_static,
            impl_self_ty,
        });
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx hir::Expr<'tcx>) {
        if !self.enabled || expr.span.from_expansion() {
            return;
        }
        if let Some((def_id, _, span)) = hir_refs::resolve_expr_def_id(cx, expr) {
            self.record_ref(expr.hir_id.owner.def_id, def_id, span);
        }
    }

    fn check_ty(&mut self, cx: &LateContext<'tcx>, ty: &'tcx hir::Ty<'tcx, hir::AmbigArg>) {
        if !self.enabled || ty.span.from_expansion() {
            return;
        }
        if let Some((def_id, _, span)) = hir_refs::resolve_ty_def_id(cx, ty) {
            self.record_ref(ty.hir_id.owner.def_id, def_id, span);
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        if !self.enabled {
            return;
        }

        let mut def_id_to_module_item: FxHashMap<LocalDefId, (LocalDefId, usize)> =
            FxHashMap::default();
        for (&module_def_id, items) in &self.modules {
            for (idx, item) in items.iter().enumerate() {
                def_id_to_module_item.insert(item.def_id, (module_def_id, idx));
            }
        }

        let mut module_refs: FxHashMap<LocalDefId, Vec<(usize, usize, Span)>> =
            FxHashMap::default();
        let mut resolve_cache: FxHashMap<LocalDefId, Option<(LocalDefId, usize)>> =
            FxHashMap::default();
        for raw_ref in &self.raw_refs {
            let source = *resolve_cache
                .entry(raw_ref.source_owner)
                .or_insert_with(|| {
                    find_module_item(cx.tcx, raw_ref.source_owner, &def_id_to_module_item)
                });
            let target = *resolve_cache.entry(raw_ref.target).or_insert_with(|| {
                find_module_item(cx.tcx, raw_ref.target, &def_id_to_module_item)
            });

            if let (Some((src_mod, src_idx)), Some((tgt_mod, tgt_idx))) = (source, target) {
                if src_mod == tgt_mod && src_idx != tgt_idx {
                    module_refs.entry(src_mod).or_default().push((
                        src_idx,
                        tgt_idx,
                        raw_ref.ref_span,
                    ));
                }
            }
        }

        for (&module_def_id, items) in &self.modules {
            if items.len() <= 1 {
                continue;
            }

            let resolved_refs = module_refs
                .get(&module_def_id)
                .map_or(&[][..], Vec::as_slice);

            let item_def_id_to_idx: FxHashMap<LocalDefId, usize> = items
                .iter()
                .enumerate()
                .map(|(idx, item)| (item.def_id, idx))
                .collect();

            let remapped_refs = remap_impl_refs(items, &item_def_id_to_idx, resolved_refs);

            let n = items.len();
            let adj = build_adj_list(&remapped_refs, n);
            let item_to_scc = compute_sccs(&adj, n);

            let ordering_violations =
                find_ordering_violations(items, &remapped_refs, &item_to_scc);

            let grouping_violations = check_impl_grouping(items, &item_def_id_to_idx);

            if ordering_violations.is_empty() && grouping_violations.is_empty() {
                continue;
            }

            emit_module_diagnostic(cx, &ordering_violations, &grouping_violations);
        }
    }
}

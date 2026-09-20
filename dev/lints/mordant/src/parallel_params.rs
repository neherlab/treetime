use std::collections::{HashMap, HashSet};

use clippy_utils::res::MaybeResPath;
use clippy_utils::source::snippet_opt;
use clippy_utils::visitors::for_each_expr;
use rustc_hir::def::{DefKind, Res};
use rustc_hir::def_id::DefId;
use rustc_hir::intravisit::FnKind;
use rustc_hir::{Body, Expr, ExprKind, FnDecl, HirId, PatKind, UnOp};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::print::with_no_trimmed_paths;
use rustc_middle::ty::{self, Mutability, Ty};
use rustc_span::def_id::LocalDefId;
use rustc_span::symbol::kw;
use rustc_span::{Span, Symbol};

use crate::MordantConfig;
use crate::baseline::{emit_with_note, join};
use crate::hir_shapes::{Callee, callee_of, owns_signature};

rustc_session::declare_lint! {
    /// Flags two or more plain-data parameters that `parallel-params-min-fns`
    /// or more crate-private functions declare alike (same names and types)
    /// and pass between them unchanged: every counted function passes the
    /// group in one call to another of them, or receives it that way. The
    /// group is one value travelling as loose parameters, so nothing keeps a
    /// caller from mixing them up or passing half of it. A struct with those
    /// fields does.
    ///
    /// Plain data is scalars, `&str`, shared slices, and structs and enums
    /// made only of those. `&mut` borrows, references and pointers to
    /// structs, trait objects, type parameters and owned buffers never count:
    /// they are the contexts, sinks and resources functions thread through by
    /// design, and bundling those is a different refactor. Also quiet on
    /// exported functions, trait methods and their impls, non-Rust ABIs,
    /// `#[no_mangle]` items and any function also used as a value (a fn
    /// pointer's signature is fixed by its type), on `self`, on `_`-prefixed
    /// parameters, on functions that merely declare the same pair without
    /// handing it on, and whenever the callee renames or retypes what it
    /// receives — that call is a translation, not a hand-off.
    ///
    /// Runs only with `parallel-params-enabled = true` in `dylint.toml`: a
    /// buffer and a cursor into it, or a level and the flags in force at it,
    /// travel together by design, and nothing in the signatures tells those
    /// from an undeclared struct.
    pub PARALLEL_PARAMS,
    Warn,
    "parameters that several functions declare and forward as a group"
}

struct Param {
    name: Symbol,
    /// The region-erased type, printed untrimmed: the identity two
    /// signatures are compared by.
    ty: String,
    /// The type as written in this signature, for the message.
    written: String,
    hir_id: HirId,
    span: Span,
}

struct FnFacts {
    /// By body param index; `None` for `self`, `_x` and destructuring
    /// patterns, so call arguments still line up with indices.
    params: Vec<Option<Param>>,
    span: Span,
}

/// One call in `from`'s body passing two or more of `from`'s own params
/// bare to `to`, as (from's param index, to's param index) pairs.
struct Forward {
    from: DefId,
    to: DefId,
    bound: Vec<(usize, usize)>,
    span: Span,
}

/// A parameter's identity across signatures: its name and erased type.
type Slot = (Symbol, String);

pub struct ParallelParams {
    min_fns: usize,
    fns: HashMap<DefId, FnFacts>,
    forwards: Vec<Forward>,
    /// Functions referenced other than by direct call: their signature is
    /// pinned by the fn-pointer type they were coerced to.
    poisoned: HashSet<DefId>,
}

rustc_session::impl_lint_pass!(ParallelParams => [PARALLEL_PARAMS]);

impl ParallelParams {
    pub fn new(config: &MordantConfig) -> Self {
        Self {
            min_fns: config.parallel_params_min_fns,
            fns: HashMap::new(),
            forwards: Vec::new(),
            poisoned: HashSet::new(),
        }
    }
}

/// The param local an argument passes on unchanged: `x`, `&x`, `&mut *x`.
fn forwarded_local(mut arg: &Expr<'_>) -> Option<HirId> {
    while let ExprKind::AddrOf(_, _, inner)
    | ExprKind::Unary(UnOp::Deref, inner)
    | ExprKind::DropTemps(inner) = arg.kind
    {
        arg = inner;
    }
    arg.res_local_id()
}

/// Plain data all the way down: scalars, `&str`, shared slices and arrays
/// of plain data, and structs and enums whose every field is. A reference or
/// pointer to a struct, a `&mut` borrow, a trait object, a type parameter or
/// an owned buffer is a context, sink or resource a function threads through
/// by design; a group of those is not a value nobody declared.
fn plain_data<'tcx>(cx: &LateContext<'tcx>, ty: Ty<'tcx>, seen: &mut HashSet<Ty<'tcx>>) -> bool {
    match ty.kind() {
        ty::Bool | ty::Char | ty::Int(_) | ty::Uint(_) | ty::Float(_) | ty::Str => true,
        ty::Array(elem, _) | ty::Slice(elem) => plain_data(cx, *elem, seen),
        ty::Tuple(tys) => tys.iter().all(|t| plain_data(cx, t, seen)),
        ty::Ref(_, inner, Mutability::Not) => {
            matches!(inner.kind(), ty::Slice(_) | ty::Str) && plain_data(cx, *inner, seen)
        }
        ty::Adt(adt, args) => {
            // A type reached again through its own fields adds no new field
            // kinds; the first visit decides.
            if !seen.insert(ty) {
                return true;
            }
            adt.all_fields()
                .all(|f| plain_data(cx, f.ty(cx.tcx, args).skip_normalization(), seen))
        }
        _ => false,
    }
}

impl<'tcx> LateLintPass<'tcx> for ParallelParams {
    fn check_fn(
        &mut self,
        cx: &LateContext<'tcx>,
        _kind: FnKind<'tcx>,
        decl: &'tcx FnDecl<'tcx>,
        body: &'tcx Body<'tcx>,
        _span: Span,
        def_id: LocalDefId,
    ) {
        if !owns_signature(cx, def_id) {
            return;
        }
        let params: Vec<Option<Param>> = body
            .params
            .iter()
            .enumerate()
            .map(|(i, param)| {
                let PatKind::Binding(_, hir_id, ident, None) = param.pat.kind else {
                    return None;
                };
                if ident.name == kw::SelfLower || ident.as_str().starts_with('_') {
                    return None;
                }
                let ty = cx
                    .tcx
                    .erase_and_anonymize_regions(cx.typeck_results().pat_ty(param.pat));
                if !plain_data(cx, ty, &mut HashSet::new()) {
                    return None;
                }
                let written = decl
                    .inputs
                    .get(i)
                    .and_then(|t| snippet_opt(cx, t.span))
                    .unwrap_or_else(|| with_no_trimmed_paths!(ty.to_string()));
                Some(Param {
                    name: ident.name,
                    ty: with_no_trimmed_paths!(ty.to_string()),
                    written,
                    hir_id,
                    span: param.span,
                })
            })
            .collect();
        if params.iter().flatten().count() < 2 {
            return;
        }
        let from = def_id.to_def_id();
        for_each_expr(cx, body.value, |e: &Expr<'tcx>| {
            if let Some(callee) = callee_of(cx, e)
                && let to = callee.def()
                && to != from
                && to.is_local()
                && matches!(cx.tcx.def_kind(to), DefKind::Fn | DefKind::AssocFn)
            {
                // A method's receiver is body param 0, so explicit args start
                // at 1 there; `T::m(x, ..)` through a path lines up as is.
                let args: Vec<(usize, &Expr<'tcx>)> = match callee {
                    Callee::Path { args, .. } => args.iter().enumerate().collect(),
                    Callee::Method { args, .. } => {
                        args.iter().enumerate().map(|(i, a)| (i + 1, a)).collect()
                    }
                };
                let mut bound: Vec<(usize, usize)> = Vec::new();
                for (to_idx, arg) in args {
                    if let Some(local) = forwarded_local(arg)
                        && let Some(from_idx) = params
                            .iter()
                            .position(|p| p.as_ref().is_some_and(|p| p.hir_id == local))
                        && !bound.iter().any(|(f, _)| *f == from_idx)
                    {
                        bound.push((from_idx, to_idx));
                    }
                }
                if bound.len() >= 2 {
                    self.forwards.push(Forward {
                        from,
                        to,
                        bound,
                        span: e.span,
                    });
                }
            }
            std::ops::ControlFlow::<()>::Continue(())
        });
        self.fns.insert(
            from,
            FnFacts {
                params,
                span: cx.tcx.def_span(from),
            },
        );
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        // A bare reference to a local fn (fn pointer, higher-order use). The
        // callee position of a direct call is not one.
        let ExprKind::Path(qpath) = &expr.kind else {
            return;
        };
        if matches!(
            clippy_utils::get_parent_expr(cx, expr),
            Some(Expr { kind: ExprKind::Call(callee, _), .. }) if callee.hir_id == expr.hir_id
        ) {
            return;
        }
        if let Res::Def(DefKind::Fn | DefKind::AssocFn, def) = cx.qpath_res(qpath, expr.hir_id)
            && def.is_local()
        {
            self.poisoned.insert(def);
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let facts = |def: &DefId| {
            if self.poisoned.contains(def) {
                None
            } else {
                self.fns.get(def)
            }
        };
        let slot = |p: &Param| -> Slot { (p.name, p.ty.clone()) };
        // Unordered slot pair -> the hand-offs that carry both.
        let mut links: HashMap<(Slot, Slot), Vec<&Forward>> = HashMap::new();
        for fwd in &self.forwards {
            let (Some(from), Some(to)) = (facts(&fwd.from), facts(&fwd.to)) else {
                continue;
            };
            // Forwarded params the callee declares under the same name and
            // type: those travel; anything renamed or converted does not.
            let kept: Vec<&Param> = fwd
                .bound
                .iter()
                .filter_map(|&(fi, ti)| {
                    let p = from.params.get(fi)?.as_ref()?;
                    let q = to.params.get(ti)?.as_ref()?;
                    (p.name == q.name && p.ty == q.ty).then_some(p)
                })
                .collect();
            for (i, a) in kept.iter().enumerate() {
                for b in &kept[i + 1..] {
                    let (sa, sb) = (slot(a), slot(b));
                    let key = if sa <= sb { (sa, sb) } else { (sb, sa) };
                    links.entry(key).or_default().push(fwd);
                }
            }
        }
        // The functions each linked pair passes between, in source order;
        // pairs sharing a function set are one group.
        let mut groups: HashMap<Vec<DefId>, (Vec<Slot>, Vec<&Forward>)> = HashMap::new();
        let mut keys: Vec<&(Slot, Slot)> = links.keys().collect();
        keys.sort();
        for key @ (a, b) in keys {
            let mut sig: Vec<DefId> = Vec::new();
            for def in links[key].iter().flat_map(|f| [f.from, f.to]) {
                if !sig.contains(&def) {
                    sig.push(def);
                }
            }
            if sig.len() < self.min_fns {
                continue;
            }
            sig.sort_by_key(|def| (self.fns[def].span.lo(), def.index));
            let (slots, witnesses) = groups.entry(sig).or_default();
            for s in [a, b] {
                if !slots.contains(s) {
                    slots.push(s.clone());
                }
            }
            witnesses.extend(links[key].iter().copied());
        }
        let mut findings: Vec<(Span, String, Span, String, String)> = Vec::new();
        for (sig, (mut slots, mut witnesses)) in groups {
            let anchor = &self.fns[&sig[0]];
            // Slots in the anchor's declaration order, printed as written there.
            let pos = |s: &Slot| {
                anchor
                    .params
                    .iter()
                    .position(|p| p.as_ref().is_some_and(|p| slot(p) == *s))
            };
            slots.sort_by_key(pos);
            let shown: Vec<String> = slots
                .iter()
                .filter_map(|s| {
                    let p = anchor.params[pos(s)?].as_ref()?;
                    Some(format!("`{}: {}`", p.name, p.written))
                })
                .collect();
            let Some(first) = slots.first().and_then(pos) else {
                continue;
            };
            let Some(at) = anchor.params[first].as_ref().map(|p| p.span) else {
                continue;
            };
            witnesses.sort_by_key(|f| f.span.lo());
            let witness = witnesses[0];
            // Two methods of different types often share a bare name; those
            // get their path so the list names each function once.
            let bare = |def: &DefId| cx.tcx.item_name(*def);
            let name = |def: &DefId| {
                if sig.iter().filter(|d| bare(d) == bare(def)).count() > 1 {
                    format!("`{}`", cx.tcx.def_path_str(*def))
                } else {
                    format!("`{}`", bare(def))
                }
            };
            let names: Vec<String> = sig.iter().map(name).collect();
            let fields: Vec<String> = slots.iter().map(|(n, _)| format!("`{n}`")).collect();
            let fns = listed(&names);
            findings.push((
                at,
                format!(
                    "{} are declared the same way by {fns} and passed between them unchanged. They are one value travelling as {} loose parameters",
                    join(&shown, "and"),
                    slots.len(),
                ),
                witness.span,
                format!(
                    "{} hands them to {} here",
                    name(&witness.from),
                    name(&witness.to)
                ),
                format!("put {} in a struct and pass that", join(&fields, "and")),
            ));
        }
        findings.sort_by_key(|(span, ..)| span.lo());
        for (span, msg, witness, note, help) in findings {
            emit_with_note(cx, PARALLEL_PARAMS, span, msg, witness, note, help);
        }
    }
}

/// Up to four names, then a count of the rest.
fn listed(names: &[String]) -> String {
    const SHOWN: usize = 4;
    if names.len() <= SHOWN {
        join(names, "and")
    } else {
        format!(
            "{} and {} more functions",
            names[..SHOWN].join(", "),
            names.len() - SHOWN
        )
    }
}

#[cfg(test)]
mod tests {
    use super::listed;
    use crate::baseline::join;

    fn owned(xs: &[&str]) -> Vec<String> {
        xs.iter().map(|s| (*s).to_owned()).collect()
    }

    #[test]
    fn join_reads_as_prose() {
        assert_eq!(join(&owned(&["a"]), "and"), "a");
        assert_eq!(join(&owned(&["a", "b"]), "and"), "a and b");
        assert_eq!(join(&owned(&["a", "b", "c"]), "or"), "a, b or c");
    }

    #[test]
    fn listed_caps_the_names_it_prints() {
        assert_eq!(listed(&owned(&["a", "b", "c", "d"])), "a, b, c and d");
        assert_eq!(
            listed(&owned(&["a", "b", "c", "d", "e", "f"])),
            "a, b, c, d and 2 more functions"
        );
    }
}

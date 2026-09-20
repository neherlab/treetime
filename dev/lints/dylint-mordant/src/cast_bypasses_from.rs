use std::collections::HashMap;

use crate::adt_facts::in_own_code_of;
use crate::baseline::{emit, emit_with_note};
use crate::ctor_flow;
use crate::hir_shapes::callee_of;
use rustc_hir::def_id::DefId;
use rustc_hir::{Expr, ExprKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::fast_reject::DeepRejectCtxt;
use rustc_middle::ty::{self, Ty, TypeVisitableExt};
use rustc_span::{Symbol, sym};

rustc_session::declare_lint! {
    /// Flags a `mem::transmute` / `mem::transmute_copy`, or a pointer cast
    /// (`p as *const T`, `p.cast::<T>()`) between different pointee types,
    /// that turns some other type into a struct, enum or union when a
    /// conversion from that same source type to that same target type
    /// already exists: an `impl From<A> for B`, an `impl TryFrom<A> for B`,
    /// or a safe receiver-less associated function of `B` taking one `A` (or
    /// `&A`) and returning `B`, `Option<B>` or `Result<B, _>`, by value or by
    /// reference. That function is where the crate says how an `A` becomes a
    /// `B`. The cast accepts any bit pattern, including the ones the
    /// conversion would reject or remap, and for an enum an unlisted
    /// discriminant is undefined behaviour on the spot.
    ///
    /// For a transmute of a value, every integer type counts as the same
    /// source: `transmute::<u16, E>(n as u16)` had to pick the repr width,
    /// and `E::from_raw(n: u32)` is still the conversion it skipped. For a
    /// pointer cast the pointee must be exactly the conversion's input type,
    /// so a byte buffer viewed as a header is not matched against
    /// `Header::new(u32)`, and a constructor returning a reference
    /// (`NameRef::new(&[u8]) -> Option<&NameRef>`) is the checked form of a
    /// pointer cast from that pointee only. A conversion for one
    /// instantiation of a generic target (`From<u32> for Id<Pkg>`) says
    /// nothing about another (`Id<Dep>`).
    ///
    /// Silent on: anything in the target type's own module or in any impl of
    /// it, trait impls included, since that code is the conversion or sits
    /// beside it; a transmute that only changes lifetimes or is otherwise
    /// between the same type; a target that is not an ADT (fn pointers,
    /// integers, type parameters); a pointer cast into a type with interior
    /// mutability, which views the pointee in place (`usize` as
    /// `AtomicUsize`) where a by-value `From` would make a new cell, unless
    /// a constructor returns that view itself; a type nothing converts
    /// into, which has no check to bypass; `unsafe fn` and unstable
    /// constructors, which promise no check or cannot be called; conversions
    /// whose input is the target type itself or generic; and a
    /// `mem::transmute` into a struct `unchecked_construction` reports, which
    /// names the check the value skipped.
    pub CAST_BYPASSES_FROM,
    Warn,
    "bits reinterpreted as a type that has a conversion from the same source, which the site skips"
}

/// One way the target type's crate (or a trait impl anywhere) turns `from`
/// into the target.
struct Conversion<'tcx> {
    /// The input as written, references included: what the message names.
    from: Ty<'tcx>,
    def: DefId,
    /// `Level::try_from`, `Code::from_raw`: how the message names it.
    name: String,
    /// Returns `Option`/`Result`: it can refuse a value, not just remap one.
    fallible: bool,
    /// Returns a reference to the target: a view of the input in place, the
    /// checked form of a pointer cast even into a type with interior
    /// mutability.
    views: bool,
}

/// How the site got from one type to the other; decides whether integer
/// widths are interchangeable when matching a conversion's input.
#[derive(Clone, Copy, PartialEq)]
enum Shape {
    Value,
    Pointer,
}

pub struct CastBypassesFrom {
    resource_errors: Vec<String>,
    /// Local struct -> whether `unchecked_construction` knows a validating
    /// constructor for it, so a `mem::transmute` into it is that lint's.
    validated: HashMap<DefId, bool>,
}

rustc_session::impl_lint_pass!(CastBypassesFrom => [CAST_BYPASSES_FROM]);

const TRANSMUTE: &str = "`mem::transmute`";
const TRANSMUTE_COPY: &str = "`mem::transmute_copy`";
const POINTER_CAST: &str = "a pointer cast";

/// `mem::transmute` / `mem::transmute_copy`, named the way the message shows
/// them. `transmute_copy` has no diagnostic item, hence the path test.
fn transmuter(cx: &LateContext<'_>, def: DefId) -> Option<&'static str> {
    if cx.tcx.is_diagnostic_item(sym::transmute, def) {
        return Some(TRANSMUTE);
    }
    (cx.tcx.crate_name(def.krate) == sym::core
        && cx.tcx.item_name(def).as_str() == "transmute_copy"
        && cx
            .tcx
            .opt_parent(def)
            .is_some_and(|m| cx.tcx.opt_item_name(m) == Some(sym::mem)))
    .then_some(TRANSMUTE_COPY)
}

/// Strips one reference or raw pointer layer off both types at once, as
/// long as both have one: `&A -> &B` and `*const A -> *mut B` compare their
/// pointees, `usize -> *const B` compares nothing.
fn peel_pointer_pair<'tcx>(mut a: Ty<'tcx>, mut b: Ty<'tcx>) -> (Ty<'tcx>, Ty<'tcx>, Shape) {
    let mut shape = Shape::Value;
    while let (Some(pa), Some(pb)) = (pointee(a), pointee(b)) {
        a = pa;
        b = pb;
        shape = Shape::Pointer;
    }
    (a, b, shape)
}

fn pointee(t: Ty<'_>) -> Option<Ty<'_>> {
    match *t.kind() {
        ty::Ref(_, inner, _) | ty::RawPtr(inner, _) => Some(inner),
        _ => None,
    }
}

impl<'tcx> Conversion<'tcx> {
    /// Whether this conversion is the checked form of a site turning `from`
    /// into the target by `shape`. A conversion taking `&A` still decides
    /// how an `A` becomes the target; one returning a reference is a view of
    /// its pointer input in place, which is what a pointer cast from that
    /// same pointee does unchecked and nothing a value transmute does; and a
    /// by-value conversion into a type with interior mutability makes a new
    /// cell where a pointer cast views the old one, so it stands in for
    /// value sites and for pointer casts into `Freeze` types only.
    fn checks(&self, from: Ty<'tcx>, shape: Shape, target_is_freeze: bool) -> bool {
        match (shape, self.views) {
            (Shape::Pointer, true) => pointee(self.from) == Some(from),
            (Shape::Value, true) => false,
            (Shape::Pointer, false) if !target_is_freeze => false,
            (_, false) => self.from == from || self.from.peel_refs() == from,
        }
    }
}

impl CastBypassesFrom {
    pub fn new(config: &crate::MordantConfig) -> Self {
        Self {
            resource_errors: config.validator_resource_errors.clone(),
            validated: HashMap::new(),
        }
    }

    /// `unchecked_construction`'s precondition for reporting a conjured value of
    /// `adt`: a local struct with a receiver-less inherent fn returning
    /// `Option`/`Result<Self>` whose body checks a field it stores.
    fn has_validator<'tcx>(&mut self, cx: &LateContext<'tcx>, adt: ty::AdtDef<'tcx>) -> bool {
        let did = adt.did();
        if !did.is_local() || !adt.is_struct() {
            return false;
        }
        if let Some(&known) = self.validated.get(&did) {
            return known;
        }
        let tcx = cx.tcx;
        let found = tcx.inherent_impls(did).iter().any(|&imp| {
            tcx.associated_items(imp)
                .in_definition_order()
                .filter(|item| item.is_fn() && !item.is_method())
                .any(|item| {
                    let output = tcx
                        .fn_sig(item.def_id)
                        .instantiate_identity()
                        .skip_normalization()
                        .skip_binder()
                        .output();
                    let ty::Adt(o, args) = *output.kind() else {
                        return false;
                    };
                    (tcx.is_diagnostic_item(sym::Option, o.did())
                        || tcx.is_diagnostic_item(sym::Result, o.did()))
                        && matches!(*args.type_at(0).kind(), ty::Adt(inner, _) if inner.did() == did)
                        && item.def_id.as_local().is_some_and(|ctor| {
                            !ctor_flow::checked_fields(cx, ctor, did, &self.resource_errors)
                                .is_empty()
                        })
                })
        });
        self.validated.insert(did, found);
        found
    }

    /// Reports `expr`, which turns a `from` into a `to` by reinterpretation
    /// (`how` names the means), when a conversion between the same pair
    /// exists and the site is not the target type's own code.
    fn check_reinterpretation<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        expr: &'tcx Expr<'tcx>,
        from: Ty<'tcx>,
        to: Ty<'tcx>,
        how: &'static str,
    ) {
        if expr.span.in_external_macro(cx.tcx.sess.source_map()) {
            return;
        }
        let tcx = cx.tcx;
        let (from, to, shape) = peel_pointer_pair(
            tcx.erase_and_anonymize_regions(from),
            tcx.erase_and_anonymize_regions(to),
        );
        if from == to {
            return;
        }
        let ty::Adt(adt, _) = *to.kind() else {
            return;
        };
        if from.ty_adt_def().is_some_and(|a| a.did() == adt.did())
            || in_own_code_of(cx, expr.hir_id, adt.did())
        {
            return;
        }
        let is_freeze = to.is_freeze(tcx, cx.typing_env());
        let conversions = collect_conversions(cx, adt, to);
        let exact = |c: &&Conversion<'tcx>| c.checks(from, shape, is_freeze);
        let same_integer = |c: &&Conversion<'tcx>| {
            shape == Shape::Value
                && !c.views
                && c.from.peel_refs().is_integral()
                && from.is_integral()
        };
        let Some(conv) = conversions
            .iter()
            .find(exact)
            .or_else(|| conversions.iter().find(same_integer))
        else {
            return;
        };
        if how == TRANSMUTE && shape == Shape::Value && self.has_validator(cx, adt) {
            return;
        }
        let (name, input) = (&conv.name, conv.from);
        let (existing, note, skipped) = if conv.fallible {
            (
                format!(
                    "`{name}` is the checked conversion from `{input}` and would reject some of them"
                ),
                format!("`{name}`, the checked conversion"),
                "the check",
            )
        } else {
            (
                format!("`{name}` is where the crate says how `{input}` becomes `{to}`"),
                format!("`{name}`, the conversion this skips"),
                "it",
            )
        };
        let msg = format!(
            "{how} turns `{from}` into `{to}` here and accepts any bit pattern. {existing}"
        );
        let help = format!(
            "call `{name}` instead. If this path must skip {skipped}, add an unchecked constructor next to it so both live beside the layout"
        );
        if conv.def.is_local() {
            emit_with_note(
                cx,
                CAST_BYPASSES_FROM,
                expr.span,
                msg,
                tcx.def_span(conv.def),
                note,
                help,
            );
        } else {
            emit(cx, CAST_BYPASSES_FROM, expr.span, msg, help);
        }
    }
}

/// Every `From`/`TryFrom` impl whose self type covers `to` and every safe,
/// stable, receiver-less inherent fn of the ADT that takes one value and
/// returns something covering `to`, bare, behind a reference, or either of
/// those in `Option`/`Result`. Inputs that are the ADT itself or mention a
/// type parameter convert nothing a reinterpreting site could have held.
fn collect_conversions<'tcx>(
    cx: &LateContext<'tcx>,
    adt: ty::AdtDef<'tcx>,
    to: Ty<'tcx>,
) -> Vec<Conversion<'tcx>> {
    let tcx = cx.tcx;
    let target = adt.did();
    let self_ty = tcx
        .type_of(target)
        .instantiate_identity()
        .skip_normalization();
    // `impl<T> From<u32> for Id<T>` converts into any `Id`; `for Id<Pkg>`
    // only into that one.
    let covers = |candidate: Ty<'tcx>| {
        DeepRejectCtxt::relate_rigid_infer(tcx)
            .types_may_unify(to, tcx.erase_and_anonymize_regions(candidate))
    };
    let is_target = |t: Ty<'tcx>| {
        t.peel_refs()
            .ty_adt_def()
            .is_some_and(|a| a.did() == target)
    };
    let usable_input = |t: Ty<'tcx>| !is_target(t) && !t.has_param();
    let callable = |def: DefId| !tcx.lookup_stability(def).is_some_and(|s| s.is_unstable());
    let name_of = |f: Symbol| format!("{}::{f}", tcx.item_name(target));
    let mut out = Vec::new();
    for (trait_sym, method, fallible) in
        [(sym::From, "from", false), (sym::TryFrom, "try_from", true)]
    {
        let Some(trait_did) = tcx.get_diagnostic_item(trait_sym) else {
            continue;
        };
        for imp in tcx.non_blanket_impls_for_ty(trait_did, self_ty) {
            let args = tcx
                .impl_trait_ref(imp)
                .instantiate_identity()
                .skip_normalization()
                .args;
            if !covers(args.type_at(0)) {
                continue;
            }
            let from = tcx.erase_and_anonymize_regions(args.type_at(1));
            if !usable_input(from) {
                continue;
            }
            let def = tcx
                .associated_items(imp)
                .in_definition_order()
                .find(|i| i.is_fn())
                .map_or(imp, |i| i.def_id);
            if !callable(imp) || !callable(def) {
                continue;
            }
            out.push(Conversion {
                from,
                def,
                name: name_of(Symbol::intern(method)),
                fallible,
                views: false,
            });
        }
    }
    for &imp in tcx.inherent_impls(target) {
        for item in tcx.associated_items(imp).in_definition_order() {
            if !item.is_fn() || item.is_method() || !callable(item.def_id) {
                continue;
            }
            let sig = tcx.instantiate_bound_regions_with_erased(
                tcx.fn_sig(item.def_id)
                    .instantiate_identity()
                    .skip_normalization(),
            );
            if sig.safety().is_unsafe() {
                continue;
            }
            let [input] = sig.inputs() else {
                continue;
            };
            let from = tcx.erase_and_anonymize_regions(*input);
            if !usable_input(from) {
                continue;
            }
            let output = sig.output();
            let (payload, fallible) = match *output.kind() {
                ty::Adt(o, args)
                    if o.did() != target
                        && (tcx.is_diagnostic_item(sym::Option, o.did())
                            || tcx.is_diagnostic_item(sym::Result, o.did())) =>
                {
                    (args.type_at(0), true)
                }
                _ => (output, false),
            };
            // `Option<&Self>` from a `&[u8]`: the checked view of those bytes.
            let (returned, views) = match *payload.kind() {
                ty::Ref(_, inner, _) => (inner, true),
                _ => (payload, false),
            };
            if !(matches!(*returned.kind(), ty::Adt(o, _) if o.did() == target) && covers(returned))
            {
                continue;
            }
            out.push(Conversion {
                from,
                def: item.def_id,
                name: name_of(item.name()),
                fallible,
                views,
            });
        }
    }
    out
}

impl<'tcx> LateLintPass<'tcx> for CastBypassesFrom {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let results = cx.typeck_results();
        match expr.kind {
            ExprKind::Call(callee, [_]) => {
                let Some(def) = callee_of(cx, expr).map(|c| c.def()) else {
                    return;
                };
                let Some(how) = transmuter(cx, def) else {
                    return;
                };
                // `transmute::<Src, Dst>` / `transmute_copy::<Src, Dst>`: the
                // types the intrinsic reinterprets between, whatever the
                // argument coerced from.
                let args = results.node_args(callee.hir_id);
                if args.len() < 2 {
                    return;
                }
                self.check_reinterpretation(cx, expr, args.type_at(0), args.type_at(1), how);
            }
            ExprKind::Cast(inner, _) if results.expr_ty(expr).is_raw_ptr() => {
                self.check_reinterpretation(
                    cx,
                    expr,
                    results.expr_ty_adjusted(inner),
                    results.expr_ty(expr),
                    POINTER_CAST,
                );
            }
            ExprKind::MethodCall(_, recv, [], _) => {
                let Some(def) = results.type_dependent_def_id(expr.hir_id) else {
                    return;
                };
                if !(cx.tcx.is_diagnostic_item(sym::ptr_cast, def)
                    || cx.tcx.is_diagnostic_item(sym::const_ptr_cast, def))
                {
                    return;
                }
                self.check_reinterpretation(
                    cx,
                    expr,
                    results.expr_ty_adjusted(recv),
                    results.expr_ty(expr),
                    POINTER_CAST,
                );
            }
            _ => {}
        }
    }
}

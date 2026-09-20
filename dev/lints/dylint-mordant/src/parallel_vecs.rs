use std::collections::{HashMap, HashSet};

use clippy_utils::res::MaybeResPath;
use clippy_utils::{SpanlessEq, get_parent_expr, hash_expr, is_default_equivalent};
use rustc_hir::def_id::{DefId, LocalDefId};
use rustc_hir::{
    BindingMode, ByRef, Expr, ExprKind, HirId, Mutability, Node, Pat, PatKind, QPath, UnOp,
};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::{self, Ty};
use rustc_span::hygiene::{ExpnKind, MacroKind};
use rustc_span::{DesugaringKind, Ident, Span, Symbol, sym};

use crate::adt_facts::field_ty;
use crate::baseline::{emit_with_note, join};
use crate::hir_shapes::{assigned_field, field_method_call, indexed_field};

rustc_session::declare_lint! {
    /// Flags two or more growable-sequence fields of one struct (`Vec`,
    /// `VecDeque`, or any type with `push`/`append`, `len` and indexing)
    /// whose lengths the crate only ever changes side by side -- every
    /// method call that takes one as `&mut`, reassignment, `&mut` borrow or
    /// `ref mut` destructuring of one sits in the same block as one of the
    /// other, on the same value, and every struct literal starts both empty
    /// -- and that some function reads at one index (`s.a[i]` with `s.b[i]`,
    /// `s.a.get(i)` with `s.b.get(i)`) or zips together. Element `i` of each
    /// is one record split across several vecs, and only the discipline of
    /// every writer keeps `a[i]` describing the same thing as `b[i]`. One
    /// `Vec` of a struct with a field for each has one length and one push.
    ///
    /// Only fields nothing outside the crate can write are considered: the
    /// struct is private to the crate, or the field is. One mutable access
    /// to either field without the other beside it disproves the pairing
    /// and the lint stays quiet, whatever the method is called, short of a
    /// few std methods that cannot change a length (`reserve`, `sort`,
    /// `iter_mut`, ..); so does a field built from anything but an empty
    /// constructor, a pair grown together but never read in step, and
    /// pushes to two different values of the type.
    pub PARALLEL_VECS,
    Warn,
    "sequence fields only ever grown together and read at one index"
}

/// Methods that take a sequence as `&mut` yet cannot change its length.
/// Anything else handed the field mutably may.
const KEEPS_LEN: &[&str] = &[
    "as_mut",
    "as_mut_ptr",
    "as_mut_slice",
    "as_mut_slices",
    "back_mut",
    "borrow_mut",
    "deref_mut",
    "fill",
    "fill_with",
    "first_mut",
    "front_mut",
    "get_mut",
    "get_unchecked_mut",
    "index_mut",
    "iter_mut",
    "last_mut",
    "make_contiguous",
    "range_mut",
    "reserve",
    "reserve_exact",
    "reverse",
    "rotate_left",
    "rotate_right",
    "shrink_to",
    "shrink_to_fit",
    "sort",
    "sort_by",
    "sort_by_key",
    "sort_unstable",
    "sort_unstable_by",
    "sort_unstable_by_key",
    "spare_capacity_mut",
    "swap",
    "try_reserve",
    "try_reserve_exact",
];

/// Constructors that yield an empty sequence whatever their arguments; a
/// bare `new()` or `default()` only with none.
const EMPTY_CTORS: &[&str] = &["new_in", "with_capacity", "with_capacity_in"];

/// Positional reads with one index argument.
const GET_OPS: &[&str] = &["get", "get_mut", "get_unchecked", "get_unchecked_mut"];

/// Adapters between a sequence and the `zip` that pairs it: `s.a.iter()`.
const ITER_ADAPTERS: &[&str] = &[
    "iter",
    "iter_mut",
    "into_iter",
    "copied",
    "cloned",
    "by_ref",
    "drain",
];

/// The value a sequence field belongs to.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Place {
    /// Read off an expression, `s` or `self.tape`: by shape, and by
    /// identity too when it is a local binding or a projection of one,
    /// since every local hashes alike.
    Expr { local: Option<HirId>, shape: u64 },
    /// Taken apart by this struct pattern, or one field's initialiser in a
    /// struct literal, which sets that length alone.
    Whole(HirId),
}

/// Where a length change happens: the body, the nearest block or match arm,
/// and the value it is applied to.
type Site = (DefId, HirId, Place);

/// An indexed read of a sequence field, kept until its body ends.
struct Read {
    owner: LocalDefId,
    adt: DefId,
    field: Symbol,
    place: Place,
    index: HirId,
}

#[derive(Default)]
pub struct ParallelVecs {
    /// Struct -> its sequence fields the crate alone can write; cached, and
    /// empty for structs that do not qualify.
    candidates: HashMap<DefId, Vec<Symbol>>,
    /// (struct, field) -> the sites that change its length.
    writes: HashMap<(DefId, Symbol), HashSet<Site>>,
    /// (struct, field, field), names ordered -> the first place the two are
    /// read at one index.
    in_step: HashMap<(DefId, Symbol, Symbol), Span>,
    reads: Vec<Read>,
}

rustc_session::impl_lint_pass!(ParallelVecs => [PARALLEL_VECS]);

fn has_inherent_method(cx: &LateContext<'_>, did: DefId, names: &[&str]) -> bool {
    cx.tcx.inherent_impls(did).iter().any(|imp| {
        names.iter().any(|n| {
            cx.tcx
                .associated_items(*imp)
                .filter_by_name_unhygienic(Symbol::intern(n))
                .next()
                .is_some()
        })
    })
}

/// `Vec`, `VecDeque`, or a type that grows (`push`/`append`), reports a
/// `len`, and is read by position (`Index` or `get`).
fn is_sequence<'tcx>(cx: &LateContext<'tcx>, ty: Ty<'tcx>) -> bool {
    let ty::Adt(adt, _) = ty.kind() else {
        return false;
    };
    let did = adt.did();
    if cx
        .tcx
        .get_diagnostic_name(did)
        .is_some_and(|n| matches!(n.as_str(), "Vec" | "VecDeque"))
    {
        return true;
    }
    let indexed = cx
        .tcx
        .lang_items()
        .index_trait()
        .is_some_and(|t| cx.tcx.non_blanket_impls_for_ty(t, ty).next().is_some())
        || has_inherent_method(cx, did, &["get", "at"]);
    indexed
        && has_inherent_method(cx, did, &["push", "push_back", "append"])
        && has_inherent_method(cx, did, &["len"])
}

impl ParallelVecs {
    /// The sequence fields of the local struct behind `ty` that only this
    /// crate can write, when there are at least two of them.
    fn candidates<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        ty: Ty<'tcx>,
    ) -> Option<(DefId, &[Symbol])> {
        let ty::Adt(adt, _) = ty.peel_refs().kind() else {
            return None;
        };
        if !adt.is_struct() || !adt.did().is_local() {
            return None;
        }
        let did = adt.did();
        let fields = self.candidates.entry(did).or_insert_with(|| {
            let exported = cx.effective_visibilities.is_exported(did.expect_local());
            let fields: Vec<Symbol> = adt
                .non_enum_variant()
                .fields
                .iter()
                .filter(|f| (!exported || !f.vis.is_public()) && is_sequence(cx, field_ty(cx, f)))
                .map(|f| f.name)
                .collect();
            if fields.len() >= 2 {
                fields
            } else {
                Vec::new()
            }
        });
        (!fields.is_empty()).then_some((did, fields.as_slice()))
    }

    /// `base.field` as a candidate sequence field: its struct, its name, and
    /// the value it is read off.
    fn sequence_field<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        base: &'tcx Expr<'tcx>,
        field: Ident,
    ) -> Option<(DefId, Symbol, Place)> {
        let (adt, fields) = self.candidates(cx, cx.typeck_results().expr_ty_adjusted(base))?;
        if !fields.contains(&field.name) {
            return None;
        }
        let place = Place::Expr {
            local: local_root(base),
            shape: hash_expr(cx, base),
        };
        Some((adt, field.name, place))
    }

    /// A length change of `field` at `at`, on `place`.
    fn note_write(
        &mut self,
        cx: &LateContext<'_>,
        adt: DefId,
        field: Symbol,
        at: HirId,
        place: Place,
    ) {
        let body = cx.tcx.hir_enclosing_body_owner(at).to_def_id();
        self.writes
            .entry((adt, field))
            .or_default()
            .insert((body, step_scope(cx, at), place));
    }

    fn record_write<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        at: &'tcx Expr<'tcx>,
        base: &'tcx Expr<'tcx>,
        field: Ident,
    ) {
        if let Some((adt, field, place)) = self.sequence_field(cx, base, field) {
            self.note_write(cx, adt, field, at.hir_id, place);
        }
    }

    /// A struct literal sets each sequence field's length on its own: every
    /// field not built empty is written there, alone.
    fn record_literal<'tcx>(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let ExprKind::Struct(_, inits, _) = expr.kind else {
            return;
        };
        // A derived `Clone` copies both lengths off one value: it keeps
        // whatever the other sites establish and proves nothing itself.
        if let ExpnKind::Macro(MacroKind::Derive, name) = expr.span.ctxt().outer_expn_data().kind
            && name == sym::Clone
        {
            return;
        }
        let Some((adt, fields)) = self.candidates(cx, cx.typeck_results().expr_ty(expr)) else {
            return;
        };
        let written: Vec<(Symbol, HirId)> = inits
            .iter()
            .filter(|i| fields.contains(&i.ident.name) && !is_empty_ctor(cx, i.expr))
            .map(|i| (i.ident.name, i.expr.hir_id))
            .collect();
        for (field, init) in written {
            self.note_write(cx, adt, field, expr.hir_id, Place::Whole(init));
        }
    }

    fn record_index<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        at: &'tcx Expr<'tcx>,
        base: &'tcx Expr<'tcx>,
        field: Ident,
        index: &'tcx Expr<'tcx>,
    ) {
        let Some((adt, field, place)) = self.sequence_field(cx, base, field) else {
            return;
        };
        let owner = cx.tcx.hir_enclosing_body_owner(at.hir_id);
        for earlier in &self.reads {
            if earlier.owner != owner
                || earlier.adt != adt
                || earlier.field == field
                || earlier.place != place
            {
                continue;
            }
            let other = cx.tcx.hir_expect_expr(earlier.index);
            if SpanlessEq::new(cx).eq_expr(at.span.ctxt(), other, index) {
                self.in_step
                    .entry(ordered(adt, earlier.field, field))
                    .or_insert(at.span);
            }
        }
        self.reads.push(Read {
            owner,
            adt,
            field,
            place,
            index: index.hir_id,
        });
    }

    /// `zip` over two sequence fields of one value reads them in step.
    fn record_zip<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        at: Span,
        l: &'tcx Expr<'tcx>,
        r: &'tcx Expr<'tcx>,
    ) {
        let (Some((lb, lf)), Some((rb, rf))) = (zipped_field(l), zipped_field(r)) else {
            return;
        };
        let (Some((adt, lf, lp)), Some((radt, rf, rp))) = (
            self.sequence_field(cx, lb, lf),
            self.sequence_field(cx, rb, rf),
        ) else {
            return;
        };
        if adt == radt && lf != rf && lp == rp {
            self.in_step.entry(ordered(adt, lf, rf)).or_insert(at);
        }
    }
}

/// The local binding a place expression projects from, through fields,
/// indexing, `&` and `*`.
fn local_root(mut e: &Expr<'_>) -> Option<HirId> {
    loop {
        match e.kind {
            ExprKind::Field(inner, _)
            | ExprKind::Index(inner, _, _)
            | ExprKind::AddrOf(_, _, inner)
            | ExprKind::Unary(UnOp::Deref, inner)
            | ExprKind::DropTemps(inner) => e = inner,
            _ => return e.res_local_id(),
        }
    }
}

fn ordered(adt: DefId, a: Symbol, b: Symbol) -> (DefId, Symbol, Symbol) {
    if a.as_str() <= b.as_str() {
        (adt, a, b)
    } else {
        (adt, b, a)
    }
}

/// The nearest block or match arm around `hir_id`: two length changes are
/// side by side when they share it.
fn step_scope(cx: &LateContext<'_>, hir_id: HirId) -> HirId {
    for (id, node) in cx.tcx.hir_parent_iter(hir_id) {
        match node {
            Node::Block(_) | Node::Arm(_) => return id,
            Node::Item(_) | Node::ImplItem(_) | Node::TraitItem(_) => break,
            _ => {}
        }
    }
    hir_id
}

/// The `base.field` a `zip` operand iterates: `&s.a`, `s.a.iter()`,
/// `s.a.iter().copied()`.
fn zipped_field<'h>(mut e: &'h Expr<'h>) -> Option<(&'h Expr<'h>, Ident)> {
    loop {
        match e.kind {
            ExprKind::AddrOf(_, _, inner) | ExprKind::DropTemps(inner) => e = inner,
            ExprKind::MethodCall(seg, recv, _, _)
                if ITER_ADAPTERS.contains(&seg.ident.name.as_str()) =>
            {
                e = recv;
            }
            ExprKind::Field(base, ident) => return Some((base, ident)),
            _ => return None,
        }
    }
}

/// `&mut s.a` in a `for` head iterates; anywhere else it is handed to code
/// that may change the length.
fn is_for_loop_head(cx: &LateContext<'_>, e: &Expr<'_>) -> bool {
    get_parent_expr(cx, e).is_some_and(|p| {
        matches!(p.kind, ExprKind::Call(..)) && p.span.is_desugaring(DesugaringKind::ForLoop)
    })
}

/// The method call auto-borrows `recv` as `&mut` (or `*mut`) to its own
/// type: the callee is handed the sequence itself mutably, not the slice
/// it derefs to, and so may change its length.
fn borrows_receiver_mut(cx: &LateContext<'_>, recv: &Expr<'_>) -> bool {
    let typeck = cx.typeck_results();
    let ty = typeck.expr_ty(recv);
    match *typeck.expr_ty_adjusted(recv).kind() {
        ty::Ref(_, inner, Mutability::Mut) | ty::RawPtr(inner, Mutability::Mut) => inner == ty,
        _ => false,
    }
}

/// An initialiser that yields an empty sequence: whatever `Default` gives,
/// or a constructor call named for building one.
fn is_empty_ctor(cx: &LateContext<'_>, e: &Expr<'_>) -> bool {
    if is_default_equivalent(cx, e) {
        return true;
    }
    let ExprKind::Call(callee, args) = e.kind else {
        return false;
    };
    let name = match callee.kind {
        ExprKind::Path(QPath::TypeRelative(_, seg)) => seg.ident.name,
        ExprKind::Path(QPath::Resolved(_, path)) => match path.segments.last() {
            Some(seg) => seg.ident.name,
            None => return false,
        },
        _ => return false,
    };
    let name = name.as_str();
    EMPTY_CTORS.contains(&name) || (args.is_empty() && matches!(name, "new" | "default"))
}

impl<'tcx> LateLintPass<'tcx> for ParallelVecs {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        match expr.kind {
            ExprKind::MethodCall(seg, recv, [arg], _) if seg.ident.name.as_str() == "zip" => {
                self.record_zip(cx, expr.span, recv, arg);
            }
            ExprKind::Call(callee, [l, r])
                if matches!(callee.kind, ExprKind::Path(ref qp)
                    if cx.qpath_res(qp, callee.hir_id).opt_def_id()
                        .is_some_and(|d| cx.tcx.opt_item_name(d)
                            .is_some_and(|n| n.as_str() == "zip"))) =>
            {
                self.record_zip(cx, expr.span, l, r);
            }
            ExprKind::MethodCall(_, recv, ..) => {
                let Some(call) = field_method_call(expr) else {
                    return;
                };
                let name = call.method.name.as_str();
                if GET_OPS.contains(&name)
                    && let [index] = call.args
                {
                    self.record_index(cx, expr, call.base, call.field, index);
                } else if !KEEPS_LEN.contains(&name) && borrows_receiver_mut(cx, recv) {
                    self.record_write(cx, expr, call.base, call.field);
                }
            }
            ExprKind::Index(..) => {
                if let Some(read) = indexed_field(expr) {
                    self.record_index(cx, expr, read.base, read.field, read.index);
                }
            }
            ExprKind::Assign(place, _, _) | ExprKind::AssignOp(_, place, _) => {
                if let Some((base, field, _)) = assigned_field(place) {
                    self.record_write(cx, expr, base, field);
                }
            }
            ExprKind::AddrOf(_, Mutability::Mut, inner) => {
                if let Some((base, field, _)) = assigned_field(inner)
                    && !is_for_loop_head(cx, expr)
                {
                    self.record_write(cx, expr, base, field);
                }
            }
            ExprKind::Struct(..) => self.record_literal(cx, expr),
            _ => {}
        }
    }

    // `let S { a, .. } = self` on `&mut self`, or `S { ref mut a, .. }`: the
    // binding is a `&mut` to the field, through which its length changes.
    fn check_pat(&mut self, cx: &LateContext<'tcx>, pat: &'tcx Pat<'tcx>) {
        let PatKind::Struct(_, bindings, _) = pat.kind else {
            return;
        };
        let Some(typeck) = cx.maybe_typeck_results() else {
            return;
        };
        let Some((adt, fields)) = self.candidates(cx, typeck.pat_ty(pat)) else {
            return;
        };
        let written: Vec<Symbol> = bindings
            .iter()
            .filter(|b| fields.contains(&b.ident.name))
            .filter(|b| {
                let mut by_mut_ref = false;
                b.pat.walk(|p| {
                    if let PatKind::Binding(..) = p.kind
                        && let Some(BindingMode(ByRef::Yes(_, Mutability::Mut), _)) =
                            typeck.pat_binding_modes().get(p.hir_id)
                    {
                        by_mut_ref = true;
                    }
                    !by_mut_ref
                });
                by_mut_ref
            })
            .map(|b| b.ident.name)
            .collect();
        for field in written {
            self.note_write(cx, adt, field, pat.hir_id, Place::Whole(pat.hir_id));
        }
    }

    fn check_body_post(&mut self, cx: &LateContext<'tcx>, body: &rustc_hir::Body<'tcx>) {
        let owner = cx.tcx.hir_body_owner_def_id(body.id());
        self.reads.retain(|r| r.owner != owner);
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let mut per_struct: HashMap<DefId, Vec<(Symbol, HashSet<Site>)>> = HashMap::new();
        for ((adt, field), sites) in self.writes.drain() {
            per_struct.entry(adt).or_default().push((field, sites));
        }
        let mut findings: Vec<(DefId, Vec<Symbol>, usize, Span)> = Vec::new();
        for (adt, mut fields) in per_struct {
            fields.sort_by_key(|(f, _)| f.as_str().to_owned());
            let mut groups: Vec<(&HashSet<Site>, Vec<Symbol>)> = Vec::new();
            for (field, sites) in &fields {
                match groups.iter_mut().find(|(s, _)| *s == sites) {
                    Some((_, members)) => members.push(*field),
                    None => groups.push((sites, vec![*field])),
                }
            }
            for (sites, members) in groups {
                let read_in_step = members
                    .iter()
                    .flat_map(|a| {
                        members
                            .iter()
                            .filter_map(|b| self.in_step.get(&ordered(adt, *a, *b)).copied())
                    })
                    .min_by_key(|span| span.lo());
                if members.len() >= 2
                    && let Some(read) = read_in_step
                {
                    findings.push((adt, members, sites.len(), read));
                }
            }
        }
        findings.sort_by_key(|(adt, ..)| cx.tcx.def_span(*adt).lo());
        for (adt, members, sites, read) in findings {
            let names: Vec<String> = members.iter().map(|m| format!("`{m}`")).collect();
            let names = join(&names, "and");
            let n = members.len();
            emit_with_note(
                cx,
                PARALLEL_VECS,
                cx.tcx.def_span(adt),
                format!(
                    "{names} of `{}` only change length together ({sites} {}) and are read at the same index. Element `i` of each is one record split across {n} vecs",
                    cx.tcx.def_path_str(adt),
                    if sites == 1 { "block" } else { "blocks" },
                ),
                read,
                "one of the reads at a shared index",
                "use one `Vec` of a struct with a field for each. One length, one push",
            );
        }
    }
}

use std::collections::{HashMap, HashSet};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::ops::ControlFlow;

use crate::adt_facts::{field_ty, has_fixed_repr, is_option_ty, struct_field};
use crate::baseline::emit_with_note;
use crate::hir_shapes::{assigned_field, field_chain, peel_blocks_unsafe};
use clippy_utils::visitors::for_each_expr_without_closures;
use clippy_utils::{as_some_expr, get_parent_expr, hash_expr, is_default_equivalent, is_none_expr};
use rustc_ast::LitKind;
use rustc_hir::def::Res;
use rustc_hir::def_id::DefId;
use rustc_hir::{
    BindingMode, BorrowKind, ByRef, Expr, ExprKind, HirId, Mutability, Node, Pat, PatKind, QPath,
    Stmt,
};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::adjustment::{Adjust, AutoBorrow, AutoBorrowMutability};
use rustc_middle::ty::{self, AdtDef};
use rustc_span::hygiene::{ExpnKind, MacroKind};
use rustc_span::{Span, Symbol, sym};

rustc_session::declare_lint! {
    /// Flags a bool field that is only ever written next to a sibling
    /// `Option` field, `true` with `Some` and `false` with `None` (or always
    /// the reverse): every write to either, struct literal or field
    /// assignment, sits beside a write to the other in the same struct
    /// expression or straight-line block. The flag is that field's
    /// `is_some()` stored twice.
    ///
    /// Only fires on named fields nothing outside the crate can name, when
    /// every write to the pair is a literal `true`/`false` and
    /// `Some(..)`/`None` (or `Default::default()`, which is `false`/`None`)
    /// and both polarities occur. A lone write to either field, one in a
    /// match arm or closure without a block of its own, a block that writes
    /// both polarities, a computed value, a compound assignment, or a `&mut`
    /// to either field -- `.take()`, `mem::replace`, `&mut s.opt`, a
    /// `ref mut` pattern binding, spelled out or through a `&mut` scrutinee
    /// -- is unprovable and silences it; `as_mut` and `as_deref_mut` cannot
    /// change the discriminant and do not. Tuple structs are built through
    /// their constructor fn and are not followed; a derived `Clone` copies
    /// the pair as it stands and is not a write of its own.
    pub BOOL_BESIDE_OPTION,
    Warn,
    "bool field that repeats whether a sibling Option field is Some"
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Bool,
    Opt,
}

/// Where a write happens: the struct expression, or for a field assignment
/// the innermost block around it, the stretch of that block no early exit
/// splits, and the place it assigns through, so `a.opt = ..` and
/// `b.flag = ..` in one block, or two writes a `?` can part, are not beside
/// each other.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct Site {
    region: HirId,
    stretch: usize,
    place: u64,
}

struct FieldWrites {
    kind: Kind,
    /// Some write to the field is not a literal of its kind, one site writes
    /// both polarities, or the field is mutably borrowed: its value is no
    /// longer a function of the sites.
    unprovable: bool,
    /// (site, set): where the field is written, and whether it is written
    /// `true` / `Some(..)` there.
    sites: HashSet<(Site, bool)>,
    /// The first write recorded, for the note.
    first: Option<Span>,
}

#[derive(Default)]
pub struct BoolBesideOption {
    writes: HashMap<(DefId, Symbol), FieldWrites>,
}

rustc_session::impl_lint_pass!(BoolBesideOption => [BOOL_BESIDE_OPTION]);

/// A bool or `Option` named field of a local struct that nothing outside the
/// crate can write: the struct with the field's kind, or None for anything
/// else.
fn tracked_field<'tcx>(
    cx: &LateContext<'tcx>,
    ty: ty::Ty<'tcx>,
    name: Symbol,
) -> Option<(AdtDef<'tcx>, Kind)> {
    let ty::Adt(adt, _) = ty.peel_refs().kind() else {
        return None;
    };
    if !adt.is_struct()
        || !adt.did().is_local()
        || has_fixed_repr(*adt)
        || adt.non_enum_variant().ctor.is_some()
    {
        return None;
    }
    let f = struct_field(*adt, name)?;
    if cx.effective_visibilities.is_exported(f.did.expect_local()) {
        return None;
    }
    let fty = field_ty(cx, f);
    let kind = if fty.is_bool() {
        Kind::Bool
    } else if is_option_ty(cx, fty) {
        Kind::Opt
    } else {
        return None;
    };
    Some((*adt, kind))
}

/// Whether `value` writes `true`/`Some(..)` (Some(true)) or `false`/`None`
/// (Some(false)); None for anything computed. `Default::default()` of either
/// kind is `false`/`None`, which is what a derived `Default` writes.
fn polarity(cx: &LateContext<'_>, kind: Kind, value: &Expr<'_>) -> Option<bool> {
    let value = peel_blocks_unsafe(value);
    let set = match kind {
        Kind::Bool => match value.kind {
            ExprKind::Lit(lit) => match lit.node {
                LitKind::Bool(b) => Some(b),
                _ => None,
            },
            _ => None,
        },
        Kind::Opt => {
            if as_some_expr(cx, value).is_some() {
                Some(true)
            } else {
                is_none_expr(cx, value).then_some(false)
            }
        }
    };
    set.or_else(|| is_default_equivalent(cx, value).then_some(false))
}

/// The straight-line region a field assignment runs in: the innermost block
/// around it, and which stretch of that block, counting the statements
/// before the assignment that can leave early (`return`, `?`, `break`,
/// `continue`), since a write after one of those need not follow a write
/// before it. None when a match arm, a closure or the body itself holds the
/// assignment bare, where nothing else runs beside it.
fn assignment_region<'tcx>(cx: &LateContext<'tcx>, id: HirId) -> Option<(HirId, usize)> {
    let mut from = id;
    for (parent_id, node) in cx.tcx.hir_parent_iter(id) {
        match node {
            Node::Block(block) => {
                let leaves_early = |s: &'tcx Stmt<'tcx>| {
                    // Loops opened inside the statement: a `break` or
                    // `continue` aimed at one of them stays inside it.
                    let mut own_loops = HashSet::new();
                    for_each_expr_without_closures(s, |e: &'tcx Expr<'tcx>| {
                        match e.kind {
                            ExprKind::Ret(_) => return ControlFlow::Break(()),
                            ExprKind::Loop(..) => {
                                own_loops.insert(e.hir_id);
                            }
                            ExprKind::Break(dest, _) | ExprKind::Continue(dest)
                                if !dest.target_id.is_ok_and(|l| own_loops.contains(&l)) =>
                            {
                                return ControlFlow::Break(());
                            }
                            _ => {}
                        }
                        ControlFlow::Continue(())
                    })
                    .is_some()
                };
                let stretch = block
                    .stmts
                    .iter()
                    .take_while(|s| s.hir_id != from)
                    .filter(|&s| leaves_early(s))
                    .count();
                return Some((parent_id, stretch));
            }
            Node::Expr(e) if matches!(e.kind, ExprKind::Closure(_)) => return None,
            Node::Expr(_) | Node::ExprField(_) | Node::Stmt(_) | Node::LetStmt(_) => {
                from = parent_id;
            }
            _ => return None,
        }
    }
    None
}

/// The place an assignment writes through, `base` of `base.field = ..`: the
/// local at its root and the fields walked from it. `SpanlessHash` hashes
/// every local alike, which would put `a.f` beside `b.g`.
fn place_key(cx: &LateContext<'_>, base: &Expr<'_>) -> u64 {
    let chain = field_chain(base);
    let mut h = DefaultHasher::new();
    if let ExprKind::Path(QPath::Resolved(None, path)) = chain.root.kind
        && let Res::Local(local) = path.res
    {
        local.hash(&mut h);
    } else {
        hash_expr(cx, chain.root).hash(&mut h);
    }
    chain.fields.hash(&mut h);
    h.finish()
}

impl BoolBesideOption {
    fn facts(&mut self, adt: AdtDef<'_>, name: Symbol, kind: Kind) -> &mut FieldWrites {
        self.writes
            .entry((adt.did(), name))
            .or_insert_with(|| FieldWrites {
                kind,
                unprovable: false,
                sites: HashSet::new(),
                first: None,
            })
    }

    fn record<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        base_ty: ty::Ty<'tcx>,
        name: Symbol,
        site: Site,
        at: Span,
        value: &Expr<'_>,
    ) {
        let Some((adt, kind)) = tracked_field(cx, base_ty, name) else {
            return;
        };
        let facts = self.facts(adt, name, kind);
        match polarity(cx, kind, value) {
            // A region that writes both polarities no longer says which of
            // its sibling writes each one sits beside.
            Some(set) if facts.sites.contains(&(site, !set)) => facts.unprovable = true,
            Some(set) => {
                facts.sites.insert((site, set));
                // A write the user spelled, over one a derive expanded to.
                if facts
                    .first
                    .is_none_or(|f| f.from_expansion() && !at.from_expansion())
                {
                    facts.first = Some(at);
                }
            }
            None => facts.unprovable = true,
        }
    }

    fn poison<'tcx>(&mut self, cx: &LateContext<'tcx>, base_ty: ty::Ty<'tcx>, name: Symbol) {
        if let Some((adt, kind)) = tracked_field(cx, base_ty, name) {
            self.facts(adt, name, kind).unprovable = true;
        }
    }
}

/// Methods that borrow an `Option` or bool field mutably without being able
/// to replace it whole.
const PROJECTING_MUT_METHODS: &[&str] = &["as_mut", "as_deref_mut", "as_pin_mut", "iter_mut"];

impl<'tcx> LateLintPass<'tcx> for BoolBesideOption {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        match expr.kind {
            ExprKind::Struct(_, fields, _) => {
                // A derived `Clone` writes every field from the same field of
                // the value it copies: it keeps whatever the other sites
                // establish and proves nothing itself. Every other derive
                // builds a value of its own and counts like a hand-written one.
                if let ExpnKind::Macro(MacroKind::Derive, name) =
                    expr.span.ctxt().outer_expn_data().kind
                    && name == sym::Clone
                {
                    return;
                }
                let ty = cx.typeck_results().expr_ty(expr);
                let site = Site {
                    region: expr.hir_id,
                    stretch: 0,
                    place: 0,
                };
                for field in fields {
                    self.record(cx, ty, field.ident.name, site, field.span, field.expr);
                }
            }
            ExprKind::Assign(place, value, _) => {
                let Some((base, ident, _)) = assigned_field(place) else {
                    return;
                };
                let base_ty = cx.typeck_results().expr_ty_adjusted(base);
                if tracked_field(cx, base_ty, ident.name).is_none() {
                    return;
                }
                match assignment_region(cx, expr.hir_id) {
                    Some((region, stretch)) => {
                        let site = Site {
                            region,
                            stretch,
                            place: place_key(cx, base),
                        };
                        self.record(cx, base_ty, ident.name, site, expr.span, value);
                    }
                    None => self.poison(cx, base_ty, ident.name),
                }
            }
            ExprKind::AssignOp(_, place, _) => {
                if let Some((base, ident, _)) = assigned_field(place) {
                    self.poison(cx, cx.typeck_results().expr_ty_adjusted(base), ident.name);
                }
            }
            ExprKind::AddrOf(BorrowKind::Ref | BorrowKind::Raw, Mutability::Mut, inner) => {
                if let ExprKind::Field(base, ident) = peel_blocks_unsafe(inner).kind {
                    self.poison(cx, cx.typeck_results().expr_ty_adjusted(base), ident.name);
                }
            }
            // The auto-`&mut` a mutating method call takes on its receiver.
            ExprKind::Field(base, ident) => {
                let mutably_borrowed = cx.typeck_results().expr_adjustments(expr).iter().any(|a| {
                    matches!(
                        a.kind,
                        Adjust::Borrow(AutoBorrow::Ref(AutoBorrowMutability::Mut { .. }))
                            | Adjust::Borrow(AutoBorrow::RawPtr(Mutability::Mut))
                    )
                });
                if !mutably_borrowed {
                    return;
                }
                let projecting = matches!(
                    get_parent_expr(cx, expr).map(|p| &p.kind),
                    Some(ExprKind::MethodCall(seg, ..))
                        if PROJECTING_MUT_METHODS.contains(&seg.ident.name.as_str())
                );
                if !projecting {
                    self.poison(cx, cx.typeck_results().expr_ty_adjusted(base), ident.name);
                }
            }
            _ => {}
        }
    }

    // `let S { opt, .. } = self` on `&mut self`, or `S { ref mut opt, .. }`:
    // the binding is a `&mut` to the field, through which it is rewritten
    // alone.
    fn check_pat(&mut self, cx: &LateContext<'tcx>, pat: &'tcx Pat<'tcx>) {
        let PatKind::Struct(_, fields, _) = pat.kind else {
            return;
        };
        let Some(typeck) = cx.maybe_typeck_results() else {
            return;
        };
        let ty = typeck.pat_ty(pat);
        for field in fields {
            let mut by_mut_ref = false;
            field.pat.walk(|p| {
                if let PatKind::Binding(..) = p.kind
                    && let Some(BindingMode(ByRef::Yes(_, Mutability::Mut), _)) =
                        typeck.pat_binding_modes().get(p.hir_id)
                {
                    by_mut_ref = true;
                }
                !by_mut_ref
            });
            if by_mut_ref {
                self.poison(cx, ty, field.ident.name);
            }
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let mut per_struct: HashMap<DefId, Vec<(Symbol, FieldWrites)>> = HashMap::new();
        for ((did, name), facts) in self.writes.drain() {
            if !facts.unprovable && !facts.sites.is_empty() {
                per_struct.entry(did).or_default().push((name, facts));
            }
        }
        let mut findings: Vec<(Span, String, Span, String)> = Vec::new();
        for (did, mut fields) in per_struct {
            fields.sort_by_key(|(name, _)| name.as_str().to_owned());
            let adt = cx.tcx.adt_def(did);
            for (flag, flag_facts) in fields.iter().filter(|(_, f)| f.kind == Kind::Bool) {
                // A flag only ever `false` (or only ever `true`) is a constant,
                // not a copy of anything.
                if !flag_facts.sites.iter().any(|&(_, set)| set)
                    || !flag_facts.sites.iter().any(|&(_, set)| !set)
                {
                    continue;
                }
                let Some(first) = flag_facts.first else {
                    continue;
                };
                for (opt, opt_facts) in fields.iter().filter(|(_, f)| f.kind == Kind::Opt) {
                    let inverted: HashSet<(Site, bool)> =
                        opt_facts.sites.iter().map(|&(s, set)| (s, !set)).collect();
                    let reading = if flag_facts.sites == opt_facts.sites {
                        "is_some"
                    } else if flag_facts.sites == inverted {
                        "is_none"
                    } else {
                        continue;
                    };
                    let (with_true, with_false) = match reading {
                        "is_some" => ("Some", "None"),
                        _ => ("None", "Some"),
                    };
                    let Some(field) = struct_field(adt, *flag) else {
                        continue;
                    };
                    findings.push((
                        cx.tcx.def_span(field.did),
                        format!(
                            "`{}.{flag}` is only ever written next to `{opt}`, `true` with `{with_true}` and `false` with `{with_false}`, in all {} places. It is `{opt}.{reading}()` stored twice",
                            cx.tcx.def_path_str(did),
                            flag_facts.sites.len(),
                        ),
                        first,
                        format!(
                            "delete `{flag}` and call `{opt}.{reading}()` where it was read"
                        ),
                    ));
                    break;
                }
            }
        }
        findings.sort_by_key(|(span, ..)| span.lo());
        for (span, msg, first, help) in findings {
            emit_with_note(
                cx,
                BOOL_BESIDE_OPTION,
                span,
                msg,
                first,
                "one of the paired writes",
                help,
            );
        }
    }
}

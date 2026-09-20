use std::collections::HashMap;

use crate::adt_facts::impl_self_adt;
use crate::baseline::emit_with_note;
use crate::ctor_flow::{self, FieldCheck};
use crate::hir_shapes::{assigned_adt_field, callee_of};
use clippy_utils::ty::ty_from_hir_ty;
use rustc_abi::FieldIdx;
use rustc_hir::def::{CtorKind, CtorOf, DefKind};
use rustc_hir::def_id::DefId;
use rustc_hir::{Expr, ExprKind, FnRetTy, HirId, ImplItem, ImplItemKind, ItemKind, StructTailExpr};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::{Span, Symbol, sym};

rustc_session::declare_lint! {
    /// Flags a value of a validated type built or changed by hand, skipping
    /// the check its validating constructor makes. A validating constructor
    /// is a receiver-less inherent function returning `Result<Self, _>` or
    /// `Option<Self>` whose body rejects some value it then stores in a field
    /// (see `ctor_flow`). Outside the type's own module and impls, each of
    /// these skips that check:
    ///
    /// * a struct literal, `S(x)` on a tuple struct included -- one with a
    ///   `..base` tail only when it names a field some validator checks,
    ///   since the fields it takes from `base` were checked when `base` was
    ///   made,
    /// * an assignment, plain or compound, to a checked field,
    /// * `mem::zeroed`, `mem::transmute` or `MaybeUninit::assume_init`
    ///   producing the type.
    ///
    /// Silent on: anything in the type's module or in any impl of it, trait
    /// impls included, so a written or derived `Default`, `From` or `Clone`
    /// is the type's own business; any site in a crate other than the one
    /// defining the type, since validators are only known for local impls,
    /// so a `pub` checked field of an exported type is not covered; writes to
    /// fields no constructor checks; `T::default()` and every other call
    /// (whatever literal it ends in is wherever the callee is), and a tuple
    /// constructor passed as a value rather than called; a `transmute` that
    /// only changes the lifetimes of a value already of the type; and
    /// `&mut s.f` handed to something else, which is not an assignment here.
    /// Constructors that only fail because their input did not parse or a
    /// resource ran out check nothing about the fields and are not
    /// validators. A constructor whose OWN return type is not `Self` is not
    /// seen at all, even when it checks and stores an argument's fields
    /// unchanged: `Thing::new(c: Cfg) -> Result<Thing, _>` registers `Thing`,
    /// never `Cfg`, because only a receiver-less fn returning `Self` is
    /// looked at in the first place, and `Cfg` is never `Self` from inside
    /// `impl Thing`. `Cfg`'s own fields and outside literals go unflagged
    /// even though `Thing::new` is the only thing standing between them and
    /// whatever `Cfg` was meant to guarantee.
    pub UNCHECKED_CONSTRUCTION,
    Warn,
    "value of a validated type made or changed without its validating constructor"
}

struct Validator {
    ctor: Symbol,
    /// A validator checks at least one field; a finding about a whole value
    /// points at this one's check.
    first: FieldCheck,
    rest: Vec<FieldCheck>,
}

impl Validator {
    fn checks(&self) -> impl Iterator<Item = &FieldCheck> {
        std::iter::once(&self.first).chain(&self.rest)
    }
}

/// Which fields a literal supplies itself.
enum Literal {
    All,
    /// A `..base` (or `..` default-fields) literal: only these are new values.
    Only(Vec<FieldIdx>),
}

enum SiteKind {
    Literal(Literal),
    Write(FieldIdx),
    /// `mem::zeroed` and friends, named the way the message shows them.
    Conjured(&'static str),
}

/// Something done to a crate-local struct outside its own module and impls.
/// Resolved against `validators` at the end of the crate: the impl holding
/// the constructor may be visited after the code that goes around it.
struct Site {
    adt: DefId,
    span: Span,
    kind: SiteKind,
}

pub struct UncheckedConstruction {
    extra_resource_errors: Vec<String>,
    /// struct -> constructors that check at least one stored field.
    validators: HashMap<DefId, Vec<Validator>>,
    sites: Vec<Site>,
}

rustc_session::impl_lint_pass!(UncheckedConstruction => [UNCHECKED_CONSTRUCTION]);

impl UncheckedConstruction {
    pub fn new(config: &crate::MordantConfig) -> Self {
        Self {
            extra_resource_errors: config.validator_resource_errors.clone(),
            validators: HashMap::new(),
            sites: Vec::new(),
        }
    }

    fn note_site(
        &mut self,
        cx: &LateContext<'_>,
        expr: &Expr<'_>,
        adt: ty::AdtDef<'_>,
        kind: SiteKind,
    ) {
        if !adt.is_struct()
            || !adt.did().is_local()
            || is_types_own_code(cx, expr.hir_id, adt.did())
        {
            return;
        }
        self.sites.push(Site {
            adt: adt.did(),
            span: expr.span,
            kind,
        });
    }
}

/// The struct an impl block is for, when it is a crate-local struct.
fn impl_self_struct(cx: &LateContext<'_>, impl_did: DefId) -> Option<DefId> {
    let adt = impl_self_adt(cx, impl_did)?;
    (adt.is_struct() && adt.did().is_local()).then(|| adt.did())
}

/// The type's own module can build and write it whatever the field
/// visibility (a static table its `Option<Self>` lookup searches, a sibling
/// helper), and so can any impl of it from any module (constructors,
/// `Default`, builders): that code is the implementation, not a caller going
/// around it.
fn is_types_own_code(cx: &LateContext<'_>, at: HirId, struct_did: DefId) -> bool {
    if cx.tcx.parent_module(at) == cx.tcx.parent_module_from_def_id(struct_did.expect_local()) {
        return true;
    }
    let mut cur = cx.tcx.hir_enclosing_body_owner(at).to_def_id();
    while let Some(parent) = cx.tcx.opt_parent(cur) {
        if matches!(cx.tcx.def_kind(parent), DefKind::Impl { .. })
            && impl_self_struct(cx, parent) == Some(struct_did)
        {
            return true;
        }
        cur = parent;
    }
    false
}

/// The name a finding gives a callee that produces a value out of nothing.
fn conjurer(cx: &LateContext<'_>, callee: DefId) -> Option<&'static str> {
    let name = cx.tcx.get_diagnostic_name(callee)?;
    if name == sym::mem_zeroed {
        Some("mem::zeroed")
    } else if name == sym::transmute {
        Some("mem::transmute")
    } else if name == sym::assume_init {
        Some("MaybeUninit::assume_init")
    } else {
        None
    }
}

impl<'tcx> LateLintPass<'tcx> for UncheckedConstruction {
    fn check_impl_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx ImplItem<'tcx>) {
        // Only inherent-impl functions count as validators; a trait impl's
        // signature is the trait's idea, not this type's promise.
        let parent = cx.tcx.local_parent(item.owner_id.def_id);
        let ItemKind::Impl(imp) = cx.tcx.hir_expect_item(parent).kind else {
            return;
        };
        if imp.of_trait.is_some() {
            return;
        }
        let Some(struct_did) = impl_self_struct(cx, parent.to_def_id()) else {
            return;
        };
        let ImplItemKind::Fn(sig, _) = &item.kind else {
            return;
        };
        // A constructor has no receiver. `fn parent(&self) -> Option<&Self>`
        // and `fn clone(&self) -> Result<Self, _>` navigate or copy a value
        // that already passed whatever check exists; they establish nothing.
        if cx.tcx.associated_item(item.owner_id).is_method() {
            return;
        }
        let FnRetTy::Return(ret_hir_ty) = sig.decl.output else {
            return;
        };
        let output = ty_from_hir_ty(cx, ret_hir_ty);
        let ty::Adt(adt, args) = output.kind() else {
            return;
        };
        // `Self` by value only: `Option<&Self>` from a receiver-less fn is a
        // lookup into a table of existing values, not construction.
        let wraps_self = (cx.tcx.is_diagnostic_item(sym::Result, adt.did())
            || cx.tcx.is_diagnostic_item(sym::Option, adt.did()))
            && matches!(args.type_at(0).kind(), ty::Adt(inner, _) if inner.did() == struct_did);
        if !wraps_self {
            return;
        }
        let mut checks = ctor_flow::checked_fields(
            cx,
            item.owner_id.def_id,
            struct_did,
            &self.extra_resource_errors,
        )
        .into_iter();
        let Some(first) = checks.next() else {
            return;
        };
        self.validators
            .entry(struct_did)
            .or_default()
            .push(Validator {
                ctor: item.ident.name,
                first,
                rest: checks.collect(),
            });
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let results = cx.typeck_results();
        match expr.kind {
            ExprKind::Struct(_, fields, tail) => {
                let Some(adt) = results.expr_ty(expr).ty_adt_def() else {
                    return;
                };
                let literal = match tail {
                    StructTailExpr::None => Literal::All,
                    _ => Literal::Only(
                        fields
                            .iter()
                            .filter_map(|f| results.opt_field_index(f.hir_id))
                            .collect(),
                    ),
                };
                self.note_site(cx, expr, adt, SiteKind::Literal(literal));
            }
            ExprKind::Assign(place, ..) | ExprKind::AssignOp(_, place, _) => {
                let Some((adt, _, place)) = assigned_adt_field(cx, place) else {
                    return;
                };
                let Some(field) = results.opt_field_index(place.hir_id) else {
                    return;
                };
                self.note_site(cx, expr, adt, SiteKind::Write(field));
            }
            ExprKind::Call(..) | ExprKind::MethodCall(..) => {
                let Some(adt) = results.expr_ty(expr).ty_adt_def() else {
                    return;
                };
                let Some(def) = callee_of(cx, expr).map(|c| c.def()) else {
                    return;
                };
                // `S(x)` has no callee body: the call is the literal.
                if matches!(
                    cx.tcx.def_kind(def),
                    DefKind::Ctor(CtorOf::Struct, CtorKind::Fn)
                ) {
                    self.note_site(cx, expr, adt, SiteKind::Literal(Literal::All));
                    return;
                }
                let Some(what) = conjurer(cx, def) else {
                    return;
                };
                // `transmute::<S<'a>, S<'static>>(s)` re-types a value that
                // already went through the constructor; it makes nothing.
                // Changing a type parameter (`S<Unchecked>` to `S<Checked>`)
                // is what a validator returning the latter exists to gate.
                if let ExprKind::Call(_, [arg]) = expr.kind
                    && cx.tcx.erase_and_anonymize_regions(results.expr_ty(arg))
                        == cx.tcx.erase_and_anonymize_regions(results.expr_ty(expr))
                {
                    return;
                }
                self.note_site(cx, expr, adt, SiteKind::Conjured(what));
            }
            _ => {}
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        for site in &self.sites {
            let Some(ctors) = self.validators.get(&site.adt) else {
                continue;
            };
            let fields = &cx.tcx.adt_def(site.adt).non_enum_variant().fields;
            let checked = |v: &Validator| {
                v.checks()
                    .map(|c| format!("`{}`", fields[c.field].name))
                    .collect::<Vec<_>>()
                    .join(", ")
            };
            // Which constructor a site is charged to: for a whole value the
            // first one, otherwise whichever checks a field the site touches,
            // which is not necessarily the first.
            let field_checker = |is_hit: &dyn Fn(FieldIdx) -> bool| {
                ctors
                    .iter()
                    .find_map(|v| v.checks().find(|c| is_hit(c.field)).map(|c| (v, c)))
            };
            let path = cx.tcx.def_path_str(site.adt);
            let ty = cx.tcx.item_name(site.adt);
            let (msg, check, note, help) = match &site.kind {
                SiteKind::Literal(literal) => {
                    let by = match literal {
                        Literal::All => ctors.first().map(|v| (v, &v.first)),
                        Literal::Only(supplied) => field_checker(&|f| supplied.contains(&f)),
                    };
                    let Some((by, check)) = by else {
                        continue;
                    };
                    let (ctor, checked) = (by.ctor, checked(by));
                    (
                        format!(
                            "`{path}` is built by hand here, skipping the check on {checked} that `{ty}::{ctor}` makes"
                        ),
                        check.check,
                        format!("the check in `{ty}::{ctor}`"),
                        format!(
                            "build it with `{ty}::{ctor}(..)`, or make {checked} private so a literal like this only compiles inside `{ty}`'s module"
                        ),
                    )
                }
                SiteKind::Conjured(what) => {
                    let Some(by) = ctors.first() else {
                        continue;
                    };
                    let (ctor, checked) = (by.ctor, checked(by));
                    (
                        format!(
                            "`{path}` is produced with `{what}` here, skipping the check on {checked} that `{ty}::{ctor}` makes"
                        ),
                        by.first.check,
                        format!("the check in `{ty}::{ctor}`"),
                        format!(
                            "build it with `{ty}::{ctor}(..)`. If this path must skip the check, add an unchecked constructor next to `{ty}::{ctor}` and call that"
                        ),
                    )
                }
                SiteKind::Write(field) => {
                    let Some((by, check)) = field_checker(&|f| f == *field) else {
                        continue;
                    };
                    let name = fields[*field].name;
                    let ctor = by.ctor;
                    (
                        format!(
                            "`{path}::{name}` is assigned directly here, skipping the check on `{name}` that `{ty}::{ctor}` makes"
                        ),
                        check.check,
                        format!("the check in `{ty}::{ctor}`"),
                        format!(
                            "change `{name}` through a method of `{ty}` that repeats the check, or make `{name}` private so this write only compiles inside `{ty}`'s module"
                        ),
                    )
                }
            };
            emit_with_note(
                cx,
                UNCHECKED_CONSTRUCTION,
                site.span,
                msg,
                check,
                note,
                help,
            );
        }
    }
}

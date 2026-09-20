use std::collections::{BTreeSet, HashMap};

use clippy_utils::last_path_segment;
use rustc_ast::LitKind;
use rustc_hir::def::Res;
use rustc_hir::def_id::DefId;
use rustc_hir::{
    Arm, BinOpKind, BindingMode, BorrowKind, ByRef, Expr, ExprKind, HirId, LetStmt, MatchSource,
    Mutability, Node, Pat, PatExpr, PatExprKind, PatKind, QPath, StructTailExpr, UnOp,
};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::adjustment::{Adjust, AutoBorrow, AutoBorrowMutability};
use rustc_middle::ty::{self, Ty, TypeckResults};
use rustc_span::{Span, Symbol, sym};

use crate::adt_facts::{field_ty, has_fixed_repr, has_positional_fields, struct_field};
use crate::baseline::{emit_with_note, join};
use crate::hir_shapes::{Callee, callee_of, peel_blocks_unsafe};

rustc_session::declare_lint! {
    /// Flags a string or byte-string field or local that only ever holds one
    /// of a fixed set of literals and is read by comparing against literals.
    /// It is an enum written as a string, so a typo on either side still
    /// compiles, a new state needs every comparison found by hand, and
    /// nothing says which strings are possible.
    ///
    /// Only fires where every store is visible: a local, or a field no other
    /// crate can name (the struct or the field is private to this one). Every
    /// store must write a literal (or an `if`/`match` choosing between
    /// literals) or copy the same place from another value, and at least two
    /// distinct literals must occur. A single non-literal store, a `..base`
    /// construction, a `&mut` borrow or `ref mut` binding, a write into part
    /// of the value (`x.f[i] = b`), an explicit `repr`, or a value that is
    /// only ever formatted or written out and never compared keeps it silent.
    pub STRINGLY_STATE,
    Warn,
    "a string only ever holding one of a closed set of literals, then compared against them"
}

/// A place whose every store this crate can see.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Slot {
    Field(DefId, Symbol),
    Local(HirId),
}

#[derive(Default)]
struct Facts {
    /// A store the lint cannot read as a literal, or a mutable borrow.
    open: bool,
    stores: usize,
    values: BTreeSet<String>,
    /// The first comparison against a literal.
    compared: Option<Span>,
}

#[derive(Default)]
pub struct StringlyState {
    slots: HashMap<Slot, Facts>,
    /// Where and under what name to report a local.
    locals: HashMap<HirId, (Span, Symbol)>,
}

rustc_session::impl_lint_pass!(StringlyState => [STRINGLY_STATE]);

/// `&str`, `String`, `Box<str>`, `&[u8]`, `Box<[u8]>`, `Vec<u8>`.
fn is_stringy(cx: &LateContext<'_>, ty: Ty<'_>) -> bool {
    fn is_text(ty: Ty<'_>) -> bool {
        match ty.kind() {
            ty::Str => true,
            ty::Slice(elem) => matches!(elem.kind(), ty::Uint(ty::UintTy::U8)),
            _ => false,
        }
    }
    match ty.kind() {
        ty::Ref(_, inner, _) => is_text(*inner),
        ty::Adt(adt, args) => {
            if cx.tcx.is_lang_item(adt.did(), rustc_hir::LangItem::String) {
                true
            } else if adt.is_box() {
                is_text(args.type_at(0))
            } else if cx.tcx.is_diagnostic_item(sym::Vec, adt.did()) {
                matches!(args.type_at(0).kind(), ty::Uint(ty::UintTy::U8))
            } else {
                false
            }
        }
        _ => false,
    }
}

/// The field this lint tracks behind `ty`.`name`, if any: a string field of
/// a struct this crate defines that no other crate can name, whether because
/// the struct is private or the field is. Either way every store to it is in
/// this crate.
fn tracked_field<'tcx>(cx: &LateContext<'tcx>, ty: Ty<'tcx>, name: Symbol) -> Option<Slot> {
    let ty::Adt(adt, _) = ty.peel_refs().kind() else {
        return None;
    };
    if !adt.is_struct() || !adt.did().is_local() {
        return None;
    }
    if has_fixed_repr(*adt) || has_positional_fields(adt.non_enum_variant()) {
        return None;
    }
    let f = struct_field(*adt, name)?;
    if cx.effective_visibilities.is_exported(f.did.expect_local()) {
        return None;
    }
    is_stringy(cx, field_ty(cx, f)).then_some(Slot::Field(adt.did(), name))
}

fn local_of(cx: &LateContext<'_>, e: &Expr<'_>) -> Option<HirId> {
    if let ExprKind::Path(QPath::Resolved(None, path)) = e.kind
        && let Res::Local(id) = path.res
        && is_stringy(cx, cx.typeck_results().node_type(id))
    {
        Some(id)
    } else {
        None
    }
}

fn lit_text(lit: &LitKind) -> Option<String> {
    match lit {
        LitKind::Str(s, _) => Some(format!("{:?}", s.as_str())),
        LitKind::ByteStr(b, _) => Some(format!("\"{}\"", b.as_byte_str().escape_ascii())),
        _ => None,
    }
}

/// The literals a stored expression spells, through the conversions that
/// turn a literal into an owned string without changing its text and through
/// an `if` or `match` each of whose arms spells one. False when any part is
/// something else.
fn stored_literals(e: &Expr<'_>, out: &mut BTreeSet<String>) -> bool {
    let e = peel_blocks_unsafe(e);
    match e.kind {
        ExprKind::Lit(lit) => match lit_text(&lit.node) {
            Some(text) => {
                out.insert(text);
                true
            }
            None => false,
        },
        ExprKind::AddrOf(BorrowKind::Ref, Mutability::Not, inner)
        | ExprKind::DropTemps(inner)
        | ExprKind::Cast(inner, _) => stored_literals(inner, out),
        ExprKind::MethodCall(seg, recv, [], _)
            if matches!(
                seg.ident.name.as_str(),
                "to_string"
                    | "to_owned"
                    | "into"
                    | "to_vec"
                    | "as_bytes"
                    | "into_bytes"
                    | "into_boxed_str"
                    | "into_boxed_slice"
                    | "clone"
            ) =>
        {
            stored_literals(recv, out)
        }
        ExprKind::Call(callee, [arg])
            if matches!(callee.kind, ExprKind::Path(ref qp)
                if last_path_segment(qp).ident.name == sym::from) =>
        {
            stored_literals(arg, out)
        }
        ExprKind::If(_, then, Some(els)) => stored_literals(then, out) && stored_literals(els, out),
        ExprKind::Match(_, arms, _) if !arms.is_empty() => {
            arms.iter().all(|a| stored_literals(a.body, out))
        }
        _ => false,
    }
}

fn is_literal_expr(e: &Expr<'_>) -> bool {
    let e = peel_blocks_unsafe(e);
    match e.kind {
        ExprKind::Lit(lit) => lit_text(&lit.node).is_some(),
        ExprKind::AddrOf(BorrowKind::Ref, Mutability::Not, inner) | ExprKind::DropTemps(inner) => {
            is_literal_expr(inner)
        }
        _ => false,
    }
}

fn pat_has_literal(pat: &Pat<'_>) -> bool {
    match pat.kind {
        PatKind::Expr(PatExpr {
            kind: PatExprKind::Lit { lit, .. },
            ..
        }) => lit_text(&lit.node).is_some(),
        PatKind::Or(pats) => pats.iter().any(|p| pat_has_literal(p)),
        PatKind::Ref(inner, ..) | PatKind::Box(inner) | PatKind::Deref(inner) => {
            pat_has_literal(inner)
        }
        _ => false,
    }
}

/// The tracked place `e` reads, through the views that leave its text as it
/// is: `&`, `*`, `x.as_str()`, `x.as_bytes()`, `&x[..]`.
fn read_slot<'tcx>(cx: &LateContext<'tcx>, e: &'tcx Expr<'tcx>) -> Option<Slot> {
    let e = peel_blocks_unsafe(e);
    match e.kind {
        ExprKind::Field(base, ident) => {
            let ty = cx.typeck_results().expr_ty_adjusted(base);
            tracked_field(cx, ty, ident.name)
        }
        ExprKind::Path(_) => local_of(cx, e).map(Slot::Local),
        ExprKind::AddrOf(BorrowKind::Ref, Mutability::Not, inner)
        | ExprKind::Unary(UnOp::Deref, inner)
        | ExprKind::DropTemps(inner) => read_slot(cx, inner),
        ExprKind::Index(inner, idx, _)
            if matches!(
                peel_blocks_unsafe(idx).kind,
                ExprKind::Struct(_, [], StructTailExpr::None)
            ) =>
        {
            read_slot(cx, inner)
        }
        ExprKind::MethodCall(seg, recv, [], _)
            if matches!(
                seg.ident.name.as_str(),
                "as_str" | "as_bytes" | "as_ref" | "as_slice" | "as_deref" | "deref" | "borrow"
            ) =>
        {
            read_slot(cx, recv)
        }
        _ => None,
    }
}

/// The tracked place a stored value copies its text from unchanged: `x.f`,
/// `x.f.clone()`, `Clone::clone(&x.f)` (what a derived `Clone` writes),
/// `s.to_owned()`.
fn copied_slot<'tcx>(cx: &LateContext<'tcx>, e: &'tcx Expr<'tcx>) -> Option<Slot> {
    let e = peel_blocks_unsafe(e);
    match e.kind {
        ExprKind::MethodCall(seg, recv, [], _)
            if matches!(
                seg.ident.name.as_str(),
                "clone" | "to_owned" | "to_string" | "to_vec" | "into"
            ) =>
        {
            read_slot(cx, recv)
        }
        ExprKind::Call(callee, [arg])
            if matches!(callee.kind, ExprKind::Path(ref qp)
                if last_path_segment(qp).ident.name == sym::clone) =>
        {
            read_slot(cx, arg)
        }
        _ => read_slot(cx, e),
    }
}

/// The place an assignment target or a `&mut` operand names.
struct Written {
    slot: Slot,
    /// It names the place itself, not an element or range of it (`x.f[i]`):
    /// only then does an assignment replace the whole value.
    whole: bool,
}

fn written_slot<'tcx>(cx: &LateContext<'tcx>, place: &'tcx Expr<'tcx>) -> Option<Written> {
    let mut place = peel_blocks_unsafe(place);
    let mut whole = true;
    while let ExprKind::Index(inner, ..)
    | ExprKind::Unary(UnOp::Deref, inner)
    | ExprKind::DropTemps(inner) = place.kind
    {
        whole &= !matches!(place.kind, ExprKind::Index(..));
        place = inner;
    }
    let slot = match place.kind {
        // The adjusted type, so a write through a `Box`, a guard or any
        // other `Deref` container reaches the struct behind it.
        ExprKind::Field(base, ident) => {
            tracked_field(cx, cx.typeck_results().expr_ty_adjusted(base), ident.name)
        }
        _ => local_of(cx, place).map(Slot::Local),
    }?;
    Some(Written { slot, whole })
}

fn binds_ref_mut(typeck: &TypeckResults<'_>, id: HirId) -> bool {
    typeck
        .pat_binding_modes()
        .get(id)
        .is_some_and(|m| matches!(m.0, ByRef::Yes(_, Mutability::Mut)))
}

/// The value a top-level binding pattern is matched against: the `let`
/// initializer, the `if let` operand or the `match` scrutinee.
fn bound_value<'tcx>(cx: &LateContext<'tcx>, pat: &Pat<'tcx>) -> Option<&'tcx Expr<'tcx>> {
    match cx.tcx.parent_hir_node(pat.hir_id) {
        Node::LetStmt(l) => l.init,
        Node::Expr(e) => match e.kind {
            ExprKind::Let(l) => Some(l.init),
            _ => None,
        },
        Node::Arm(arm) => match cx.tcx.parent_hir_node(arm.hir_id) {
            Node::Expr(&Expr {
                kind: ExprKind::Match(scrut, ..),
                ..
            }) => Some(scrut),
            _ => None,
        },
        _ => None,
    }
}

impl StringlyState {
    fn facts(&mut self, slot: Slot) -> &mut Facts {
        self.slots.entry(slot).or_default()
    }

    fn store<'tcx>(&mut self, cx: &LateContext<'tcx>, slot: Slot, value: &'tcx Expr<'tcx>) {
        // Copying the same place from another value adds nothing to the set.
        if copied_slot(cx, value) == Some(slot) {
            return;
        }
        let facts = self.facts(slot);
        if stored_literals(value, &mut facts.values) {
            facts.stores += 1;
        } else {
            facts.open = true;
        }
    }

    fn compared(&mut self, slot: Slot, at: Span) {
        self.facts(slot).compared.get_or_insert(at);
    }

    /// `a == "lit"`, `a != b"lit"`, either way round.
    fn note_comparison<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        at: Span,
        l: &'tcx Expr<'tcx>,
        r: &'tcx Expr<'tcx>,
    ) {
        for (place, other) in [(l, r), (r, l)] {
            if is_literal_expr(other)
                && let Some(slot) = read_slot(cx, place)
            {
                self.compared(slot, at);
            }
        }
    }

    /// `eql(x.f, b"lit")`, `x.eq_ignore_ascii_case("lit")` and kin: a call
    /// whose name says it compares, with the place and a literal among its
    /// operands.
    fn note_comparing_call<'tcx>(&mut self, cx: &LateContext<'tcx>, e: &'tcx Expr<'tcx>) {
        let Some(callee) = callee_of(cx, e) else {
            return;
        };
        let Some(name) = cx.tcx.opt_item_name(callee.def()) else {
            return;
        };
        let name = name.as_str();
        let compares = name.starts_with("eql")
            || name.starts_with("eq_ignore_ascii_case")
            || name == "eq"
            || name == "ne"
            || name == "starts_with"
            || name == "ends_with"
            || name.starts_with("has_prefix")
            || name.starts_with("has_suffix");
        if !compares {
            return;
        }
        let operands: Vec<&'tcx Expr<'tcx>> = match callee {
            Callee::Path { args, .. } => args.iter().collect(),
            Callee::Method { recv, args, .. } => std::iter::once(recv).chain(args.iter()).collect(),
        };
        if !operands.iter().any(|o| is_literal_expr(o)) {
            return;
        }
        for o in operands {
            if let Some(slot) = read_slot(cx, o) {
                self.compared(slot, e.span);
            }
        }
    }

    fn note_match<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        at: Span,
        scrut: &'tcx Expr<'tcx>,
        arms: &'tcx [Arm<'tcx>],
    ) {
        if let Some(slot) = read_slot(cx, scrut)
            && arms.iter().any(|a| pat_has_literal(a.pat))
        {
            self.compared(slot, at);
        }
    }

    /// Where to report, the message, and the help.
    fn message(
        &self,
        cx: &LateContext<'_>,
        slot: Slot,
        facts: &Facts,
    ) -> Option<(Span, String, String)> {
        let (span, name) = match slot {
            Slot::Field(did, name) => {
                let field = struct_field(cx.tcx.adt_def(did), name)?;
                (cx.tcx.def_span(field.did), name)
            }
            Slot::Local(id) => *self.locals.get(&id)?,
        };
        const SHOWN: usize = 6;
        let mut list: Vec<String> = facts
            .values
            .iter()
            .take(SHOWN)
            .map(|v| format!("`{v}`"))
            .collect();
        let more = facts.values.len().saturating_sub(SHOWN);
        if more > 0 {
            list.push(format!("{more} more"));
        }
        let values = join(&list, "or");
        Some((
            span,
            format!(
                "`{name}` only ever holds {values} ({} {}) and is read by comparing against \
                 literals. It is an enum written as a string, and a typo still compiles",
                facts.stores,
                if facts.stores == 1 { "store" } else { "stores" },
            ),
            format!(
                "declare an enum with a variant per string and store that in `{name}`. Keep the \
                 text in an `as_str` method for printing"
            ),
        ))
    }
}

impl<'tcx> LateLintPass<'tcx> for StringlyState {
    fn check_local(&mut self, cx: &LateContext<'tcx>, l: &'tcx LetStmt<'tcx>) {
        let PatKind::Binding(BindingMode(ByRef::No, _), id, ident, None) = l.pat.kind else {
            return;
        };
        if l.span.from_expansion() || !is_stringy(cx, cx.typeck_results().node_type(id)) {
            return;
        }
        self.locals.insert(id, (l.pat.span, ident.name));
        if let Some(init) = l.init {
            self.store(cx, Slot::Local(id), init);
        }
    }

    /// A `ref mut` binding, spelt out or implied by matching through a
    /// `&mut`, is a write this lint cannot read.
    fn check_pat(&mut self, cx: &LateContext<'tcx>, pat: &'tcx Pat<'tcx>) {
        let Some(typeck) = cx.maybe_typeck_results() else {
            return;
        };
        match pat.kind {
            PatKind::Struct(_, fields, _) => {
                let ty = typeck.pat_ty(pat);
                for f in fields {
                    let mut by_ref_mut = false;
                    f.pat
                        .each_binding(|_, id, _, _| by_ref_mut |= binds_ref_mut(typeck, id));
                    if by_ref_mut && let Some(slot) = tracked_field(cx, ty, f.ident.name) {
                        self.facts(slot).open = true;
                    }
                }
            }
            PatKind::Binding(..) if binds_ref_mut(typeck, pat.hir_id) => {
                if let Some(value) = bound_value(cx, pat)
                    && let Some(slot) = read_slot(cx, value)
                {
                    self.facts(slot).open = true;
                }
            }
            _ => {}
        }
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        match expr.kind {
            ExprKind::Struct(_, fields, tail) => {
                let ty = cx.typeck_results().expr_ty(expr);
                let Some(adt) = ty.ty_adt_def().filter(|a| a.is_struct()) else {
                    return;
                };
                for f in fields {
                    if let Some(slot) = tracked_field(cx, ty, f.ident.name) {
                        self.store(cx, slot, f.expr);
                    }
                }
                if matches!(tail, StructTailExpr::None) {
                    return;
                }
                // `..base` fills every unlisted field from a value this site
                // does not spell.
                for f in &adt.non_enum_variant().fields {
                    if !fields.iter().any(|l| l.ident.name == f.name)
                        && let Some(slot) = tracked_field(cx, ty, f.name)
                    {
                        self.facts(slot).open = true;
                    }
                }
            }
            ExprKind::Assign(place, value, _) => match written_slot(cx, place) {
                Some(Written { slot, whole: true }) => self.store(cx, slot, value),
                Some(Written { slot, whole: false }) => self.facts(slot).open = true,
                None => {}
            },
            ExprKind::AssignOp(_, place, _)
            | ExprKind::AddrOf(BorrowKind::Ref | BorrowKind::Raw, Mutability::Mut, place) => {
                if let Some(w) = written_slot(cx, place) {
                    self.facts(w.slot).open = true;
                }
            }
            // The auto-`&mut` a mutating method call takes: the place is
            // written through something this lint does not read.
            ExprKind::Field(..) | ExprKind::Path(..) | ExprKind::Index(..) => {
                let mutably_borrowed = cx.typeck_results().expr_adjustments(expr).iter().any(|a| {
                    matches!(
                        a.kind,
                        Adjust::Borrow(AutoBorrow::Ref(AutoBorrowMutability::Mut { .. }))
                            | Adjust::Borrow(AutoBorrow::RawPtr(Mutability::Mut))
                    )
                });
                if mutably_borrowed && let Some(w) = written_slot(cx, expr) {
                    self.facts(w.slot).open = true;
                }
            }
            ExprKind::Binary(op, l, r) if matches!(op.node, BinOpKind::Eq | BinOpKind::Ne) => {
                self.note_comparison(cx, expr.span, l, r);
            }
            ExprKind::Match(scrut, arms, MatchSource::Normal | MatchSource::Postfix) => {
                self.note_match(cx, scrut.span, scrut, arms);
            }
            ExprKind::Call(..) | ExprKind::MethodCall(..) => self.note_comparing_call(cx, expr),
            _ => {}
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let mut findings: Vec<(Span, String, Span, String)> = self
            .slots
            .iter()
            .filter(|(_, f)| !f.open && f.values.len() >= 2)
            .filter_map(|(slot, f)| {
                let compared = f.compared?;
                let (span, msg, help) = self.message(cx, *slot, f)?;
                Some((span, msg, compared, help))
            })
            .collect();
        findings.sort_by_key(|(span, ..)| span.lo());
        for (span, msg, compared, help) in findings {
            emit_with_note(
                cx,
                STRINGLY_STATE,
                span,
                msg,
                compared,
                "one of the comparisons",
                help,
            );
        }
    }
}

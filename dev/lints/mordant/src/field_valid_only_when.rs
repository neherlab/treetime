use std::collections::{HashMap, HashSet};

use crate::adt_facts::{has_fixed_repr, has_positional_fields, private_local_struct, struct_field};
use crate::baseline::emit_with_note;
use crate::enum_facts::{arm_variant, ctor_literal_variant};
use clippy_utils::eq_expr_value;
use clippy_utils::{get_enclosing_block, in_automatically_derived, is_default_equivalent};
use rustc_ast::LitKind;
use rustc_hir::def_id::DefId;
use rustc_hir::{
    Arm, BinOpKind, Block, BorrowKind, Expr, ExprKind, HirId, Mutability, Node, Pat, PatExpr,
    PatExprKind, PatKind, StmtKind, StructTailExpr, UnOp,
};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_middle::ty::adjustment::{Adjust, AutoBorrow, AutoBorrowMutability};
use rustc_span::{Span, Symbol, SyntaxContext, sym};

rustc_session::declare_lint! {
    /// Flags a struct field that only means something when a sibling field
    /// has one particular value: every read in the crate comes after an `if`,
    /// `match`, `let .. else` or diverging guard that tests the sibling
    /// against the same variant or literal, and every construction that
    /// gives the sibling another value fills the field with a placeholder
    /// (`None`, `0`, `""`, `false`, `Default::default()`, a null pointer).
    /// The struct carries the field in every case anyway, so any new reader
    /// can use it without the test. Making that case of the sibling an enum
    /// variant that carries the field leaves the other cases nothing to fill
    /// in or misread.
    ///
    /// Only fires on structs private to the crate with named fields and no
    /// explicit `repr`, and only on proof: one read the lint cannot place
    /// under such a test — an accessor, a destructuring pattern, a `Debug`
    /// written by hand, a test on a copy of the sibling rather than the
    /// sibling itself — keeps it quiet. So does any write of a real value
    /// outside that case: a construction that gives the field one beside
    /// another value of the sibling (a field taken from `..other` counts as
    /// real, one from `..Default::default()` as a placeholder), or an
    /// assignment to it neither under the test nor in a block that also
    /// assigns the sibling that value. A
    /// sibling never given the tested value by any literal or assignment
    /// keeps it quiet too (the field is then dead, not dependent). Reads in
    /// derived impls are not counted.
    pub FIELD_VALID_ONLY_WHEN,
    Warn,
    "a field that only means something when a sibling field has one value"
}

/// A value a field can be tested against and built with, compared exactly.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Value {
    Variant(DefId),
    Bool(bool),
    Int(u128),
}

impl Value {
    /// The one other value the field can hold, when there is exactly one.
    fn complement(self) -> Option<Value> {
        match self {
            Value::Bool(b) => Some(Value::Bool(!b)),
            Value::Variant(_) | Value::Int(_) => None,
        }
    }
}

/// `sibling == value`, known to hold where a read happens.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct Test {
    sibling: Symbol,
    value: Value,
}

/// What one construction site wrote into one field.
#[derive(Clone, Copy)]
struct Init {
    value: Option<Value>,
    placeholder: bool,
}

/// One read of a field: where it is and the sibling tests that dominate it.
/// An empty set is a read nothing is known at.
struct Read {
    at: Span,
    tests: HashSet<Test>,
}

#[derive(Default)]
pub struct FieldValidOnlyWhen {
    /// (struct, field) -> its reads.
    reads: HashMap<(DefId, Symbol), Vec<Read>>,
    /// struct -> per literal construction site, what each named field got.
    sites: HashMap<DefId, Vec<HashMap<Symbol, Init>>>,
    /// (struct, field) -> the values assigned to it after construction
    /// (`None` for one not spelled out), which the sites do not bound.
    assigned: HashMap<(DefId, Symbol), HashSet<Option<Value>>>,
    /// (struct, field) -> per assignment of a real value, the sibling tests
    /// dominating it plus the sibling values assigned beside it in the same
    /// block: the cases that write can belong to.
    writes: HashMap<(DefId, Symbol), Vec<HashSet<Test>>>,
}

rustc_session::impl_lint_pass!(FieldValidOnlyWhen => [FIELD_VALID_ONLY_WHEN]);

/// The struct behind `ty` when every read and construction of it is this
/// crate's to see and its fields have names a message can use.
fn relevant<'tcx>(cx: &LateContext<'tcx>, ty: ty::Ty<'tcx>) -> Option<ty::AdtDef<'tcx>> {
    let adt = private_local_struct(cx, ty)?;
    let v = adt.non_enum_variant();
    (!has_fixed_repr(adt) && !has_positional_fields(v) && v.fields.len() >= 2).then_some(adt)
}

/// `e` under HIR temporaries and shared borrows: `&x.g` names what `x.g`
/// does for a comparison or a `match`.
fn strip_ref<'h>(mut e: &'h Expr<'h>) -> &'h Expr<'h> {
    while let ExprKind::DropTemps(inner)
    | ExprKind::AddrOf(BorrowKind::Ref, Mutability::Not, inner) = e.kind
    {
        e = inner;
    }
    e
}

/// `e` with `!`s removed, and whether they flip it (their parity, which
/// `hir_shapes::peel_not` does not keep).
fn peel_parity<'h>(mut e: &'h Expr<'h>) -> (&'h Expr<'h>, bool) {
    let mut flipped = false;
    loop {
        match e.kind {
            ExprKind::Unary(UnOp::Not, inner) => {
                flipped = !flipped;
                e = inner;
            }
            ExprKind::DropTemps(inner) => e = inner,
            _ => return (e, flipped),
        }
    }
}

fn lit_value(lit: &LitKind) -> Option<Value> {
    match *lit {
        LitKind::Bool(b) => Some(Value::Bool(b)),
        LitKind::Int(n, _) => Some(Value::Int(n.get())),
        _ => None,
    }
}

/// The value `e` spells out: a variant literal, a bool or an unsigned
/// integer literal. Named constants are not values here: two names can be
/// one value.
fn expr_value(cx: &LateContext<'_>, e: &Expr<'_>) -> Option<Value> {
    let e = strip_ref(e);
    if let Some(v) = ctor_literal_variant(cx, e) {
        return Some(Value::Variant(v));
    }
    match e.kind {
        ExprKind::Lit(lit) => lit_value(&lit.node),
        _ => None,
    }
}

/// The one value a pattern admits at its head; or-patterns, ranges,
/// bindings without a sub-pattern and wildcards admit more than one.
fn pat_value(cx: &LateContext<'_>, pat: &Pat<'_>) -> Option<Value> {
    match pat.kind {
        PatKind::Binding(.., Some(sub))
        | PatKind::Ref(sub, ..)
        | PatKind::Deref(sub)
        | PatKind::Box(sub) => pat_value(cx, sub),
        PatKind::Expr(PatExpr {
            kind:
                PatExprKind::Lit {
                    lit,
                    negated: false,
                },
            ..
        }) => lit_value(&lit.node),
        _ => arm_variant(cx, pat).map(Value::Variant),
    }
}

/// Whether variant `v` has no fields, so `!= E::V` rules the variant out and
/// not just one payload of it.
fn negation_rules_out(cx: &LateContext<'_>, v: Value) -> bool {
    match v {
        Value::Bool(_) | Value::Int(_) => true,
        Value::Variant(did) => cx
            .tcx
            .adt_def(cx.tcx.parent(did))
            .variant_with_id(did)
            .fields
            .is_empty(),
    }
}

/// A don't-care initialiser: whatever `Default` would give, or a null
/// pointer.
fn is_placeholder(cx: &LateContext<'_>, e: &Expr<'_>) -> bool {
    if is_default_equivalent(cx, e) {
        return true;
    }
    if let ExprKind::Call(f, []) = e.kind
        && let ExprKind::Path(qpath) = &f.kind
        && let Some(def) = cx.qpath_res(qpath, f.hir_id).opt_def_id()
    {
        return matches!(
            cx.tcx.get_diagnostic_name(def),
            Some(sym::ptr_null | sym::ptr_null_mut)
        );
    }
    false
}

/// Which sibling of the read field `e` names: `e` must be `<base>.g` over
/// the very place expression the read was `<base>.f` over.
fn sibling_of(cx: &LateContext<'_>, base: &Expr<'_>, e: &Expr<'_>) -> Option<Symbol> {
    let ExprKind::Field(other, ident) = strip_ref(e).kind else {
        return None;
    };
    eq_expr_value(cx, SyntaxContext::root(), base, other).then_some(ident.name)
}

/// Every `sibling == value` fact that follows from `cond` evaluating to
/// `holds`, for siblings read off `base`.
fn cond_tests(
    cx: &LateContext<'_>,
    base: &Expr<'_>,
    cond: &Expr<'_>,
    holds: bool,
    out: &mut HashSet<Test>,
) {
    let (inner, flipped) = peel_parity(cond);
    let holds = holds != flipped;
    match inner.kind {
        ExprKind::Binary(op, l, r) if op.node == BinOpKind::And && holds => {
            cond_tests(cx, base, l, true, out);
            cond_tests(cx, base, r, true, out);
        }
        ExprKind::Binary(op, l, r) if op.node == BinOpKind::Or && !holds => {
            cond_tests(cx, base, l, false, out);
            cond_tests(cx, base, r, false, out);
        }
        ExprKind::Binary(op, l, r) if matches!(op.node, BinOpKind::Eq | BinOpKind::Ne) => {
            let equal = (op.node == BinOpKind::Eq) == holds;
            let sides = sibling_of(cx, base, l)
                .map(|s| (s, r))
                .or_else(|| sibling_of(cx, base, r).map(|s| (s, l)));
            let Some((sibling, ve)) = sides else {
                return;
            };
            let Some(value) = expr_value(cx, ve) else {
                return;
            };
            known(cx, sibling, value, equal, out);
        }
        ExprKind::Let(l) if holds => {
            if let Some(sibling) = sibling_of(cx, base, l.init)
                && let Some(value) = pat_value(cx, l.pat)
            {
                known(cx, sibling, value, true, out);
            }
        }
        // `matches!(base.g, P)`.
        ExprKind::Match(scrut, [yes, no], _)
            if is_bool_lit(yes.body, true)
                && is_bool_lit(no.body, false)
                && matches!(no.pat.kind, PatKind::Wild)
                && (holds || yes.guard.is_none()) =>
        {
            if let Some(sibling) = sibling_of(cx, base, scrut)
                && let Some(value) = pat_value(cx, yes.pat)
            {
                known(cx, sibling, value, holds, out);
            }
        }
        _ => {
            if let Some(sibling) = sibling_of(cx, base, inner)
                && cx.typeck_results().expr_ty(inner).is_bool()
            {
                known(cx, sibling, Value::Bool(holds), true, out);
            }
        }
    }
}

fn is_bool_lit(e: &Expr<'_>, b: bool) -> bool {
    matches!(strip_ref(e).kind, ExprKind::Lit(lit) if lit.node == LitKind::Bool(b))
}

/// Record `sibling == value` (or, for `equal == false`, what `!=` pins down,
/// which is something only when the value has exactly one complement).
fn known(
    cx: &LateContext<'_>,
    sibling: Symbol,
    value: Value,
    equal: bool,
    out: &mut HashSet<Test>,
) {
    if equal {
        out.insert(Test { sibling, value });
    } else if negation_rules_out(cx, value)
        && let Some(value) = value.complement()
    {
        out.insert(Test { sibling, value });
    }
}

fn is_never(cx: &LateContext<'_>, e: &Expr<'_>) -> bool {
    cx.typeck_results().expr_ty(e).is_never()
}

/// Facts established by the statements of `block` that run before `child`:
/// a guard whose failing branch diverges, a `let .. else`, a `match` whose
/// every other arm diverges.
fn preceding_guards(
    cx: &LateContext<'_>,
    base: &Expr<'_>,
    block: &Block<'_>,
    child: HirId,
    out: &mut HashSet<Test>,
) {
    let pos = block
        .stmts
        .iter()
        .position(|s| s.hir_id == child || matches!(s.kind, StmtKind::Let(l) if l.hir_id == child))
        .unwrap_or_else(|| {
            if block.expr.is_some_and(|e| e.hir_id == child) {
                block.stmts.len()
            } else {
                0
            }
        });
    for stmt in &block.stmts[..pos] {
        match stmt.kind {
            StmtKind::Let(l) => {
                if l.els.is_some()
                    && let Some(init) = l.init
                    && let Some(sibling) = sibling_of(cx, base, init)
                    && let Some(value) = pat_value(cx, l.pat)
                {
                    known(cx, sibling, value, true, out);
                }
            }
            StmtKind::Expr(e) | StmtKind::Semi(e) => match e.kind {
                ExprKind::If(cond, then, els) => {
                    let then_diverges = is_never(cx, then);
                    let else_diverges = els.is_some_and(|x| is_never(cx, x));
                    if then_diverges && !else_diverges {
                        cond_tests(cx, base, cond, false, out);
                    } else if else_diverges && !then_diverges {
                        cond_tests(cx, base, cond, true, out);
                    }
                }
                ExprKind::Match(scrut, arms, _) => {
                    let mut live = arms.iter().filter(|a| !is_never(cx, a.body));
                    if let (Some(only), None) = (live.next(), live.next())
                        && only.guard.is_none()
                        && let Some(sibling) = sibling_of(cx, base, scrut)
                        && let Some(value) = pat_value(cx, only.pat)
                    {
                        known(cx, sibling, value, true, out);
                    }
                }
                _ => {}
            },
            StmtKind::Item(_) => {}
        }
    }
}

/// Every sibling test that dominates `read` (a `<base>.f` expression):
/// enclosing `if` conditions and `match` arms on `<base>.g`, and guards
/// that ran earlier in an enclosing block.
fn dominating_tests(cx: &LateContext<'_>, read: &Expr<'_>, base: &Expr<'_>) -> HashSet<Test> {
    let mut out = HashSet::new();
    let mut child = read.hir_id;
    let mut via_arm: Option<&Arm<'_>> = None;
    for (id, node) in cx.tcx.hir_parent_iter(read.hir_id) {
        match node {
            Node::Expr(e) => match e.kind {
                ExprKind::If(cond, then, els) => {
                    if then.hir_id == child {
                        cond_tests(cx, base, cond, true, &mut out);
                    } else if els.is_some_and(|x| x.hir_id == child) {
                        cond_tests(cx, base, cond, false, &mut out);
                    }
                }
                ExprKind::Match(scrut, ..) => {
                    if let Some(arm) = via_arm.take()
                        && let Some(sibling) = sibling_of(cx, base, scrut)
                        && let Some(value) = pat_value(cx, arm.pat)
                    {
                        known(cx, sibling, value, true, &mut out);
                    }
                }
                _ => {}
            },
            Node::Arm(arm) => via_arm = Some(arm),
            Node::Block(b) => preceding_guards(cx, base, b, child, &mut out),
            Node::Item(_) | Node::ImplItem(_) | Node::TraitItem(_) | Node::ForeignItem(_) => break,
            _ => {}
        }
        child = id;
    }
    out
}

/// The `<base>.g = <value>` assignments among the statements of the block
/// enclosing `at`: writes that put the struct into a case at the same step.
fn assigned_beside(cx: &LateContext<'_>, base: &Expr<'_>, at: HirId, out: &mut HashSet<Test>) {
    let Some(block) = get_enclosing_block(cx, at) else {
        return;
    };
    for stmt in block.stmts {
        if let StmtKind::Semi(e) | StmtKind::Expr(e) = stmt.kind
            && let ExprKind::Assign(lhs, rhs, _) = e.kind
            && let Some(sibling) = sibling_of(cx, base, lhs)
            && let Some(value) = expr_value(cx, rhs)
        {
            out.insert(Test { sibling, value });
        }
    }
}

/// `x.g += ..`, `&mut x.g`, or the auto-`&mut` a mutating method call takes:
/// `g` changes to something not spelled out.
fn mutated_in_place(cx: &LateContext<'_>, place: &Expr<'_>, parent: Node<'_>) -> bool {
    if let Node::Expr(p) = parent
        && let ExprKind::AssignOp(_, lhs, _)
        | ExprKind::AddrOf(BorrowKind::Ref | BorrowKind::Raw, Mutability::Mut, lhs) = p.kind
        && lhs.hir_id == place.hir_id
    {
        return true;
    }
    cx.typeck_results().expr_adjustments(place).iter().any(|a| {
        matches!(
            a.kind,
            Adjust::Borrow(AutoBorrow::Ref(AutoBorrowMutability::Mut { .. }))
                | Adjust::Borrow(AutoBorrow::RawPtr(Mutability::Mut))
        )
    })
}

impl FieldValidOnlyWhen {
    fn record_read(&mut self, adt: DefId, field: Symbol, at: Span, tests: HashSet<Test>) {
        self.reads
            .entry((adt, field))
            .or_default()
            .push(Read { at, tests });
    }
}

impl<'tcx> LateLintPass<'tcx> for FieldValidOnlyWhen {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        match expr.kind {
            ExprKind::Struct(_, fields, tail) => {
                let Some(adt) = relevant(cx, cx.typeck_results().expr_ty(expr)) else {
                    return;
                };
                // What `..` fills the unlisted fields with: don't-cares from
                // `..Default::default()` or declared defaults, somebody's real
                // values from `..other`.
                let rest = match tail {
                    StructTailExpr::None => None,
                    StructTailExpr::Base(b) => Some(is_placeholder(cx, b)),
                    StructTailExpr::DefaultFields(_) => Some(true),
                    StructTailExpr::NoneWithError(_) => return,
                };
                let mut site: HashMap<Symbol, Init> = fields
                    .iter()
                    .map(|f| {
                        let init = Init {
                            value: expr_value(cx, f.expr),
                            placeholder: is_placeholder(cx, f.expr),
                        };
                        (f.ident.name, init)
                    })
                    .collect();
                if let Some(placeholder) = rest {
                    for f in &adt.non_enum_variant().fields {
                        site.entry(f.name).or_insert(Init {
                            value: None,
                            placeholder,
                        });
                    }
                }
                self.sites.entry(adt.did()).or_default().push(site);
            }
            ExprKind::Field(base, ident) => {
                let Some(adt) = relevant(cx, cx.typeck_results().expr_ty_adjusted(base)) else {
                    return;
                };
                if struct_field(adt, ident.name).is_none()
                    || in_automatically_derived(cx.tcx, expr.hir_id)
                {
                    return;
                }
                // `base.f = ..` writes; everything else, `base.f += ..` and
                // `&mut base.f` included, reads (and, through the `&mut`,
                // may leave behind a value no site spells out).
                let parent = cx.tcx.parent_hir_node(expr.hir_id);
                if let Node::Expr(parent) = parent
                    && let ExprKind::Assign(lhs, rhs, _) = parent.kind
                    && lhs.hir_id == expr.hir_id
                {
                    self.assigned
                        .entry((adt.did(), ident.name))
                        .or_default()
                        .insert(expr_value(cx, rhs));
                    if !is_placeholder(cx, rhs) {
                        let mut cases = dominating_tests(cx, expr, base);
                        assigned_beside(cx, base, parent.hir_id, &mut cases);
                        self.writes
                            .entry((adt.did(), ident.name))
                            .or_default()
                            .push(cases);
                    }
                    return;
                }
                if mutated_in_place(cx, expr, parent) {
                    self.assigned
                        .entry((adt.did(), ident.name))
                        .or_default()
                        .insert(None);
                }
                let tests = dominating_tests(cx, expr, base);
                self.record_read(adt.did(), ident.name, expr.span, tests);
            }
            _ => {}
        }
    }

    fn check_pat(&mut self, cx: &LateContext<'tcx>, pat: &'tcx Pat<'tcx>) {
        let PatKind::Struct(_, fields, _) = pat.kind else {
            return;
        };
        let Some(typeck) = cx.maybe_typeck_results() else {
            return;
        };
        let Some(adt) = relevant(cx, typeck.pat_ty(pat)) else {
            return;
        };
        if in_automatically_derived(cx.tcx, pat.hir_id) {
            return;
        }
        for f in fields {
            self.record_read(adt.did(), f.ident.name, f.span, HashSet::new());
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let mut findings: Vec<(Span, String, Span, String)> = Vec::new();
        for ((did, field), reads) in &self.reads {
            // The tests every read agrees on; one bare read empties it.
            let mut common: Option<HashSet<Test>> = None;
            for Read { tests, .. } in reads {
                common = Some(match common {
                    None => tests.clone(),
                    Some(c) => c.intersection(tests).copied().collect(),
                });
            }
            let Some(common) = common else { continue };
            let adt = cx.tcx.adt_def(*did);
            let order: Vec<Symbol> = adt
                .non_enum_variant()
                .fields
                .iter()
                .map(|f| f.name)
                .collect();
            let mut candidates: Vec<Test> =
                common.into_iter().filter(|t| t.sibling != *field).collect();
            candidates.sort_by_key(|t| order.iter().position(|n| *n == t.sibling));
            let sites = self.sites.get(did).map_or(&[][..], Vec::as_slice);
            for test in candidates {
                // The tested case must be one the crate makes: a literal site
                // or a later assignment giving the sibling that value or one
                // not spelled out. Otherwise the field is never read at all,
                // which is not this lint's claim.
                let admits = |v: Option<Value>| v.is_none_or(|v| v == test.value);
                let reached = self
                    .assigned
                    .get(&(*did, test.sibling))
                    .is_some_and(|vs| vs.iter().copied().any(admits))
                    || sites
                        .iter()
                        .any(|s| s.get(&test.sibling).is_some_and(|i| admits(i.value)));
                if !reached {
                    continue;
                }
                // Every real value assigned later goes in under the test or
                // together with the sibling being set to that case.
                let writes = self
                    .writes
                    .get(&(*did, *field))
                    .map_or(&[][..], Vec::as_slice);
                if !writes.iter().all(|cases| cases.contains(&test)) {
                    continue;
                }
                // Sites that give the sibling some other spelled-out value.
                let elsewhere = sites.iter().filter(|s| {
                    s.get(&test.sibling)
                        .and_then(|i| i.value)
                        .is_some_and(|v| v != test.value)
                });
                let (mut placeholders, mut real) = (0usize, false);
                for site in elsewhere {
                    match site.get(field) {
                        Some(init) if init.placeholder => placeholders += 1,
                        Some(_) => real = true,
                        None => {}
                    }
                }
                if placeholders == 0 || real {
                    continue;
                }
                let Some(fdef) = struct_field(adt, *field) else {
                    continue;
                };
                let sibling = test.sibling;
                let (holds, case) = match test.value {
                    Value::Variant(v) => {
                        let variant = cx.tcx.item_name(v);
                        (
                            format!(
                                "{sibling} == {}::{variant}",
                                cx.tcx.item_name(cx.tcx.parent(v))
                            ),
                            format!("`{variant}`"),
                        )
                    }
                    Value::Bool(true) => (sibling.to_string(), "`true`".to_owned()),
                    Value::Bool(false) => (format!("!{sibling}"), "`false`".to_owned()),
                    Value::Int(n) => (format!("{sibling} == {n}"), format!("`{n}`")),
                };
                let first_read = reads.iter().map(|r| r.at).min_by_key(|at| at.lo());
                let Some(first_read) = first_read else {
                    continue;
                };
                findings.push((
                    cx.tcx.def_span(fdef.did),
                    format!(
                        "`{field}` is only read after checking `{holds}` ({} read{}), and every other `{}` fills it with a placeholder ({placeholders} place{}). The field only means something in that one case",
                        reads.len(),
                        if reads.len() == 1 { "" } else { "s" },
                        cx.tcx.item_name(*did),
                        if placeholders == 1 { "" } else { "s" },
                    ),
                    first_read,
                    format!(
                        "make the {case} case of `{sibling}` an enum variant that carries `{field}`, and drop the flat field"
                    ),
                ));
                break;
            }
        }
        findings.sort_by_key(|(span, ..)| span.lo());
        for (span, msg, read, help) in findings {
            emit_with_note(
                cx,
                FIELD_VALID_ONLY_WHEN,
                span,
                msg,
                read,
                "one of the guarded reads",
                help,
            );
        }
    }
}

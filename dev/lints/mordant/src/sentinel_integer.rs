use std::collections::{HashMap, HashSet};

use clippy_utils::higher::Range;
use clippy_utils::macros::{find_assert_eq_args, root_macro_call_first_node};
use clippy_utils::res::MaybeResPath;
use rustc_ast::LitKind;
use rustc_hir::def::{DefKind, Res};
use rustc_hir::def_id::DefId;
use rustc_hir::{
    BinOpKind, Expr, ExprKind, HirId, LetStmt, Pat, PatExpr, PatExprKind, PatKind, QPath, UnOp,
};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::{self, Ty};
use rustc_span::{Span, Symbol, sym};

use crate::adt_facts::{field_ty, is_option_ty, struct_field};
use crate::baseline::emit_with_note;
use crate::hir_shapes::{assigned_field, callee_of, peel_blocks_unsafe};

rustc_session::declare_lint! {
    /// Flags an integer struct field that can hold a sentinel (some
    /// function compares it `==`/`!=` against `T::MAX`, `-1`, or a constant
    /// named `INVALID`/`NONE`/`SENTINEL`, treating that value as "no value")
    /// and that another function indexes with (`v[x.f as usize]`, a slice
    /// range end, `buf[off..off + len]`, `get_unchecked`) or offsets a pointer
    /// by without checking. One reader knows the magic value means "none".
    /// The other turns it into an out-of-range index or a wild offset. The
    /// type is `u32` when the value set is `Option<u32>`, and only convention
    /// tells the readers apart.
    ///
    /// Reported on the unchecked reader. A function counts as checking the
    /// field if it compares the field, or a local read off it, against
    /// anything (`==`, `!=`, `assert_ne!`, an ordering test against a length,
    /// a `match` arm or `matches!` with a literal or constant pattern), clamps
    /// it (`min`, `checked_add`, ..), looks it up with `.get(i)`, or directly
    /// calls a predicate (a `bool`-returning function) that does; a function
    /// all of whose visible callers check is their unchecked half and stays
    /// quiet too. `.get(i)` and keyed lookups (`map.remove(&x.f)`, `map[&x.f]`)
    /// already answer for a value that is not there and are not uses.
    /// `wrapping_*` arithmetic opts out of any range decision and is neither
    /// a check nor, on an integer, a use. Plain arithmetic on the field is not
    /// a use: positions and lengths are summed everywhere and the sum is only
    /// wrong where it meets memory. A field only ever *assigned* `MAX` and
    /// never compared to it is a bound, not a missing value, and is left
    /// alone.
    pub SENTINEL_INTEGER,
    Warn,
    "an integer field compared to a sentinel by one reader and indexed with unchecked by another"
}

/// A struct (local or not) and one of its integer fields.
type Field = (DefId, Symbol);

#[derive(Default)]
struct Evidence {
    /// How the first comparison seen spells the sentinel.
    spelling: String,
    compared: usize,
    /// Where that first comparison is.
    at: Option<Span>,
}

#[derive(Clone, Copy)]
enum Use {
    Index,
    Offset,
}

struct Read {
    body: DefId,
    span: Span,
    how: Use,
}

#[derive(Default)]
pub struct SentinelInteger {
    evidence: HashMap<Field, Evidence>,
    reads: HashMap<Field, Vec<Read>>,
    /// The function tests the field, or a value read off it, somewhere.
    checked: HashSet<(DefId, Field)>,
    /// `let i = x.f as usize + y.g`: locals that carry fields' values.
    locals: HashMap<HirId, Vec<Field>>,
    /// Function -> local `bool`-returning functions it calls directly.
    calls: HashMap<DefId, HashSet<DefId>>,
    /// Local function -> functions that call it directly.
    callers: HashMap<DefId, HashSet<DefId>>,
    /// Local functions referenced other than by a direct call: their caller
    /// set is unknowable.
    poisoned: HashSet<DefId>,
}

rustc_session::impl_lint_pass!(SentinelInteger => [SENTINEL_INTEGER]);

enum Sentinel {
    Max,
    MinusOne,
    Named(DefId),
}

/// `INVALID_ID`, `Slot::NONE`, `NOT_SET_SENTINEL`: a constant whose name says
/// the value stands for no value.
fn names_absence(name: &str) -> bool {
    name.split('_')
        .any(|w| matches!(w, "INVALID" | "NONE" | "SENTINEL"))
}

fn is_minus_one(lit: LitKind) -> bool {
    matches!(lit, LitKind::Int(v, _) if v.get() == 1)
}

/// A constant that spells a sentinel: `core`'s `MAX`, or one named for
/// absence.
fn sentinel_const(cx: &LateContext<'_>, res: Res) -> Option<Sentinel> {
    let Res::Def(DefKind::Const { .. } | DefKind::AssocConst { .. }, did) = res else {
        return None;
    };
    let name = cx.tcx.item_name(did);
    if name.as_str() == "MAX" && cx.tcx.crate_name(did.krate).as_str() == "core" {
        Some(Sentinel::Max)
    } else if names_absence(name.as_str()) {
        Some(Sentinel::Named(did))
    } else {
        None
    }
}

/// The sentinel an expression spells, if it is one of the three forms.
fn sentinel_of<'tcx>(cx: &LateContext<'tcx>, e: &'tcx Expr<'tcx>) -> Option<Sentinel> {
    let e = peel_blocks_unsafe(e);
    if !cx.typeck_results().expr_ty(e).is_integral() {
        return None;
    }
    match e.kind {
        ExprKind::Unary(UnOp::Neg, inner) => match peel_blocks_unsafe(inner).kind {
            ExprKind::Lit(lit) if is_minus_one(lit.node) => Some(Sentinel::MinusOne),
            _ => None,
        },
        ExprKind::Path(ref qpath) => sentinel_const(cx, cx.qpath_res(qpath, e.hir_id)),
        _ => None,
    }
}

/// The sentinel a `match` arm pattern spells: `-1` or a constant path.
fn pat_sentinel(cx: &LateContext<'_>, pe: &PatExpr<'_>) -> Option<Sentinel> {
    match pe.kind {
        PatExprKind::Lit { lit, negated: true } if is_minus_one(lit.node) => {
            Some(Sentinel::MinusOne)
        }
        PatExprKind::Path(ref qpath) => sentinel_const(cx, cx.qpath_res(qpath, pe.hir_id)),
        PatExprKind::Lit { .. } => None,
    }
}

/// The arm patterns that test the scrutinee's value — a literal, constant or
/// range — through `|`, `&` and guards.
fn value_pats<'a>(pat: &'a Pat<'a>, out: &mut Vec<&'a Pat<'a>>) {
    match pat.kind {
        PatKind::Expr(_) | PatKind::Range(..) => out.push(pat),
        PatKind::Or(alts) => {
            for alt in alts {
                value_pats(alt, out);
            }
        }
        PatKind::Ref(inner, _, _)
        | PatKind::Deref(inner)
        | PatKind::Box(inner)
        | PatKind::Guard(inner, _) => value_pats(inner, out),
        _ => {}
    }
}

fn spelling<'tcx>(cx: &LateContext<'tcx>, s: &Sentinel, ty: Ty<'tcx>) -> String {
    match s {
        Sentinel::Max => format!("{ty}::MAX"),
        Sentinel::MinusOne => "-1".to_owned(),
        Sentinel::Named(did) => cx.tcx.item_name(*did).to_string(),
    }
}

/// `base.name` as a struct and its integer field.
fn field_key(cx: &LateContext<'_>, base: &Expr<'_>, name: Symbol) -> Option<Field> {
    let adt = cx
        .typeck_results()
        .expr_ty_adjusted(base)
        .peel_refs()
        .ty_adt_def()?;
    if !adt.is_struct() {
        return None;
    }
    let f = struct_field(adt, name)?;
    field_ty(cx, f).is_integral().then_some((adt.did(), name))
}

/// The function an expression belongs to, with closures folded into the
/// function that wrote them: a check before a `.map(|..| v[x.f])` covers it.
fn owner_fn(cx: &LateContext<'_>, hir_id: HirId) -> DefId {
    let mut did = cx.tcx.hir_enclosing_body_owner(hir_id).to_def_id();
    while cx.tcx.is_closure_like(did) || matches!(cx.tcx.def_kind(did), DefKind::InlineConst) {
        did = cx.tcx.parent(did);
    }
    did
}

/// Memory a `usize` positions into, where the INDEXERS calls panic or are UB
/// past the end.
fn is_contiguous<'tcx>(cx: &LateContext<'tcx>, ty: Ty<'tcx>) -> bool {
    match ty.kind() {
        ty::Slice(_) | ty::Array(..) | ty::Str => true,
        ty::Adt(adt, _) => {
            cx.tcx.is_diagnostic_item(sym::Vec, adt.did())
                || cx.tcx.is_diagnostic_item(sym::String, adt.did())
        }
        _ => false,
    }
}

/// An INDEXERS call that does not answer for an index that is not there:
/// on contiguous memory, or on anything else unless its result is the
/// `Option`/`bool` a keyed collection (`map.remove`, `deque.remove`) hands
/// back for an absent key.
fn indexes_positionally<'tcx>(
    cx: &LateContext<'tcx>,
    call: &'tcx Expr<'tcx>,
    recv: &'tcx Expr<'tcx>,
) -> bool {
    if is_contiguous(cx, cx.typeck_results().expr_ty_adjusted(recv).peel_refs()) {
        return true;
    }
    let out = cx.typeck_results().expr_ty(call);
    !(out.is_bool() || out.is_unit() || is_option_ty(cx, out))
}

/// Value-preserving wrappers a field read is still visible through.
const ADAPTERS: &[&str] = &[
    "clone",
    "into",
    "try_into",
    "unwrap",
    "expect",
    "cast_signed",
    "cast_unsigned",
];

/// Calls that index contiguous memory by their one argument and do not
/// answer for an index that is not there.
const INDEXERS: &[&str] = &[
    "get_unchecked",
    "get_unchecked_mut",
    "split_at",
    "split_at_mut",
    "split_off",
    "remove",
    "swap_remove",
];

const OFFSETS: &[&str] = &[
    "add",
    "sub",
    "offset",
    "byte_add",
    "byte_sub",
    "byte_offset",
    "wrapping_add",
    "wrapping_sub",
    "wrapping_offset",
    "wrapping_byte_add",
    "wrapping_byte_sub",
    "wrapping_byte_offset",
];

impl SentinelInteger {
    /// Every field whose value `e` carries: the field itself through casts,
    /// borrows, derefs and value-preserving adapters, a local bound to one,
    /// or both operands of a sum.
    fn reads_of<'tcx>(&self, cx: &LateContext<'tcx>, e: &'tcx Expr<'tcx>) -> Vec<Field> {
        let mut out = Vec::new();
        self.collect_reads(cx, e, &mut out);
        out
    }

    fn collect_reads<'tcx>(
        &self,
        cx: &LateContext<'tcx>,
        mut e: &'tcx Expr<'tcx>,
        out: &mut Vec<Field>,
    ) {
        loop {
            e = peel_blocks_unsafe(e);
            match e.kind {
                ExprKind::Cast(inner, _)
                | ExprKind::AddrOf(_, _, inner)
                | ExprKind::Unary(UnOp::Deref, inner) => e = inner,
                ExprKind::MethodCall(seg, recv, args, _)
                    if args.len() <= 1 && ADAPTERS.contains(&seg.ident.as_str()) =>
                {
                    e = recv;
                }
                // `usize::from(x.f)`, `u32::try_from(x.f)`.
                ExprKind::Call(callee, [arg])
                    if matches!(callee.kind, ExprKind::Path(QPath::TypeRelative(_, seg))
                        if matches!(seg.ident.as_str(), "from" | "try_from")) =>
                {
                    e = arg;
                }
                // `off + len`, `idx - 1`: the sum still carries the sentinel
                // of either side.
                ExprKind::Binary(op, l, r)
                    if matches!(op.node, BinOpKind::Add | BinOpKind::Sub | BinOpKind::Mul) =>
                {
                    self.collect_reads(cx, l, out);
                    e = r;
                }
                ExprKind::Field(base, ident) => {
                    out.extend(field_key(cx, base, ident.name));
                    return;
                }
                ExprKind::Path(_) => {
                    if let Some(id) = e.res_local_id()
                        && let Some(fields) = self.locals.get(&id)
                    {
                        out.extend_from_slice(fields);
                    }
                    return;
                }
                _ => return,
            }
        }
    }

    fn compared<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        field: Field,
        s: &Sentinel,
        ty: Ty<'tcx>,
        at: Span,
    ) {
        let ev = self.evidence.entry(field).or_default();
        if ev.compared == 0 {
            ev.spelling = spelling(cx, s, ty);
            ev.at = Some(at);
        }
        ev.compared += 1;
    }

    /// The enclosing function tests every field `e` carries; those fields.
    fn tested<'tcx>(&mut self, cx: &LateContext<'tcx>, e: &'tcx Expr<'tcx>) -> Vec<Field> {
        let fields = self.reads_of(cx, e);
        if !fields.is_empty() {
            let body = owner_fn(cx, e.hir_id);
            for f in &fields {
                self.checked.insert((body, *f));
            }
        }
        fields
    }

    /// `l == r` / `l != r`, spelled as an operator or an `assert_eq!`.
    fn equated<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        at: Span,
        l: &'tcx Expr<'tcx>,
        r: &'tcx Expr<'tcx>,
    ) {
        for (a, b) in [(l, r), (r, l)] {
            let fields = self.tested(cx, a);
            if let Some(s) = sentinel_of(cx, b) {
                let ty = cx.typeck_results().expr_ty(peel_blocks_unsafe(b));
                for f in fields {
                    self.compared(cx, f, &s, ty, at);
                }
            }
        }
    }

    /// `match scrut { CONST => .., 0..=9 => .. }`, `matches!(scrut, -1)`,
    /// `if let CONST = scrut`: an arm that names a value tests the scrutinee.
    fn matched<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        scrut: &'tcx Expr<'tcx>,
        pats: impl Iterator<Item = &'tcx Pat<'tcx>>,
    ) {
        let mut leaves = Vec::new();
        for pat in pats {
            value_pats(pat, &mut leaves);
        }
        if leaves.is_empty() {
            return;
        }
        let fields = self.tested(cx, scrut);
        if fields.is_empty() {
            return;
        }
        for leaf in leaves {
            if let PatKind::Expr(pe) = leaf.kind
                && let Some(s) = pat_sentinel(cx, pe)
            {
                let ty = cx.typeck_results().node_type(leaf.hir_id);
                for f in &fields {
                    self.compared(cx, *f, &s, ty, leaf.span);
                }
            }
        }
    }

    fn read<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        operand: &'tcx Expr<'tcx>,
        span: Span,
        how: Use,
    ) {
        // A position is an integer by value; `&x.f` is a key being looked up.
        if !cx.typeck_results().expr_ty(operand).is_integral() {
            return;
        }
        let body = owner_fn(cx, operand.hir_id);
        for field in self.reads_of(cx, operand) {
            self.reads
                .entry(field)
                .or_default()
                .push(Read { body, span, how });
        }
    }

    /// An index operand: the value itself, or either end of a range.
    fn indexed<'tcx>(&mut self, cx: &LateContext<'tcx>, idx: &'tcx Expr<'tcx>, span: Span) {
        match Range::hir(cx, idx) {
            Some(range) => {
                for end in [range.start, range.end].into_iter().flatten() {
                    self.read(cx, end, span, Use::Index);
                }
            }
            None => self.read(cx, idx, span, Use::Index),
        }
    }

    fn checks(&self, body: DefId, field: Field) -> bool {
        self.checked.contains(&(body, field))
            || self
                .calls
                .get(&body)
                .is_some_and(|cs| cs.iter().any(|c| self.checked.contains(&(*c, field))))
    }

    /// Every visible caller of `body` checks the field first: `body` is the
    /// unchecked half of a checked pair, not an unchecked reader.
    fn callers_check(&self, body: DefId, field: Field) -> bool {
        if self.poisoned.contains(&body) {
            return false;
        }
        self.callers
            .get(&body)
            .is_some_and(|cs| !cs.is_empty() && cs.iter().all(|c| self.checks(*c, field)))
    }

    fn record_call<'tcx>(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        match callee_of(cx, expr) {
            Some(callee) => {
                let def = callee.def();
                if def.is_local() && matches!(cx.tcx.def_kind(def), DefKind::Fn | DefKind::AssocFn)
                {
                    let body = owner_fn(cx, expr.hir_id);
                    self.callers.entry(def).or_default().insert(body);
                    // Only a predicate (`is_root()`, `has_parent()`) stands
                    // in for a comparison; a call that happens to compare
                    // inside says nothing about the caller's own reads.
                    let returns_bool = cx
                        .tcx
                        .fn_sig(def)
                        .instantiate_identity()
                        .skip_normalization()
                        .output()
                        .skip_binder()
                        .is_bool();
                    if returns_bool {
                        self.calls.entry(body).or_default().insert(def);
                    }
                }
            }
            None => {
                let ExprKind::Path(qpath) = &expr.kind else {
                    return;
                };
                if matches!(
                    clippy_utils::get_parent_expr(cx, expr),
                    Some(Expr { kind: ExprKind::Call(callee, _), .. }) if callee.hir_id == expr.hir_id
                ) {
                    return;
                }
                if let Res::Def(DefKind::Fn | DefKind::AssocFn, def) =
                    cx.qpath_res(qpath, expr.hir_id)
                    && def.is_local()
                {
                    self.poisoned.insert(def);
                }
            }
        }
    }

    fn assert_eq_args<'tcx>(
        cx: &LateContext<'tcx>,
        expr: &'tcx Expr<'tcx>,
    ) -> Option<(&'tcx Expr<'tcx>, &'tcx Expr<'tcx>)> {
        let mac = root_macro_call_first_node(cx, expr)?;
        if !matches!(
            cx.tcx.get_diagnostic_name(mac.def_id),
            Some(
                sym::assert_eq_macro
                    | sym::assert_ne_macro
                    | sym::debug_assert_eq_macro
                    | sym::debug_assert_ne_macro
            )
        ) {
            return None;
        }
        find_assert_eq_args(cx, expr, mac.expn).map(|(l, r, _)| (l, r))
    }
}

impl<'tcx> LateLintPass<'tcx> for SentinelInteger {
    fn check_local(&mut self, cx: &LateContext<'tcx>, local: &'tcx LetStmt<'tcx>) {
        if let Some(init) = local.init
            && let PatKind::Binding(_, id, _, None) = local.pat.kind
        {
            let fields = self.reads_of(cx, init);
            if !fields.is_empty() {
                self.locals.insert(id, fields);
            }
        }
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        self.record_call(cx, expr);
        // `assert_eq!(a, b)` compares through match-arm bindings the operator
        // arm below cannot see back through; read its arguments directly.
        if let Some((l, r)) = Self::assert_eq_args(cx, expr) {
            self.equated(cx, expr.span, l, r);
        }
        match expr.kind {
            ExprKind::Struct(_, fields, _) => {
                let Some(adt) = cx.typeck_results().expr_ty(expr).ty_adt_def() else {
                    return;
                };
                if !adt.is_struct() {
                    return;
                }
                for init in fields {
                    if let Some(s) = sentinel_of(cx, init.expr)
                        && let Some(f) = struct_field(adt, init.ident.name)
                        && field_ty(cx, f).is_integral()
                    {
                        // A literal that writes the sentinel is where the
                        // spelling is clearest, but it is not a check.
                        let ev = self
                            .evidence
                            .entry((adt.did(), init.ident.name))
                            .or_default();
                        if ev.spelling.is_empty() {
                            ev.spelling = spelling(cx, &s, cx.typeck_results().expr_ty(init.expr));
                        }
                    }
                }
            }
            ExprKind::Assign(place, val, _) => {
                if let Some((base, ident, _)) = assigned_field(place)
                    && let Some(field) = field_key(cx, base, ident.name)
                    && let Some(s) = sentinel_of(cx, val)
                {
                    let ev = self.evidence.entry(field).or_default();
                    if ev.spelling.is_empty() {
                        ev.spelling = spelling(cx, &s, cx.typeck_results().expr_ty(val));
                    }
                }
            }
            ExprKind::Binary(op, l, r) => match op.node {
                BinOpKind::Eq | BinOpKind::Ne => self.equated(cx, expr.span, l, r),
                BinOpKind::Lt | BinOpKind::Le | BinOpKind::Gt | BinOpKind::Ge => {
                    self.tested(cx, l);
                    self.tested(cx, r);
                }
                _ => {}
            },
            ExprKind::Match(scrut, arms, _) => self.matched(cx, scrut, arms.iter().map(|a| a.pat)),
            ExprKind::Let(l) => self.matched(cx, l.init, std::iter::once(l.pat)),
            ExprKind::Index(_, idx, _) => self.indexed(cx, idx, expr.span),
            ExprKind::MethodCall(seg, recv, args, _) => {
                let name = seg.ident.as_str();
                // The author deciding what an out-of-range value does:
                // overflow-aware arithmetic on it, a clamp of it, or a lookup
                // that answers for absence. `wrapping_*` decides nothing.
                let bounded = ["checked_", "saturating_", "overflowing_"]
                    .iter()
                    .any(|p| name.starts_with(p))
                    || matches!(
                        name,
                        "min" | "max" | "clamp" | "get" | "get_mut" | "contains" | "contains_key"
                    );
                if bounded {
                    for operand in std::iter::once(recv).chain(args.iter()) {
                        self.tested(cx, operand);
                    }
                }
                if let [arg] = args {
                    if INDEXERS.contains(&name) && indexes_positionally(cx, expr, recv) {
                        self.indexed(cx, arg, expr.span);
                    } else if OFFSETS.contains(&name)
                        && cx.typeck_results().expr_ty_adjusted(recv).is_raw_ptr()
                    {
                        self.read(cx, arg, expr.span, Use::Offset);
                    }
                }
            }
            _ => {}
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let mut findings: Vec<(Span, String, Span, String, String)> = Vec::new();
        for (field, reads) in &self.reads {
            let Some(ev) = self.evidence.get(field) else {
                continue;
            };
            let Some(evidence_at) = ev.at else {
                continue;
            };
            // One report per unchecked function, at its first use.
            let mut first: HashMap<DefId, &Read> = HashMap::new();
            for read in reads {
                first
                    .entry(read.body)
                    .and_modify(|r| {
                        if read.span.lo() < r.span.lo() {
                            *r = read;
                        }
                    })
                    .or_insert(read);
            }
            for (body, read) in first {
                if self.checks(body, *field) || self.callers_check(body, *field) {
                    continue;
                }
                let (how, becomes, use_it) = match read.how {
                    Use::Index => ("indexes with", "an out-of-range index", "index"),
                    Use::Offset => ("offsets a pointer by", "a wild offset", "offset"),
                };
                let reader = match cx.tcx.opt_item_name(body) {
                    Some(name) => format!("`{name}`"),
                    None => "This body".to_owned(),
                };
                let (owner, name, sentinel) = (cx.tcx.item_name(field.0), field.1, &ev.spelling);
                findings.push((
                    read.span,
                    format!(
                        "`{owner}.{name}` can be `{sentinel}`, which {} other place{} treat{} as \"no value\". {reader} {how} it here without checking, so the sentinel becomes {becomes}",
                        ev.compared,
                        if ev.compared == 1 { "" } else { "s" },
                        if ev.compared == 1 { "s" } else { "" },
                    ),
                    evidence_at,
                    format!("one of the places that checks for `{sentinel}`"),
                    format!(
                        "store `{name}` as an `Option` (over `NonZero` or `NonMax` to keep the size). {reader} then has to handle the empty case before it can {use_it}"
                    ),
                ));
            }
        }
        findings.sort_by_key(|(span, ..)| span.lo());
        findings.dedup_by_key(|(span, ..)| *span);
        for (span, msg, evidence_at, note, help) in findings {
            emit_with_note(cx, SENTINEL_INTEGER, span, msg, evidence_at, note, help);
        }
    }
}

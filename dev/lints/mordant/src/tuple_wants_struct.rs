use std::collections::{HashMap, HashSet};
use std::ops::ControlFlow;

use clippy_utils::visitors::for_each_expr_without_closures;
use rustc_hir::def::{DefKind, Res};
use rustc_hir::def_id::DefId;
use rustc_hir::intravisit::FnKind;
use rustc_hir::{
    Body, Expr, ExprKind, FnDecl, HirId, LangItem, LetStmt, LocalSource, MatchSource, Node, Pat,
    PatKind, QPath, StmtKind,
};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::{self, AssocContainer, Ty, TyCtxt};
use rustc_span::def_id::LocalDefId;
use rustc_span::{Ident, Span, Symbol, sym};

use crate::baseline::{emit_with_note, join};
use crate::hir_shapes::{callee_of, peel_blocks_unsafe};

rustc_session::declare_lint! {
    /// Flags a crate-private function returning a tuple (bare, or inside an
    /// `Option`/`Result`) with two members of one type, when every call in
    /// the crate unpacks it on the spot (`let (a, b) = f()`, an
    /// `if let`/`match` arm, or `(a, self.b) = f()`) and all of them, across
    /// at least two functions, give the members the same names. The names
    /// exist everywhere except in the type, and since the members share a
    /// type a swapped pair compiles. When the function's own body already
    /// returns them under those names in another order, the finding says
    /// where. A struct with those fields puts the names in the signature.
    ///
    /// Silent when the member types all differ (the compiler already rejects
    /// a transposition), when any caller keeps the tuple whole (stores,
    /// returns, passes or indexes it with `.0`), binds a member as `_` or
    /// under a different name, when fewer than two functions call it, and
    /// when the function is exported, referenced as a value, or has its
    /// signature fixed by a trait.
    pub TUPLE_WANTS_STRUCT,
    Warn,
    "a tuple with same-typed members that every caller destructures under the same names"
}

/// A `(a, b)` / `Some((a, b))` in return position whose elements are all
/// plain locals or fields, as their names.
struct Returned {
    names: Vec<Symbol>,
    span: Span,
}

/// A function returning a tuple with two members of one type.
struct Candidate {
    ret_span: Span,
    name: Symbol,
    /// "`.0` and `.2` are both `f32`".
    same: String,
    returns: Vec<Returned>,
}

#[derive(Default)]
struct Sites {
    /// The member names every destructuring site so far agrees on.
    names: Option<Vec<Symbol>>,
    /// The items the destructuring sites sit in.
    owners: HashSet<LocalDefId>,
    count: usize,
    /// The first destructuring call, for the note.
    first: Option<Span>,
    /// A site kept the tuple whole, left a member unnamed, or disagreed.
    opaque: bool,
}

#[derive(Default)]
pub struct TupleWantsStruct {
    candidates: HashMap<DefId, Candidate>,
    sites: HashMap<DefId, Sites>,
    /// Referenced other than by a direct call: those callers are invisible.
    poisoned: HashSet<DefId>,
}

rustc_session::impl_lint_pass!(TupleWantsStruct => [TUPLE_WANTS_STRUCT]);

/// What one call site does with the returned tuple.
enum Use {
    /// Destructures it, naming every member.
    Names(Vec<Symbol>),
    /// Never reaches the members: drops the value, tests `is_some()`, takes
    /// the `None`/`Err` arm.
    Neutral,
    /// Anything else: the tuple survives whole or a member goes unnamed.
    Opaque,
}

/// A declared return type that is a tuple of two or more members.
struct TupleReturn<'tcx> {
    /// An `Option` or `Result` sits around the tuple.
    wrapped: bool,
    members: &'tcx [Ty<'tcx>],
}

/// Read off the signature, not a call's type, so a generic `T` that happens
/// to be a tuple at one call is not one.
fn tuple_return(tcx: TyCtxt<'_>, def: DefId) -> Option<TupleReturn<'_>> {
    let mut out = tcx
        .fn_sig(def)
        .instantiate_identity()
        .skip_normalization()
        .output()
        .skip_binder();
    let mut wrapped = false;
    if let ty::Adt(adt, args) = out.kind()
        && (tcx.is_diagnostic_item(sym::Option, adt.did())
            || tcx.is_diagnostic_item(sym::Result, adt.did()))
    {
        out = args.type_at(0);
        wrapped = true;
    }
    match out.kind() {
        ty::Tuple(tys) if tys.len() >= 2 => Some(TupleReturn {
            wrapped,
            members: tys.as_slice(),
        }),
        _ => None,
    }
}

/// The first two members of one type, as "both members are `T`" for a pair
/// and "`.i` and `.j` are both `T`" for anything longer.
fn same_typed_pair(members: &[Ty<'_>]) -> Option<String> {
    members.iter().enumerate().find_map(|(i, a)| {
        let j = i + 1 + members[i + 1..].iter().position(|b| a == b)?;
        Some(if members.len() == 2 {
            format!("both members are `{a}`")
        } else {
            format!("`.{i}` and `.{j}` are both `{a}`")
        })
    })
}

fn is_lang_ctor(cx: &LateContext<'_>, variant: Option<DefId>, items: &[LangItem]) -> bool {
    variant
        .and_then(|v| cx.tcx.lang_items().from_def_id(v))
        .is_some_and(|item| items.contains(&item))
}

/// The expressions a body evaluates to: every `return e` and every tail,
/// through blocks, `if` and `match`. A closure or `async` block's `return`
/// is that closure's value, not the body's.
fn return_values<'tcx>(body: &'tcx Body<'tcx>) -> Vec<&'tcx Expr<'tcx>> {
    fn tails<'h>(e: &'h Expr<'h>, out: &mut Vec<&'h Expr<'h>>) {
        match e.kind {
            ExprKind::Block(b, _) => {
                if let Some(t) = b.expr {
                    tails(t, out);
                }
            }
            ExprKind::DropTemps(inner) => tails(inner, out),
            ExprKind::If(_, then, els) => {
                tails(then, out);
                if let Some(els) = els {
                    tails(els, out);
                }
            }
            ExprKind::Match(_, arms, _) => {
                for arm in arms {
                    tails(arm.body, out);
                }
            }
            ExprKind::Ret(_) => {}
            _ => out.push(e),
        }
    }
    let mut out = Vec::new();
    for_each_expr_without_closures(body.value, |e| {
        if let ExprKind::Ret(Some(v)) = e.kind {
            out.push(v);
        }
        ControlFlow::<()>::Continue(())
    });
    tails(body.value, &mut out);
    out
}

/// `(a, self.b)` or, behind a wrapper, `Some((a, self.b))` / `Ok(..)`: the
/// names the body itself gives the members at this return, when every
/// element is a bare local or field.
fn returned_names(cx: &LateContext<'_>, e: &Expr<'_>, wrapped: bool) -> Option<Returned> {
    let mut e = peel_blocks_unsafe(e);
    if wrapped {
        let ExprKind::Call(_, [arg]) = e.kind else {
            return None;
        };
        let variant = crate::enum_facts::ctor_literal_variant(cx, e);
        if !is_lang_ctor(cx, variant, &[LangItem::OptionSome, LangItem::ResultOk]) {
            return None;
        }
        e = peel_blocks_unsafe(arg);
    }
    let ExprKind::Tup(elems) = e.kind else {
        return None;
    };
    let names = elems
        .iter()
        .map(|el| place_name(peel_blocks_unsafe(el)))
        .collect::<Option<Vec<_>>>()?;
    Some(Returned {
        names,
        span: e.span,
    })
}

/// A bare local `a` or a field `self.b`, as that name.
fn place_name(e: &Expr<'_>) -> Option<Symbol> {
    match e.kind {
        ExprKind::Path(QPath::Resolved(None, p)) if matches!(p.res, Res::Local(_)) => {
            Some(p.segments[0].ident.name)
        }
        ExprKind::Field(_, ident) => Some(ident.name),
        _ => None,
    }
}

/// Where a pattern's bindings get their names. rustc lowers
/// `(a, self.b) = f()` to `let (lhs, lhs) = f(); a = lhs; self.b = lhs;`:
/// there each binding is called `lhs` and its name is what it is assigned
/// on to.
#[derive(Clone, Copy)]
enum Naming<'a> {
    Bound,
    Assigned(&'a HashMap<HirId, Symbol>),
}

impl Naming<'_> {
    fn name(self, binding: HirId, ident: Ident) -> Option<Symbol> {
        match self {
            Naming::Bound => Some(ident.name),
            Naming::Assigned(targets) => targets.get(&binding).copied(),
        }
    }
}

/// The left-hand sides of a lowered `(a, self.b) = f()`, keyed by the `lhs`
/// binding each is assigned from; `None` when one of them is more than a
/// local or a field.
fn assign_targets(cx: &LateContext<'_>, desugared: &LetStmt<'_>) -> Option<HashMap<HirId, Symbol>> {
    let block = cx
        .tcx
        .hir_parent_iter(desugared.hir_id)
        .find_map(|(_, node)| match node {
            Node::Block(b) => Some(b),
            _ => None,
        })?;
    let mut targets = HashMap::new();
    for stmt in block.stmts {
        let (StmtKind::Expr(e) | StmtKind::Semi(e)) = stmt.kind else {
            continue;
        };
        let ExprKind::Assign(target, value, _) = e.kind else {
            continue;
        };
        let ExprKind::Path(QPath::Resolved(None, p)) = value.kind else {
            continue;
        };
        let Res::Local(binding) = p.res else {
            continue;
        };
        targets.insert(binding, place_name(target)?);
    }
    Some(targets)
}

fn pat_use(cx: &LateContext<'_>, pat: &Pat<'_>, wrapped: bool, naming: Naming<'_>) -> Use {
    match pat.kind {
        PatKind::Wild => Use::Neutral,
        PatKind::Ref(inner, _, _) => pat_use(cx, inner, wrapped, naming),
        PatKind::Tuple(pats, dotdot) if !wrapped => {
            if dotdot.as_opt_usize().is_some() {
                return Use::Opaque;
            }
            let mut names = Vec::with_capacity(pats.len());
            for p in pats {
                let PatKind::Binding(_, id, ident, None) = p.kind else {
                    return Use::Opaque;
                };
                let Some(name) = naming.name(id, ident) else {
                    return Use::Opaque;
                };
                if name.as_str().starts_with('_') {
                    return Use::Opaque;
                }
                names.push(name);
            }
            Use::Names(names)
        }
        PatKind::TupleStruct(_, [inner], _) if wrapped => {
            let variant = crate::enum_facts::arm_variant(cx, pat);
            if is_lang_ctor(cx, variant, &[LangItem::OptionSome, LangItem::ResultOk]) {
                pat_use(cx, inner, false, naming)
            } else if is_lang_ctor(cx, variant, &[LangItem::ResultErr]) {
                Use::Neutral
            } else {
                Use::Opaque
            }
        }
        PatKind::Expr(_)
            if wrapped
                && is_lang_ctor(
                    cx,
                    crate::enum_facts::arm_variant(cx, pat),
                    &[LangItem::OptionNone],
                ) =>
        {
            Use::Neutral
        }
        _ => Use::Opaque,
    }
}

/// Folds the uses of several arms over one value: any opaque arm makes the
/// whole opaque, and arms that name the members must agree.
fn merge(uses: impl Iterator<Item = Use>) -> Use {
    let mut acc = Use::Neutral;
    for u in uses {
        acc = match (acc, u) {
            (Use::Opaque, _) | (_, Use::Opaque) => return Use::Opaque,
            (Use::Neutral, u) | (u, Use::Neutral) => u,
            (Use::Names(a), Use::Names(b)) if a == b => Use::Names(a),
            (Use::Names(_), Use::Names(_)) => return Use::Opaque,
        };
    }
    acc
}

/// What the code around `call` does with its value: up through `?`,
/// `.unwrap()`, the `Option`/`Result` adapters that keep the tuple intact and
/// value-carrying blocks, to the pattern that receives it.
fn use_of_call<'tcx>(cx: &LateContext<'tcx>, call: &'tcx Expr<'tcx>, mut wrapped: bool) -> Use {
    let mut current = call.hir_id;
    for (parent_id, node) in cx.tcx.hir_parent_iter(call.hir_id) {
        match node {
            Node::LetStmt(l) if l.init.is_some_and(|i| i.hir_id == current) => {
                return match l.source {
                    LocalSource::Normal => pat_use(cx, l.pat, wrapped, Naming::Bound),
                    LocalSource::AssignDesugar => match assign_targets(cx, l) {
                        Some(targets) => pat_use(cx, l.pat, wrapped, Naming::Assigned(&targets)),
                        None => Use::Opaque,
                    },
                    _ => Use::Opaque,
                };
            }
            // `f();`: the value is dropped unread.
            Node::Stmt(s) if matches!(s.kind, StmtKind::Semi(e) if e.hir_id == current) => {
                return Use::Neutral;
            }
            Node::Block(b) if b.expr.is_some_and(|e| e.hir_id == current) => {}
            Node::Expr(parent) => match parent.kind {
                ExprKind::DropTemps(_) | ExprKind::Block(..) => {}
                // `call?`: `match Try::branch(call) { Continue(val) => val, .. }`.
                ExprKind::Call(_, [arg])
                    if wrapped
                        && arg.hir_id == current
                        && matches!(
                            cx.tcx.parent_hir_node(parent_id),
                            Node::Expr(Expr {
                                kind: ExprKind::Match(_, _, MatchSource::TryDesugar(_)),
                                ..
                            })
                        ) => {}
                ExprKind::Match(_, _, MatchSource::TryDesugar(_)) => wrapped = false,
                ExprKind::MethodCall(seg, recv, _, _) if wrapped && recv.hir_id == current => {
                    match seg.ident.as_str() {
                        "unwrap" | "expect" | "unwrap_unchecked" => wrapped = false,
                        // Still an `Option`/`Result` around the same tuple.
                        "ok" | "ok_or" | "ok_or_else" | "map_err" => {}
                        "is_some" | "is_none" | "is_ok" | "is_err" => return Use::Neutral,
                        _ => return Use::Opaque,
                    }
                }
                ExprKind::Match(scrut, arms, MatchSource::Normal | MatchSource::Postfix)
                    if scrut.hir_id == current =>
                {
                    return merge(
                        arms.iter()
                            .map(|arm| pat_use(cx, arm.pat, wrapped, Naming::Bound)),
                    );
                }
                ExprKind::Let(l) if l.init.hir_id == current => {
                    return pat_use(cx, l.pat, wrapped, Naming::Bound);
                }
                _ => return Use::Opaque,
            },
            _ => return Use::Opaque,
        }
        current = parent_id;
    }
    Use::Opaque
}

fn tuple_of(names: &[Symbol]) -> String {
    let names: Vec<&str> = names.iter().map(Symbol::as_str).collect();
    format!("({})", names.join(", "))
}

impl TupleWantsStruct {
    fn record(&mut self, cx: &LateContext<'_>, def: DefId, site: &Expr<'_>, used: Use) {
        let sites = self.sites.entry(def).or_default();
        match used {
            Use::Neutral => {}
            Use::Opaque => sites.opaque = true,
            Use::Names(names) => {
                sites.count += 1;
                sites.first.get_or_insert(site.span);
                sites
                    .owners
                    .insert(cx.tcx.hir_get_parent_item(site.hir_id).def_id);
                match &sites.names {
                    None => sites.names = Some(names),
                    Some(agreed) if *agreed == names => {}
                    Some(_) => sites.opaque = true,
                }
            }
        }
    }
}

impl<'tcx> LateLintPass<'tcx> for TupleWantsStruct {
    fn check_fn(
        &mut self,
        cx: &LateContext<'tcx>,
        kind: FnKind<'tcx>,
        decl: &'tcx FnDecl<'tcx>,
        body: &'tcx Body<'tcx>,
        span: Span,
        def_id: LocalDefId,
    ) {
        if matches!(kind, FnKind::Closure)
            || span.from_expansion()
            || cx.effective_visibilities.is_exported(def_id)
        {
            return;
        }
        let def = def_id.to_def_id();
        // A trait method's return type is the trait's to change, not this
        // function's.
        if let Some(assoc) = cx.tcx.opt_associated_item(def)
            && !matches!(assoc.container, AssocContainer::InherentImpl)
        {
            return;
        }
        let Some(ret) = tuple_return(cx.tcx, def) else {
            return;
        };
        let Some(same) = same_typed_pair(ret.members) else {
            return;
        };
        let returns = return_values(body)
            .into_iter()
            .filter_map(|e| returned_names(cx, e, ret.wrapped))
            .collect();
        self.candidates.insert(
            def,
            Candidate {
                ret_span: decl.output.span(),
                name: cx.tcx.item_name(def),
                same,
                returns,
            },
        );
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let Some(callee) = callee_of(cx, expr) else {
            // A bare reference to a local fn (fn pointer, higher-order use):
            // whatever destructures its result is out of sight.
            if let ExprKind::Path(qpath) = &expr.kind
                && !matches!(
                    clippy_utils::get_parent_expr(cx, expr),
                    Some(Expr { kind: ExprKind::Call(callee, _), .. }) if callee.hir_id == expr.hir_id
                )
                && let Res::Def(DefKind::Fn | DefKind::AssocFn, def) =
                    cx.qpath_res(qpath, expr.hir_id)
                && def.is_local()
            {
                self.poisoned.insert(def);
            }
            return;
        };
        let def = callee.def();
        if !def.is_local() || !matches!(cx.tcx.def_kind(def), DefKind::Fn | DefKind::AssocFn) {
            return;
        }
        let Some(ret) = tuple_return(cx.tcx, def) else {
            return;
        };
        let used = if expr.span.from_expansion() {
            Use::Opaque
        } else {
            use_of_call(cx, expr, ret.wrapped)
        };
        self.record(cx, def, expr, used);
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let mut findings: Vec<(Span, String, Span, &str, String)> = Vec::new();
        for (def, cand) in &self.candidates {
            if self.poisoned.contains(def) {
                continue;
            }
            let Some(sites) = self.sites.get(def) else {
                continue;
            };
            let (Some(names), Some(first)) = (&sites.names, sites.first) else {
                continue;
            };
            if sites.opaque || sites.count < 2 || sites.owners.len() < 2 {
                continue;
            }
            let lead = format!(
                "{} {} calls to `{}` unpack the result as `{}`",
                if sites.count == 2 {
                    "both of the"
                } else {
                    "all"
                },
                sites.count,
                cand.name,
                tuple_of(names),
            );
            // The body spelling the same names in another order is the
            // transposition itself, on one side or the other.
            let transposed = cand.returns.iter().find(|r| {
                r.names != *names && {
                    let (mut a, mut b) = (r.names.clone(), names.clone());
                    a.sort();
                    b.sort();
                    a == b
                }
            });
            let (msg, at, note) = match transposed {
                Some(r) => (
                    format!(
                        "{lead}, but the body returns them as `{}`. Since {} both orders compile, and one of them is wrong",
                        tuple_of(&r.names),
                        cand.same,
                    ),
                    r.span,
                    "the body returns them in the other order here",
                ),
                None => (
                    format!(
                        "{lead}. The names exist everywhere except in the type, and since {} a swapped pair compiles",
                        cand.same,
                    ),
                    first,
                    "one of the calls",
                ),
            };
            let fields: Vec<String> = names.iter().map(|n| format!("`{n}`")).collect();
            let help = format!("return a struct with fields {}", join(&fields, "and"));
            findings.push((cand.ret_span, msg, at, note, help));
        }
        // `candidates` is a HashMap; report in source order.
        findings.sort_by_key(|(span, ..)| span.lo());
        for (span, msg, at, note, help) in findings {
            emit_with_note(cx, TUPLE_WANTS_STRUCT, span, msg, at, note, help);
        }
    }
}

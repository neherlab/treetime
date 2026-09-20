use std::collections::HashMap;

use rustc_hir::def::{CtorKind, DefKind, Res};
use rustc_hir::def_id::{DefId, LocalDefId};
use rustc_hir::{Block, Expr, ExprKind, LetStmt, Stmt, StmtKind, StructTailExpr, UnOp};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::{Span, sym};

use crate::adt_facts::{matches_config_path, result_err_ty};
use crate::baseline::{emit, emit_with_note};
use crate::enum_facts::{arm_variant, ctor_literal_variant};
use crate::hir_shapes::{callee_of, peel_blocks_unsafe, sole_expr};

rustc_session::declare_lint! {
    /// Flags a call to a function that can reject its input, where the
    /// rejection is replaced with a fixed value and never looked at:
    /// `f(x).unwrap_or(0)`, `.unwrap_or_default()`,
    /// `.unwrap_or_else(|_| CONST)`, `.ok()` feeding one of those, `let
    /// Ok(v) = f(x) else { return Ok(()) }` or `let Some(v) = f(x).ok() else {
    /// return }`, when `f` is a `Result`-returning function of this crate
    /// whose body rejects some of what it is handed (see
    /// `ctor_flow::argument_decided_failure`). Rejected input carries on as
    /// if it were fine. Callees the
    /// analysis cannot see into (other crates, `Option` returners,
    /// combinator-built failures) are covered only when listed in
    /// `defaulted-failure-callees`.
    ///
    /// Silent on `Option` callees (absent-with-a-default is what `Option` is
    /// for; a `Result` names a failure); on failures decided only by the
    /// callee's receiver (its own state, not the caller's value), including
    /// ones an argument test merely stands in front of; on error types in
    /// `validator-resource-errors` or `defaulted-failure-ignored-errors` (the
    /// environment refused, or the failure is already recorded somewhere
    /// else); on a `bool` answer defaulted to `false` by any of the three
    /// forms (the ordinary way to fold "could not tell" into "no"); on
    /// failures handled by a `match`, `if let` or `?`; on a computed fallback
    /// (`unwrap_or(prev)`); on a default in statement position (`f(x)
    /// .unwrap_or_default();` carries nothing on -- like `.ok();`, it is a
    /// discard, `discarded_error`'s shape); and on an else block returning
    /// `None`/`false`, which still reports the failure.
    pub DEFAULTED_FAILURE,
    Warn,
    "a callee's rejection of its argument is replaced by a default"
}

pub struct DefaultedFailure {
    /// Error types whose failure exits do not count: the configured resource
    /// errors plus the ones this lint is told are recorded elsewhere.
    ignored_errors: Vec<String>,
    /// Callees taken as failing on their argument without reading their body.
    listed: Vec<String>,
    /// callee -> the check it fails on, once per callee.
    facts: HashMap<LocalDefId, Option<Span>>,
}

rustc_session::impl_lint_pass!(DefaultedFailure => [DEFAULTED_FAILURE]);

impl DefaultedFailure {
    pub fn new(config: &crate::MordantConfig) -> Self {
        let mut ignored_errors = config.validator_resource_errors.clone();
        ignored_errors.extend(config.defaulted_failure_ignored_errors.iter().cloned());
        Self {
            ignored_errors,
            listed: config.defaulted_failure_callees.clone(),
            facts: HashMap::new(),
        }
    }

    fn fact(&mut self, cx: &LateContext<'_>, callee: LocalDefId) -> Option<Span> {
        if let Some(f) = self.facts.get(&callee) {
            return *f;
        }
        let f = crate::ctor_flow::argument_decided_failure(cx.tcx, callee, &self.ignored_errors);
        self.facts.insert(callee, f);
        f
    }

    fn is_listed(&self, cx: &LateContext<'_>, callee: DefId) -> bool {
        !self.listed.is_empty()
            && matches_config_path(cx.tcx, callee, self.listed.iter().map(String::as_str))
    }

    /// `call` produces the `Result`/`Option` whose failure the caller
    /// replaces; `replaced` says how. Reported when the callee is one this
    /// lint knows rejects its argument.
    fn report<'tcx>(
        &mut self,
        cx: &LateContext<'tcx>,
        call: &Expr<'tcx>,
        at: Span,
        replaced: Replaced,
    ) {
        let Some(FallibleCall { callee, wrapper }) = fallible_call(cx, call) else {
            return;
        };
        let name = cx.tcx.def_path_str(callee);
        let consequence = match replaced {
            Replaced::Default(method) => format!(
                "`.{method}` here replaces the rejection with a fixed value. Rejected input carries on as if it were fine"
            ),
            Replaced::ElseSuccess => "the `else` here returns success when that happens. Rejected input is reported as handled".to_string(),
        };
        let help =
            "pass the error on with `?`, or match on it here while the reason is still attached";
        if wrapper == Wrapper::Result
            && let Some(local) = callee.as_local()
            && let Some(check) = self.fact(cx, local)
        {
            emit_with_note(
                cx,
                DEFAULTED_FAILURE,
                at,
                format!("`{name}` can reject its input, and {consequence}"),
                check,
                format!("the check in `{name}` that gets overridden"),
                help,
            );
        } else if self.is_listed(cx, callee) {
            emit(
                cx,
                DEFAULTED_FAILURE,
                at,
                format!(
                    "`{name}` is listed in `defaulted-failure-callees` as able to reject its input, and {consequence}"
                ),
                help,
            );
        }
    }
}

/// How the caller disposed of the failure.
#[derive(Clone, Copy)]
enum Replaced {
    /// `.unwrap_or(..)` and kin, by method name.
    Default(rustc_span::Symbol),
    /// `let Ok(v) = .. else { return <success> }`.
    ElseSuccess,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Wrapper {
    Result,
    Option,
}

struct FallibleCall {
    callee: DefId,
    wrapper: Wrapper,
}

/// The function `call` invokes, when `call` produces a `Result` or `Option`
/// and is not merely a constructor of one.
fn fallible_call<'tcx>(cx: &LateContext<'tcx>, call: &Expr<'tcx>) -> Option<FallibleCall> {
    let ty::Adt(adt, _) = cx.typeck_results().expr_ty(call).kind() else {
        return None;
    };
    let wrapper = if cx.tcx.is_diagnostic_item(sym::Result, adt.did()) {
        Wrapper::Result
    } else if cx.tcx.is_diagnostic_item(sym::Option, adt.did()) {
        Wrapper::Option
    } else {
        return None;
    };
    let callee = callee_of(cx, call)?.def();
    matches!(cx.tcx.def_kind(callee), DefKind::Fn | DefKind::AssocFn)
        .then_some(FallibleCall { callee, wrapper })
}

/// `recv`, or the `Result` under a value-position `.ok()` on it.
fn through_ok<'h>(cx: &LateContext<'_>, recv: &'h Expr<'h>) -> &'h Expr<'h> {
    let recv = peel_blocks_unsafe(recv);
    if let ExprKind::MethodCall(seg, inner, [], _) = recv.kind
        && seg.ident.as_str() == "ok"
        && result_err_ty(cx.tcx, cx.typeck_results().expr_ty_adjusted(inner)).is_some()
    {
        return peel_blocks_unsafe(inner);
    }
    recv
}

/// A fallback that is the same whatever failed: a literal, a constant, a
/// unit value, `T::default()` / `T::new()`, or those assembled into a tuple,
/// array or struct. A fallback computed from the surroundings is a decision
/// this lint cannot judge.
fn is_fixed_value<'tcx>(cx: &LateContext<'tcx>, e: &'tcx Expr<'tcx>) -> bool {
    let e = peel_blocks_unsafe(e);
    match e.kind {
        ExprKind::Lit(_) => true,
        ExprKind::Unary(UnOp::Neg, inner)
        | ExprKind::Cast(inner, _)
        | ExprKind::AddrOf(_, _, inner) => is_fixed_value(cx, inner),
        ExprKind::Tup(items) | ExprKind::Array(items) => {
            items.iter().all(|i| is_fixed_value(cx, i))
        }
        ExprKind::Path(ref qpath) => is_constant_res(cx.qpath_res(qpath, e.hir_id)),
        ExprKind::Call(_, args) => {
            callee_of(cx, e).is_some_and(|c| match cx.tcx.def_kind(c.def()) {
                DefKind::Fn | DefKind::AssocFn => {
                    args.is_empty() && is_nullary_ctor_name(cx, c.def())
                }
                DefKind::Ctor(..) => args.iter().all(|a| is_fixed_value(cx, a)),
                _ => false,
            })
        }
        ExprKind::Struct(_, fields, tail) => {
            fields.iter().all(|f| is_fixed_value(cx, f.expr))
                && match tail {
                    StructTailExpr::None | StructTailExpr::DefaultFields(_) => true,
                    StructTailExpr::Base(base) => is_fixed_value(cx, base),
                    _ => false,
                }
        }
        _ => false,
    }
}

fn is_constant_res(res: Res) -> bool {
    matches!(
        res,
        Res::Def(
            DefKind::Const { .. } | DefKind::AssocConst { .. } | DefKind::Ctor(_, CtorKind::Const),
            _
        )
    )
}

fn is_nullary_ctor_name(cx: &LateContext<'_>, def: DefId) -> bool {
    matches!(cx.tcx.item_name(def).as_str(), "default" | "new")
}

/// What a caller puts in place of the failure, when it is the same whatever
/// failed.
enum Fallback<'h> {
    /// The type's own `default()` / `new()`: `unwrap_or_default()`, or
    /// `unwrap_or_else(T::default)`.
    Default,
    /// The expression written out: `unwrap_or(v)`, or the body of the
    /// closure in `unwrap_or_else(|_| v)`.
    Value(&'h Expr<'h>),
}

/// The fallback of a `unwrap_or*` call named `name` with `args`, when it is
/// fixed. A fallback computed from the surroundings (`unwrap_or(prev)`, a
/// closure reading the error it is handed) is None.
fn fixed_fallback<'tcx>(
    cx: &LateContext<'tcx>,
    name: &str,
    args: &'tcx [Expr<'tcx>],
) -> Option<Fallback<'tcx>> {
    match (name, args) {
        ("unwrap_or_default", []) => Some(Fallback::Default),
        ("unwrap_or", [value]) => is_fixed_value(cx, value).then_some(Fallback::Value(value)),
        ("unwrap_or_else", [thunk]) => {
            let thunk = peel_blocks_unsafe(thunk);
            match thunk.kind {
                ExprKind::Closure(c) => {
                    let value = cx.tcx.hir_body(c.body).value;
                    is_fixed_value(cx, value).then_some(Fallback::Value(value))
                }
                ExprKind::Path(ref qpath) => match cx.qpath_res(qpath, thunk.hir_id) {
                    Res::Def(DefKind::Fn | DefKind::AssocFn, def)
                        if is_nullary_ctor_name(cx, def) =>
                    {
                        Some(Fallback::Default)
                    }
                    _ => None,
                },
                _ => None,
            }
        }
        _ => None,
    }
}

/// The else block of a `let .. else` reports success or nothing at all:
/// `return`, `return ()` or `return Ok(())`, and nothing else in the block.
fn else_reports_success(cx: &LateContext<'_>, els: &Block<'_>) -> bool {
    let Some(only) = sole_expr(els) else {
        return false;
    };
    let ExprKind::Ret(value) = peel_blocks_unsafe(only).kind else {
        return false;
    };
    match value {
        None => true,
        Some(v) => {
            let v = peel_blocks_unsafe(v);
            match v.kind {
                ExprKind::Tup([]) => true,
                ExprKind::Call(_, [payload]) => {
                    ctor_literal_variant(cx, v)
                        .is_some_and(|variant| cx.tcx.item_name(variant) == sym::Ok)
                        && matches!(peel_blocks_unsafe(payload).kind, ExprKind::Tup([]))
                }
                _ => false,
            }
        }
    }
}

/// A `bool` answer defaulted to `false`, whichever way it is spelled
/// (`unwrap_or(false)`, `unwrap_or_default()`, `unwrap_or_else(|_| false)`,
/// `|_| Default::default()`): on a `Result<bool, _>` this folds "could not
/// tell" into "no", which is the answer's own vocabulary rather than a value
/// smuggled past a check. `unwrap_or(true)` is not this: it answers "yes"
/// unasked.
fn folds_to_no<'tcx>(
    cx: &LateContext<'tcx>,
    unwrapped: &Expr<'_>,
    fallback: &Fallback<'tcx>,
) -> bool {
    if !cx.typeck_results().expr_ty(unwrapped).is_bool() {
        return false;
    }
    let value = match fallback {
        Fallback::Default => return true,
        Fallback::Value(value) => peel_blocks_unsafe(value),
    };
    match value.kind {
        ExprKind::Lit(lit) => lit.node == rustc_ast::LitKind::Bool(false),
        ExprKind::Call(_, []) => callee_of(cx, value).is_some_and(|c| {
            matches!(cx.tcx.def_kind(c.def()), DefKind::Fn | DefKind::AssocFn)
                && cx.tcx.item_name(c.def()).as_str() == "default"
        }),
        _ => false,
    }
}

/// `expr` is a whole statement (`f(x).unwrap_or_default();`): the value is
/// dropped, so nothing carries on with it. That is a discard, `.ok();`'s
/// shape, not a default.
fn is_discarded(cx: &LateContext<'_>, expr: &Expr<'_>) -> bool {
    matches!(
        cx.tcx.parent_hir_node(expr.hir_id),
        rustc_hir::Node::Stmt(Stmt {
            kind: StmtKind::Semi(_),
            ..
        })
    )
}

impl<'tcx> LateLintPass<'tcx> for DefaultedFailure {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let ExprKind::MethodCall(seg, recv, args, _) = expr.kind else {
            return;
        };
        let name = seg.ident.name;
        let Some(fallback) = fixed_fallback(cx, name.as_str(), args) else {
            return;
        };
        if expr.span.from_expansion() || folds_to_no(cx, expr, &fallback) || is_discarded(cx, expr)
        {
            return;
        }
        let call = through_ok(cx, recv);
        self.report(cx, call, expr.span, Replaced::Default(name));
    }

    /// `let Ok(v) = f(x) else { return Ok(()) };`, or the same with the
    /// `Result` turned into an `Option` first: `let Some(v) = f(x).ok() else
    /// { .. }`. A `Some` head over anything else is not taken: an `Option`
    /// callee's `None` under a bare `return` is the ordinary "not this shape,
    /// nothing to do", and a lint cannot tell that from a swallowed rejection.
    fn check_stmt(&mut self, cx: &LateContext<'tcx>, stmt: &'tcx Stmt<'tcx>) {
        let StmtKind::Let(LetStmt {
            pat,
            init: Some(init),
            els: Some(els),
            span,
            ..
        }) = stmt.kind
        else {
            return;
        };
        if span.from_expansion() {
            return;
        }
        let init = peel_blocks_unsafe(init);
        let call = through_ok(cx, init);
        let head_is_failure_of_call = arm_variant(cx, pat).is_some_and(|variant| {
            let head = cx.tcx.item_name(variant);
            head == sym::Ok || (head == sym::Some && !std::ptr::eq(call, init))
        });
        if !head_is_failure_of_call || !else_reports_success(cx, els) {
            return;
        }
        self.report(cx, call, *span, Replaced::ElseSuccess);
    }
}

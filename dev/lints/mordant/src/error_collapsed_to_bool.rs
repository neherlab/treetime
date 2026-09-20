use std::collections::{HashMap, HashSet};
use std::ops::ControlFlow;

use clippy_utils::visitors::for_each_expr_without_closures;
use clippy_utils::{
    get_expr_use_or_unification_node, is_def_id_trait_method, is_in_test, is_refutable,
};
use rustc_abi::ExternAbi;
use rustc_hir::def::DefKind;
use rustc_hir::def_id::{DefId, LocalDefId};
use rustc_hir::intravisit::FnKind;
use rustc_hir::{
    Block, Body, Expr, ExprKind, FnDecl, HirId, LangItem, LetStmt, MatchSource, Node, Pat, PatKind,
    StmtKind,
};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::Ty;
use rustc_middle::ty::layout::LayoutOf;
use rustc_middle::ty::print::with_no_trimmed_paths;
use rustc_span::Span;

use crate::adt_facts::{is_option_ty, result_err_ty};
use crate::baseline::emit_hir_then;
use crate::enum_facts::{arm_variant, ctor_literal_variant};
use crate::hir_shapes::{callee_of, peel_blocks_unsafe, peel_not, sole_expr};

rustc_session::declare_lint! {
    /// Flags a call that ignores the `bool` or `Option` a function of this
    /// crate turned its typed error into, so a failure inside that function
    /// is treated as success. The callee returns `bool` (or
    /// `Option<T>`) and somewhere in its body the `Err` of a `Result` it held
    /// becomes a bare `false` (or `None`) and nothing else -- `match r {
    /// Err(_) => return false, .. }`, `if r.is_err() { return false }`, `let
    /// Ok(v) = r else { return None }`, `r.ok()?`, or `r.is_ok()` / `r.ok()`
    /// as the value returned -- so the error's kind is already gone from its
    /// signature. The reported call then ignores even that bit (`f(x);`,
    /// `let _ = f(x);`). A `Result` return would have made this caller decide.
    ///
    /// Silent when the `Err` arm or `else` block does anything besides exit
    /// (logs, stores or converts the error: it was looked at), or an `Err`
    /// pattern of that `match`/`if let` names a kind (`Err(Errno::NOENT) =>
    /// false` answers a question, and the sibling arm saw the rest); when the
    /// error type had no kind to lose -- zero-sized (`()`, `AllocError`,
    /// `TryFromIntError`, a lone unit variant), where `false` says as much as
    /// the error did, or a bare primitive (`binary_search`'s `Err(idx)` is an
    /// answer, not a failure); on trait methods and non-Rust-ABI functions
    /// (the signature is not the function's to choose); on collapses inside
    /// closures; on callees in other crates; on calls in tests or produced
    /// by macros; on every call whose value is read (`if`, `&&`, `let x =`,
    /// `?`, a tail); and on `unwrap_or(false)` over a `Result<bool, _>`,
    /// which is `defaulted_failure`'s vocabulary. `discarded_error` owns
    /// `.ok();` on the `Result` itself: there the collapse and the drop are
    /// one expression; here they are a signature apart.
    pub ERROR_COLLAPSED_TO_BOOL,
    Warn,
    "a call drops the bool or Option a callee collapsed a typed error into"
}

/// What the collapsing function hands back in place of the error.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Exit {
    False,
    None,
}

impl Exit {
    fn as_str(self) -> &'static str {
        match self {
            Exit::False => "false",
            Exit::None => "None",
        }
    }
}

/// One function's collapse: the first site by position, how many more
/// there are, and the error type lost at that first site.
struct Collapse {
    exit: Exit,
    at: Span,
    err: String,
    more: usize,
}

/// A call to a local `bool`/`Option` function whose value nothing reads.
struct Dropped {
    callee: LocalDefId,
    hir_id: HirId,
    span: Span,
}

#[derive(Default)]
pub struct ErrorCollapsedToBool {
    collapses: HashMap<LocalDefId, Collapse>,
    /// In visit order, so reports come out in source order per file.
    drops: Vec<Dropped>,
}

rustc_session::impl_lint_pass!(ErrorCollapsedToBool => [ERROR_COLLAPSED_TO_BOOL]);

/// The error type of `e`'s `Result`, when it has a kind that `false` loses:
/// not a bare primitive, and not zero-sized (`()`, `!`, a unit struct, a
/// lone unit variant), which distinguishes nothing a `bool` does not. A type
/// whose layout is unknown here (generic, unsized) is given the benefit of
/// the doubt.
fn err_of<'tcx>(cx: &LateContext<'tcx>, e: &Expr<'tcx>) -> Option<Ty<'tcx>> {
    let ty = cx.typeck_results().expr_ty(e).peel_refs();
    let err = result_err_ty(cx.tcx, ty)?;
    let bare = err.peel_refs();
    if bare.is_primitive() || cx.layout_of(bare).is_ok_and(|l| l.is_zst()) {
        return None;
    }
    Some(err)
}

/// `pat` takes its variant's payload whole (`Err(_)`, `Err(e)`, `Err(..)`)
/// rather than naming a kind of it (`Err(Errno::NOENT)`).
fn takes_whole(cx: &LateContext<'_>, pat: &Pat<'_>) -> bool {
    match pat.kind {
        PatKind::TupleStruct(_, subs, _) => !subs.iter().any(|p| is_refutable(cx, p)),
        PatKind::Struct(_, fields, _) => !fields.iter().any(|f| is_refutable(cx, f.pat)),
        _ => false,
    }
}

/// `v` is the bare exit value: the literal `false`, or the path `None`.
fn is_exit(cx: &LateContext<'_>, v: &Expr<'_>, exit: Exit) -> bool {
    let v = peel_blocks_unsafe(v);
    match exit {
        Exit::False => {
            matches!(v.kind, ExprKind::Lit(lit) if lit.node == rustc_ast::LitKind::Bool(false))
        }
        Exit::None => ctor_literal_variant(cx, v)
            .is_some_and(|d| Some(d) == cx.tcx.lang_items().option_none_variant()),
    }
}

/// `e` does nothing but leave the function with the exit value: `return
/// false`, `{ return false; }`, or -- when `e` is itself the function's
/// value (`value_position`) -- the bare `false`. Returns that leaving
/// expression.
fn pure_exit<'h>(
    cx: &LateContext<'_>,
    e: &'h Expr<'h>,
    exit: Exit,
    value_position: bool,
) -> Option<&'h Expr<'h>> {
    let mut e = peel_blocks_unsafe(e);
    if let ExprKind::Block(b, _) = e.kind {
        e = peel_blocks_unsafe(sole_expr(b)?);
    }
    let exits = match e.kind {
        ExprKind::Ret(Some(v)) => is_exit(cx, v, exit),
        ExprKind::Ret(None) => false,
        _ => value_position && is_exit(cx, e, exit),
    };
    exits.then_some(e)
}

/// The expressions whose value is the function's value: the body's tail
/// and every `return` operand, followed down through blocks, `if` and
/// `match`. `branches` are the `if`/`match` nodes passed through (their
/// arms are in value position); `leaves` are where the descent stopped.
#[derive(Default)]
struct Tails<'tcx> {
    branches: HashSet<HirId>,
    leaves: Vec<&'tcx Expr<'tcx>>,
}

impl<'tcx> Tails<'tcx> {
    fn of(body: &'tcx Body<'tcx>) -> Self {
        let mut tails = Self::default();
        tails.descend(body.value);
        for_each_expr_without_closures(body.value, |e: &'tcx Expr<'tcx>| {
            if let ExprKind::Ret(Some(v)) = e.kind {
                tails.descend(v);
            }
            ControlFlow::<()>::Continue(())
        });
        tails
    }

    fn descend(&mut self, e: &'tcx Expr<'tcx>) {
        match e.kind {
            ExprKind::DropTemps(inner) => self.descend(inner),
            ExprKind::Block(b, _) => {
                if let Some(tail) = b.expr {
                    self.descend(tail);
                }
            }
            ExprKind::If(_, then, els) => {
                self.branches.insert(e.hir_id);
                self.descend(then);
                if let Some(els) = els {
                    self.descend(els);
                }
            }
            ExprKind::Match(_, arms, MatchSource::Normal) => {
                self.branches.insert(e.hir_id);
                for arm in arms {
                    self.descend(arm.body);
                }
            }
            _ => self.leaves.push(e),
        }
    }
}

/// Every site in `body` where a `Result`'s error becomes the bare `exit`
/// value and nothing else happens to it.
fn collapses<'tcx>(
    cx: &LateContext<'tcx>,
    body: &'tcx Body<'tcx>,
    exit: Exit,
) -> Vec<(Span, Ty<'tcx>)> {
    let li = cx.tcx.lang_items();
    let (err_variant, ok_variant) = (li.result_err_variant(), li.result_ok_variant());
    let tails = Tails::of(body);
    let mut sites = Vec::new();

    for leaf in &tails.leaves {
        if let ExprKind::MethodCall(seg, recv, [], _) = leaf.kind
            && let Some(err) = err_of(cx, recv)
            && match exit {
                Exit::False => seg.ident.as_str() == "is_ok",
                Exit::None => seg.ident.as_str() == "ok",
            }
        {
            sites.push((leaf.span, err));
        }
    }

    for_each_expr_without_closures(body.value, |e: &'tcx Expr<'tcx>| {
        let value_position = tails.branches.contains(&e.hir_id);
        match e.kind {
            // One `Err` arm that names a kind, is guarded, or does more than
            // exit means the error was looked at in this `match`.
            ExprKind::Match(scrut, arms, MatchSource::Normal) => {
                if let Some(err) = err_of(cx, scrut) {
                    let err_arms: Vec<_> = arms
                        .iter()
                        .filter(|arm| {
                            arm_variant(cx, arm.pat).is_some_and(|v| Some(v) == err_variant)
                        })
                        .collect();
                    if !err_arms.is_empty()
                        && err_arms.iter().all(|arm| {
                            arm.guard.is_none()
                                && takes_whole(cx, arm.pat)
                                && pure_exit(cx, arm.body, exit, value_position).is_some()
                        })
                    {
                        sites.extend(err_arms.iter().map(|arm| (arm.span, err)));
                    }
                }
            }
            ExprKind::Match(scrut, _, MatchSource::TryDesugar(_)) if exit == Exit::None => {
                if let ExprKind::Call(branch, [operand]) = scrut.kind
                    && let ExprKind::Path(qpath) = branch.kind
                    && cx.tcx.qpath_is_lang_item(qpath, LangItem::TryTraitBranch)
                    && let ExprKind::MethodCall(seg, recv, [], _) = peel_blocks_unsafe(operand).kind
                    && seg.ident.as_str() == "ok"
                    && let Some(err) = err_of(cx, recv)
                {
                    sites.push((operand.span, err));
                }
            }
            ExprKind::If(cond, then, els) => {
                let (cond, negated) = peel_not(cond);
                let taken = match cond.kind {
                    ExprKind::MethodCall(seg, recv, [], _) => err_of(cx, recv).and_then(|err| {
                        let on_err = match (seg.ident.as_str(), negated) {
                            ("is_err", false) | ("is_ok", true) => Some(then),
                            ("is_ok", false) | ("is_err", true) => els,
                            _ => None,
                        };
                        on_err.map(|b| (b, err))
                    }),
                    ExprKind::Let(l) => err_of(cx, l.init).and_then(|err| {
                        let head = arm_variant(cx, l.pat)?;
                        let on_err = if Some(head) == err_variant {
                            takes_whole(cx, l.pat).then_some(then)
                        } else if Some(head) == ok_variant {
                            els
                        } else {
                            None
                        };
                        on_err.map(|b| (b, err))
                    }),
                    _ => None,
                };
                if let Some((on_err, err)) = taken
                    && let Some(leave) = pure_exit(cx, on_err, exit, value_position)
                {
                    sites.push((leave.span, err));
                }
            }
            // A `loop` body is the one block that is not also an expression.
            ExprKind::Block(block, _) | ExprKind::Loop(block, ..) => {
                let_else_collapses(cx, block, exit, ok_variant, &mut sites)
            }
            _ => {}
        }
        ControlFlow::<()>::Continue(())
    });
    sites
}

/// `let Ok(v) = r else { return false };` statements of `block`.
fn let_else_collapses<'tcx>(
    cx: &LateContext<'tcx>,
    block: &'tcx Block<'tcx>,
    exit: Exit,
    ok_variant: Option<DefId>,
    sites: &mut Vec<(Span, Ty<'tcx>)>,
) {
    for stmt in block.stmts {
        if let StmtKind::Let(LetStmt {
            pat,
            init: Some(init),
            els: Some(els),
            ..
        }) = stmt.kind
            && let Some(err) = err_of(cx, init)
            && arm_variant(cx, pat).is_some_and(|v| Some(v) == ok_variant)
            && let Some(only) = sole_expr(els)
            && let Some(leave) = pure_exit(cx, only, exit, false)
        {
            sites.push((leave.span, err));
        }
    }
}

/// Nothing reads the call's value: it is a statement of its own (through
/// any blocks or one-armed matches that merely pass it up), or the
/// initializer of `let _ =`.
fn is_dropped(cx: &LateContext<'_>, call: &Expr<'_>) -> bool {
    match get_expr_use_or_unification_node(cx.tcx, call) {
        None => true,
        Some((Node::Stmt(stmt), _)) => matches!(stmt.kind, StmtKind::Semi(_) | StmtKind::Expr(_)),
        Some((Node::LetStmt(l), _)) => matches!(l.pat.kind, PatKind::Wild),
        Some(_) => false,
    }
}

impl<'tcx> LateLintPass<'tcx> for ErrorCollapsedToBool {
    fn check_fn(
        &mut self,
        cx: &LateContext<'tcx>,
        kind: FnKind<'tcx>,
        _decl: &'tcx FnDecl<'tcx>,
        body: &'tcx Body<'tcx>,
        span: Span,
        def_id: LocalDefId,
    ) {
        if matches!(kind, FnKind::Closure)
            || span.from_expansion()
            || kind.header().is_some_and(|h| h.abi != ExternAbi::Rust)
            || is_def_id_trait_method(cx, def_id)
            || cx.tcx.trait_of_assoc(def_id.to_def_id()).is_some()
        {
            return;
        }
        let ret = cx
            .tcx
            .fn_sig(def_id)
            .instantiate_identity()
            .skip_normalization()
            .output()
            .skip_binder();
        let exit = if ret.is_bool() {
            Exit::False
        } else if is_option_ty(cx, ret) {
            Exit::None
        } else {
            return;
        };
        let mut sites = collapses(cx, body, exit);
        sites.sort_by_key(|(span, _)| span.lo());
        sites.dedup_by_key(|(span, _)| *span);
        // Printed now because `Ty` cannot outlive the pass; untrimmed because
        // most of these are never reported, and a trimmed print with no
        // diagnostic after it is a delayed ICE.
        if let Some(&(at, err)) = sites.first() {
            self.collapses.insert(
                def_id,
                Collapse {
                    exit,
                    at,
                    err: with_no_trimmed_paths!(err.to_string()),
                    more: sites.len() - 1,
                },
            );
        }
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        if !matches!(expr.kind, ExprKind::Call(..) | ExprKind::MethodCall(..))
            || expr.span.from_expansion()
        {
            return;
        }
        let Some(callee) = callee_of(cx, expr) else {
            return;
        };
        let Some(local) = callee.def().as_local() else {
            return;
        };
        if !matches!(cx.tcx.def_kind(local), DefKind::Fn | DefKind::AssocFn) {
            return;
        }
        let ty = cx.typeck_results().expr_ty(expr);
        if !(ty.is_bool() || is_option_ty(cx, ty))
            || !is_dropped(cx, expr)
            || is_in_test(cx.tcx, expr.hir_id)
        {
            return;
        }
        self.drops.push(Dropped {
            callee: local,
            hir_id: expr.hir_id,
            span: expr.span,
        });
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        for dropped in &self.drops {
            let Some(collapse) = self.collapses.get(&dropped.callee) else {
                continue;
            };
            let name = with_no_trimmed_paths!(cx.tcx.def_path_str(dropped.callee));
            let exit = collapse.exit.as_str();
            let err = &collapse.err;
            emit_hir_then(
                cx,
                ERROR_COLLAPSED_TO_BOOL,
                dropped.hir_id,
                dropped.span,
                format!(
                    "`{name}` turns its `{err}` into `{exit}`, and this call ignores the `{exit}`. A failure inside `{name}` is treated as success"
                ),
                |diag| {
                    let more = match collapse.more {
                        0 => String::new(),
                        n => format!(" (and in {n} more places)"),
                    };
                    diag.span_note(
                        collapse.at,
                        format!("the error becomes `{exit}` here{more}"),
                    );
                    diag.help(format!(
                        "return the `Result` from `{name}` and decide here what failure means. At minimum mark it `#[must_use]`"
                    ));
                },
            );
        }
    }
}

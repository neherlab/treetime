use std::collections::HashMap;

use rustc_abi::ExternAbi;
use rustc_ast::LitKind;
use rustc_hir::def::DefKind;
use rustc_hir::def_id::DefId;
use rustc_hir::{Expr, ExprKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::Span;

use crate::baseline::emit_with_note;
use crate::hir_shapes::{Callee, callee_of};

rustc_session::declare_lint! {
    /// Flags a crate-private function that takes two or more bools when a
    /// call passes bare `true` / `false` for at least two of them. At
    /// `f(x, true, false)` nothing says which is which, and the swapped call
    /// compiles too. A two-variant enum per flag, or an options struct,
    /// makes every call name what it sets.
    ///
    /// Stays quiet on a single `bool` parameter (nothing to confuse it
    /// with), on exported, `extern` and trait functions (their signature is
    /// not the crate's to change), when every call passes named values or
    /// at most one bare literal, and on a literal a comment beside it names
    /// (`/* force */ true`, or `true, // force` before a line break).
    pub BARE_BOOL_ARGS,
    Warn,
    "several bool parameters filled with bare literals at a call site"
}

#[derive(Default)]
struct Calls {
    sites: usize,
    /// Calls with two or more unnamed bool literals, as written (`f(.., true, false)`).
    bare: Vec<(Span, String)>,
}

#[derive(Default)]
pub struct BareBoolArgs {
    calls: HashMap<DefId, Calls>,
}

rustc_session::impl_lint_pass!(BareBoolArgs => [BARE_BOOL_ARGS]);

/// Signature indices of `def`'s `bool` parameters, when it is a crate-private
/// Rust-ABI fn of this crate whose signature is its own (not a trait's).
fn bool_params(cx: &LateContext<'_>, def: DefId) -> Option<Vec<usize>> {
    let local = def.as_local()?;
    if !matches!(cx.tcx.def_kind(def), DefKind::Fn | DefKind::AssocFn)
        || cx.effective_visibilities.is_exported(local)
        || cx.tcx.trait_of_assoc(def).is_some()
        || cx.tcx.trait_impl_of_assoc(def).is_some()
        || cx.tcx.def_span(def).from_expansion()
    {
        return None;
    }
    let sig = cx
        .tcx
        .fn_sig(def)
        .instantiate_identity()
        .skip_normalization();
    if sig.abi() != ExternAbi::Rust {
        return None;
    }
    let bools: Vec<usize> = sig
        .inputs()
        .skip_binder()
        .iter()
        .enumerate()
        .filter(|(_, ty)| ty.is_bool())
        .map(|(i, _)| i)
        .collect();
    (bools.len() >= 2).then_some(bools)
}

fn bare_bool(e: &Expr<'_>) -> Option<bool> {
    match e.kind {
        ExprKind::Lit(lit) if !e.span.from_expansion() => match lit.node {
            LitKind::Bool(b) => Some(b),
            _ => None,
        },
        _ => None,
    }
}

/// A comment beside the argument names it: `/* force */ true` or a comment
/// line above it (leading), `true, // force` before a line break (trailing).
/// The first line of a gap between two arguments belongs to the earlier one,
/// the rest to the later one, so one comment never names both neighbours.
fn commented(cx: &LateContext<'_>, before: Span, arg: Span, after: Span) -> bool {
    let sm = cx.tcx.sess.source_map();
    let has_comment = |s: &str| s.contains("//") || s.contains("/*");
    let lead = sm
        .span_to_snippet(before.between(arg))
        .is_ok_and(|s| match s.split_once('\n') {
            Some((_, rest)) => has_comment(rest),
            None => has_comment(&s),
        });
    let trail = sm.span_to_snippet(arg.between(after)).is_ok_and(|s| {
        s.split_once('\n')
            .is_some_and(|(line, _)| has_comment(line))
    });
    lead || trail
}

impl<'tcx> LateLintPass<'tcx> for BareBoolArgs {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        if expr.span.from_expansion() {
            return;
        }
        // (callee, anchor before the first argument, arguments, signature
        // index of argument 0: a method's receiver is input 0).
        let (def, anchor, args, offset) = match (callee_of(cx, expr), expr.kind) {
            (Some(Callee::Path { def, args }), ExprKind::Call(callee, _)) => {
                (def, callee.span, args, 0)
            }
            (Some(Callee::Method { def, args, .. }), ExprKind::MethodCall(seg, ..)) => {
                (def, seg.ident.span, args, 1)
            }
            _ => return,
        };
        let Some(bools) = bool_params(cx, def) else {
            return;
        };
        let calls = self.calls.entry(def).or_default();
        calls.sites += 1;
        let close = expr.span.shrink_to_hi();
        let mut unnamed = 0;
        let mut shape: Vec<&str> = Vec::new();
        for (i, arg) in args.iter().enumerate() {
            let before = if i == 0 { anchor } else { args[i - 1].span };
            let after = args.get(i + 1).map_or(close, |a| a.span);
            let word = match bare_bool(arg) {
                Some(b)
                    if bools.contains(&(i + offset)) && !commented(cx, before, arg.span, after) =>
                {
                    unnamed += 1;
                    if b { "true" } else { "false" }
                }
                _ => "..",
            };
            if word != ".." || shape.last() != Some(&"..") {
                shape.push(word);
            }
        }
        if unnamed >= 2 {
            let name = cx.tcx.item_name(def);
            calls
                .bare
                .push((expr.span, format!("{name}({})", shape.join(", "))));
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let mut findings: Vec<(Span, String, Span, String)> = Vec::new();
        for (def, calls) in &mut self.calls {
            calls.bare.sort_by_key(|(span, _)| span.lo());
            let Some((first, written)) = calls.bare.first() else {
                continue;
            };
            let Some(bools) = bool_params(cx, *def) else {
                continue;
            };
            let idents = cx.tcx.fn_arg_idents(*def);
            let mut names: Vec<String> = bools
                .iter()
                .map(|&i| match idents.get(i).copied().flatten() {
                    Some(ident) => format!("`{}`", ident.name),
                    None => "`_`".to_string(),
                })
                .collect();
            let last = names.pop().unwrap_or_default();
            let which = if names.len() == 1 {
                "both"
            } else {
                "at least two of them"
            };
            let (k, n) = (calls.bare.len(), calls.sites);
            let pass = if k == 1 { "passes" } else { "pass" };
            let span = cx
                .tcx
                .def_ident_span(*def)
                .unwrap_or_else(|| cx.tcx.def_span(*def));
            let params = format!("{} and {last}", names.join(", "));
            findings.push((
                span,
                format!(
                    "`{}` takes bools {params}, and {k} of its {n} calls {pass} bare `true`/`false` for {which}. At `{written}` nothing says which is which",
                    cx.tcx.item_name(*def),
                ),
                *first,
                format!(
                    "give {params} a two-variant enum each, or an options struct, so every call names what it sets"
                ),
            ));
        }
        // `calls` is a HashMap; report in source order.
        findings.sort_by_key(|(span, ..)| span.lo());
        for (span, msg, call, help) in findings {
            emit_with_note(
                cx,
                BARE_BOOL_ARGS,
                span,
                msg,
                call,
                "one of those calls",
                help,
            );
        }
    }
}

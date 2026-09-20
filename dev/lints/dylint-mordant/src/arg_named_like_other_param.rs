use std::ops::ControlFlow;

use clippy_utils::get_parent_expr;
use clippy_utils::visitors::for_each_expr;
use rustc_hir::def::{DefKind, Res};
use rustc_hir::def_id::DefId;
use rustc_hir::{BinOpKind, Expr, ExprKind, QPath, UnOp};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::Ty;
use rustc_span::{Span, Symbol};

use crate::baseline::{emit, emit_with_note};
use crate::hir_shapes::{Callee, callee_of};

rustc_session::declare_lint! {
    /// Flags a call argument passed where the callee expects a parameter of
    /// another name, when the callee also has a parameter of the argument's
    /// own name and the same type: `resize(height, width)` against
    /// `fn resize(width: u32, height: u32)`, or `spawn(opts.inherit_stderr,
    /// opts.inherit_stdout)`. It looks like an argument in the wrong
    /// position, and since the two share a type the swap compiles. Distinct
    /// types per parameter (a newtype per quantity, an enum per flag) would
    /// make it a compile error.
    ///
    /// The argument's name is a local or parameter it is spelled as, or the
    /// last field of a field access. Silent when the names agree after
    /// dropping a leading `is_`/`has_`/`_` and a trailing `_`, when the
    /// parameter the name points at already receives an argument of that
    /// name, when the bound parameter's name contains the argument's as a
    /// word (`from_index` receiving `index`), when the same condition also
    /// calls the callee with the two values the other way round (`sub(a, b)
    /// && sub(b, a)` is a symmetric use), for `self`/`this` on either side
    /// (a receiver slot is a grammatical position, not a role), for
    /// one-character names, and for calls through closures or fn pointers.
    /// A literal or constant in the namesake slot does not excuse the call:
    /// `resize(height, 0)` is the transposition with one side spelled out.
    pub ARG_NAMED_LIKE_OTHER_PARAM,
    Warn,
    "argument named as a different same-typed parameter of the callee"
}

rustc_session::declare_lint_pass!(ArgNamedLikeOtherParam => [ARG_NAMED_LIKE_OTHER_PARAM]);

/// The name an argument is spelled with: a local `x`, or the last field of
/// `x.a.b`, through `&`, `*` and casts. Anything computed has no name to
/// cross.
fn arg_name(mut e: &Expr<'_>) -> Option<Symbol> {
    loop {
        match &e.kind {
            ExprKind::Cast(inner, _)
            | ExprKind::AddrOf(_, _, inner)
            | ExprKind::Unary(UnOp::Deref, inner)
            | ExprKind::DropTemps(inner) => e = inner,
            ExprKind::Field(_, ident) => return Some(ident.name),
            ExprKind::Path(QPath::Resolved(None, p))
                if p.segments.len() == 1 && matches!(p.res, Res::Local(_)) =>
            {
                return Some(p.segments[0].ident.name);
            }
            _ => return None,
        }
    }
}

/// The role a name claims: `is_open_` and `open` claim the same one;
/// `self`, `self_` and `this` name a position, not a role, so they claim
/// none, and neither does a single character.
fn role(name: Symbol) -> Option<Symbol> {
    let name = name.as_str();
    let name = name.trim_start_matches('_').trim_end_matches('_');
    let name = name
        .strip_prefix("is_")
        .or_else(|| name.strip_prefix("has_"))
        .unwrap_or(name);
    (name.len() >= 2 && name != "self" && name != "this").then(|| Symbol::intern(name))
}

/// The role the argument in signature slot `slot` of a call is named for.
fn received(args: &[Expr<'_>], offset: usize, slot: usize) -> Option<Symbol> {
    let arg = args.get(slot.checked_sub(offset)?)?;
    role(arg_name(arg)?)
}

fn args_of<'tcx>(callee: &Callee<'tcx>) -> (&'tcx [Expr<'tcx>], usize) {
    match *callee {
        Callee::Path { args, .. } => (args, 0),
        // The receiver is bound to `self`; explicit args start one slot on.
        Callee::Method { args, .. } => (args, 1),
    }
}

/// `f(a, b) && f(b, a)`: the condition `expr` sits in also calls `def`
/// with the same two names in slots `slot` and `other` the other way round,
/// so the crossed call is the deliberate second half of a symmetric test.
fn mirrored_in_condition<'tcx>(
    cx: &LateContext<'tcx>,
    expr: &'tcx Expr<'tcx>,
    def: DefId,
    (slot, slot_role): (usize, Symbol),
    (other, other_role): (usize, Symbol),
) -> bool {
    let mut root = expr;
    while let Some(parent) = get_parent_expr(cx, root)
        && match parent.kind {
            ExprKind::Unary(UnOp::Not, _) | ExprKind::DropTemps(_) => true,
            ExprKind::Binary(op, ..) => matches!(op.node, BinOpKind::And | BinOpKind::Or),
            _ => false,
        }
    {
        root = parent;
    }
    if root.hir_id == expr.hir_id {
        return false;
    }
    for_each_expr(cx, root, |e: &'tcx Expr<'tcx>| {
        if e.hir_id != expr.hir_id
            && let Some(c) = callee_of(cx, e)
            && c.def() == def
        {
            let (args, offset) = args_of(&c);
            if received(args, offset, slot) == Some(slot_role)
                && received(args, offset, other) == Some(other_role)
            {
                return ControlFlow::Break(());
            }
        }
        ControlFlow::Continue(())
    })
    .is_some()
}

/// `from_file_index` qualifies `file_index` rather than naming another slot.
fn qualifies(param: &str, arg: &str) -> bool {
    param.len() > arg.len()
        && ((param.ends_with(arg) && param[..param.len() - arg.len()].ends_with('_'))
            || (param.starts_with(arg) && param[arg.len()..].starts_with('_')))
}

struct Crossing<'tcx> {
    /// Index into the call's argument list.
    arg: usize,
    name: Symbol,
    bound_to: Symbol,
    /// Signature index of the parameter the name points at.
    names_param: usize,
    ty: Ty<'tcx>,
}

impl<'tcx> LateLintPass<'tcx> for ArgNamedLikeOtherParam {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        if expr.span.from_expansion() {
            return;
        }
        let Some(callee) = callee_of(cx, expr) else {
            return;
        };
        let def = callee.def();
        if !matches!(cx.tcx.def_kind(def), DefKind::Fn | DefKind::AssocFn) {
            return;
        }
        let (args, offset) = args_of(&callee);
        let idents = cx.tcx.fn_arg_idents(def);
        let sig = cx.tcx.erase_and_anonymize_regions(
            cx.tcx.instantiate_bound_regions_with_erased(
                cx.tcx
                    .fn_sig(def)
                    .instantiate_identity()
                    .skip_normalization(),
            ),
        );
        let inputs = sig.inputs();
        if idents.len() != inputs.len() || args.len() + offset > inputs.len() {
            return;
        }
        let param_name = |i: usize| idents[i].map(|id| id.name);
        let mut crossings: Vec<Crossing<'tcx>> = Vec::new();
        for (i, arg) in args.iter().enumerate() {
            if arg.span.from_expansion() {
                continue;
            }
            let slot = i + offset;
            let (Some(name), Some(bound_to)) = (arg_name(arg), param_name(slot)) else {
                continue;
            };
            let (Some(an), Some(bn)) = (role(name), role(bound_to)) else {
                continue;
            };
            if an == bn || qualifies(bn.as_str(), an.as_str()) {
                continue;
            }
            let Some(other) = (0..inputs.len()).find(|&q| {
                q != slot && inputs[q] == inputs[slot] && param_name(q).and_then(role) == Some(an)
            }) else {
                continue;
            };
            // `f(name, name)`: the namesake parameter already gets its name,
            // so nothing is transposed, one value fills two roles.
            let at_other = received(args, offset, other);
            if at_other == Some(an)
                || at_other.is_some_and(|there| {
                    mirrored_in_condition(cx, expr, def, (slot, there), (other, an))
                })
            {
                continue;
            }
            crossings.push(Crossing {
                arg: i,
                name,
                bound_to,
                names_param: other,
                ty: inputs[slot],
            });
        }
        let fn_name = cx.tcx.item_name(def);
        // The callee's signature is the other half of the crossing, when this
        // crate has it to show.
        let report = |span: Span, msg: String, help: String| {
            if def.is_local() {
                emit_with_note(
                    cx,
                    ARG_NAMED_LIKE_OTHER_PARAM,
                    span,
                    msg,
                    cx.tcx.def_span(def),
                    format!("`{fn_name}` declares them in this order"),
                    help,
                );
            } else {
                emit(cx, ARG_NAMED_LIKE_OTHER_PARAM, span, msg, help);
            }
        };
        let mut reported = vec![false; crossings.len()];
        for (k, c) in crossings.iter().enumerate() {
            if reported[k] {
                continue;
            }
            reported[k] = true;
            let partner = (k + 1..crossings.len()).find(|&j| {
                let d = &crossings[j];
                d.arg + offset == c.names_param && c.arg + offset == d.names_param
            });
            let other = param_name(c.names_param).unwrap_or(c.name);
            let types = format!(
                "Distinct types for `{}` and `{other}` would make the next swap a compile error",
                c.bound_to,
            );
            if let Some(j) = partner {
                reported[j] = true;
                let d = &crossings[j];
                report(
                    expr.span,
                    format!(
                        "`{}` and `{}` are passed where `{fn_name}` expects `{}` then `{}`. All \
                         are `{}`, so the swap compiles",
                        c.name, d.name, c.bound_to, d.bound_to, c.ty,
                    ),
                    format!(
                        "swap the arguments. If the call is right, rename the values to match. \
                         {types}"
                    ),
                );
            } else {
                report(
                    args[c.arg].span,
                    format!(
                        "`{}` is passed as `{fn_name}`'s `{}` parameter, but `{fn_name}` also \
                         has a parameter named `{other}` and both are `{}`. This looks like an \
                         argument in the wrong position",
                        c.name, c.bound_to, c.ty,
                    ),
                    format!(
                        "move `{}` to the `{other}` position. If it really is the `{}` here, \
                         rename it. {types}",
                        c.name, c.bound_to,
                    ),
                );
            }
        }
    }
}

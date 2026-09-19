use clippy_utils::diagnostics::span_lint_and_then;
use clippy_utils::is_trait_impl_item;
use rustc_hir::intravisit::{self, Visitor};
use rustc_hir::{Expr, ExprKind, ImplItem, ImplItemKind, LangItem, Node};
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::Span;

use super::hir_refs;

rustc_session::declare_lint! {
    /// Warns when a `Drop::drop` implementation contains operations that can
    /// panic, since panicking during unwinding causes an immediate process abort.
    pub PANIC_IN_DROP,
    Warn,
    "panic-able expression in `Drop` impl -- this will abort during unwinding"
}

/// If the expression is `std::thread::panicking()` or
/// `!std::thread::panicking()`, returns `Some(negated)` where `negated` is
/// `true` for the `!...` form. Returns `None` otherwise.
///
/// Uses `qpath_res` rather than `clippy_utils::fn_def_id` because the latter
/// calls `cx.typeck_results()`, which panics in `check_impl_item` contexts
/// that run outside an active body.
fn is_panicking_guard<'tcx>(cx: &LateContext<'tcx>, cond: &Expr<'tcx>) -> Option<bool> {
    let (negated, inner) = if let ExprKind::Unary(rustc_hir::UnOp::Not, inner) = &cond.kind {
        (true, *inner)
    } else {
        (false, cond)
    };
    // Cheap name prefilter -- only then pay for `def_path_str`.
    if let ExprKind::Call(callee, _) = &inner.kind
        && let ExprKind::Path(qpath) = &callee.kind
        && let last = match qpath {
            rustc_hir::QPath::Resolved(_, p) => p.segments.last().map(|s| s.ident),
            rustc_hir::QPath::TypeRelative(_, seg) => Some(seg.ident),
        }
        && last.is_some_and(|i| i.as_str() == "panicking")
        && let Some(def_id) = cx.qpath_res(qpath, callee.hir_id).opt_def_id()
        && cx.tcx.def_path_str(def_id) == "std::thread::panicking"
    {
        return Some(negated);
    }
    None
}

struct DropPanicFinder<'a, 'tcx> {
    cx: &'a LateContext<'tcx>,
    typeck: &'a rustc_middle::ty::TypeckResults<'tcx>,
    findings: Vec<(Span, &'static str)>,
}

// No NestedFilter -- stored/passed closures don't run during drop; only
// immediately-invoked ones (IIFEs) do, and those are handled explicitly.
impl<'tcx> Visitor<'tcx> for DropPanicFinder<'_, 'tcx> {
    fn visit_expr(&mut self, expr: &'tcx Expr<'tcx>) {
        if matches!(expr.kind, ExprKind::Closure(_)) {
            if let Some(body) = hir_refs::iife_closure_body(self.cx.tcx, expr) {
                intravisit::walk_body(self, body);
            }
            return;
        }

        // `if [!]std::thread::panicking() { ... } else { ... }` -- exactly one branch
        // runs while unwinding. Skip the safe (not-unwinding) branch, but still
        // visit the unwinding branch: a panic there is the double-panic-abort.
        if let ExprKind::If(cond, then, opt_else) = &expr.kind
            && let Some(negated) = is_panicking_guard(self.cx, cond)
        {
            // `!panicking()` -> `then` is safe, `else` unwinds.
            // `panicking()`  -> `then` unwinds, `else` is safe.
            let unwinding = if negated { *opt_else } else { Some(*then) };
            if let Some(branch) = unwinding {
                intravisit::walk_expr(self, branch);
            }
            return;
        }

        // Check panic macros before method calls to avoid reporting macro internals.
        if expr.span.from_expansion()
            && let Some((call_site, kind)) = hir_refs::find_panic_macro(expr.span)
        {
            self.findings.push((call_site, kind.desc()));
            return;
        }

        if let Some(finding) = hir_refs::panicking_unwrap_or_expect(self.cx, self.typeck, expr) {
            self.findings.push(finding);
        }

        intravisit::walk_expr(self, expr);
    }
}

/// Returns `true` if the parent `impl` block is `impl Drop for T`.
fn is_drop_impl<'tcx>(cx: &LateContext<'tcx>, impl_item: &'tcx ImplItem<'tcx>) -> bool {
    // Chose trait DefId comparison over string matching because it works
    // regardless of imports or re-exports.
    let parent_id = cx.tcx.parent_hir_id(impl_item.hir_id());
    let Node::Item(item) = cx.tcx.hir_node(parent_id) else {
        return false;
    };
    let rustc_hir::ItemKind::Impl(impl_block) = &item.kind else {
        return false;
    };
    let Some(trait_header) = impl_block.of_trait else {
        return false;
    };
    let Some(trait_def_id) = trait_header.trait_ref.trait_def_id() else {
        return false;
    };
    cx.tcx.is_lang_item(trait_def_id, LangItem::Drop)
}

pub struct PanicInDrop;

impl PanicInDrop {
    pub const fn new() -> Self {
        Self
    }
}

rustc_session::impl_lint_pass!(PanicInDrop => [PANIC_IN_DROP]);

impl<'tcx> LateLintPass<'tcx> for PanicInDrop {
    fn check_impl_item(&mut self, cx: &LateContext<'tcx>, impl_item: &'tcx ImplItem<'tcx>) {
        let ImplItemKind::Fn(_sig, body_id) = &impl_item.kind else {
            return;
        };

        if impl_item.span.from_expansion()
            || impl_item.ident.as_str() != "drop"
            // Cheap pre-check before `is_drop_impl`'s HIR parent walk.
            || !is_trait_impl_item(cx, impl_item.hir_id())
            || !is_drop_impl(cx, impl_item)
        {
            return;
        }

        let body = cx.tcx.hir_body(*body_id);
        let typeck = cx.tcx.typeck(impl_item.owner_id.def_id);
        let mut finder = DropPanicFinder {
            cx,
            typeck,
            findings: Vec::new(),
        };
        intravisit::walk_body(&mut finder, body);

        if finder.findings.is_empty() {
            return;
        }

        span_lint_and_then(
            cx,
            PANIC_IN_DROP,
            impl_item.span,
            "panic-able expression in `Drop` impl -- this will abort during unwinding",
            |diag| {
                for (span, desc) in &finder.findings {
                    diag.span_note(
                        *span,
                        format!(
                            "`{desc}` can panic -- handle the error or use \
                             `let _ = ...` to ignore it"
                        ),
                    );
                }
                diag.help(
                    "panicking in `drop()` while already unwinding causes an \
                     immediate process abort",
                );
            },
        );
    }
}

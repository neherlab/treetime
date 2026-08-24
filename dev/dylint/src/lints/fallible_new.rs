use clippy_utils::diagnostics::span_lint_and_then;
use clippy_utils::is_trait_impl_item;
use rustc_hir::intravisit::{self, Visitor};
use rustc_hir::{Expr, ExprKind, ImplItem, ImplItemKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::{Span, sym};

use super::hir_refs;
use crate::config::FallibleNewConfig;

rustc_session::declare_lint! {
    /// Warns when a `fn new()` constructor contains operations that can panic,
    /// suggesting it return `Result` or be renamed to convey fallibility.
    pub FALLIBLE_NEW,
    Deny,
    "constructor `new` can panic -- consider returning `Result` or renaming to `try_new`"
}

/// Returns `true` if the function's return type resolves to `Result<_, _>`
/// (after resolving type aliases).
fn returns_result<'tcx>(cx: &LateContext<'tcx>, impl_item: &'tcx ImplItem<'tcx>) -> bool {
    let def_id = impl_item.owner_id.to_def_id();
    let fn_sig = cx.tcx.fn_sig(def_id).instantiate_identity().skip_norm_wip();
    let ret_ty = fn_sig.output().skip_binder();
    if let ty::Adt(adt, _) = ret_ty.kind() {
        return cx.tcx.is_diagnostic_item(sym::Result, adt.did());
    }
    false
}

struct PanicFinder<'a, 'tcx> {
    cx: &'a LateContext<'tcx>,
    typeck: &'a rustc_middle::ty::TypeckResults<'tcx>,
    findings: Vec<(Span, &'static str)>,
}

// No NestedFilter -- stored/returned closures don't run during construction;
// only immediately-invoked ones (IIFEs) do, and those are handled explicitly.
impl<'tcx> Visitor<'tcx> for PanicFinder<'_, 'tcx> {
    fn visit_expr(&mut self, expr: &'tcx Expr<'tcx>) {
        if matches!(expr.kind, ExprKind::Closure(_)) {
            if let Some(body) = hir_refs::iife_closure_body(self.cx.tcx, expr) {
                intravisit::walk_body(self, body);
            }
            return;
        }

        if let Some(finding) = hir_refs::panicking_unwrap_or_expect(self.cx, self.typeck, expr) {
            self.findings.push(finding);
        }

        // For assert macros, walk into the expansion so inner unwrap()/expect()
        // are still caught but the assert itself is not reported.
        if expr.span.from_expansion()
            && let Some((call_site, kind)) = hir_refs::find_panic_macro(expr.span)
            && matches!(kind, hir_refs::PanicMacro::Panic | hir_refs::PanicMacro::Unreachable)
        {
            self.findings.push((call_site, kind.desc()));
            return;
        }

        intravisit::walk_expr(self, expr);
    }
}

pub struct FallibleNew {
    check_new_variants: bool,
}

impl FallibleNew {
    pub fn new() -> Self {
        let config: FallibleNewConfig = dylint_linting::config_or_default("fallible_new");
        Self {
            check_new_variants: config.check_new_variants,
        }
    }
}

rustc_session::impl_lint_pass!(FallibleNew => [FALLIBLE_NEW]);

impl<'tcx> LateLintPass<'tcx> for FallibleNew {
    fn check_impl_item(&mut self, cx: &LateContext<'tcx>, impl_item: &'tcx ImplItem<'tcx>) {
        let ImplItemKind::Fn(_sig, body_id) = &impl_item.kind else {
            return;
        };

        if impl_item.span.from_expansion() {
            return;
        }

        let name = impl_item.ident.as_str();

        if name != "new" && !(self.check_new_variants && name.starts_with("new_")) {
            return;
        }

        // Skip trait impls (signature dictated by trait).
        if is_trait_impl_item(cx, impl_item.hir_id()) || returns_result(cx, impl_item) {
            return;
        }

        let body = cx.tcx.hir_body(*body_id);
        let typeck = cx.tcx.typeck(impl_item.owner_id.def_id);
        let mut finder = PanicFinder {
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
            FALLIBLE_NEW,
            impl_item.span,
            format!(
                "constructor `{name}` can panic -- consider returning `Result` or renaming to `try_{name}`"
            ),
            |diag| {
                for (span, desc) in &finder.findings {
                    diag.span_note(
                        *span,
                        format!("`{desc}` can panic -- use `?` with a `Result` return type instead"),
                    );
                }
            },
        );
    }
}

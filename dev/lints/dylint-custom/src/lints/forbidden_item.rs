//! Flags uses of items that `dylint.toml` forbids in the crate being linted.
//!
//! Each `[[forbidden_item]]` entry names items by absolute path, the crates in
//! which they are forbidden, and the reason that the diagnostic shows. A use is
//! any path that resolves to a forbidden item: a type in any position
//! (signature, field, `let` annotation, generic argument, `impl` target), an
//! import, a call, a method call, or a function reference. The paths resolve at
//! lint time in the dependency graph of the linted crate, so a path that does
//! not resolve, for example to an item that its crate does not define yet,
//! forbids nothing. Test code (`#[test]` functions and items under
//! `#[cfg(test)]`) and code that a macro expands are not checked.

use clippy_utils::diagnostics::span_lint_and_then;
use clippy_utils::is_in_test;
use clippy_utils::paths::{PathNS, lookup_path_str};
use rustc_data_structures::fx::{FxHashMap, FxHashSet};
use rustc_hir::def::Res;
use rustc_hir::def_id::{DefId, LOCAL_CRATE};
use rustc_hir::{Expr, ExprKind, HirId, Path, QPath};
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::Span;

use crate::config::ForbiddenItemRule;
use crate::lints::hir_refs::owning_item;

rustc_session::declare_lint! {
    /// Flags a use of an item that a `[[forbidden_item]]` entry of `dylint.toml`
    /// forbids in the crate being linted, outside test code, and shows the
    /// reason of the entry.
    pub FORBIDDEN_ITEM,
    Warn,
    "use of an item that `dylint.toml` forbids in this crate"
}

pub struct ForbiddenItem {
    rules: Vec<ForbiddenItemRule>,
    forbidden: FxHashMap<DefId, Forbidden>,
    reported: FxHashSet<Span>,
}

impl ForbiddenItem {
    pub fn new() -> Self {
        Self {
            rules: dylint_linting::config_or_default("forbidden_item"),
            forbidden: FxHashMap::default(),
            reported: FxHashSet::default(),
        }
    }

    fn report(&mut self, cx: &LateContext<'_>, def_id: DefId, span: Span, hir_id: HirId) {
        if span.from_expansion() || is_in_test(cx.tcx, hir_id) {
            return;
        }
        let Some(forbidden) = self.forbidden.get(&owning_item(cx, def_id)) else {
            return;
        };
        if !self.reported.insert(span) {
            return;
        }
        let crate_name = cx.tcx.crate_name(LOCAL_CRATE);
        span_lint_and_then(
            cx,
            FORBIDDEN_ITEM,
            span,
            format!("`{}` is forbidden in crate `{crate_name}`", forbidden.path),
            |diag| {
                diag.note(forbidden.reason.clone());
            },
        );
    }
}

rustc_session::impl_lint_pass!(ForbiddenItem => [FORBIDDEN_ITEM]);

impl<'tcx> LateLintPass<'tcx> for ForbiddenItem {
    fn check_crate(&mut self, cx: &LateContext<'tcx>) {
        let crate_name = cx.tcx.crate_name(LOCAL_CRATE);
        for rule in &self.rules {
            if !rule.crates.iter().any(|name| name.as_str() == crate_name.as_str()) {
                continue;
            }
            for path in &rule.items {
                for ns in [PathNS::Type, PathNS::Value] {
                    for def_id in lookup_path_str(cx.tcx, ns, path) {
                        self.forbidden.entry(def_id).or_insert_with(|| Forbidden {
                            path: path.clone(),
                            reason: rule.reason.clone(),
                        });
                    }
                }
            }
        }
    }

    fn check_path(&mut self, cx: &LateContext<'tcx>, path: &Path<'tcx>, hir_id: HirId) {
        if self.forbidden.is_empty() {
            return;
        }
        if let Res::Def(_, def_id) = path.res {
            self.report(cx, def_id, path.span, hir_id);
        }
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        if self.forbidden.is_empty() {
            return;
        }
        let resolved = match expr.kind {
            ExprKind::MethodCall(segment, ..) => cx
                .typeck_results()
                .type_dependent_def_id(expr.hir_id)
                .map(|def_id| (def_id, segment.ident.span.with_hi(expr.span.hi()))),
            ExprKind::Path(ref qpath @ QPath::TypeRelative(..)) => cx
                .qpath_res(qpath, expr.hir_id)
                .opt_def_id()
                .map(|def_id| (def_id, expr.span)),
            _ => None,
        };
        if let Some((def_id, span)) = resolved {
            self.report(cx, def_id, span, expr.hir_id);
        }
    }
}

struct Forbidden {
    path: String,
    reason: String,
}

#[cfg(test)]
mod tests {
    use dylint_testing::ui::Test;

    const DYLINT_TOML: &str = r#"
[[forbidden_item]]
items = [
  "graph::schedule::TreeSchedule",
  "graph::schedule::GraphSchedule::into_tree",
  "graph::schedule::NotDefinedYet",
]
crates = ["core_algorithms"]
reason = "core algorithms take a `GraphSchedule`; the tree handle serves tree-only output formats"

[no_comments]
allowed_comment_prefixes = ["aux-build:", "compile-flags:"]
"#;

    #[test]
    fn forbidden_item_ui() {
        Test::src_base(env!("CARGO_PKG_NAME"), "ui/forbidden_item")
            .rustc_flags(["--edition=2024"])
            .dylint_toml(DYLINT_TOML)
            .run();
    }
}

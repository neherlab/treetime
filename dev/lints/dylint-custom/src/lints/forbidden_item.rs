//! Flags uses of items that `dylint.toml` forbids in the crate being linted.
//!
//! Each `[[forbidden_item]]` entry names items by absolute path, the crates in
//! which they are forbidden, and the reason that the diagnostic shows. A use is
//! any path that resolves to a forbidden item: a type in any position
//! (signature, field, `let` annotation, generic argument, `impl` target), an
//! import, a call, a method call, or a function reference. A use that a macro
//! expands to, from this crate or another, is reported at the outermost macro
//! call in the linted crate; for a macro of the linted crate, a note points at
//! the use in the macro definition.
//! Test code (`#[test]` functions and items under `#[cfg(test)]`) is not
//! checked.
//!
//! The paths resolve at lint time in the dependency graph of the linted crate.
//! A path that resolves to no item, after a rename, a move, or a typo, is
//! reported at the crate root of every crate that the entry names, because the
//! entry then forbids nothing there.

use clippy_utils::diagnostics::{span_lint_and_then, span_lint_hir_and_then};
use clippy_utils::is_in_test;
use clippy_utils::paths::{PathNS, lookup_path_str};
use rustc_data_structures::fx::{FxHashMap, FxHashSet};
use rustc_hir::def::Res;
use rustc_hir::def_id::{CRATE_DEF_ID, DefId, LOCAL_CRATE};
use rustc_hir::{CRATE_HIR_ID, Expr, ExprKind, HirId, Path, QPath};
use rustc_lint::{LateContext, LateLintPass, LintContext as _};
use rustc_span::Span;

use crate::config::ForbiddenItemRule;
use crate::lints::hir_refs::owning_item;

rustc_session::declare_lint! {
    /// Flags a use of an item that a `[[forbidden_item]]` entry of `dylint.toml`
    /// forbids in the crate being linted, outside test code, and shows the
    /// reason of the entry. Also flags an entry path that resolves to no item.
    pub FORBIDDEN_ITEM,
    Warn,
    "use of an item that `dylint.toml` forbids in this crate"
}

pub struct ForbiddenItem {
    rules: Vec<ForbiddenItemRule>,
    forbidden: FxHashMap<DefId, Forbidden>,
    reported: FxHashSet<(Span, DefId)>,
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
        if is_in_test(cx.tcx, hir_id) {
            return;
        }
        let item = owning_item(cx, def_id);
        let Some(forbidden) = self.forbidden.get(&item) else {
            return;
        };
        let call_site = span.source_callsite();
        if !self.reported.insert((call_site, item)) {
            return;
        }
        let crate_name = cx.tcx.crate_name(LOCAL_CRATE);
        span_lint_and_then(
            cx,
            FORBIDDEN_ITEM,
            call_site,
            format!("`{}` is forbidden in crate `{crate_name}`", forbidden.path),
            |diag| {
                if call_site != span {
                    if span.in_external_macro(cx.sess().source_map()) {
                        diag.note("the macro called here, from another crate, expands to the use");
                    } else {
                        diag.span_note(span, "the macro expands to this use");
                    }
                }
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
                let def_ids = [PathNS::Type, PathNS::Value]
                    .into_iter()
                    .flat_map(|ns| lookup_path_str(cx.tcx, ns, path))
                    .collect::<Vec<_>>();
                if def_ids.is_empty() {
                    report_unresolved(cx, path);
                }
                for def_id in def_ids {
                    self.forbidden.entry(def_id).or_insert_with(|| Forbidden {
                        path: path.clone(),
                        reason: rule.reason.clone(),
                    });
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

fn report_unresolved(cx: &LateContext<'_>, path: &str) {
    let crate_name = cx.tcx.crate_name(LOCAL_CRATE);
    span_lint_hir_and_then(
        cx,
        FORBIDDEN_ITEM,
        CRATE_HIR_ID,
        cx.tcx.def_span(CRATE_DEF_ID).shrink_to_lo(),
        format!("`{path}` resolves to no item, so its `[[forbidden_item]]` entry forbids nothing in crate `{crate_name}`"),
        |diag| {
            diag.help("correct the path in `dylint.toml` after a rename or move, or remove it from the entry");
        },
    );
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
  "graph::schedule::RenamedTreeSchedule",
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

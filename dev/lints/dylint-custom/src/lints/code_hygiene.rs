//! Structural hygiene lints on items: `use super::`, `use` inside a function,
//! nested `fn` items, and versioned/aged identifier names.

use rustc_hir::{Item, ItemKind, Node};
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::kw;

use clippy_utils::diagnostics::span_lint_and_help;

use crate::lints::suppression::is_hir_in_test_zone;

rustc_session::declare_lint! {
    /// Forbids `use super::...` imports in production code. Paths should start
    /// from `crate::` so a module can be read and moved without tracing relative
    /// ancestors. Test modules may import the module under test with `super::`.
    pub SUPER_IMPORT,
    Warn,
    "`use super::` import -- import from `crate::` instead"
}

rustc_session::declare_lint! {
    /// Forbids `use` declarations inside a function body. Imports belong at the
    /// top of the file so every dependency of a module is visible in one place.
    pub LOCAL_USE,
    Warn,
    "`use` inside a function -- move the import to the top of the file"
}

rustc_session::declare_lint! {
    /// Forbids `fn` items nested inside another function body. Extract the helper
    /// to module level; a closure is the tool for capturing local state.
    pub NESTED_FUNCTION,
    Warn,
    "nested `fn` item -- extract it to module level"
}

rustc_session::declare_lint! {
    /// Forbids version- or age-suffixed identifiers such as `_v2`, `_old`, and
    /// `_new`. These encode edit history in the name; rename to describe the role.
    pub VERSIONED_NAME,
    Warn,
    "version- or age-suffixed name -- rename to describe the role, not the revision"
}

pub struct CodeHygiene;

impl CodeHygiene {
    pub const fn new() -> Self {
        Self
    }
}

rustc_session::impl_lint_pass!(CodeHygiene => [SUPER_IMPORT, LOCAL_USE, NESTED_FUNCTION, VERSIONED_NAME]);

impl<'tcx> LateLintPass<'tcx> for CodeHygiene {
    fn check_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx Item<'tcx>) {
        if item.span.from_expansion() {
            return;
        }

        let in_block = matches!(cx.tcx.parent_hir_node(item.hir_id()), Node::Block(_));

        match &item.kind {
            ItemKind::Use(path, _) => {
                if path
                    .segments
                    .first()
                    .is_some_and(|seg| seg.ident.name == kw::Super)
                    && !is_hir_in_test_zone(cx, item.hir_id())
                {
                    span_lint_and_help(
                        cx,
                        SUPER_IMPORT,
                        item.span,
                        "`use super::` import is banned",
                        None,
                        "import the item from `crate::` with an absolute path",
                    );
                }
                if in_block {
                    span_lint_and_help(
                        cx,
                        LOCAL_USE,
                        item.span,
                        "`use` declaration inside a function body",
                        None,
                        "move the import to the top of the file, after the `mod` declarations",
                    );
                }
            }
            ItemKind::Fn { .. } if in_block => {
                span_lint_and_help(
                    cx,
                    NESTED_FUNCTION,
                    item.span,
                    "`fn` item nested inside another function",
                    None,
                    "extract the helper to module level; use a closure to capture local state",
                );
            }
            _ => {}
        }

        if let Some(name) = cx.tcx.opt_item_name(item.owner_id.to_def_id())
            && is_versioned_name(name.as_str())
        {
            span_lint_and_help(
                cx,
                VERSIONED_NAME,
                item.span,
                format!("`{name}` carries a version or age suffix"),
                None,
                "rename to describe the item's role rather than its revision",
            );
        }
    }
}

/// Returns `true` if `name` ends in a revision/age marker: `_v<N>`, `_old`,
/// `_new`, `_fixed`, or `_tmp` (case-insensitive on the suffix word).
fn is_versioned_name(name: &str) -> bool {
    let lower = name.to_ascii_lowercase();
    for suffix in ["_old", "_new", "_fixed", "_tmp"] {
        if lower == suffix.trim_start_matches('_') || lower.ends_with(suffix) {
            return true;
        }
    }
    if let Some(rest) = lower.rsplit_once("_v").map(|(_, digits)| digits) {
        return !rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit());
    }
    false
}

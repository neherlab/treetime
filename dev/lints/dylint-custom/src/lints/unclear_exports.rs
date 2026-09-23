use clippy_utils::diagnostics::span_lint_and_help;
use rustc_hir::{Item, ItemKind, UsePath, UseKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::kw;

use crate::lints::suppression::is_hir_in_test_zone;

rustc_session::declare_lint! {
    /// Forbids glob imports (`use foo::*`) and renamed imports (`use foo::Bar as Baz`).
    /// Every imported name must be listed explicitly under its original name so the
    /// module's API surface is intentional, auditable, and traceable. As in clippy's
    /// `wildcard_imports`, a `prelude` module and `use super::*` in a test module
    /// may be glob-imported: preludes exist to be imported whole, and a test module
    /// tests its parent's items.
    pub UNCLEAR_EXPORTS,
    Warn,
    "unclear exports -- glob imports and renamed imports are banned"
}

const GLOB_MSG: &str =
    "glob imports (`use foo::*`) are banned -- list each imported name explicitly";
const GLOB_HELP: &str = "replace `use foo::*` with an explicit list: `use foo::{Bar, Baz}`";

const RENAME_MSG: &str =
    "renamed imports (`use foo::Bar as Baz`) are banned -- use the original name";
const RENAME_HELP: &str =
    "import the item under its original name, or create a type alias if a new name is truly needed";

pub struct UnclearExports;

impl UnclearExports {
    pub const fn new() -> Self {
        Self
    }
}

rustc_session::impl_lint_pass!(UnclearExports => [UNCLEAR_EXPORTS]);

impl<'tcx> LateLintPass<'tcx> for UnclearExports {
    fn check_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx Item<'tcx>) {
        if item.span.from_expansion() {
            return;
        }

        let ItemKind::Use(path, kind) = &item.kind else {
            return;
        };

        if *kind == UseKind::Glob {
            let exempt = is_prelude_glob(path) || (is_super_glob(path) && is_hir_in_test_zone(cx, item.hir_id()));
            if !exempt {
                span_lint_and_help(cx, UNCLEAR_EXPORTS, item.span, GLOB_MSG, None, GLOB_HELP);
            }
            return;
        }

        // `Bar as _` is the idiomatic way to import a trait for its methods
        // without binding the name -- not a rename.
        if let UseKind::Single(bound_ident) = *kind {
            if let Some(last_seg) = path.segments.last() {
                let original = last_seg.ident.name;
                let bound = bound_ident.name;
                if original != bound && bound != kw::Underscore {
                    span_lint_and_help(
                        cx,
                        UNCLEAR_EXPORTS,
                        item.span,
                        RENAME_MSG,
                        None,
                        RENAME_HELP,
                    );
                }
            }
        }
    }
}

/// Returns `true` for `use some::path::prelude::*`.
fn is_prelude_glob(path: &UsePath<'_>) -> bool {
    path.segments.last().is_some_and(|segment| segment.ident.as_str() == "prelude")
}

/// Returns `true` for `use super::*`.
fn is_super_glob(path: &UsePath<'_>) -> bool {
    matches!(path.segments, [segment] if segment.ident.name == kw::Super)
}

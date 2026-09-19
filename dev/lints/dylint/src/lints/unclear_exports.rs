use clippy_utils::diagnostics::span_lint_and_help;
use rustc_hir::{Item, ItemKind, UseKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::kw;

rustc_session::declare_lint! {
    /// Forbids glob imports (`use foo::*`) and renamed imports (`use foo::Bar as Baz`).
    /// Every imported name must be listed explicitly under its original name so the
    /// module's API surface is intentional, auditable, and traceable.
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
            span_lint_and_help(cx, UNCLEAR_EXPORTS, item.span, GLOB_MSG, None, GLOB_HELP);
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

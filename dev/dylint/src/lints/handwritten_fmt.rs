//! Flags hand-written `Display` and `Debug` implementations, which should be
//! derived (`Debug`) or replaced by a dedicated method rather than written by
//! hand.

use rustc_hir::def::Res;
use rustc_hir::{Item, ItemKind};
use rustc_lint::{LateContext, LateLintPass};

use clippy_utils::diagnostics::span_lint_and_help;

rustc_session::declare_lint! {
    /// Flags a hand-written `impl std::fmt::Display` or `impl std::fmt::Debug`.
    /// `Debug` should be derived; a `Display`-like rendering belongs in a named
    /// method so its behavior is explicit and testable.
    pub HANDWRITTEN_FMT_IMPL,
    Warn,
    "hand-written Display/Debug impl -- derive Debug or use a named method"
}

pub struct HandwrittenFmt;

impl HandwrittenFmt {
    pub const fn new() -> Self {
        Self
    }
}

rustc_session::impl_lint_pass!(HandwrittenFmt => [HANDWRITTEN_FMT_IMPL]);

impl<'tcx> LateLintPass<'tcx> for HandwrittenFmt {
    fn check_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx Item<'tcx>) {
        if item.span.from_expansion() {
            return;
        }
        let ItemKind::Impl(impl_block) = &item.kind else {
            return;
        };
        let Some(of_trait) = &impl_block.of_trait else {
            return;
        };
        let Res::Def(_, def_id) = of_trait.trait_ref.path.res else {
            return;
        };
        let path = cx.tcx.def_path_str(def_id);
        let which = if path.ends_with("fmt::Display") {
            "Display"
        } else if path.ends_with("fmt::Debug") {
            "Debug"
        } else {
            return;
        };
        span_lint_and_help(
            cx,
            HANDWRITTEN_FMT_IMPL,
            item.span,
            format!("hand-written `{which}` implementation"),
            None,
            "derive `Debug`, or move a `Display`-style rendering into a named method",
        );
    }
}

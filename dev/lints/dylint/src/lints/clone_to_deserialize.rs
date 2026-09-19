//! Flags a value cloned only to be handed to a deserializer, which borrows its
//! input; the clone is a wasted allocation.

use rustc_hir::{Expr, ExprKind};
use rustc_lint::{LateContext, LateLintPass};

use clippy_utils::diagnostics::span_lint_and_help;

rustc_session::declare_lint! {
    /// Flags `T::from_str(&value.clone())` and the other serde entry points
    /// (`from_slice`, `from_reader`, `from_value`) applied to a freshly cloned
    /// value. Deserializers borrow their input, so the clone is unnecessary.
    pub VALUE_CLONED_TO_DESERIALIZE,
    Warn,
    "value cloned only to deserialize -- pass a reference to the original"
}

pub struct CloneToDeserialize;

impl CloneToDeserialize {
    pub const fn new() -> Self {
        Self
    }
}

rustc_session::impl_lint_pass!(CloneToDeserialize => [VALUE_CLONED_TO_DESERIALIZE]);

impl<'tcx> LateLintPass<'tcx> for CloneToDeserialize {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        if expr.span.from_expansion() {
            return;
        }
        let ExprKind::Call(callee, args) = &expr.kind else {
            return;
        };
        if !is_deserialize_entry_point(callee) {
            return;
        }
        for arg in *args {
            if let Some(clone_span) = clone_call_span(arg) {
                span_lint_and_help(
                    cx,
                    VALUE_CLONED_TO_DESERIALIZE,
                    clone_span,
                    "value is cloned only to be deserialized",
                    None,
                    "a deserializer borrows its input; pass a reference to the original value",
                );
            }
        }
    }
}

/// Returns `true` if `callee` is a serde deserialization free function.
fn is_deserialize_entry_point(callee: &Expr<'_>) -> bool {
    let ExprKind::Path(qpath) = &callee.kind else {
        return false;
    };
    let rustc_hir::QPath::Resolved(_, path) = qpath else {
        return false;
    };
    path.segments.last().is_some_and(|seg| {
        matches!(
            seg.ident.as_str(),
            "from_str" | "from_slice" | "from_reader" | "from_value"
        )
    })
}

/// If `arg` is `x.clone()` or `&x.clone()`, returns the span of the clone call.
fn clone_call_span(arg: &Expr<'_>) -> Option<rustc_span::Span> {
    let inner = match &arg.kind {
        ExprKind::AddrOf(_, _, inner) => inner,
        _ => arg,
    };
    if let ExprKind::MethodCall(method, _, _, span) = &inner.kind
        && method.ident.as_str() == "clone"
    {
        return Some(*span);
    }
    None
}

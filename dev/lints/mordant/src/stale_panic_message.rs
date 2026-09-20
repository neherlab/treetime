use rustc_ast::LitKind;
use rustc_hir::{Expr, ExprKind};
use rustc_lint::{LateContext, LateLintPass};

use crate::baseline::emit;
use crate::claims::{self, DefNames};

rustc_session::declare_lint! {
    /// Flags a panic-family message (`panic!`, `unreachable!`, `assert!`,
    /// `expect`) that mentions a backticked identifier when nothing with
    /// that name exists in the file, this crate, or any crate it links.
    /// Whoever hits the panic is sent looking for something that is gone.
    pub STALE_PANIC_MESSAGE,
    Warn,
    "panic or assert message names an identifier that no longer exists"
}

#[derive(Default)]
pub struct StalePanicMessage {
    defs: Option<DefNames>,
}

rustc_session::impl_lint_pass!(StalePanicMessage => [STALE_PANIC_MESSAGE]);

const PANIC_MACROS: &[&str] = &[
    "panic",
    "unreachable",
    "todo",
    "unimplemented",
    "assert",
    "assert_eq",
    "assert_ne",
    "debug_assert",
    "debug_assert_eq",
    "debug_assert_ne",
];

impl StalePanicMessage {
    fn check_message(&self, cx: &LateContext<'_>, text: &str, at: rustc_span::Span) {
        let Some(defs) = &self.defs else { return };
        let code = claims::file_code_only(cx, at);
        // The message text is itself part of the file's code lines; blank
        // every copy of it so the message cannot vouch for its own names.
        let code = code.replace(text, "");
        for ident in claims::backticked_idents(text) {
            if !claims::word_in(&code, &ident) && !defs.contains(cx, &ident) {
                emit(
                    cx,
                    STALE_PANIC_MESSAGE,
                    at,
                    format!(
                        "this message mentions `{ident}`, but nothing with that name exists in this file, this crate, or any crate it links. Whoever hits the panic is sent looking for something that is gone"
                    ),
                    format!("rewrite the message in terms of what replaced `{ident}`"),
                );
            }
        }
    }
}

impl<'tcx> LateLintPass<'tcx> for StalePanicMessage {
    fn check_crate(&mut self, cx: &LateContext<'tcx>) {
        self.defs = Some(DefNames::collect(cx));
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        match &expr.kind {
            // `expect("...")` takes the message as a plain argument.
            ExprKind::MethodCall(seg, _, [arg], _) if seg.ident.as_str() == "expect" => {
                if let ExprKind::Lit(lit) = &arg.kind
                    && let LitKind::Str(s, _) = lit.node
                    && !arg.span.from_expansion()
                {
                    self.check_message(cx, s.as_str(), arg.span);
                }
            }
            // Panic-family macros: the user-written message literal keeps its
            // call-site span, so the macro shows up on the spans of the
            // EXPANDED ancestors it sits inside, not on the literal itself.
            ExprKind::Lit(lit) => {
                let LitKind::Str(s, _) = lit.node else {
                    return;
                };
                if expr.span.from_expansion() {
                    return;
                }
                let mut from_panic = false;
                for (_, node) in cx.tcx.hir_parent_iter(expr.hir_id).take(8) {
                    let rustc_hir::Node::Expr(parent) = node else {
                        break;
                    };
                    if parent.span.from_expansion()
                        && clippy_utils::macros::macro_backtrace(parent.span).any(|mac| {
                            let name = cx.tcx.item_name(mac.def_id);
                            PANIC_MACROS.contains(&name.as_str())
                        })
                    {
                        from_panic = true;
                        break;
                    }
                }
                if from_panic {
                    self.check_message(cx, s.as_str(), expr.span);
                }
            }
            _ => {}
        }
    }
}

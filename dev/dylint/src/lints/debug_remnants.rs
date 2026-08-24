use clippy_utils::diagnostics::span_lint_and_help;
use rustc_data_structures::fx::FxHashSet;
use rustc_hir::Expr;
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::{ExpnKind, Span};

use crate::config::{DebugRemnantsConfig, LogFramework};
use crate::lints::suppression::is_in_test_zone;

rustc_session::declare_lint! {
    /// Flags debugging macros (`println!`, `print!`, `eprintln!`, `dbg!`) and
    /// suggests structured logging replacements (`tracing` or `log`).
    pub DEBUG_REMNANTS,
    Warn,
    "debug macro in committed code -- replace with structured logging"
}

pub struct DebugRemnants {
    framework: LogFramework,
    /// Dedup: one diagnostic per macro call site, not per expanded HIR node.
    seen_callsites: FxHashSet<Span>,
}

impl DebugRemnants {
    pub fn new() -> Self {
        let config: DebugRemnantsConfig = dylint_linting::config_or_default("debug_remnants");
        Self {
            framework: config.suggested_framework,
            seen_callsites: FxHashSet::default(),
        }
    }
}

rustc_session::impl_lint_pass!(DebugRemnants => [DEBUG_REMNANTS]);

impl<'tcx> LateLintPass<'tcx> for DebugRemnants {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        if !expr.span.from_expansion() {
            return;
        }

        let expn_data = expr.span.ctxt().outer_expn_data();
        let ExpnKind::Macro(_, macro_name) = expn_data.kind else {
            return;
        };

        // `eprint!` intentionally excluded: typically used for progress
        // indicators and prompts, not leftover debugging.
        let level = match macro_name.as_str() {
            "println" | "print" => "info",
            "eprintln" => "warn",
            "dbg" => "debug",
            _ => return,
        };

        let call_site = expn_data.call_site;
        if !self.seen_callsites.insert(call_site) {
            return;
        }

        // No `fn main()` exemption -- keeps the lint visible for CLI tools
        // unless explicitly `#[allow]`ed.
        if is_in_test_zone(cx, expr) {
            return;
        }

        let framework = self.framework.as_str();

        span_lint_and_help(
            cx,
            DEBUG_REMNANTS,
            call_site,
            format!("debug remnant: `{macro_name}!` in committed code"),
            None,
            format!("replace with `{framework}::{level}!(...)` for structured logging"),
        );
    }
}

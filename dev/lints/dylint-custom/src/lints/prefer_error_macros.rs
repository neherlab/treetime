//! Flags raw `eyre!` / `bail!` invocations, steering them to the project's
//! error helper macros.
//!
//! A project that wraps `eyre::eyre!` in its own helper macros, for consistent
//! formatting and internal-error annotations, lists them in the
//! `[prefer_error_macros]` config: `helper_macros` names every wrapper,
//! `report_macro` the replacement for `eyre!` (same arguments, returns a
//! `Report`), and `error_macro` the replacement for `bail!` (returns an
//! `Err`). Writing `eyre!` or `bail!` directly bypasses that layer.
//!
//! The helper macros themselves expand to `eyre::eyre!`, and `bail!` expands to
//! `Err(eyre!(...))`, so this lint walks the macro expansion backtrace: an
//! `eyre!` produced by a helper macro is left alone, and a `bail!` is reported
//! as `bail!` rather than as its inner `eyre!`.

use clippy_utils::diagnostics::{span_lint_and_help, span_lint_and_sugg};
use clippy_utils::is_in_test;
use clippy_utils::macros::expn_backtrace;
use clippy_utils::source::snippet_with_applicability;
use rustc_data_structures::fx::FxHashSet;
use rustc_errors::Applicability;
use rustc_hir::Expr;
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::{ExpnKind, Span, Symbol};

use crate::config::PreferErrorMacrosConfig;

rustc_session::declare_lint! {
    /// Flags raw `eyre!` / `bail!` and recommends the project error macros
    /// configured in `[prefer_error_macros]`.
    pub PREFER_ERROR_MACROS,
    Warn,
    "raw `eyre!` / `bail!` -- use the project error helper macros"
}

pub struct PreferErrorMacros {
    config: PreferErrorMacrosConfig,
    /// Dedup: one diagnostic per macro call site, not per expanded HIR node.
    seen_callsites: FxHashSet<Span>,
    sym_eyre: Symbol,
}

impl PreferErrorMacros {
    pub fn new() -> Self {
        Self {
            config: dylint_linting::config_or_default("prefer_error_macros"),
            seen_callsites: FxHashSet::default(),
            sym_eyre: Symbol::intern("eyre"),
        }
    }
}

rustc_session::impl_lint_pass!(PreferErrorMacros => [PREFER_ERROR_MACROS]);

impl<'tcx> LateLintPass<'tcx> for PreferErrorMacros {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        if !expr.span.from_expansion() {
            return;
        }

        // Walk the expansion chain (innermost to outermost). Skip when a helper
        // macro is present (generated code), otherwise attribute to the
        // outermost `eyre!` / `bail!` that comes from the `eyre` crate -- the
        // macro the user actually wrote.
        let mut candidate: Option<(&'static str, Span)> = None;
        for (_, data) in expn_backtrace(expr.span) {
            let ExpnKind::Macro(_, name) = data.kind else {
                continue;
            };
            let name = name.as_str();
            if self.config.helper_macros.iter().any(|helper| helper == name) {
                return;
            }
            let raw = match name {
                "eyre" => "eyre",
                "bail" => "bail",
                _ => continue,
            };
            let from_eyre = data
                .macro_def_id
                .is_some_and(|did| cx.tcx.crate_name(did.krate) == self.sym_eyre);
            if from_eyre {
                candidate = Some((raw, data.call_site));
            }
        }

        let Some((raw, call_site)) = candidate else {
            return;
        };

        if is_in_test(cx.tcx, expr.hir_id) {
            return;
        }

        if !self.seen_callsites.insert(call_site) {
            return;
        }

        // Build the replacement from the call-site source, swapping only the
        // macro path (everything up to the argument delimiter) and keeping the
        // original arguments verbatim. This handles both `eyre!(...)` and the
        // path-qualified `eyre::eyre!(...)` form.
        //
        // The report macro takes the arguments of `eyre!`, so `eyre!` -> report
        // macro is equivalent and safe to auto-apply. An error macro may route
        // its arguments through `format!`, so `bail!` -> `return` error macro can
        // change behavior for a single non-format-string argument; offer it as
        // display-only.
        let replacement = match raw {
            "bail" => self.config.error_macro.as_deref().map(|error| {
                (
                    format!("return `{error}!(...)` instead of `bail!(...)`"),
                    format!("return {error}!"),
                    Applicability::MaybeIncorrect,
                )
            }),
            _ => self.config.report_macro.as_deref().map(|report| {
                let help = match self.config.error_macro.as_deref() {
                    Some(error) => format!(
                        "use `{report}!(...)` for a `Report`, or `{error}!(...)` to build an `Err(...)`"
                    ),
                    None => format!("use `{report}!(...)` for a `Report`"),
                };
                (help, format!("{report}!"), Applicability::MachineApplicable)
            }),
        };

        let msg = format!("raw `{raw}!` -- use the project error helper macros");
        let Some((help, prefix, mut applicability)) = replacement else {
            span_lint_and_help(cx, PREFER_ERROR_MACROS, call_site, msg, None, "use the project error helper macros");
            return;
        };
        let snippet = snippet_with_applicability(cx, call_site, "", &mut applicability);
        match snippet.find(|c| matches!(c, '(' | '[' | '{')) {
            Some(delim) => span_lint_and_sugg(
                cx,
                PREFER_ERROR_MACROS,
                call_site,
                msg,
                help,
                format!("{prefix}{}", &snippet[delim..]),
                applicability,
            ),
            // No source snippet available: keep the warning, drop the suggestion.
            None => span_lint_and_help(cx, PREFER_ERROR_MACROS, call_site, msg, None, help),
        }
    }
}

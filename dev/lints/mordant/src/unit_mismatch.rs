use rustc_hir::{BinOpKind, Expr, ExprKind};
use rustc_lint::{LateContext, LateLintPass};

use crate::baseline::emit;
use crate::hir_shapes::value_name;

rustc_session::declare_lint! {
    /// Flags addition, subtraction, and comparison between values whose names
    /// say they are in different units: `timeout_ms + deadline_ns` mixes
    /// them, compiles, and is always wrong.
    /// Multiplication and division stay silent, since they are how units
    /// legitimately convert.
    pub UNIT_MISMATCH,
    Warn,
    "arithmetic between values whose names claim different units"
}

rustc_session::declare_lint_pass!(UnitMismatch => [UNIT_MISMATCH]);

/// Unit classes by name suffix. Aliases share a class; a mismatch is two
/// operands from different classes.
fn unit_class(suffix: &str) -> Option<&'static str> {
    Some(match suffix {
        "ns" | "nanos" => "nanoseconds",
        "us" | "micros" => "microseconds",
        "ms" | "millis" => "milliseconds",
        "sec" | "secs" | "seconds" => "seconds",
        "mins" | "minutes" => "minutes",
        "bytes" => "bytes",
        "kb" | "kib" => "kilobytes",
        "mb" | "mib" => "megabytes",
        "gb" | "gib" => "gigabytes",
        "tenths" => "tenths",
        _ => return None,
    })
}

/// The unit an expression's name claims: the suffix after the last `_` of
/// its `value_name`.
fn claimed_unit(e: &Expr<'_>) -> Option<(&'static str, String)> {
    let name = value_name(e)?;
    let name = name.name.as_str();
    let (_, suffix) = name.rsplit_once('_')?;
    unit_class(suffix).map(|c| (c, name.to_string()))
}

impl<'tcx> LateLintPass<'tcx> for UnitMismatch {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let ExprKind::Binary(op, lhs, rhs) = expr.kind else {
            return;
        };
        if !matches!(
            op.node,
            BinOpKind::Add
                | BinOpKind::Sub
                | BinOpKind::Lt
                | BinOpKind::Le
                | BinOpKind::Gt
                | BinOpKind::Ge
                | BinOpKind::Eq
                | BinOpKind::Ne
        ) {
            return;
        }
        if expr.span.from_expansion() {
            return;
        }
        let (Some((lc, ln)), Some((rc, rn))) = (claimed_unit(lhs), claimed_unit(rhs)) else {
            return;
        };
        if lc != rc {
            let op = if matches!(op.node, BinOpKind::Add | BinOpKind::Sub) {
                "arithmetic"
            } else {
                "comparison"
            };
            emit(
                cx,
                UNIT_MISMATCH,
                expr.span,
                format!("`{ln}` says {lc} and `{rn}` says {rc}, and this {op} mixes them"),
                "convert one side first. If a name is wrong about its unit, rename it. A `Duration` or a unit newtype would catch the next one",
            );
        }
    }
}

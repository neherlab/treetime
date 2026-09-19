//! Flags typographic characters in string and character literals: curly quotes,
//! em and en dashes, and emoji. Source text should use plain ASCII punctuation.

use rustc_ast::LitKind;
use rustc_hir::{Expr, ExprKind};
use rustc_lint::{LateContext, LateLintPass, LintContext as _};

use clippy_utils::diagnostics::span_lint_and_help;

rustc_session::declare_lint! {
    /// Flags curly quotes, em/en dashes, and emoji inside string and character
    /// literals. These characters are easy to introduce by copy-paste and hard
    /// to tell apart from their ASCII counterparts when reading source.
    pub TYPOGRAPHIC_CHARACTERS,
    Warn,
    "typographic character in a literal -- use plain ASCII punctuation"
}

pub struct Typographic;

impl Typographic {
    pub const fn new() -> Self {
        Self
    }
}

rustc_session::impl_lint_pass!(Typographic => [TYPOGRAPHIC_CHARACTERS]);

impl<'tcx> LateLintPass<'tcx> for Typographic {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        if expr.span.from_expansion() {
            return;
        }
        let ExprKind::Lit(lit) = &expr.kind else {
            return;
        };
        if !matches!(lit.node, LitKind::Str(..) | LitKind::Char(_)) {
            return;
        }
        let Ok(snippet) = cx.sess().source_map().span_to_snippet(expr.span) else {
            return;
        };
        if let Some(category) = snippet.chars().find_map(banned_category) {
            span_lint_and_help(
                cx,
                TYPOGRAPHIC_CHARACTERS,
                expr.span,
                format!("literal contains {category}"),
                None,
                "replace it with plain ASCII punctuation (straight quotes, `-`, or a word)",
            );
        }
    }
}

/// Returns the category name if `c` is a banned typographic character.
fn banned_category(c: char) -> Option<&'static str> {
    match c {
        '\u{2018}' | '\u{2019}' | '\u{201C}' | '\u{201D}' => Some("a curly quote"),
        '\u{2013}' | '\u{2014}' | '\u{2015}' => Some("an em or en dash"),
        _ if is_emoji(c) => Some("an emoji"),
        _ => None,
    }
}

/// Returns `true` if `c` falls in a common emoji or pictograph code-point range.
fn is_emoji(c: char) -> bool {
    matches!(u32::from(c),
        0x1F000..=0x1FAFF | 0x2600..=0x27BF | 0x2B00..=0x2BFF | 0xFE00..=0xFE0F | 0x200D)
}

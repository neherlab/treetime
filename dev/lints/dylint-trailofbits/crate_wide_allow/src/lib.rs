#![feature(rustc_private)]
#![warn(unused_extern_crates)]

extern crate rustc_ast;
extern crate rustc_span;

use clippy_utils::diagnostics::span_lint_and_help;
use rustc_ast::{AttrStyle, Crate, MetaItem, MetaItemKind};
use rustc_lint::{EarlyContext, EarlyLintPass};
use rustc_span::sym;

dylint_linting::declare_early_lint! {
    /// ### What it does
    ///
    /// Checks for use of `#![allow(...)]` at the crate level.
    ///
    /// ### Why is this bad?
    ///
    /// Such uses cannot be overridden with `--warn` or `--deny` from the command line. They _can_
    /// be overridden with `--force-warn` or `--forbid`, but one must know the `#![allow(...)]`
    /// are present to use these unconventional options.
    ///
    /// ### Example
    ///
    /// ```rust
    /// #![allow(clippy::assertions_on_constants)] // in code
    /// ```
    ///
    /// Use instead:
    ///
    /// ```rust
    /// // Allow `clippy::assertions-on-constants` in Cargo.toml. See:
    /// // - https://doc.rust-lang.org/cargo/reference/manifest.html#the-lints-section
    /// // - https://doc.rust-lang.org/clippy/configuration.html#lints-section-in-cargotoml
    /// ```
    pub CRATE_WIDE_ALLOW,
    Warn,
    "use of `#![allow(...)]` at the crate level"
}

impl EarlyLintPass for CrateWideAllow {
    fn check_crate(&mut self, cx: &EarlyContext, krate: &Crate) {
        for attr in &krate.attrs {
            assert_eq!(AttrStyle::Inner, attr.style);
            if attr.has_name(sym::allow)
                && let Some([arg]) = attr.meta_item_list().as_deref()
                && let Some(MetaItem {
                    path,
                    kind: MetaItemKind::Word,
                    ..
                }) = arg.meta_item()
            {
                let path = path
                    .segments
                    .iter()
                    .map(|segment| segment.ident.as_str())
                    .collect::<Vec<_>>()
                    .join("::")
                    .replace('_', "-");
                span_lint_and_help(
                    cx,
                    CRATE_WIDE_ALLOW,
                    attr.span,
                    format!("silently overrides `--warn {path}` and `--deny {path}`"),
                    None,
                    format!("allow `{path}` in Cargo.toml"),
                );
            }
        }
    }
}

//! Flags source files longer than a configurable line budget, complementing the
//! clippy `too_many_lines` per-function cap with a per-file cap.

use std::collections::HashSet;

use rustc_hir::Item;
use rustc_lint::{LateContext, LateLintPass, LintContext as _};
use rustc_span::FileName;

use clippy_utils::diagnostics::span_lint_and_help;

use crate::config::FileLengthConfig;

rustc_session::declare_lint! {
    /// Flags a source file whose line count exceeds the configured budget. A very
    /// long file usually mixes several responsibilities; split it by concern.
    pub FILE_TOO_LONG,
    Warn,
    "source file exceeds the configured line budget -- split it by responsibility"
}

pub struct FileLength {
    threshold: usize,
    reported: HashSet<String>,
}

impl FileLength {
    pub fn new() -> Self {
        let config: FileLengthConfig = dylint_linting::config_or_default("file_too_long");
        Self {
            threshold: config.threshold,
            reported: HashSet::new(),
        }
    }
}

rustc_session::impl_lint_pass!(FileLength => [FILE_TOO_LONG]);

impl<'tcx> LateLintPass<'tcx> for FileLength {
    fn check_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx Item<'tcx>) {
        if item.span.from_expansion() {
            return;
        }
        let sm = cx.sess().source_map();
        let file = sm.lookup_source_file(item.span.lo());
        let FileName::Real(real) = &file.name else {
            return;
        };
        let Some(path) = real.local_path() else {
            return;
        };
        // Only workspace sources, not dependencies pulled from a registry cache.
        let display = path.to_string_lossy().into_owned();
        if display.contains("/.cargo/") || display.contains("/registry/") {
            return;
        }
        let lines = file.count_lines();
        if lines <= self.threshold || !self.reported.insert(display.clone()) {
            return;
        }
        span_lint_and_help(
            cx,
            FILE_TOO_LONG,
            item.span,
            format!("`{display}` is {lines} lines, over the {} line budget", self.threshold),
            None,
            "split the file by responsibility so each module has one reason to change",
        );
    }
}

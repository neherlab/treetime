#![feature(rustc_private)]
#![cfg_attr(
    not(dylint_lib = "treetime_lints"),
    allow(unknown_lints, reason = "enable dylint lint annotations")
)]

extern crate rustc_ast;
extern crate rustc_data_structures;
extern crate rustc_hir;
extern crate rustc_lint;
extern crate rustc_middle;
extern crate rustc_session;
extern crate rustc_span;

mod config;
mod lints;

use rustc_lint::LintStore;
use rustc_session::Session;

dylint_linting::dylint_library!();

#[allow(
    clippy::no_mangle_with_rust_abi,
    reason = "dylint requires extern fn signature"
)]
#[allow(
    unsafe_code,
    reason = "dylint requires #[no_mangle] for plugin registration"
)]
#[unsafe(no_mangle)]
pub fn register_lints(sess: &Session, lint_store: &mut LintStore) {
    dylint_linting::init_config(sess);
    lint_store.register_lints(&[
        lints::debug_remnants::DEBUG_REMNANTS,
        lints::unclear_exports::UNCLEAR_EXPORTS,
        lints::result_result::RESULT_RESULT,
        lints::fallible_new::FALLIBLE_NEW,
        lints::panic_in_drop::PANIC_IN_DROP,
        lints::topological_ordering::TOPOLOGICAL_ORDERING,
        lints::proper_error_type::PROPER_ERROR_TYPE,
        lints::prefer_error_macros::PREFER_ERROR_MACROS,
        lints::suggest_builder::SUGGEST_BUILDER,
        lints::needless_builder::NEEDLESS_BUILDER,
        lints::bon_builder_collector::BON_BUILDER_COLLECTOR,
    ]);
    lint_store.register_pre_expansion_pass(|| {
        Box::new(lints::bon_builder_collector::BonBuilderCollector)
    });
    lint_store.register_late_pass(|_| Box::new(lints::debug_remnants::DebugRemnants::new()));
    lint_store.register_late_pass(|_| Box::new(lints::unclear_exports::UnclearExports::new()));
    lint_store.register_late_pass(|_| Box::new(lints::result_result::ResultResult::new()));
    lint_store.register_late_pass(|_| Box::new(lints::fallible_new::FallibleNew::new()));
    lint_store.register_late_pass(|_| Box::new(lints::panic_in_drop::PanicInDrop::new()));
    lint_store
        .register_late_pass(|_| Box::new(lints::topological_ordering::TopologicalOrdering::new()));
    lint_store
        .register_late_pass(|_| Box::new(lints::proper_error_type::ProperErrorType::default()));
    lint_store
        .register_late_pass(|_| Box::new(lints::prefer_error_macros::PreferErrorMacros::new()));
    lint_store.register_late_pass(|_| Box::new(lints::suggest_builder::SuggestBuilder::new()));
    lint_store.register_late_pass(|_| Box::new(lints::needless_builder::NeedlessBuilder::new()));
}

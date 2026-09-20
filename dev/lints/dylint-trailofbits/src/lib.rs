#![feature(rustc_private)]
#![warn(unused_extern_crates)]

dylint_linting::dylint_library!();

extern crate rustc_lint;
extern crate rustc_session;

#[expect(clippy::no_mangle_with_rust_abi)]
#[unsafe(no_mangle)]
pub fn register_lints(sess: &rustc_session::Session, lint_store: &mut rustc_lint::LintStore) {
    assert_eq_arg_misordering::register_lints(sess, lint_store);
    await_holding_span_guard::register_lints(sess, lint_store);
    crate_wide_allow::register_lints(sess, lint_store);
    env_literal::register_lints(sess, lint_store);
    try_io_result::register_lints(sess, lint_store);
    unnamed_constant::register_lints(sess, lint_store);
    wrong_serialize_struct_arg::register_lints(sess, lint_store);
}

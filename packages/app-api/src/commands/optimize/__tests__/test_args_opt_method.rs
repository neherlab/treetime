#[cfg(test)]
mod tests {
  use crate::commands::optimize::args::BranchOptMethodCli;
  use clap::Parser;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::optimize::params::BranchOptMethod;

  /// Minimal parser exercising the `clap::ValueEnum` derive on the `BranchOptMethodCli` mirror. The
  /// derive lives on the adapter enum; pinning the flag mapping here pins the `--opt-method` contract.
  #[derive(Parser)]
  struct OptMethodArgs {
    #[arg(long, value_enum, default_value_t = BranchOptMethodCli::default())]
    opt_method: BranchOptMethodCli,
  }

  /// Each kebab-case label on `--opt-method` parses to the corresponding `BranchOptMethodCli` variant.
  /// Pins the `#[serde(rename_all = "kebab-case")]`-derived `ValueEnum` mapping so a future variant
  /// rename or deletion fails this test.
  #[rustfmt::skip]
  #[rstest]
  #[case::brent(      "brent",       BranchOptMethodCli::Brent)]
  #[case::brent_sqrt( "brent-sqrt",  BranchOptMethodCli::BrentSqrt)]
  #[case::brent_log(  "brent-log",   BranchOptMethodCli::BrentLog)]
  #[case::newton(     "newton",      BranchOptMethodCli::Newton)]
  #[case::newton_sqrt("newton-sqrt", BranchOptMethodCli::NewtonSqrt)]
  #[case::newton_log( "newton-log",  BranchOptMethodCli::NewtonLog)]
  #[trace]
  fn test_args_opt_method_kebab_case_parses(#[case] flag: &str, #[case] expected: BranchOptMethodCli) {
    let args = OptMethodArgs::try_parse_from(["treetime", &format!("--opt-method={flag}")]).unwrap();
    assert_eq!(expected, args.opt_method);
  }

  /// Omitting `--opt-method` selects `BrentSqrt`. This is the v0-matching default carried by the
  /// enum's `#[default]` annotation; a regression that moved it would silently route runs to the
  /// wrong optimizer.
  #[test]
  fn test_args_opt_method_default_is_brent_sqrt() {
    let args = OptMethodArgs::try_parse_from(["treetime"]).unwrap();
    assert_eq!(BranchOptMethodCli::BrentSqrt, args.opt_method);
  }

  /// An unknown `--opt-method` value is rejected at parse time. Pins the `value_enum` constraint and
  /// prevents a typo from silently routing to a fallback variant.
  #[test]
  fn test_args_opt_method_rejects_unknown() {
    let result = OptMethodArgs::try_parse_from(["treetime", "--opt-method=brent-foo"]);
    assert!(result.is_err(), "expected parse error for unknown --opt-method value");
  }

  /// Each mirror variant converts to the matching core `BranchOptMethod`. Pins the adapter-to-core
  /// seam that replaced the core `ValueEnum` derive.
  #[rustfmt::skip]
  #[rstest]
  #[case::brent(      BranchOptMethodCli::Brent,      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethodCli::BrentSqrt,  BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethodCli::BrentLog,   BranchOptMethod::BrentLog)]
  #[case::newton(     BranchOptMethodCli::Newton,     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethodCli::NewtonSqrt, BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethodCli::NewtonLog,  BranchOptMethod::NewtonLog)]
  #[trace]
  fn test_args_opt_method_converts_to_core(#[case] cli: BranchOptMethodCli, #[case] expected: BranchOptMethod) {
    assert_eq!(expected, BranchOptMethod::from(cli));
  }
}

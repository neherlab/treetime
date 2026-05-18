#[cfg(all(test, feature = "clap"))]
mod tests {
  use crate::commands::optimize::args::TreetimeOptimizeArgs;
  use crate::optimize::params::BranchOptMethod;
  use clap::Parser;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  /// Each kebab-case label on `--opt-method` parses to the corresponding
  /// `BranchOptMethod` variant. Pins the `clap::ValueEnum` `rename_all = "kebab-case"`
  /// mapping so a future variant rename or deletion fails this test.
  #[rustfmt::skip]
  #[rstest]
  #[case::brent(     "brent",      BranchOptMethod::Brent)]
  #[case::brent_sqrt("brent-sqrt", BranchOptMethod::BrentSqrt)]
  #[case::brent_log( "brent-log",  BranchOptMethod::BrentLog)]
  #[case::newton(    "newton",     BranchOptMethod::Newton)]
  #[case::newton_sqrt("newton-sqrt", BranchOptMethod::NewtonSqrt)]
  #[case::newton_log("newton-log", BranchOptMethod::NewtonLog)]
  #[trace]
  fn test_args_opt_method_kebab_case_parses(#[case] flag: &str, #[case] expected: BranchOptMethod) {
    let args = TreetimeOptimizeArgs::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--outdir=/dev/null",
      &format!("--opt-method={flag}"),
    ])
    .unwrap();
    assert_eq!(expected, args.opt_method);
  }

  /// Omitting `--opt-method` selects `BrentSqrt`. This is the v0-matching
  /// default; a regression that moved the `#[default]` annotation or changed
  /// the `default_value_t` would silently route CLI runs to the wrong optimizer.
  #[test]
  fn test_args_opt_method_default_is_brent_sqrt() {
    let args = TreetimeOptimizeArgs::try_parse_from(["treetime", "--tree=/dev/null", "--outdir=/dev/null"]).unwrap();
    assert_eq!(BranchOptMethod::BrentSqrt, args.opt_method);
  }

  /// An unknown `--opt-method` value is rejected at parse time. Pins the
  /// `value_enum` constraint and prevents a typo from silently routing to a
  /// fallback variant.
  #[test]
  fn test_args_opt_method_rejects_unknown() {
    let result = TreetimeOptimizeArgs::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--outdir=/dev/null",
      "--opt-method=brent-foo",
    ]);
    assert!(result.is_err(), "expected parse error for unknown --opt-method value");
  }
}

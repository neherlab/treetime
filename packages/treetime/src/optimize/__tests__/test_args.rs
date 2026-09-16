#[cfg(test)]
mod tests {
  use crate::optimize::params::BranchOptMethod;
  use clap::Parser;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  /// Minimal parser exercising the `clap::ValueEnum` derive on `BranchOptMethod`. The derive lives on
  /// the core enum behind the `clap` feature; the CLI adapter reuses it verbatim, so pinning the
  /// mapping here pins it for every consumer.
  #[derive(Parser)]
  struct OptMethodArgs {
    #[arg(long, value_enum, default_value_t = BranchOptMethod::default())]
    opt_method: BranchOptMethod,
  }

  /// Each kebab-case label on `--opt-method` parses to the corresponding `BranchOptMethod` variant.
  /// Pins the `#[serde(rename_all = "kebab-case")]`-derived `ValueEnum` mapping so a future variant
  /// rename or deletion fails this test.
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
    let args = OptMethodArgs::try_parse_from(["treetime", &format!("--opt-method={flag}")]).unwrap();
    assert_eq!(expected, args.opt_method);
  }

  /// Omitting `--opt-method` selects `BrentSqrt`. This is the v0-matching default carried by the
  /// enum's `#[default]` annotation; a regression that moved it would silently route runs to the
  /// wrong optimizer.
  #[test]
  fn test_args_opt_method_default_is_brent_sqrt() {
    let args = OptMethodArgs::try_parse_from(["treetime"]).unwrap();
    assert_eq!(BranchOptMethod::BrentSqrt, args.opt_method);
  }

  /// An unknown `--opt-method` value is rejected at parse time. Pins the `value_enum` constraint and
  /// prevents a typo from silently routing to a fallback variant.
  #[test]
  fn test_args_opt_method_rejects_unknown() {
    let result = OptMethodArgs::try_parse_from(["treetime", "--opt-method=brent-foo"]);
    assert!(result.is_err(), "expected parse error for unknown --opt-method value");
  }
}

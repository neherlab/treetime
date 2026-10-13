#[cfg(test)]
mod tests {
  use crate::commands::optimize::args::BranchOptMethodCli;
  use clap::Parser;
  use clap::error::ErrorKind;
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

  #[rustfmt::skip]
  #[rstest]
  #[case::brent_sqrt( "brent-sqrt",  BranchOptMethodCli::BrentSqrt)]
  #[case::newton_sqrt("newton-sqrt", BranchOptMethodCli::NewtonSqrt)]
  #[trace]
  fn test_args_opt_method_kebab_case_parses(#[case] flag: &str, #[case] expected: BranchOptMethodCli) {
    let args = OptMethodArgs::try_parse_from(["treetime", &format!("--opt-method={flag}")]).unwrap();
    assert_eq!(expected, args.opt_method);
  }

  #[test]
  fn test_args_opt_method_default_is_brent_sqrt() {
    let args = OptMethodArgs::try_parse_from(["treetime"]).unwrap();
    assert_eq!(BranchOptMethodCli::BrentSqrt, args.opt_method);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::unknown(   "brent-foo")]
  #[case::brent(     "brent")]
  #[case::brent_log( "brent-log")]
  #[case::newton(    "newton")]
  #[case::newton_log("newton-log")]
  #[trace]
  fn test_args_opt_method_rejects_unsupported(#[case] flag: &str) {
    let Err(error) = OptMethodArgs::try_parse_from(["treetime", &format!("--opt-method={flag}")]) else {
      panic!("expected a parse error for --opt-method={flag}");
    };
    assert_eq!(ErrorKind::InvalidValue, error.kind());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::brent_sqrt( BranchOptMethodCli::BrentSqrt,  BranchOptMethod::BrentSqrt)]
  #[case::newton_sqrt(BranchOptMethodCli::NewtonSqrt, BranchOptMethod::NewtonSqrt)]
  #[trace]
  fn test_args_opt_method_converts_to_core(#[case] cli: BranchOptMethodCli, #[case] expected: BranchOptMethod) {
    assert_eq!(expected, BranchOptMethod::from(cli));
  }
}

#[cfg(test)]
mod tests {
  use crate::commands::optimize::args::InitialGuessModeCli;
  use clap::Parser;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::optimize::params::InitialGuessMode;

  /// Minimal parser exercising the `clap::ValueEnum` derive on the `InitialGuessModeCli` mirror.
  /// Pinning the flag mapping here pins the `--initial-branch-lengths` contract.
  #[derive(Parser)]
  struct InitialGuessArgs {
    #[arg(long, value_enum, default_value_t = InitialGuessModeCli::default())]
    initial_branch_lengths: InitialGuessModeCli,
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::auto(  "auto",   InitialGuessModeCli::Auto)]
  #[case::always("always", InitialGuessModeCli::Always)]
  #[case::never( "never",  InitialGuessModeCli::Never)]
  #[trace]
  fn test_args_initial_guess_kebab_case_parses(#[case] flag: &str, #[case] expected: InitialGuessModeCli) {
    let args =
      InitialGuessArgs::try_parse_from(["treetime", &format!("--initial-branch-lengths={flag}")]).unwrap();
    assert_eq!(expected, args.initial_branch_lengths);
  }

  #[test]
  fn test_args_initial_guess_default_is_auto() {
    let args = InitialGuessArgs::try_parse_from(["treetime"]).unwrap();
    assert_eq!(InitialGuessModeCli::Auto, args.initial_branch_lengths);
  }

  #[test]
  fn test_args_initial_guess_rejects_unknown() {
    let result = InitialGuessArgs::try_parse_from(["treetime", "--initial-branch-lengths=sometimes"]);
    assert!(
      result.is_err(),
      "expected parse error for unknown --initial-branch-lengths value"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::auto(  InitialGuessModeCli::Auto,   InitialGuessMode::Auto)]
  #[case::always(InitialGuessModeCli::Always, InitialGuessMode::Always)]
  #[case::never( InitialGuessModeCli::Never,  InitialGuessMode::Never)]
  #[trace]
  fn test_args_initial_guess_converts_to_core(#[case] cli: InitialGuessModeCli, #[case] expected: InitialGuessMode) {
    assert_eq!(expected, InitialGuessMode::from(cli));
  }
}

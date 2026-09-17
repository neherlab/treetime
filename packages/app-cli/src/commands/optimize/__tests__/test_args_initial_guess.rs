#[cfg(test)]
mod tests {
  use crate::commands::optimize::args::InitialGuessModeCli;
  use clap::Parser;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::optimize::params::InitialGuessMode;

  /// Minimal parser exercising the `clap::ValueEnum` derive on the `InitialGuessModeCli` mirror.
  /// Pinning the flag mapping here pins the `--branch-length-initial-guess` contract.
  #[derive(Parser)]
  struct InitialGuessArgs {
    #[arg(long, value_enum, default_value_t = InitialGuessModeCli::default())]
    branch_length_initial_guess: InitialGuessModeCli,
  }

  /// Each kebab-case label parses to the corresponding `InitialGuessModeCli` variant. Pins the
  /// `#[serde(rename_all = "kebab-case")]`-derived `ValueEnum` mapping.
  #[rustfmt::skip]
  #[rstest]
  #[case::auto(  "auto",   InitialGuessModeCli::Auto)]
  #[case::always("always", InitialGuessModeCli::Always)]
  #[case::never( "never",  InitialGuessModeCli::Never)]
  #[trace]
  fn test_args_initial_guess_kebab_case_parses(#[case] flag: &str, #[case] expected: InitialGuessModeCli) {
    let args =
      InitialGuessArgs::try_parse_from(["treetime", &format!("--branch-length-initial-guess={flag}")]).unwrap();
    assert_eq!(expected, args.branch_length_initial_guess);
  }

  /// Omitting the flag selects `Auto`, the mode that estimates only edges with missing or invalid
  /// branch lengths and preserves valid input values.
  #[test]
  fn test_args_initial_guess_default_is_auto() {
    let args = InitialGuessArgs::try_parse_from(["treetime"]).unwrap();
    assert_eq!(InitialGuessModeCli::Auto, args.branch_length_initial_guess);
  }

  /// An unknown value is rejected at parse time. Pins the `value_enum` constraint.
  #[test]
  fn test_args_initial_guess_rejects_unknown() {
    let result = InitialGuessArgs::try_parse_from(["treetime", "--branch-length-initial-guess=sometimes"]);
    assert!(
      result.is_err(),
      "expected parse error for unknown --branch-length-initial-guess value"
    );
  }

  /// Each mirror variant converts to the matching core `InitialGuessMode`. Pins the adapter-to-core
  /// seam that replaced the core `ValueEnum` derive.
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

use clap::{Args, ColorChoice, ValueEnum};
use deser::ser::Chunk;
use deser::{Atom, Error, Serialize, State};
use std::ffi::OsString;
use treetime_utils::init::global::stderr_has_colors;

pub(crate) fn color_choice_of_args(argv: &[OsString]) -> ColorChoice {
  let mut choice = ColorChoice::Auto;
  let mut args = argv.iter().skip(1).take_while(|arg| *arg != "--");
  while let Some(arg) = args.next() {
    let value = match arg.to_str() {
      Some("--color") => args.next().and_then(|value| value.to_str()),
      Some(arg) => arg.strip_prefix("--color="),
      None => None,
    };
    if let Some(given) = value.and_then(parse_color_choice) {
      choice = given;
    }
  }
  choice
}

#[derive(Args, Debug, Clone)]
pub(crate) struct ColorArgs {
  /// When to color error messages, help, and the log: `auto` colors them when the output is a terminal and the
  /// `NO_COLOR` environment variable is not set
  #[clap(
    long,
    global = true,
    value_enum,
    value_name = "WHEN",
    default_value_t = ColorChoice::Auto,
    hide_short_help = true,
    display_order = 100
  )]
  pub color: ColorChoice,
}

impl ColorArgs {
  pub(crate) fn stderr_has_colors(&self) -> bool {
    match self.color {
      ColorChoice::Always => true,
      ColorChoice::Never => false,
      ColorChoice::Auto => stderr_has_colors(),
    }
  }
}

impl Serialize for ColorArgs {
  fn serialize(&self, _state: &mut State) -> Result<Chunk<'_>, Error> {
    Ok(Chunk::Atom(Atom::Str(self.color.to_string().into())))
  }
}

fn parse_color_choice(value: &str) -> Option<ColorChoice> {
  ColorChoice::value_variants().iter().copied().find(|choice| {
    choice
      .to_possible_value()
      .is_some_and(|possible| possible.matches(value, false))
  })
}

#[cfg(test)]
mod tests {
  use crate::cli::color::color_choice_of_args;
  use clap::ColorChoice;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::ffi::OsString;

  #[rustfmt::skip]
  #[rstest]
  #[case::absent(           &["treetime", "clock"],                                    ColorChoice::Auto)]
  #[case::separate_value(   &["treetime", "--color", "never", "clock"],                ColorChoice::Never)]
  #[case::joined_value(     &["treetime", "clock", "--color=always"],                  ColorChoice::Always)]
  #[case::last_wins(        &["treetime", "--color=never", "clock", "--color=always"], ColorChoice::Always)]
  #[case::invalid_is_auto(  &["treetime", "clock", "--color=sometimes"],               ColorChoice::Auto)]
  #[case::after_separator(  &["treetime", "clock", "--", "--color=never"],             ColorChoice::Auto)]
  #[case::program_name(     &["--color=never"],                                        ColorChoice::Auto)]
  #[trace]
  fn test_color_choice_of_args(#[case] argv: &[&str], #[case] expected: ColorChoice) {
    let argv: Vec<OsString> = argv.iter().map(OsString::from).collect();
    assert_eq!(expected, color_choice_of_args(&argv));
  }
}

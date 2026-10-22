use crate::commands::shared::argument_error::argument_error;
use eyre::Report;

pub fn check_applies(
  setting: &str,
  given: bool,
  other: &str,
  applies: bool,
  rule: impl FnOnce() -> String,
) -> Result<(), Report> {
  if given && !applies {
    return Err(argument_error(
      &[setting, other],
      format!("{} {}", setting_flag(setting), rule()),
    ));
  }
  Ok(())
}

pub fn setting_flag(key: &str) -> String {
  format!("--{}", key.replace('_', "-"))
}

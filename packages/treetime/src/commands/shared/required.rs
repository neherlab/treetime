use crate::make_report;
use eyre::Report;

/// Build the clap-style "required arguments were not provided" error for a set of missing arg ids.
///
/// The flag spelling is read from the command's clap definition, so renaming a flag cannot desync the
/// message from what the parser accepts. The `--config` overlay stores required inputs as `Option`
/// (clap must not reject them at parse time, because the config file may supply them), so presence is
/// enforced here when the raw args are converted to their validated form.
///
/// Without the `clap` feature (non-CLI builds such as the Node bindings) the spelling is reconstructed
/// from the id with clap's default conventions (`--{kebab} <{UPPER}>`), which matches the derived flag
/// for every required argument.
#[cfg(feature = "clap")]
pub fn missing_required_args<C: clap::CommandFactory>(missing_ids: &[&str]) -> Report {
  let command = C::command();
  let list = missing_ids
    .iter()
    .map(|id| required_flag(&command, id))
    .collect::<Vec<_>>()
    .join("\n  ");
  make_report!("the following required arguments were not provided:\n  {list}")
}

/// Format one missing argument as `--{long} <{VALUE}>` from its clap definition.
#[cfg(feature = "clap")]
fn required_flag(command: &clap::Command, id: &str) -> String {
  let arg = command.get_arguments().find(|arg| arg.get_id() == id);
  let long = arg
    .and_then(clap::Arg::get_long)
    .map_or_else(|| id.replace('_', "-"), str::to_owned);
  let value = arg
    .and_then(|arg| arg.get_value_names().and_then(|names| names.first()))
    .map_or_else(|| id.to_uppercase(), ToString::to_string);
  format!("--{long} <{value}>")
}

/// Non-CLI fallback: reconstruct the flag spelling from the id with clap's default conventions.
#[cfg(not(feature = "clap"))]
#[allow(clippy::extra_unused_type_parameters)] // keep one call signature across the `clap` cfg split
pub fn missing_required_args<C>(missing_ids: &[&str]) -> Report {
  let list = missing_ids
    .iter()
    .map(|id| format!("--{} <{}>", id.replace('_', "-"), id.to_uppercase()))
    .collect::<Vec<_>>()
    .join("\n  ");
  make_report!("the following required arguments were not provided:\n  {list}")
}

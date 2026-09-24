use eyre::Report;
use treetime::make_report;

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

#[cfg(not(feature = "clap"))]
#[expect(
  clippy::extra_unused_type_parameters,
  reason = "keeps the signature of the clap variant, which reads argument definitions from the type"
)]
pub fn missing_required_args<C>(missing_ids: &[&str]) -> Report {
  let list = missing_ids
    .iter()
    .map(|id| format!("--{} <{}>", id.replace('_', "-"), id.to_uppercase()))
    .collect::<Vec<_>>()
    .join("\n  ");
  make_report!("the following required arguments were not provided:\n  {list}")
}

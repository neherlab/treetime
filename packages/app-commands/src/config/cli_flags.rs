use crate::config::properties::{CLI_FLAG_KEY, leaf_properties};
use clap::Command;
use eyre::Report;
use schemars::Schema;
use serde_json::Value;
use treetime_utils::{make_error, make_report};

pub fn annotate_cli_flags(schema: &mut Schema, command: &Command) -> Result<(), Report> {
  let value = schema.as_value().clone();
  let leaves = leaf_properties(&value)?;
  let object = schema.ensure_object();
  let mut annotated = Value::Object(object.clone());
  for leaf in &leaves {
    let flag = cli_flag(command, leaf.key())
      .ok_or_else(|| make_report!("config key `{}` has no command-line flag", leaf.key_path.join(".")))?;
    let Some(Value::Object(property)) = annotated.pointer_mut(&leaf.schema_pointer) else {
      return make_error!("schema has no property at `{}`", leaf.schema_pointer);
    };
    property.insert(CLI_FLAG_KEY.to_owned(), Value::String(flag));
  }
  let Value::Object(annotated) = annotated else {
    return make_error!("a command schema must be a JSON object");
  };
  *object = annotated;
  Ok(())
}

pub fn cli_flag(command: &Command, id: &str) -> Option<String> {
  command
    .get_arguments()
    .find(|arg| arg.get_id() == id)
    .and_then(clap::Arg::get_long)
    .map(|long| format!("--{long}"))
}

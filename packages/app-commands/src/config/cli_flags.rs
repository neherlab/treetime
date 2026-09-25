use crate::command::AppCommand;
use crate::config::catalog::setting_arg;
use crate::config::properties::{CLI_FLAG_KEY, leaf_properties};
use clap::Command;
use eyre::Report;
use schemars::Schema;
use serde_json::Value;
use treetime_utils::{make_error, make_report};

pub fn annotated_config_schema(command: AppCommand) -> Result<Schema, Report> {
  let mut schema = command.config_schema();
  annotate_cli_flags(&mut schema, &command.cli_command())?;
  Ok(schema)
}

pub fn annotate_cli_flags(schema: &mut Schema, command: &Command) -> Result<(), Report> {
  let mut command = command.clone();
  command.build();
  let leaves = leaf_properties(schema.as_value())?;
  for leaf in &leaves {
    let arg = setting_arg(&command, leaf)?;
    let long = arg
      .get_long()
      .ok_or_else(|| make_report!("config key `{}` has no command-line flag", leaf.key_path.join(".")))?;
    let Some(Value::Object(property)) = schema.pointer_mut(&leaf.schema_pointer) else {
      return make_error!("schema has no property at `{}`", leaf.schema_pointer);
    };
    property.insert(CLI_FLAG_KEY.to_owned(), Value::String(format!("--{long}")));
  }
  Ok(())
}

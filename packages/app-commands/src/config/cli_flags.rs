use crate::command::AppCommand;
use crate::config::properties::{CLI_FLAG_KEY, SettingProperty, setting_properties};
use clap::{Arg, Command};
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
  for setting in setting_properties(schema.as_value())? {
    let arg = setting_arg(&command, &setting)?;
    let long = arg
      .get_long()
      .ok_or_else(|| make_report!("config key `{}` has no command-line flag", setting.key))?;
    let Some(Value::Object(property)) = schema.pointer_mut(&setting.schema_pointer) else {
      return make_error!("schema has no property at `{}`", setting.schema_pointer);
    };
    property.insert(CLI_FLAG_KEY.to_owned(), Value::String(format!("--{long}")));
  }
  Ok(())
}

pub fn setting_arg<'a>(cli: &'a Command, setting: &SettingProperty) -> Result<&'a Arg, Report> {
  cli
    .get_arguments()
    .find(|arg| arg.get_id() == setting.key.as_str() && arg.get_long().is_some())
    .ok_or_else(|| make_report!("config key `{}` has no command-line flag", setting.key))
}

use app_commands::command_config::CommandConfig;
use app_commands::commands::timetree::branch_length_mode::BranchLengthModeCli;
use app_commands::config::resolve_paths::resolve_config_paths;
use app_commands::yaml::yaml_document;
use eyre::{Report, WrapErr};
use serde_json::Value;
use std::env;
use std::fmt::Write as _;
use std::io::{self, Write};
use std::path::Path;

const AUTO_MODE_NOTE: &str = "# branch_length_mode: auto resolves to marginal or input from the input tree when the run \
                              starts; the settings that need re-estimated branch lengths are checked then";

pub(crate) fn print_dry_run(config: &CommandConfig) -> Result<(), Report> {
  let workdir = env::current_dir().wrap_err("When resolving the current directory")?;
  let text = dry_run_text(config, &workdir)?;
  io::stdout()
    .lock()
    .write_all(text.as_bytes())
    .wrap_err("When writing the dry run to standard output")
}

pub(crate) fn dry_run_text(config: &CommandConfig, workdir: &Path) -> Result<String, Report> {
  let command = config.command();
  let mut settings = Value::Object(config.settings()?);
  resolve_config_paths(&mut settings, command.config_schema().as_value(), workdir)?;
  let config = command.command_config(&settings)?;
  let args = config.args()?;
  let outputs = args.resolve_outputs()?;
  let planned = args.planned_paths(&outputs)?;

  let mut text = yaml_document(&settings)?;
  if is_auto_branch_length_mode(&config) {
    writeln!(text, "{AUTO_MODE_NOTE}")?;
  }
  writeln!(text, "# Planned output files:")?;
  let conditions = outputs.conditional_selections();
  for (selection, paths) in planned {
    for path in paths {
      let path = one_line(&path.display().to_string());
      match conditions.get(&selection) {
        Some(condition) => writeln!(text, "# {path} ({})", condition.description())?,
        None => writeln!(text, "# {path}")?,
      }
    }
  }
  Ok(text)
}

fn is_auto_branch_length_mode(config: &CommandConfig) -> bool {
  matches!(config, CommandConfig::Timetree { config } if config.branch_length_mode == BranchLengthModeCli::Auto)
}

fn one_line(text: &str) -> String {
  let mut line = String::with_capacity(text.len());
  for character in text.chars() {
    if character.is_control() {
      line.extend(character.escape_default());
    } else {
      line.push(character);
    }
  }
  line
}

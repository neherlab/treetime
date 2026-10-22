use app_commands::command_config::CommandConfig;
use app_commands::commands::shared::run_paths::{RunPath, run_inputs};
use eyre::Report;

pub(crate) fn step_inputs(command: &CommandConfig) -> Result<Vec<RunPath>, Report> {
  run_inputs(command.command(), &command.settings()?)
}

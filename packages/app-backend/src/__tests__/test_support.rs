#[cfg(test)]
pub(crate) mod tests {
  use crate::command::{AppCommandExt, PreparedCommand};
  use crate::job::JobId;
  use crate::runs::manager::{ConfigHook, RunManager, UploadedInput};
  use app_commands::command::AppCommand;
  use app_commands::json_value::SparseConfig;
  use deser::adapters::As;
  use deser_serde::Serde;
  use eyre::Report;
  use serde_json::Value;
  use std::path::PathBuf;
  use std::sync::Arc;
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  pub(crate) fn sparse(value: Value) -> SparseConfig {
    let Value::Object(settings) = value else {
      panic!("a test config is a mapping of settings");
    };
    SparseConfig(settings)
  }

  pub(crate) fn prepare_value(command: AppCommand, config: &Value) -> Result<PreparedCommand, Report> {
    let text = json_write_str(&As::<_, Serde>::new(config), JsonPretty(true))?;
    command.prepare_text("config.json", &text)
  }

  pub(crate) fn unconfined() -> ConfigHook {
    Box::new(|_config: &mut Value| Ok(()))
  }

  pub(crate) fn upload_bytes(
    runs: &Arc<RunManager>,
    id: &JobId,
    name: &str,
    replaces: &[PathBuf],
    bytes: &[u8],
    limit: usize,
  ) -> Result<UploadedInput, Report> {
    let mut upload = runs.begin_upload(id, name, replaces, limit)?;
    upload.write(bytes)?;
    upload.finish()
  }
}

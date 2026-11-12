use crate::app_settings::store::AppSettingsStore;
use app_commands::command::AppCommand;
use app_commands::config::properties::setting_properties;
use eyre::Report;
use serde_json::{Map, Value};
use std::sync::Arc;

const MAX_GRID_POINTS_KEY: &str = "max_grid_points";

pub enum RunLimits {
  Unset,
  Settings(Arc<AppSettingsStore>),
}

impl RunLimits {
  pub fn apply(&self, command: AppCommand, settings: &mut Map<String, Value>) -> Result<(), Report> {
    let Self::Settings(store) = self else {
      return Ok(());
    };
    if settings.contains_key(MAX_GRID_POINTS_KEY) || !has_max_grid_points(command)? {
      return Ok(());
    }
    if let Some(limit) = store.read()?.analysis.max_grid_points {
      settings.insert(MAX_GRID_POINTS_KEY.to_owned(), Value::from(limit.get()));
    }
    Ok(())
  }
}

fn has_max_grid_points(command: AppCommand) -> Result<bool, Report> {
  Ok(
    setting_properties(command.config_schema().as_value())?
      .iter()
      .any(|property| property.key == MAX_GRID_POINTS_KEY),
  )
}

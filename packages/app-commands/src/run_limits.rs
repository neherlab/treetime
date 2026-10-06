use crate::app_settings::store::AppSettingsStore;
use crate::command::AppCommand;
use crate::config::properties::leaf_properties;
use eyre::Report;
use serde_json::{Map, Value};
use std::sync::Arc;
use treetime_grid::MaxGridPoints;
use treetime_utils::make_error;

pub const MAX_GRID_POINTS_KEY: &str = "max_grid_points";

pub enum RunLimits {
  Unset,
  Settings(Arc<AppSettingsStore>),
  Server(MaxGridPoints),
}

impl RunLimits {
  pub fn apply(&self, command: AppCommand, settings: &mut Map<String, Value>) -> Result<(), Report> {
    if matches!(self, Self::Unset) || !has_max_grid_points(command)? {
      return Ok(());
    }
    match (self, settings.get(MAX_GRID_POINTS_KEY)) {
      (Self::Server(limit), Some(value)) => check_server_limit(*limit, value),
      (Self::Server(limit), None) => {
        settings.insert(MAX_GRID_POINTS_KEY.to_owned(), Value::from(limit.get()));
        Ok(())
      },
      (Self::Settings(store), None) => {
        if let Some(limit) = store.read()?.analysis.max_grid_points {
          settings.insert(MAX_GRID_POINTS_KEY.to_owned(), Value::from(limit.get()));
        }
        Ok(())
      },
      (Self::Settings(_) | Self::Unset, Some(_)) | (Self::Unset, None) => Ok(()),
    }
  }
}

fn has_max_grid_points(command: AppCommand) -> Result<bool, Report> {
  Ok(
    leaf_properties(command.config_schema().as_value())?
      .iter()
      .any(|leaf| leaf.key_path.len() == 1 && leaf.key() == MAX_GRID_POINTS_KEY),
  )
}

fn check_server_limit(limit: MaxGridPoints, value: &Value) -> Result<(), Report> {
  let Some(requested) = value.as_u64() else {
    return Ok(());
  };
  if usize::try_from(requested).map_or(true, |requested| requested > limit.get()) {
    return make_error!(
      "`{MAX_GRID_POINTS_KEY}` is {requested}, more than the limit of this server, {limit} \
       (treetime-server --max-grid-points); set {limit} or less, or leave the setting unset"
    );
  }
  Ok(())
}

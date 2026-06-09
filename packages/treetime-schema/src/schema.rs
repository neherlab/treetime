use crate::{ErrorResponse, ProgressEvent, VersionInfo};
use eyre::Report;
use log::info;
use schemars::JsonSchema;
use schemars::generate::SchemaSettings;
use std::io::Write;
use std::path::{Path, PathBuf};
use strum::IntoEnumIterator;
use strum_macros::EnumIter;
use treetime_utils::io::json::{JsonPretty, json_write_str};

#[cfg(feature = "clap")]
use clap::ValueEnum;

#[derive(Debug, Clone, Default, EnumIter, serde::Serialize)]
#[cfg_attr(feature = "clap", derive(ValueEnum))]
pub enum TreetimeSchemaFormat {
  #[default]
  All,
  VersionInfo,
  ProgressEvent,
  ErrorResponse,
}

impl TreetimeSchemaFormat {
  const fn default_filename(&self) -> Option<&'static str> {
    match self {
      TreetimeSchemaFormat::All => None,
      TreetimeSchemaFormat::VersionInfo => Some("version-info.schema.json"),
      TreetimeSchemaFormat::ProgressEvent => Some("progress-event.schema.json"),
      TreetimeSchemaFormat::ErrorResponse => Some("error-response.schema.json"),
    }
  }
}

fn generate_schema_for<T: JsonSchema>(output: &Path) -> Result<(), Report> {
  let settings = SchemaSettings::draft07();
  let schema = settings.into_generator().into_root_schema_for::<T>();
  let json = json_write_str(&schema, JsonPretty(true))?;

  if output == Path::new("-") {
    std::io::stdout().write_all(json.as_bytes())?;
  } else {
    if let Some(parent) = output.parent() {
      std::fs::create_dir_all(parent)?;
    }
    std::fs::write(output, json)?;
  }

  info!("Wrote JSON schema to '{}'", output.display());
  Ok(())
}

pub fn generate_schema(format: &TreetimeSchemaFormat, output: Option<&PathBuf>) -> Result<(), Report> {
  match format {
    TreetimeSchemaFormat::All => {
      let output_dir = output.map_or_else(|| PathBuf::from("."), |p| p.clone());
      for fmt in TreetimeSchemaFormat::iter() {
        if let Some(filename) = fmt.default_filename() {
          let path = output_dir.join(filename);
          generate_schema(&fmt, Some(&path))?;
        }
      }
    },
    TreetimeSchemaFormat::VersionInfo => {
      let path = output.map_or_else(|| PathBuf::from("-"), |p| p.clone());
      generate_schema_for::<VersionInfo>(&path)?;
    },
    TreetimeSchemaFormat::ProgressEvent => {
      let path = output.map_or_else(|| PathBuf::from("-"), |p| p.clone());
      generate_schema_for::<ProgressEvent>(&path)?;
    },
    TreetimeSchemaFormat::ErrorResponse => {
      let path = output.map_or_else(|| PathBuf::from("-"), |p| p.clone());
      generate_schema_for::<ErrorResponse>(&path)?;
    },
  }
  Ok(())
}

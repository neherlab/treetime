use crate::app_settings::settings::AppSettings;
use crate::atomic_write::write_atomically;
use crate::yaml::{yaml_document, yaml_read_str};
use eyre::{Report, WrapErr};
use parking_lot::Mutex;
use std::fs;
use std::io::{ErrorKind, Write};
use std::path::{Path, PathBuf};
use treetime_utils::io::fs::read_file_to_string_if_exists;
use treetime_utils::io::json::to_json_value;
use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
use treetime_utils::make_error;

pub const SETTINGS_YAML: &str = "settings.yaml";

pub const SETTINGS_JSON: &str = "settings.json";

#[derive(Debug)]
pub struct AppSettingsStore {
  path: PathBuf,
  format: SettingsFormat,
  lock: Mutex<()>,
}

impl AppSettingsStore {
  pub fn open(dir: &Path) -> Result<Self, Report> {
    let yaml = dir.join(SETTINGS_YAML);
    let json = dir.join(SETTINGS_JSON);
    let (path, format) = match (is_file(&yaml)?, is_file(&json)?) {
      (true, true) => {
        return make_error!(
          "both '{}' and '{}' exist; keep one of them",
          yaml.display(),
          json.display()
        );
      },
      (false, true) => (json, SettingsFormat::Json),
      (_, false) => (yaml, SettingsFormat::Yaml),
    };
    Ok(Self {
      path,
      format,
      lock: Mutex::new(()),
    })
  }

  pub fn path(&self) -> &Path {
    &self.path
  }

  pub fn read(&self) -> Result<AppSettings, Report> {
    let _guard = self.lock.lock();
    self.read_file()
  }

  pub fn update(&self, change: impl FnOnce(&mut AppSettings)) -> Result<AppSettings, Report> {
    let _guard = self.lock.lock();
    let mut settings = self.read_file()?;
    change(&mut settings);
    self.write_file(&settings)?;
    Ok(settings)
  }

  fn read_file(&self) -> Result<AppSettings, Report> {
    let Some(text) = read_file_to_string_if_exists(&self.path)? else {
      return Ok(AppSettings::default());
    };
    self
      .format
      .parse(&text)
      .wrap_err_with(|| format!("When reading the settings file '{}'", self.path.display()))
  }

  fn write_file(&self, settings: &AppSettings) -> Result<(), Report> {
    let text = self.format.write(settings)?;
    let dir = self.path.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(dir).wrap_err_with(|| format!("When creating the directory '{}'", dir.display()))?;
    write_atomically(&self.path, |file| Ok(file.write_all(text.as_bytes())?))
      .wrap_err_with(|| format!("When writing the settings file '{}'", self.path.display()))
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SettingsFormat {
  Yaml,
  Json,
}

impl SettingsFormat {
  fn parse(self, text: &str) -> Result<AppSettings, Report> {
    if text.trim().is_empty() {
      return Ok(AppSettings::default());
    }
    match self {
      Self::Yaml => yaml_read_str(text),
      Self::Json => json_read_str(text),
    }
  }

  fn write(self, settings: &AppSettings) -> Result<String, Report> {
    match self {
      Self::Yaml => yaml_document(&to_json_value(&settings)?),
      Self::Json => Ok(format!("{}\n", json_write_str(settings, JsonPretty(true))?)),
    }
  }
}

fn is_file(path: &Path) -> Result<bool, Report> {
  match fs::metadata(path) {
    Ok(metadata) => Ok(metadata.is_file()),
    Err(err) if err.kind() == ErrorKind::NotFound => Ok(false),
    Err(err) => Err(Report::new(err)).wrap_err_with(|| format!("When checking '{}'", path.display())),
  }
}

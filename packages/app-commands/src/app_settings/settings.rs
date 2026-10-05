use crate::command::AppCommand;
use crate::job::JobId;
use crate::json_value::SparseConfig;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_with::skip_serializing_none;
use std::collections::BTreeMap;
use std::path::PathBuf;

/// Settings of a local TreeTime installation, kept in `settings.yaml` or `settings.json` in the app folder. Every
/// setting is optional.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AppSettings {
  /// Folders of the app. Each defaults to a folder of the same name in the app folder.
  #[serde(default, skip_serializing_if = "AppPathSettings::is_unset")]
  pub paths: AppPathSettings,

  /// Preferences of the user interface.
  #[serde(default, skip_serializing_if = "UiSettings::is_unset")]
  pub ui: UiSettings,
}

/// Folders of the app. A relative path is relative to the app folder. The environment variables
/// `TREETIME_PROFILE_DIR`, `TREETIME_RUNS_DIR`, `TREETIME_LOGS_DIR`, and `TREETIME_EXAMPLES_DIR` take precedence.
#[skip_serializing_none]
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AppPathSettings {
  /// Browser profile of the desktop app: cache, local storage, and crash reports. Default: `profile`.
  #[serde(default)]
  pub profile: Option<PathBuf>,

  /// Runs, each in its own folder. Default: `runs`.
  #[serde(default)]
  pub runs: Option<PathBuf>,

  /// Logs and crash diagnostics. Default: `logs`.
  #[serde(default)]
  pub logs: Option<PathBuf>,

  /// Example datasets and configurations that the app lists. Default: `examples`.
  #[serde(default)]
  pub examples: Option<PathBuf>,
}

impl AppPathSettings {
  pub fn is_unset(&self) -> bool {
    *self == Self::default()
  }
}

/// Preferences of the user interface. Unset preferences take the defaults of the user interface.
#[skip_serializing_none]
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct UiSettings {
  /// Color theme.
  #[serde(default)]
  pub theme: Option<UiTheme>,

  /// Width of the sidebar in pixels.
  #[serde(default)]
  pub sidebar_width: Option<u32>,

  /// The unfinished analysis form.
  #[serde(default)]
  pub draft: Option<UiDraft>,
}

impl UiSettings {
  pub fn is_unset(&self) -> bool {
    *self == Self::default()
  }
}

/// Color theme of the user interface.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum UiTheme {
  /// Follow the theme of the operating system.
  System,
  Light,
  Dark,
}

/// The unfinished analysis form: the command, its settings, and how the form is shown.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct UiDraft {
  pub command: AppCommand,

  /// Settings of the command, as in a configuration file.
  pub config: SparseConfig,

  /// Where each input file came from, by setting key.
  pub sources: BTreeMap<String, UiDraftSource>,

  /// The run whose settings the form was loaded from.
  pub from_run_id: Option<JobId>,

  /// The run that holds the uploaded input files of the form.
  pub upload_run_id: Option<JobId>,

  pub view: UiSettingsView,

  /// Text of the settings search field.
  pub search: String,

  /// Whether the form shows only the settings that differ from their defaults.
  pub changed_only: bool,

  pub code_format: UiCodeFormat,
}

/// Origin of an input file of the form.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct UiDraftSource {
  /// Name of the file shown in the form.
  pub label: String,
  pub origin: UiDraftOrigin,
  /// Size of the file in bytes, when known.
  pub size: Option<u64>,
}

/// Where an input file of the form came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum UiDraftOrigin {
  /// An example dataset.
  Dataset,
  /// A file uploaded to the server.
  Upload,
  /// A file on the computer that runs TreeTime.
  Local,
  /// An input of an earlier run.
  Run,
  /// A path from a loaded configuration file.
  Config,
}

/// Which settings the form shows.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum UiSettingsView {
  /// The main settings of the command.
  Main,
  /// Every setting of the command.
  All,
}

/// Format of the code that reproduces the form.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum UiCodeFormat {
  /// A command line.
  Cli,
  /// A YAML configuration file.
  Yaml,
}

/// The runs folder of the running back end, and the folder used when the settings name none.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct Workspace {
  /// The runs folder in use.
  pub path: PathBuf,
  /// The runs folder used when the settings name none: `runs` in the app folder.
  pub default_path: PathBuf,
  /// The environment variable that sets the runs folder; the app cannot change the folder then.
  pub fixed_by: Option<String>,
  /// Why the runs folder named in the settings could not be opened. The default runs folder is in use then.
  pub error: Option<String>,
}

/// A new runs folder. It takes effect when the back end starts again.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct WorkspaceUpdate {
  /// Absolute path of the folder. Unset: `runs` in the app folder.
  pub path: Option<PathBuf>,
}

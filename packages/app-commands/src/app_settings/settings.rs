use crate::command::AppCommand;
use crate::job::JobId;
use crate::json_value::SparseConfig;
use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_grid::MaxGridPoints;
use treetime_schema::skip_serializing_optionals;

/// Settings of a local TreeTime installation, kept in `settings.yaml` or `settings.json` in the app folder. Every
/// setting is optional.
#[derive(Clone, Debug, Default, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct AppSettings {
  /// Folders of the app. Each defaults to a folder of the same name in the app folder.
  #[schemars(default, skip_serializing_if = "AppPathSettings::is_unset")]
  #[deser(default, skip_serializing_if = AppPathSettings::is_unset)]
  pub paths: AppPathSettings,

  /// Preferences of the user interface.
  #[schemars(default, skip_serializing_if = "UiSettings::is_unset")]
  #[deser(default, skip_serializing_if = UiSettings::is_unset)]
  pub ui: UiSettings,

  /// Settings of the analyses that the app runs.
  #[schemars(default, skip_serializing_if = "AnalysisSettings::is_unset")]
  #[deser(default, skip_serializing_if = AnalysisSettings::is_unset)]
  pub analysis: AnalysisSettings,
}

/// Folders of the app. A relative path is relative to the app folder. The environment variables
/// `TREETIME_PROFILE_DIR`, `TREETIME_RUNS_DIR`, `TREETIME_LOGS_DIR`, and `TREETIME_EXAMPLES_DIR` take precedence.
#[derive(Clone, Debug, Default, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct AppPathSettings {
  /// Browser profile of the desktop app: cache, local storage, and crash reports. Default: `profile`.
  #[schemars(default)]
  #[deser(default)]
  pub profile: Option<PathBuf>,

  /// Runs, each in its own folder. Default: `runs`.
  #[schemars(default)]
  #[deser(default)]
  pub runs: Option<PathBuf>,

  /// Logs and crash diagnostics. Default: `logs`.
  #[schemars(default)]
  #[deser(default)]
  pub logs: Option<PathBuf>,

  /// Example datasets and configurations that the app lists. Default: `examples`.
  #[schemars(default)]
  #[deser(default)]
  pub examples: Option<PathBuf>,
}

impl AppPathSettings {
  pub fn is_unset(&self) -> bool {
    *self == Self::default()
  }
}

/// Preferences of the user interface. Unset preferences take the defaults of the user interface.
#[derive(Clone, Debug, Default, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct UiSettings {
  /// Color theme.
  #[schemars(default)]
  #[deser(default)]
  pub theme: Option<UiTheme>,

  /// Width of the sidebar in pixels.
  #[schemars(default)]
  #[deser(default)]
  pub sidebar_width: Option<u32>,

  /// The unfinished analysis form.
  #[schemars(default)]
  #[deser(default)]
  pub draft: Option<UiDraft>,
}

impl UiSettings {
  pub fn is_unset(&self) -> bool {
    *self == Self::default()
  }
}

/// Settings of the analyses that the app runs. A setting here applies to every run whose configuration does not set
/// it.
#[derive(Clone, Debug, Default, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct AnalysisSettings {
  /// Largest number of points of one probability grid during time inference, for runs whose configuration does not
  /// set `max_grid_points`. Unset: 1000000.
  #[schemars(default)]
  #[deser(default)]
  pub max_grid_points: Option<MaxGridPoints>,
}

impl AnalysisSettings {
  pub fn is_unset(&self) -> bool {
    *self == Self::default()
  }
}

/// Color theme of the user interface.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(rename_all = "lowercase")]
#[deser(rename_all = "lowercase")]
pub enum UiTheme {
  /// Follow the theme of the operating system.
  System,
  Light,
  Dark,
}

/// The unfinished analysis form: the command, its settings, and how the form is shown.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
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
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct UiDraftSource {
  /// Name of the file shown in the form.
  pub label: String,
  pub origin: UiDraftOrigin,
  /// Size of the file in bytes, when known.
  pub size: Option<u64>,
}

/// Where an input file of the form came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(rename_all = "lowercase")]
#[deser(rename_all = "lowercase")]
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
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(rename_all = "lowercase")]
#[deser(rename_all = "lowercase")]
pub enum UiSettingsView {
  /// The main settings of the command.
  Main,
  /// Every setting of the command.
  All,
}

/// Format of the code that reproduces the form.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(rename_all = "lowercase")]
#[deser(rename_all = "lowercase")]
pub enum UiCodeFormat {
  /// A command line.
  Cli,
  /// A YAML configuration file.
  Yaml,
}

/// The runs folder of the running back end, and the folder used when the settings name none.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
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
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Deserialize)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct WorkspaceUpdate {
  /// Absolute path of the folder. Unset: `runs` in the app folder.
  pub path: Option<PathBuf>,
}

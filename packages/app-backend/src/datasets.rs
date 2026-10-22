use crate::check_inputs::InputKind;
use app_commands::command::AppCommand;
use app_commands::config::properties::{PathRole, setting_properties};
use app_commands::config::resolve_paths::resolve_config_paths_where;
use app_commands::yaml::yaml_value_read_str;
use app_datasets::{DatasetFiles, ExampleFile, TREE_FILE, discover_datasets, example_configs, group_datasets};
use deser::{Deserialize, Serialize};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use relative_path::{Component, RelativePath};
use schemars::JsonSchema;
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use strum::VariantNames;
use treetime_utils::io::compression::remove_compression_ext;
use treetime_utils::io::fs::absolute_path;
use treetime_utils::{make_error, make_report};

pub const DATA_INDEX_FILE: &str = "index.json";

const PIPELINE_COMMAND: &str = "pipeline";

const DATASET_INPUTS: [(InputKind, &[&str]); 3] = [
  (InputKind::Tree, &[TREE_FILE]),
  (InputKind::Alignment, &["aln.fasta"]),
  (InputKind::Metadata, &["metadata.tsv", "metadata.csv"]),
];

pub fn dataset_catalog(examples_dir: &Path) -> Result<DatasetCatalog, Report> {
  let discovered = discover_datasets(examples_dir, AppCommand::VARIANTS)?;
  let folder = absolute_path(examples_dir)?;
  Ok(DatasetCatalog {
    datasets: discovered
      .datasets
      .iter()
      .map(|dataset| Dataset::new(&discovered.examples_dir, dataset))
      .collect(),
    examples: discovered
      .examples
      .into_iter()
      .map(|example| ExampleConfig::new(&folder, example))
      .try_collect()?,
  })
}

/// Example datasets and example command configurations found in the examples folder.
#[derive(Clone, Debug, JsonSchema, Serialize)]
pub struct DatasetCatalog {
  /// Directories that hold a `tree.nwk`, with their files.
  pub datasets: Vec<Dataset>,
  /// Example configurations of the commands the application runs.
  pub examples: Vec<ExampleConfig>,
}

/// A directory of example input files.
#[derive(Clone, Debug, JsonSchema, Serialize)]
pub struct Dataset {
  /// Path of the directory relative to the examples folder, with `/` separators.
  pub name: String,
  /// Names of the files in the directory, sorted.
  pub files: Vec<String>,
  /// The files of the directory that fill command inputs, at most one per kind.
  pub inputs: Vec<DatasetInput>,
}

impl Dataset {
  fn new(examples_dir: &str, dataset: &DatasetFiles) -> Self {
    let inputs = dataset_inputs(dataset)
      .map(|(kind, file)| DatasetInput {
        kind,
        file: format!("{}/{file}", dataset.name),
        path: [examples_dir, dataset.name.as_str(), file.as_str()]
          .iter()
          .filter(|part| !part.is_empty())
          .copied()
          .collect::<Vec<_>>()
          .join("/"),
      })
      .collect();
    Self {
      name: dataset.name.clone(),
      files: dataset.files.clone(),
      inputs,
    }
  }
}

/// A file of an example dataset that fills a command input.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct DatasetInput {
  /// Input the file fills.
  pub kind: InputKind,
  /// Path of the file relative to the examples folder, with `/` separators.
  pub file: String,
  /// Path of the file as a run configuration names it.
  pub path: String,
}

/// An example configuration file of one command.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct ExampleConfig {
  /// Path of the file relative to the examples folder, with `/` separators.
  pub path: String,
  /// Command the configuration is for, taken from the `$schema` URL of its `yaml-language-server` directive.
  pub command: AppCommand,
  /// Title of the example: the first comment line after the directive.
  pub title: String,
  /// Absolute folder of the file, which relative paths in the file resolve from.
  pub folder: PathBuf,
  /// Text of the file.
  pub content: String,
}

impl ExampleConfig {
  fn new(examples_dir: &Path, example: ExampleFile) -> Result<Self, Report> {
    let file = examples_dir.join(&example.path);
    Ok(Self {
      command: AppCommand::from_str(&example.command)
        .map_err(|err| make_report!("example `{}` names an unknown command: {err}", example.path))?,
      folder: file
        .parent()
        .map_or_else(|| examples_dir.to_path_buf(), Path::to_path_buf),
      path: example.path,
      title: example.title,
      content: example.content,
    })
  }
}

/// Index of the example data in the `data/` folder of the repository: every tracked file, the datasets, the example
/// configurations of the commands, and the example pipelines. Every path is relative to `data/` and uses `/`
/// separators. The apps read the example data only through this index.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct DataIndex {
  /// Every tracked file of the data folder except the index itself, sorted by path.
  pub files: Vec<IndexedFile>,
  /// Folders that hold a `tree.nwk`, with the files that fill command inputs, sorted by name.
  pub datasets: Vec<IndexedDataset>,
  /// Example configurations of the commands the application runs, sorted by path.
  pub configs: Vec<IndexedConfig>,
  /// Example pipeline configurations, sorted by path.
  pub pipelines: Vec<IndexedPipeline>,
}

impl DataIndex {
  pub fn new(
    paths: Vec<String>,
    mut size: impl FnMut(&str) -> Result<u64, Report>,
    read: impl FnMut(&str) -> Result<String, Report>,
    pipeline_steps: impl Fn(&str) -> Result<Vec<IndexedStep>, Report>,
  ) -> Result<Self, Report> {
    if let Some(path) = paths.iter().find(|path| !is_data_path(path)) {
      return make_error!("the file list names '{path}', which is not a normalized path relative to the data folder");
    }
    let files: Vec<IndexedFile> = paths
      .into_iter()
      .filter(|path| path != DATA_INDEX_FILE)
      .sorted()
      .dedup()
      .map(|path| {
        let size = size(&path)?;
        Ok::<_, Report>(IndexedFile { path, size })
      })
      .try_collect()?;
    let tracked: BTreeSet<&str> = files.iter().map(|file| file.path.as_str()).collect();
    let datasets = group_datasets(tracked.iter().copied())
      .iter()
      .map(IndexedDataset::new)
      .collect();
    let commands = AppCommand::VARIANTS
      .iter()
      .copied()
      .chain([PIPELINE_COMMAND])
      .collect_vec();
    let (pipelines, configs): (Vec<_>, Vec<_>) = example_configs(tracked.iter().copied(), read, &commands)?
      .into_iter()
      .partition(|example| example.command == PIPELINE_COMMAND);
    let configs = configs
      .into_iter()
      .map(|example| IndexedConfig::new(example, &tracked))
      .try_collect()?;
    let pipelines = pipelines
      .into_iter()
      .map(|example| IndexedPipeline::new(example, &pipeline_steps))
      .try_collect()?;
    Ok(Self {
      files,
      datasets,
      configs,
      pipelines,
    })
  }
}

/// A tracked file of the data folder.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct IndexedFile {
  /// Path relative to the data folder.
  pub path: String,
  /// Size in bytes.
  pub size: u64,
}

/// A folder of the data folder that holds a `tree.nwk`.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct IndexedDataset {
  /// Path of the folder relative to the data folder.
  pub name: String,
  /// The files of the folder that fill command inputs, at most one per kind. A compressed file fills the input of its
  /// uncompressed name.
  pub inputs: Vec<IndexedInput>,
}

impl IndexedDataset {
  fn new(dataset: &DatasetFiles) -> Self {
    Self {
      name: dataset.name.clone(),
      inputs: dataset_inputs(dataset)
        .map(|(kind, file)| IndexedInput {
          kind,
          path: join_data_path(&dataset.name, file),
        })
        .collect(),
    }
  }
}

/// A file of a dataset that fills a command input.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct IndexedInput {
  /// Input the file fills.
  pub kind: InputKind,
  /// Path of the file relative to the data folder.
  pub path: String,
}

/// An example configuration of one command.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct IndexedConfig {
  /// Path of the file relative to the data folder.
  pub path: String,
  /// Command the configuration is for, taken from the `$schema` URL of its `yaml-language-server` directive.
  pub command: AppCommand,
  /// Title of the example: the first comment line after the directive.
  pub title: String,
  /// Folder of the file relative to the data folder, which relative paths in the file resolve from.
  pub folder: String,
  /// Text of the file.
  pub content: String,
  /// The files that each input setting of the configuration reads, relative to the data folder, by setting key.
  pub inputs: BTreeMap<String, Vec<String>>,
}

impl IndexedConfig {
  fn new(example: ExampleFile, tracked: &BTreeSet<&str>) -> Result<Self, Report> {
    let path = example.path;
    let command = AppCommand::from_str(&example.command)
      .map_err(|err| make_report!("example `{path}` names an unknown command: {err}"))?;
    let folder = data_folder(&path).to_owned();
    let inputs = config_inputs(command, &example.content, &folder, tracked)
      .wrap_err_with(|| format!("When indexing the input files of example configuration '{path}'"))?;
    Ok(Self {
      path,
      command,
      title: example.title,
      folder,
      content: example.content,
      inputs,
    })
  }
}

/// An example pipeline configuration.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct IndexedPipeline {
  /// Path of the file relative to the data folder.
  pub path: String,
  /// Title of the example: the first comment line after the directive.
  pub title: String,
  /// Steps of the pipeline, in the order they run.
  pub steps: Vec<IndexedStep>,
}

impl IndexedPipeline {
  fn new(
    example: ExampleFile,
    pipeline_steps: &impl Fn(&str) -> Result<Vec<IndexedStep>, Report>,
  ) -> Result<Self, Report> {
    let steps = pipeline_steps(&example.content)
      .wrap_err_with(|| format!("When reading the steps of example pipeline '{}'", example.path))?;
    Ok(Self {
      path: example.path,
      title: example.title,
      steps,
    })
  }
}

/// One step of an example pipeline.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct IndexedStep {
  /// Name of the step.
  pub name: String,
  /// Command the step runs.
  pub command: AppCommand,
}

fn dataset_inputs(dataset: &DatasetFiles) -> impl Iterator<Item = (InputKind, &String)> {
  DATASET_INPUTS.iter().filter_map(move |(kind, names)| {
    names
      .iter()
      .find_map(|name| {
        dataset
          .files
          .iter()
          .find(|file| remove_compression_ext(file.as_str()) == Path::new(name))
      })
      .map(|file| (*kind, file))
  })
}

fn is_data_path(path: &str) -> bool {
  let normalized = RelativePath::new(path).normalize();
  !path.is_empty()
    && normalized.as_str() == path
    && !matches!(normalized.components().next(), Some(Component::ParentDir))
}

fn join_data_path(folder: &str, name: &str) -> String {
  if folder.is_empty() {
    name.to_owned()
  } else {
    format!("{folder}/{name}")
  }
}

fn data_folder(path: &str) -> &str {
  path.rsplit_once('/').map_or("", |(folder, _)| folder)
}

fn config_inputs(
  command: AppCommand,
  content: &str,
  folder: &str,
  tracked: &BTreeSet<&str>,
) -> Result<BTreeMap<String, Vec<String>>, Report> {
  let schema = command.config_schema();
  let properties = setting_properties(schema.as_value())?;
  let mut document = yaml_value_read_str(content)?;
  if let Some(template) = properties
    .iter()
    .find(|property| property.path_role == Some(PathRole::InputTemplate) && document.get(&property.key).is_some())
  {
    return make_error!(
      "the input template setting `{}` names files by a placeholder, which the index cannot list",
      template.key
    );
  }
  let input_keys: BTreeSet<&str> = properties
    .iter()
    .filter(|property| property.path_role == Some(PathRole::Input))
    .map(|property| property.key.as_str())
    .collect();
  resolve_config_paths_where(&mut document, schema.as_value(), Path::new(folder), |key| {
    input_keys.contains(key)
  })?;
  input_keys
    .iter()
    .filter_map(|key| document.get(*key).map(|value| (*key, value)))
    .map(|(key, value)| {
      let files: Vec<String> = setting_paths(value)
        .and_then(|paths| paths.into_iter().map(|path| tracked_path(path, tracked)).try_collect())
        .wrap_err_with(|| format!("When reading the input setting `{key}`"))?;
      Ok::<_, Report>((key.to_owned(), files))
    })
    .try_collect()
}

fn setting_paths(value: &Value) -> Result<Vec<&str>, Report> {
  match value {
    Value::String(path) => Ok(vec![path.as_str()]),
    Value::Array(paths) => paths
      .iter()
      .map(|path| {
        path
          .as_str()
          .ok_or_else(|| make_report!("the list holds {path}, which is not a path"))
      })
      .try_collect(),
    other => make_error!("the value {other} is not a path"),
  }
}

fn tracked_path(path: &str, tracked: &BTreeSet<&str>) -> Result<String, Report> {
  if Path::new(path).has_root() {
    return make_error!("the path '{path}' is outside the data folder");
  }
  let normalized = RelativePath::new(path).normalize();
  if matches!(normalized.components().next(), Some(Component::ParentDir)) {
    return make_error!("the path '{path}' leaves the data folder");
  }
  if !tracked.contains(normalized.as_str()) {
    return make_error!("the path '{path}' names no tracked file of the data folder");
  }
  Ok(normalized.into_string())
}

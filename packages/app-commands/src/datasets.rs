use crate::check_inputs::InputKind;
use crate::command::AppCommand;
use app_datasets::{DatasetFiles, ExampleFile, TREE_FILE, discover_datasets};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use strum::VariantNames;
use treetime_utils::io::fs::absolute_path;
use treetime_utils::make_report;

const DATASET_FILES: [(InputKind, &[&str]); 3] = [
  (InputKind::Tree, &[TREE_FILE]),
  (InputKind::Alignment, &["aln.fasta.xz", "aln.fasta"]),
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
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct DatasetCatalog {
  /// Directories that hold a `tree.nwk`, with their files.
  pub datasets: Vec<Dataset>,
  /// Example configurations of the commands the application runs.
  pub examples: Vec<ExampleConfig>,
}

/// A directory of example input files.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
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
    let inputs = DATASET_FILES
      .iter()
      .filter_map(|(kind, names)| {
        let file = names
          .iter()
          .find(|name| dataset.files.iter().any(|file| file == *name))?;
        Some(DatasetInput {
          kind: *kind,
          file: format!("{}/{file}", dataset.name),
          path: [examples_dir, &dataset.name, file]
            .iter()
            .filter(|part| !part.is_empty())
            .copied()
            .collect::<Vec<_>>()
            .join("/"),
        })
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
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct DatasetInput {
  /// Input the file fills.
  pub kind: InputKind,
  /// Path of the file relative to the examples folder, with `/` separators.
  pub file: String,
  /// Path of the file as a run configuration names it.
  pub path: String,
}

/// An example configuration file of one command.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
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

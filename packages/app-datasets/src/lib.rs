#[cfg(test)]
mod __tests__;

use eyre::{Report, WrapErr};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::fs::{self, ReadDir};
use std::io::{self, ErrorKind};
use std::path::{Path, PathBuf};

const TREE_FILE: &str = "tree.nwk";
const SCHEMA_DIRECTIVE: &str = "yaml-language-server:";
const SCHEMA_KEY: &str = "$schema=";
const SCHEMA_FILE_PREFIX: &str = "input-config-";
const SCHEMA_FILE_SUFFIX: &str = ".schema.json";
const YAML_EXTENSIONS: [&str; 2] = ["yaml", "yml"];

pub fn discover_datasets(data_dir: &Path, commands: &[&str]) -> Result<DatasetCatalog, Report> {
  let mut datasets = Vec::new();
  let mut examples = Vec::new();
  match fs::read_dir(data_dir) {
    Ok(entries) => collect(data_dir, data_dir, entries, commands, &mut datasets, &mut examples)?,
    Err(error) if error.kind() == ErrorKind::NotFound => {},
    Err(error) => return Err(read_dir_report(error, data_dir)),
  }
  datasets.sort_by(|a, b| a.name.cmp(&b.name));
  examples.sort_by(|a, b| a.path.cmp(&b.path));
  Ok(DatasetCatalog {
    data_dir: data_dir.to_string_lossy().replace('\\', "/"),
    datasets,
    examples,
  })
}

/// Example datasets and example command configurations found in the data directory.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct DatasetCatalog {
  /// Data directory as a run configuration names it: file `f` of dataset `d` is `<data_dir>/<d>/<f>`, a path
  /// relative to the working directory of the process that runs the commands, as in the example configurations.
  pub data_dir: String,
  /// Directories that hold a `tree.nwk`, with their files.
  pub datasets: Vec<DatasetInfo>,
  /// Example configurations of the commands the application runs.
  pub examples: Vec<ExampleConfig>,
}

/// A directory of example input files.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct DatasetInfo {
  /// Path of the directory relative to the data directory, with `/` separators.
  pub name: String,
  /// Names of the files in the directory, sorted.
  pub files: Vec<String>,
}

/// An example configuration file of one command.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct ExampleConfig {
  /// Path of the file relative to the data directory, with `/` separators.
  pub path: String,
  /// Command the configuration is for, taken from the `$schema` URL of its `yaml-language-server` directive.
  pub command: String,
  /// Title of the example: the first comment line after the directive.
  pub title: String,
  /// Text of the file.
  pub content: String,
}

pub fn parse_example_config(path: &str, content: &str, commands: &[&str]) -> Option<ExampleConfig> {
  let mut lines = content.lines().map(str::trim).skip_while(|line| line.is_empty());
  let command = schema_command(lines.next()?)?;
  if !commands.contains(&command) {
    return None;
  }
  let title = lines
    .filter_map(|line| line.strip_prefix('#'))
    .map(str::trim)
    .find(|text| !text.is_empty())?;
  Some(ExampleConfig {
    path: path.to_owned(),
    command: command.to_owned(),
    title: title.to_owned(),
    content: content.to_owned(),
  })
}

fn schema_command(directive: &str) -> Option<&str> {
  let directive = directive
    .strip_prefix('#')?
    .trim()
    .strip_prefix(SCHEMA_DIRECTIVE)?
    .trim();
  let url = directive.strip_prefix(SCHEMA_KEY)?.trim();
  let file = url.rsplit('/').next()?;
  file.strip_prefix(SCHEMA_FILE_PREFIX)?.strip_suffix(SCHEMA_FILE_SUFFIX)
}

fn collect(
  base: &Path,
  dir: &Path,
  entries: ReadDir,
  commands: &[&str],
  datasets: &mut Vec<DatasetInfo>,
  examples: &mut Vec<ExampleConfig>,
) -> Result<(), Report> {
  let mut files: Vec<String> = Vec::new();
  let mut subdirs: Vec<PathBuf> = Vec::new();

  for entry in entries {
    let entry = entry.map_err(|error| read_dir_report(error, dir))?;
    let file_type = entry
      .file_type()
      .wrap_err_with(|| format!("When reading the file type of '{}'", entry.path().display()))?;
    let name = entry.file_name().to_string_lossy().into_owned();
    if file_type.is_dir() {
      subdirs.push(entry.path());
    } else if !name.starts_with('.') {
      if is_yaml(&entry.path()) {
        let content = fs::read_to_string(entry.path())
          .wrap_err_with(|| format!("When reading example configuration '{}'", entry.path().display()))?;
        if let Some(example) = parse_example_config(&relative(base, &entry.path()), &content, commands) {
          examples.push(example);
        }
      }
      files.push(name);
    }
  }

  if files.iter().any(|file| file == TREE_FILE) {
    files.sort();
    datasets.push(DatasetInfo {
      name: relative(base, dir),
      files,
    });
  }

  for subdir in subdirs {
    let entries = fs::read_dir(&subdir).map_err(|error| read_dir_report(error, &subdir))?;
    collect(base, &subdir, entries, commands, datasets, examples)?;
  }
  Ok(())
}

fn is_yaml(path: &Path) -> bool {
  path
    .extension()
    .and_then(|extension| extension.to_str())
    .is_some_and(|extension| YAML_EXTENSIONS.contains(&extension))
}

fn relative(base: &Path, path: &Path) -> String {
  path
    .strip_prefix(base)
    .unwrap_or(path)
    .to_string_lossy()
    .replace('\\', "/")
}

fn read_dir_report(error: io::Error, dir: &Path) -> Report {
  Report::new(error).wrap_err(format!("When listing dataset directory '{}'", dir.display()))
}

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor(unsafe)]
  fn init() {
    global_init();
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}

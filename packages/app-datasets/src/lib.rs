#[cfg(test)]
mod __tests__;

use eyre::{Report, WrapErr};
use std::fs::{self, ReadDir};
use std::io::{self, ErrorKind};
use std::path::{Path, PathBuf};

pub const TREE_FILE: &str = "tree.nwk";
const SCHEMA_DIRECTIVE: &str = "yaml-language-server:";
const SCHEMA_KEY: &str = "$schema=";
const SCHEMA_FILE_PREFIX: &str = "input-config-";
const SCHEMA_FILE_SUFFIX: &str = ".schema.json";
const SCHEMA_URL_BASE: &str = "https://raw.githubusercontent.com/neherlab/treetime/rust/packages/schemas";
const YAML_EXTENSIONS: [&str; 2] = ["yaml", "yml"];

pub fn discover_datasets(data_dir: &Path, commands: &[&str]) -> Result<DataDirectory, Report> {
  let mut datasets = Vec::new();
  let mut examples = Vec::new();
  match fs::read_dir(data_dir) {
    Ok(entries) => collect(data_dir, data_dir, entries, commands, &mut datasets, &mut examples)?,
    Err(error) if error.kind() == ErrorKind::NotFound => {},
    Err(error) => return Err(read_dir_report(error, data_dir)),
  }
  datasets.sort_by(|a, b| a.name.cmp(&b.name));
  examples.sort_by(|a, b| a.path.cmp(&b.path));
  Ok(DataDirectory {
    data_dir: data_dir.to_string_lossy().replace('\\', "/"),
    datasets,
    examples,
  })
}

#[derive(Clone, Debug)]
pub struct DataDirectory {
  pub data_dir: String,
  pub datasets: Vec<DatasetFiles>,
  pub examples: Vec<ExampleFile>,
}

#[derive(Clone, Debug)]
pub struct DatasetFiles {
  pub name: String,
  pub files: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExampleFile {
  pub path: String,
  pub command: String,
  pub title: String,
  pub content: String,
}

pub fn schema_directive(command: &str) -> String {
  format!("# {SCHEMA_DIRECTIVE} {SCHEMA_KEY}{SCHEMA_URL_BASE}/{SCHEMA_FILE_PREFIX}{command}{SCHEMA_FILE_SUFFIX}")
}

pub fn text_schema_command(content: &str) -> Option<&str> {
  schema_command(content.lines().map(str::trim).find(|line| !line.is_empty())?)
}

pub fn parse_example_config(path: &str, content: &str, commands: &[&str]) -> Option<ExampleFile> {
  let command = text_schema_command(content)?;
  if !commands.contains(&command) {
    return None;
  }
  let title = content
    .lines()
    .map(str::trim)
    .skip_while(|line| line.is_empty())
    .skip(1)
    .filter_map(|line| line.strip_prefix('#'))
    .map(str::trim)
    .find(|text| !text.is_empty())?;
  Some(ExampleFile {
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
  datasets: &mut Vec<DatasetFiles>,
  examples: &mut Vec<ExampleFile>,
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
    datasets.push(DatasetFiles {
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

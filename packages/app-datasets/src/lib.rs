#[cfg(test)]
mod __tests__;

use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::fs;
use std::io::{self, ErrorKind};
use std::path::Path;

pub const TREE_FILE: &str = "tree.nwk";
const SCHEMA_DIRECTIVE: &str = "yaml-language-server:";
const SCHEMA_KEY: &str = "$schema=";
const SCHEMA_FILE_PREFIX: &str = "input-config-";
const SCHEMA_FILE_SUFFIX: &str = ".schema.json";
const SCHEMA_URL_BASE: &str = "https://raw.githubusercontent.com/neherlab/treetime/rust/packages/schemas";
const YAML_EXTENSIONS: [&str; 2] = ["yaml", "yml"];

pub fn discover_datasets(examples_dir: &Path, commands: &[&str]) -> Result<ExamplesFolder, Report> {
  let files = list_files(examples_dir)?;
  let examples = example_configs(
    files.iter().map(String::as_str),
    |path| {
      let file = examples_dir.join(path);
      fs::read_to_string(&file).wrap_err_with(|| format!("When reading example configuration '{}'", file.display()))
    },
    commands,
  )?;
  Ok(ExamplesFolder {
    examples_dir: examples_dir.to_string_lossy().replace('\\', "/"),
    datasets: group_datasets(files.iter().map(String::as_str)),
    examples,
  })
}

#[derive(Clone, Debug)]
pub struct ExamplesFolder {
  pub examples_dir: String,
  pub datasets: Vec<DatasetFiles>,
  pub examples: Vec<ExampleFile>,
}

pub fn group_datasets<'a>(paths: impl IntoIterator<Item = &'a str>) -> Vec<DatasetFiles> {
  paths
    .into_iter()
    .filter_map(split_path)
    .filter(|(_, name)| !name.starts_with('.'))
    .into_group_map()
    .into_iter()
    .filter(|(_, files)| files.contains(&TREE_FILE))
    .map(|(folder, files)| DatasetFiles {
      name: folder.to_owned(),
      files: files.into_iter().map(str::to_owned).sorted().collect(),
    })
    .sorted_by(|a, b| a.name.cmp(&b.name))
    .collect()
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DatasetFiles {
  pub name: String,
  pub files: Vec<String>,
}

pub fn example_configs<'a>(
  paths: impl IntoIterator<Item = &'a str>,
  mut read: impl FnMut(&str) -> Result<String, Report>,
  commands: &[&str],
) -> Result<Vec<ExampleFile>, Report> {
  let mut examples = Vec::new();
  for path in paths.into_iter().sorted() {
    if !is_yaml(path) || split_path(path).is_none_or(|(_, name)| name.starts_with('.')) {
      continue;
    }
    if let Some(example) = parse_example_config(path, &read(path)?, commands) {
      examples.push(example);
    }
  }
  Ok(examples)
}

pub fn schema_directive(command: &str) -> String {
  format!("# {SCHEMA_DIRECTIVE} {SCHEMA_KEY}{SCHEMA_URL_BASE}/{SCHEMA_FILE_PREFIX}{command}{SCHEMA_FILE_SUFFIX}")
}

pub(crate) fn parse_example_config(path: &str, content: &str, commands: &[&str]) -> Option<ExampleFile> {
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExampleFile {
  pub path: String,
  pub command: String,
  pub title: String,
  pub content: String,
}

pub fn text_schema_command(content: &str) -> Option<&str> {
  schema_command(content.lines().map(str::trim).find(|line| !line.is_empty())?)
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

fn is_yaml(path: &str) -> bool {
  Path::new(path)
    .extension()
    .and_then(|extension| extension.to_str())
    .is_some_and(|extension| YAML_EXTENSIONS.contains(&extension))
}

fn split_path(path: &str) -> Option<(&str, &str)> {
  match path.rsplit_once('/') {
    Some((folder, name)) => Some((folder, name)),
    None if path.is_empty() => None,
    None => Some(("", path)),
  }
}

fn list_files(examples_dir: &Path) -> Result<Vec<String>, Report> {
  let mut files = Vec::new();
  let mut folders = vec![examples_dir.to_path_buf()];
  while let Some(dir) = folders.pop() {
    let entries = match fs::read_dir(&dir) {
      Ok(entries) => entries,
      Err(error) if error.kind() == ErrorKind::NotFound && dir == examples_dir => return Ok(files),
      Err(error) => return Err(read_dir_report(error, &dir)),
    };
    for entry in entries {
      let entry = entry.map_err(|error| read_dir_report(error, &dir))?;
      let file_type = entry
        .file_type()
        .wrap_err_with(|| format!("When reading the file type of '{}'", entry.path().display()))?;
      if file_type.is_dir() {
        folders.push(entry.path());
      } else {
        files.push(relative(examples_dir, &entry.path()));
      }
    }
  }
  files.sort();
  Ok(files)
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

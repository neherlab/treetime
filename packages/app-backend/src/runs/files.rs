use crate::command::{OutputAlignment, OutputFile};
use crate::runs::errors::{invalid, not_found};
use app_output::output_plan::OutputSelection;
use deser::Serialize;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::JsonSchema;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::Read;
use std::path::{Component, Path, PathBuf};
use treetime_schema::skip_serializing_optionals;

pub const RUN_DOWNLOAD_CHUNK_SIZE: usize = 1024 * 1024;

pub fn list_run_files(out_dir: &Path, output_files: &[OutputFile]) -> Result<Vec<RunFile>, Report> {
  let planned: BTreeMap<String, &OutputFile> = output_files
    .iter()
    .map(|file| (relative_name(&file.path), file))
    .collect();
  let mut files = vec![];
  for path in walk_files(out_dir)? {
    let relative = path.strip_prefix(out_dir)?;
    let name = relative_name(relative);
    let size = usize::try_from(
      fs::metadata(&path)
        .wrap_err_with(|| format!("When reading the size of '{}'", path.display()))?
        .len(),
    )?;
    let file = planned.get(&name);
    let kind = file.map(|file| file.kind);
    files.push(RunFile {
      description: kind.map(OutputSelection::description).unwrap_or_default().to_owned(),
      kind,
      alignment: file.and_then(|file| file.alignment.clone()),
      path: name,
      size,
    });
  }
  Ok(files.into_iter().sorted_by(|a, b| a.path.cmp(&b.path)).collect())
}

/// One file in a run's `out/` folder.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct RunFile {
  /// Path relative to the run's `out/` folder, with `/` separators.
  pub path: String,
  /// Size in bytes.
  pub size: usize,
  /// Output selection that produced the file; absent for files the output plan does not name.
  pub kind: Option<OutputSelection>,
  /// What the file holds; empty for files the output plan does not name.
  pub description: String,
  /// Alignment file the output belongs to, for an output written once per alignment file.
  pub alignment: Option<OutputAlignment>,
}

pub fn resolve_run_file(out_dir: &Path, relative: &str) -> Result<PathBuf, Report> {
  let requested = inside_path(relative)?;
  let root = out_dir
    .canonicalize()
    .wrap_err_with(|| format!("When resolving '{}'", out_dir.display()))?;
  let resolved = root
    .join(requested)
    .canonicalize()
    .map_err(|err| not_found(format!("the run has no output file `{relative}`: {err}")))?;
  if !resolved.starts_with(&root) {
    return Err(invalid(format!(
      "file path `{relative}` must name a file inside the run's output folder"
    )));
  }
  if !resolved.is_file() {
    return Err(not_found(format!("the run has no output file `{relative}`")));
  }
  Ok(resolved)
}

pub fn run_file_destination(out_dir: &Path, relative: &str) -> Result<PathBuf, Report> {
  Ok(out_dir.join(inside_path(relative)?))
}

pub struct FileChunks {
  file: File,
  path: PathBuf,
}

impl FileChunks {
  pub fn open(path: &Path) -> Result<Self, Report> {
    Ok(Self {
      file: File::open(path).wrap_err_with(|| format!("When opening '{}'", path.display()))?,
      path: path.to_owned(),
    })
  }
}

impl Iterator for FileChunks {
  type Item = Result<Vec<u8>, Report>;

  fn next(&mut self) -> Option<Self::Item> {
    let mut block = vec![0; RUN_DOWNLOAD_CHUNK_SIZE];
    match self.file.read(&mut block) {
      Ok(0) => None,
      Ok(length) => {
        block.truncate(length);
        Some(Ok(block))
      },
      Err(error) => Some(Err(error).wrap_err_with(|| format!("When reading '{}'", self.path.display()))),
    }
  }
}

fn inside_path(relative: &str) -> Result<&Path, Report> {
  let requested = Path::new(relative);
  let escapes = requested
    .components()
    .any(|component| !matches!(component, Component::Normal(_)));
  if relative.is_empty() || escapes {
    return Err(invalid(format!(
      "file path `{relative}` must name a file inside the run's output folder"
    )));
  }
  Ok(requested)
}

pub fn walk_files(dir: &Path) -> Result<Vec<PathBuf>, Report> {
  let mut files = vec![];
  if !dir.is_dir() {
    return Ok(files);
  }
  for entry in fs::read_dir(dir).wrap_err_with(|| format!("When listing '{}'", dir.display()))? {
    let entry = entry.wrap_err_with(|| format!("When listing '{}'", dir.display()))?;
    let file_type = entry
      .file_type()
      .wrap_err_with(|| format!("When reading the type of '{}'", entry.path().display()))?;
    if file_type.is_dir() {
      files.extend(walk_files(&entry.path())?);
    } else if !file_type.is_symlink() {
      files.push(entry.path());
    }
  }
  Ok(files.into_iter().sorted().collect())
}

pub fn relative_name(path: &Path) -> String {
  path
    .components()
    .map(|component| component.as_os_str().to_string_lossy().into_owned())
    .join("/")
}

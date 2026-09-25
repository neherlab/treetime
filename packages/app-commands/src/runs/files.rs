use crate::command::OutputFile;
use crate::runs::errors::{invalid, not_found};
use app_output::output_plan::OutputSelection;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, Seek, Write};
use std::path::{Component, Path, PathBuf};
use zip::CompressionMethod;
use zip::ZipWriter;
use zip::write::SimpleFileOptions;

/// One file in a run's `out/` folder.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct RunFile {
  /// Path relative to the run's `out/` folder, with `/` separators.
  pub path: String,
  /// Size in bytes.
  pub size: usize,
  /// Output selection that produced the file; absent for files the output plan does not name.
  pub kind: Option<OutputSelection>,
  /// What the file holds; empty for files the output plan does not name.
  pub description: String,
}

pub fn list_run_files(out_dir: &Path, output_files: &[OutputFile]) -> Result<Vec<RunFile>, Report> {
  let kinds: BTreeMap<String, OutputSelection> = output_files
    .iter()
    .map(|file| (relative_name(&file.path), file.kind))
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
    let kind = kinds.get(&name).copied();
    files.push(RunFile {
      description: kind.map(OutputSelection::description).unwrap_or_default().to_owned(),
      kind,
      path: name,
      size,
    });
  }
  Ok(files.into_iter().sorted_by(|a, b| a.path.cmp(&b.path)).collect())
}

pub fn resolve_run_file(out_dir: &Path, relative: &str) -> Result<PathBuf, Report> {
  let requested = Path::new(relative);
  let escapes = requested
    .components()
    .any(|component| !matches!(component, Component::Normal(_)));
  if relative.is_empty() || escapes {
    return Err(invalid(format!(
      "file path `{relative}` must name a file inside the run's output folder"
    )));
  }
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

pub fn write_run_zip(out_dir: &Path, folder_name: &str, writer: impl Write + Seek) -> Result<(), Report> {
  let mut zip = ZipWriter::new(writer);
  let options = SimpleFileOptions::default().compression_method(CompressionMethod::Deflated);
  for path in walk_files(out_dir)? {
    let name = format!("{folder_name}/{}", relative_name(path.strip_prefix(out_dir)?));
    zip
      .start_file(name, options)
      .wrap_err_with(|| format!("When adding '{}' to the archive", path.display()))?;
    let mut file = File::open(&path).wrap_err_with(|| format!("When opening '{}'", path.display()))?;
    io::copy(&mut file, &mut zip).wrap_err_with(|| format!("When archiving '{}'", path.display()))?;
  }
  zip.finish().wrap_err("When finishing the archive")?;
  Ok(())
}

fn walk_files(dir: &Path) -> Result<Vec<PathBuf>, Report> {
  let mut files = vec![];
  if !dir.is_dir() {
    return Ok(files);
  }
  for entry in fs::read_dir(dir).wrap_err_with(|| format!("When listing '{}'", dir.display()))? {
    let entry = entry?;
    let file_type = entry.file_type()?;
    if file_type.is_dir() {
      files.extend(walk_files(&entry.path())?);
    } else if !file_type.is_symlink() {
      files.push(entry.path());
    }
  }
  Ok(files.into_iter().sorted().collect())
}

fn relative_name(path: &Path) -> String {
  path
    .components()
    .map(|component| component.as_os_str().to_string_lossy().into_owned())
    .join("/")
}

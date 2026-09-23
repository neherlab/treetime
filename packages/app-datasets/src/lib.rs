use eyre::{Report, WrapErr};
use serde::Serialize;
use std::fs::{self, ReadDir};
use std::io::{self, ErrorKind};
use std::path::{Path, PathBuf};

pub fn discover_datasets(data_dir: &Path) -> Result<Vec<DatasetInfo>, Report> {
  let mut datasets = Vec::new();
  match fs::read_dir(data_dir) {
    Ok(entries) => collect_datasets(data_dir, data_dir, entries, &mut datasets)?,
    Err(error) if error.kind() == ErrorKind::NotFound => {},
    Err(error) => return Err(read_dir_report(error, data_dir)),
  }
  datasets.sort_by(|a, b| a.name.cmp(&b.name));
  Ok(datasets)
}

fn collect_datasets(base: &Path, dir: &Path, entries: ReadDir, out: &mut Vec<DatasetInfo>) -> Result<(), Report> {
  let mut files: Vec<String> = Vec::new();
  let mut subdirs: Vec<PathBuf> = Vec::new();

  for entry in entries {
    let entry = entry.map_err(|error| read_dir_report(error, dir))?;
    let ft = entry
      .file_type()
      .wrap_err_with(|| format!("When reading the file type of '{}'", entry.path().display()))?;
    let name = entry.file_name().to_string_lossy().into_owned();
    if ft.is_dir() {
      subdirs.push(entry.path());
    } else if !ft.is_dir() && !name.starts_with('.') {
      files.push(name);
    }
  }

  if files.iter().any(|f| f == "tree.nwk") {
    let rel = dir
      .strip_prefix(base)
      .unwrap_or(dir)
      .to_string_lossy()
      .replace('\\', "/");
    files.sort();
    out.push(DatasetInfo { name: rel, files });
  }

  for subdir in subdirs {
    let entries = fs::read_dir(&subdir).map_err(|error| read_dir_report(error, &subdir))?;
    collect_datasets(base, &subdir, entries, out)?;
  }
  Ok(())
}

#[derive(Debug, Serialize)]
pub struct DatasetInfo {
  name: String,
  files: Vec<String>,
}

fn read_dir_report(error: io::Error, dir: &Path) -> Report {
  Report::new(error).wrap_err(format!("When listing dataset directory '{}'", dir.display()))
}

use serde::Serialize;
use std::path::{Path, PathBuf};

#[derive(Debug, Serialize)]
pub struct DatasetInfo {
  pub name: String,
  pub files: Vec<String>,
}

pub fn discover_datasets(data_dir: &Path) -> Vec<DatasetInfo> {
  let mut datasets = Vec::new();
  collect_datasets(data_dir, data_dir, &mut datasets);
  datasets.sort_by(|a, b| a.name.cmp(&b.name));
  datasets
}

fn collect_datasets(base: &Path, dir: &Path, out: &mut Vec<DatasetInfo>) {
  let Ok(entries) = std::fs::read_dir(dir) else {
    return;
  };

  let mut files: Vec<String> = Vec::new();
  let mut subdirs: Vec<PathBuf> = Vec::new();

  for entry in entries.flatten() {
    let Ok(ft) = entry.file_type() else {
      continue;
    };
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
    collect_datasets(base, &subdir, out);
  }
}

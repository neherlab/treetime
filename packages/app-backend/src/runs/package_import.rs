use crate::runs::errors::invalid;
use crate::runs::input_store::StoredUpload;
use crate::runs::inputs::file_sha256;
use crate::runs::manager::RunManager;
use crate::runs::record::{RunInput, RunRecord};
use crate::runs::run_paths::{INPUTS_DIR, OUT_DIR, map_path_values};
use crate::runs::store::{EVENTS_FILE, RUN_FILE};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::{self, Read, Seek};
use std::path::{Component, Path};
use treetime_utils::io::json::json_read_str;
use zip::ZipArchive;

const PACKAGE_FILE: &str = "package.zip";

const STAGED_RUN_DIR: &str = "run";

pub fn import_run(manager: &RunManager, body: &mut dyn Read) -> Result<RunRecord, Report> {
  let staging = manager.store().staging_dir()?;
  let mut archive = receive(body, &staging.path().join(PACKAGE_FILE))?;
  let entries = PackageEntries::of(&mut archive)?;
  let record = read_record(&mut archive, &entries)?;
  check_record(&record)?;
  let run_dir = staging.path().join(STAGED_RUN_DIR);
  let out_dir = run_dir.join(OUT_DIR);
  fs::create_dir_all(&out_dir).wrap_err_with(|| format!("When creating '{}'", out_dir.display()))?;
  let mut uploads = vec![];
  let mut found = BTreeSet::new();
  for (name, index) in &entries.files {
    let path = Path::new(name);
    if name == EVENTS_FILE || path.starts_with(OUT_DIR) {
      extract(&mut archive, *index, name, &run_dir.join(path))?;
    } else if let Some(input) = record.inputs.iter().find(|input| input.path == path) {
      if found.insert(input.path.clone()) {
        uploads.extend(extract_input(manager, &mut archive, *index, input, &run_dir)?);
      }
    }
  }
  if let Some(missing) = record.inputs.iter().find(|input| !found.contains(&input.path)) {
    return Err(invalid(format!(
      "the run package holds no input file `{}`, which its run.json lists",
      missing.path.display()
    )));
  }
  if let Some(missing) = record
    .output_files
    .iter()
    .find(|file| !out_dir.join(&file.path).is_file())
  {
    return Err(invalid(format!(
      "the run package holds no output file `{OUT_DIR}/{}`, which its run.json lists",
      missing.path.display()
    )));
  }
  manager.publish_import(record, &run_dir, uploads)
}

fn receive(body: &mut dyn Read, path: &Path) -> Result<ZipArchive<File>, Report> {
  let mut file = File::options()
    .read(true)
    .write(true)
    .create_new(true)
    .open(path)
    .wrap_err_with(|| format!("When creating '{}'", path.display()))?;
  io::copy(body, &mut file).wrap_err("When receiving the run package")?;
  file.rewind().wrap_err("When reading the run package")?;
  ZipArchive::new(file).map_err(|err| invalid(format!("the file is not a run package: it is not a zip archive: {err}")))
}

struct PackageEntries {
  files: Vec<(String, usize)>,
}

impl PackageEntries {
  fn of(archive: &mut ZipArchive<File>) -> Result<Self, Report> {
    let mut entries = vec![];
    for index in 0..archive.len() {
      let entry = archive
        .by_index_raw(index)
        .map_err(|err| invalid(format!("the run package cannot be read: {err}")))?;
      let name = entry
        .name()
        .map_err(|err| invalid(format!("the run package cannot be read: {err}")))?
        .into_owned();
      if !is_inside(Path::new(name.trim_end_matches('/'))) {
        return Err(outside_package(&name));
      }
      if entry.is_symlink() {
        return Err(invalid(format!(
          "the run package holds `{name}`, a symbolic link; a run package holds files and folders only"
        )));
      }
      if entry.is_file() {
        entries.push((name, index));
      }
    }
    let root = if entries.iter().any(|(name, _)| name == RUN_FILE) {
      None
    } else {
      let folders = entries
        .iter()
        .filter_map(|(name, _)| name.split_once('/').map(|(folder, _)| folder))
        .unique()
        .collect_vec();
      match folders.as_slice() {
        [folder] if entries.iter().any(|(name, _)| *name == format!("{folder}/{RUN_FILE}")) => {
          Some(format!("{folder}/"))
        },
        _ => {
          return Err(invalid(format!(
            "the file is not a run package: it holds no `{RUN_FILE}`"
          )));
        },
      }
    };
    let files = entries
      .iter()
      .filter_map(|(name, index)| match &root {
        Some(root) => name.strip_prefix(root.as_str()).map(|name| (name.to_owned(), *index)),
        None => Some((name.clone(), *index)),
      })
      .collect();
    Ok(Self { files })
  }

  fn index(&self, name: &str) -> Option<usize> {
    self
      .files
      .iter()
      .find_map(|(entry, index)| (entry == name).then_some(*index))
  }
}

fn read_record(archive: &mut ZipArchive<File>, entries: &PackageEntries) -> Result<RunRecord, Report> {
  let index = entries
    .index(RUN_FILE)
    .ok_or_else(|| invalid(format!("the file is not a run package: it holds no `{RUN_FILE}`")))?;
  let mut text = String::new();
  archive
    .by_index(index)
    .map_err(|err| invalid(format!("the run package cannot be read: {err}")))?
    .read_to_string(&mut text)
    .map_err(|err| invalid(format!("`{RUN_FILE}` of the run package cannot be read: {err}")))?;
  json_read_str::<RunRecord>(&text).map_err(|report| {
    invalid(format!(
      "`{RUN_FILE}` of the run package is not a run record: {report:#}"
    ))
  })
}

fn check_record(record: &RunRecord) -> Result<(), Report> {
  if !record.status.has_ended() {
    return Err(invalid(format!(
      "the run package holds run `{}`, which has not ended (status `{}`); a run package holds a run that has ended",
      record.id.as_str(),
      record.status
    )));
  }
  if !record.uploads.is_empty() {
    return Err(invalid(
      "the run package lists uploads, which only a run that has not started holds",
    ));
  }
  let command = record.config.command();
  let mut outside = vec![];
  map_path_values(command, &mut record.config.settings()?, |path| {
    if !is_inside(path) {
      outside.push(path.to_path_buf());
    }
    path.to_path_buf()
  })?;
  let paths = outside
    .into_iter()
    .chain(record.inputs.iter().map(|input| input.path.clone()))
    .chain(
      record
        .output_files
        .iter()
        .flat_map(|file| {
          [
            Some(file.path.clone()),
            file.alignment.as_ref().map(|alignment| alignment.path.clone()),
          ]
        })
        .flatten(),
    );
  for path in paths {
    if !is_inside(&path) {
      return Err(outside_package(&path.to_string_lossy()));
    }
  }
  if let Some(input) = record.inputs.iter().find(|input| !input.path.starts_with(INPUTS_DIR)) {
    return Err(invalid(format!(
      "run.json of the run package lists input `{}`, which is not in the package's `{INPUTS_DIR}/` folder",
      input.path.display()
    )));
  }
  Ok(())
}

fn extract_input(
  manager: &RunManager,
  archive: &mut ZipArchive<File>,
  index: usize,
  input: &RunInput,
  run_dir: &Path,
) -> Result<Option<StoredUpload>, Report> {
  let mismatch = |size: usize, sha256: &str| {
    invalid(format!(
      "input `{}` of the run package does not match its run.json: the file has {size} bytes and SHA-256 {sha256}, \
       run.json lists {} bytes and SHA-256 {}",
      input.path.display(),
      input.size,
      input.sha256
    ))
  };
  if let Some(key) = input.stored() {
    let mut entry = archive
      .by_index(index)
      .map_err(|err| invalid(format!("the run package cannot be read: {err}")))?;
    let upload = manager
      .input_store()
      .store(&key.name, &mut entry, usize::MAX, usize::MAX)
      .wrap_err_with(|| format!("When reading the input `{}` of the run package", input.path.display()))?;
    if (upload.size(), upload.sha256()) != (input.size, input.sha256.as_str()) {
      return Err(mismatch(upload.size(), upload.sha256()));
    }
    return Ok(Some(upload));
  }
  let target = run_dir.join(&input.path);
  extract(archive, index, &input.path.to_string_lossy(), &target)?;
  let (size, sha256) = file_sha256(&target)?;
  if (size, &sha256) != (input.size, &input.sha256) {
    return Err(mismatch(size, &sha256));
  }
  Ok(None)
}

fn extract(archive: &mut ZipArchive<File>, index: usize, name: &str, target: &Path) -> Result<(), Report> {
  if let Some(parent) = target.parent() {
    fs::create_dir_all(parent).wrap_err_with(|| format!("When creating '{}'", parent.display()))?;
  }
  let mut entry = archive
    .by_index(index)
    .map_err(|err| invalid(format!("the run package cannot be read: {err}")))?;
  let mut file = File::create(target).wrap_err_with(|| format!("When creating '{}'", target.display()))?;
  io::copy(&mut entry, &mut file)
    .map_err(|err| invalid(format!("the file `{name}` of the run package cannot be read: {err}")))?;
  Ok(())
}

fn is_inside(path: &Path) -> bool {
  let mut components = path.components().peekable();
  components.peek().is_some() && components.all(|component| matches!(component, Component::Normal(_)))
}

fn outside_package(path: &str) -> Report {
  invalid(format!(
    "the run package names `{path}`, a path outside the package; every path of a run package is relative and inside it"
  ))
}

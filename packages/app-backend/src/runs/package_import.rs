use crate::runs::errors::{PackageTooLarge, invalid};
use crate::runs::input_store::{StoredUpload, format_size, set_sha256};
use crate::runs::inputs::hex_digest;
use crate::runs::manager::RunManager;
use crate::runs::package_reader::{PackageFile, PackageReader, is_inside, outside_package};
use crate::runs::record::{InputSet, RunInput, RunRecord};
use crate::runs::run_paths::{INPUTS_DIR, OUT_DIR, map_path_values};
use crate::runs::store::EVENTS_FILE;
use deser::Deserialize;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::JsonSchema;
use sha2::{Digest, Sha256};
use std::fs::{self, File};
use std::io::{Seek, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tempfile::TempDir;
use treetime_utils::make_internal_report;

const PACKAGE_FILE: &str = "package.zip";

const STAGED_RUN_DIR: &str = "run";

const STAGED_SETS_DIR: &str = "sets";

pub struct PackageImport {
  runs: Arc<RunManager>,
  staging: TempDir,
  file: File,
  size: usize,
  limit: usize,
}

impl PackageImport {
  pub fn begin(runs: Arc<RunManager>, limit: usize) -> Result<Self, Report> {
    let staging = runs.store().staging_dir()?;
    let path = staging.path().join(PACKAGE_FILE);
    let file = File::options()
      .read(true)
      .write(true)
      .create_new(true)
      .open(&path)
      .wrap_err_with(|| format!("When creating '{}'", path.display()))?;
    Ok(Self {
      runs,
      staging,
      file,
      size: 0,
      limit,
    })
  }

  pub fn write(&mut self, chunk: &[u8]) -> Result<(), Report> {
    self.size = self.size.saturating_add(chunk.len());
    if self.size > self.limit {
      return Err(Report::new(PackageTooLarge {
        limit: format_size(self.limit),
      }));
    }
    self.file.write_all(chunk).wrap_err("When receiving the run package")
  }

  pub fn finish(self) -> Result<RunRecord, Report> {
    let Self {
      runs,
      staging,
      mut file,
      limit,
      ..
    } = self;
    file.rewind().wrap_err("When reading the run package")?;
    import_run(&runs, &staging, file, limit)
  }
}

pub fn import_package_file(runs: &RunManager, path: &Path, limit: usize) -> Result<RunRecord, Report> {
  let file = File::open(path).wrap_err_with(|| format!("When opening '{}'", path.display()))?;
  let size = file
    .metadata()
    .wrap_err_with(|| format!("When reading '{}'", path.display()))?
    .len();
  if size > u64::try_from(limit)? {
    return Err(Report::new(PackageTooLarge {
      limit: format_size(limit),
    }));
  }
  let staging = runs.store().staging_dir()?;
  import_run(runs, &staging, file, limit)
}

/// A run package to import from a file on disk.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct ImportFileRequest {
  /// Absolute path of the run package.
  pub path: PathBuf,
}

fn import_run(manager: &RunManager, staging: &TempDir, file: File, limit: usize) -> Result<RunRecord, Report> {
  let mut package = PackageReader::open(file)?;
  package.ensure_within(limit)?;
  let record = package.record()?;
  check_record(&record)?;
  let run_dir = staging.path().join(STAGED_RUN_DIR);
  let out_dir = run_dir.join(OUT_DIR);
  fs::create_dir_all(&out_dir).wrap_err_with(|| format!("When creating '{}'", out_dir.display()))?;
  let sets_dir = staging.path().join(STAGED_SETS_DIR);
  let mut uploads = vec![];
  let mut found = vec![];
  for file in package.files().to_vec() {
    let path = Path::new(file.name());
    if file.name() == EVENTS_FILE || path.starts_with(OUT_DIR) {
      extract(&mut package, &file, &run_dir.join(path))?;
    } else if let Some(input) = record.inputs.iter().find(|input| input.path == path)
      && !found.contains(&input.path)
    {
      found.push(input.path.clone());
      if let Some(upload) = extract_input(
        manager,
        &mut package,
        &file,
        input,
        &record.input_sets,
        &sets_dir,
        limit,
      )? {
        uploads.push(upload);
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
  let sets = record
    .input_sets
    .iter()
    .enumerate()
    .map(|(number, set)| (set.sha256.clone(), sets_dir.join(number.to_string())))
    .collect_vec();
  manager.publish_import(record, &run_dir, uploads, &sets)
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
    .chain(record.input_sets.iter().map(|set| set.folder.clone()))
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
  if let Some(set) = record
    .input_sets
    .iter()
    .find(|set| !set.folder.starts_with(INPUTS_DIR) || set.folder == Path::new(INPUTS_DIR))
  {
    return Err(invalid(format!(
      "run.json of the run package lists the file set `{}`, which is not a folder in the package's `{INPUTS_DIR}/` \
       folder",
      set.folder.display()
    )));
  }
  if let Some(input) = record
    .inputs
    .iter()
    .find(|input| input.stored().is_none() && set_of(&record.input_sets, &input.path).is_none())
  {
    return Err(invalid(format!(
      "run.json of the run package lists input `{}`, which is neither an uploaded file nor a file of a set in \
       `input_sets`",
      input.path.display()
    )));
  }
  record.input_sets.iter().try_for_each(|set| check_set(record, set))
}

fn extract_input(
  manager: &RunManager,
  package: &mut PackageReader<File>,
  file: &PackageFile,
  input: &RunInput,
  sets: &[InputSet],
  sets_dir: &Path,
  limit: usize,
) -> Result<Option<StoredUpload>, Report> {
  let (size, sha256, upload) = if let Some(key) = input.stored() {
    let mut writer = manager.input_store().begin(&key.name, limit, limit)?;
    package.copy(file, |chunk| writer.write(chunk))?;
    let upload = writer.finish()?;
    (upload.size(), upload.sha256().to_owned(), Some(upload))
  } else {
    let (number, rest) = set_of(sets, &input.path)
      .ok_or_else(|| make_internal_report!("input `{}` is in no file set of the package", input.path.display()))?;
    let target = sets_dir.join(number.to_string()).join(rest);
    let (size, sha256) = extract(package, file, &target)?;
    (size, sha256, None)
  };
  if (size, sha256.as_str()) != (input.size, input.sha256.as_str()) {
    return Err(invalid(format!(
      "input `{}` of the run package does not match its run.json: the file has {size} bytes and SHA-256 {sha256}, \
       run.json lists {} bytes and SHA-256 {}",
      input.path.display(),
      input.size,
      input.sha256
    )));
  }
  Ok(upload)
}

fn check_set(record: &RunRecord, set: &InputSet) -> Result<(), Report> {
  let files = record
    .inputs
    .iter()
    .filter_map(|input| {
      input
        .path
        .strip_prefix(&set.folder)
        .ok()
        .map(|rest| (rest, input.sha256.as_str()))
    })
    .unique()
    .collect_vec();
  if files.is_empty() {
    return Err(invalid(format!(
      "the file set `{}` of the run package holds no input that run.json lists",
      set.folder.display()
    )));
  }
  let sha256 = set_sha256(files)?;
  if sha256 != set.sha256 {
    return Err(invalid(format!(
      "the file set `{}` of the run package does not match its run.json: its files give the key {sha256}, run.json \
       lists {}",
      set.folder.display(),
      set.sha256
    )));
  }
  Ok(())
}

fn set_of<'a>(sets: &'a [InputSet], path: &'a Path) -> Option<(usize, &'a Path)> {
  sets
    .iter()
    .enumerate()
    .find_map(|(number, set)| path.strip_prefix(&set.folder).ok().map(|rest| (number, rest)))
}

fn extract(package: &mut PackageReader<File>, file: &PackageFile, target: &Path) -> Result<(usize, String), Report> {
  if let Some(parent) = target.parent() {
    fs::create_dir_all(parent).wrap_err_with(|| format!("When creating '{}'", parent.display()))?;
  }
  let mut output = File::create(target).wrap_err_with(|| format!("When creating '{}'", target.display()))?;
  let mut hasher = Sha256::new();
  let mut size = 0_usize;
  package.copy(file, |chunk| {
    hasher.update(chunk);
    size += chunk.len();
    output
      .write_all(chunk)
      .wrap_err_with(|| format!("When writing '{}'", target.display()))
  })?;
  Ok((size, hex_digest(&hasher.finalize())?))
}

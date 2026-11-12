use crate::command::{OutputAlignment, OutputFile, config_from_settings};
use crate::config::code::config_code;
use crate::job::JobId;
use crate::runs::errors::{conflict, invalid};
use crate::runs::files::{FileChunks, relative_name, walk_files};
use crate::runs::input_store::validate_upload_name;
use crate::runs::inputs::file_sha256;
use crate::runs::manager::RunManager;
use crate::runs::record::{RunInput, RunRecord};
use crate::runs::run_paths::{INPUTS_DIR, OUT_DIR, map_path_values, path_strings};
use crate::runs::store::{EVENTS_FILE, RUN_FILE, record_text};
use app_commands::command::AppCommand;
use app_commands::commands::ancestral::aa_node_data::{cds_names, translation_path};
use app_commands::config::properties::{PathRole, setting_properties};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use parking_lot::Mutex;
use serde_json::{Map, Value};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fs::File;
use std::io::{self, Seek, Write};
use std::iter;
use std::mem;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use treetime_utils::make_internal_report;
use zip::CompressionMethod;
use zip::ZipWriter;
use zip::write::{SimpleFileOptions, StreamWriter};

pub const PACKAGE_SCRIPT_FILE: &str = "command.sh";

pub fn export_run(manager: &RunManager, id: &JobId) -> Result<RunPackage, Report> {
  let record = manager.get(id)?;
  if !record.status.has_ended() {
    return Err(conflict(format!(
      "run `{}` has not ended (status `{}`); a run package holds a run that has ended",
      id.as_str(),
      record.status
    )));
  }
  let run_dir = manager.store().run_dir(id);
  let command = record.config.command();
  let mut settings = record.config.settings()?;
  let paths = PackagePaths::of(command, &settings, &record.inputs)?;
  let inputs = package_inputs(manager, &run_dir, &record.inputs, &paths)?;
  map_path_values(command, &mut settings, |path| paths.package(path))?;
  let code = config_code(command, &settings)?;
  let record = RunRecord {
    config: config_from_settings(command, &Value::Object(settings))?,
    inputs: record
      .inputs
      .iter()
      .map(|input| RunInput {
        path: paths.package(&input.path),
        ..input.clone()
      })
      .collect(),
    output_files: record
      .output_files
      .iter()
      .map(|file| OutputFile {
        alignment: file.alignment.as_ref().map(|alignment| OutputAlignment {
          number: alignment.number,
          path: paths.package(&alignment.path),
        }),
        ..file.clone()
      })
      .collect(),
    ..record
  };
  let events = run_dir.join(EVENTS_FILE);
  let out_dir = run_dir.join(OUT_DIR);
  let outputs = walk_files(&out_dir)?
    .into_iter()
    .map(|path| {
      let name = format!("{OUT_DIR}/{}", relative_name(path.strip_prefix(&out_dir)?));
      Ok(PackageEntry::file(name, path))
    })
    .collect::<Result<Vec<_>, Report>>()?;
  let entries = [
    PackageEntry::bytes(RUN_FILE, record_text(&record)?),
    PackageEntry::bytes(&code.config_file, code.yaml_text.clone()),
    PackageEntry::bytes(PACKAGE_SCRIPT_FILE, package_script(&code.config_command)),
  ]
  .into_iter()
  .chain(events.is_file().then(|| PackageEntry::file(EVENTS_FILE, events)))
  .chain(inputs)
  .chain(outputs)
  .collect();
  Ok(RunPackage {
    folder: package_folder(id),
    entries,
  })
}

pub fn package_folder(id: &JobId) -> String {
  format!("treetime-{}", id.as_str())
}

pub struct RunPackage {
  folder: String,
  entries: Vec<PackageEntry>,
}

impl RunPackage {
  pub fn folder(&self) -> &str {
    &self.folder
  }

  pub fn file_name(&self) -> String {
    format!("{}.zip", self.folder)
  }

  pub fn write(self, writer: impl Write + Seek) -> Result<(), Report> {
    let mut zip = ZipWriter::new(writer);
    let options = SimpleFileOptions::default().compression_method(CompressionMethod::Deflated);
    for entry in self.entries {
      let name = format!("{}/{}", self.folder, entry.name);
      zip
        .start_file(&name, options)
        .wrap_err_with(|| format!("When adding '{name}' to the run package"))?;
      match entry.source {
        EntrySource::Bytes(bytes) => zip.write_all(&bytes).wrap_err("When writing the run package")?,
        EntrySource::File(path) => {
          let mut file = File::open(&path).wrap_err_with(|| format!("When opening '{}'", path.display()))?;
          io::copy(&mut file, &mut zip).wrap_err_with(|| format!("When packing '{}'", path.display()))?;
        },
      }
    }
    zip.finish().wrap_err("When finishing the run package")?;
    Ok(())
  }

  pub fn stream(self, mut writer: impl Write) -> Result<(), Report> {
    for chunk in self.chunks() {
      writer.write_all(&chunk?).wrap_err("When writing the run package")?;
    }
    Ok(())
  }

  pub fn chunks(self) -> RunPackageChunks {
    let buffer = SharedBuffer::default();
    RunPackageChunks {
      zip: Some(ZipWriter::new_stream(buffer.clone())),
      buffer,
      folder: self.folder,
      pending: self.entries.into(),
      current: None,
    }
  }
}

pub struct RunPackageChunks {
  zip: Option<ZipWriter<StreamWriter<SharedBuffer>>>,
  buffer: SharedBuffer,
  folder: String,
  pending: VecDeque<PackageEntry>,
  current: Option<FileChunks>,
}

impl RunPackageChunks {
  fn advance(&mut self) -> Result<Option<Vec<u8>>, Report> {
    loop {
      let written = self.buffer.take();
      if !written.is_empty() {
        return Ok(Some(written));
      }
      let Some(zip) = self.zip.as_mut() else {
        return Ok(None);
      };
      if let Some(file) = self.current.as_mut() {
        match file.next() {
          Some(block) => zip.write_all(&block?).wrap_err("When writing the run package")?,
          None => self.current = None,
        }
      } else if let Some(entry) = self.pending.pop_front() {
        let name = format!("{}/{}", self.folder, entry.name);
        zip
          .start_file(
            &name,
            SimpleFileOptions::default().compression_method(CompressionMethod::Stored),
          )
          .wrap_err_with(|| format!("When adding '{name}' to the run package"))?;
        match entry.source {
          EntrySource::Bytes(bytes) => zip.write_all(&bytes).wrap_err("When writing the run package")?,
          EntrySource::File(path) => self.current = Some(FileChunks::open(&path)?),
        }
      } else if let Some(zip) = self.zip.take() {
        zip.finish().wrap_err("When finishing the run package")?;
      }
    }
  }
}

impl Iterator for RunPackageChunks {
  type Item = Result<Vec<u8>, Report>;

  fn next(&mut self) -> Option<Self::Item> {
    let chunk = self.advance().transpose();
    if chunk.as_ref().is_some_and(Result::is_err) {
      self.zip = None;
      self.current = None;
      self.buffer.take();
    }
    chunk
  }
}

struct PackageEntry {
  name: String,
  source: EntrySource,
}

impl PackageEntry {
  fn bytes(name: &str, text: String) -> Self {
    Self {
      name: name.to_owned(),
      source: EntrySource::Bytes(text.into_bytes()),
    }
  }

  fn file(name: impl Into<String>, path: PathBuf) -> Self {
    Self {
      name: name.into(),
      source: EntrySource::File(path),
    }
  }
}

enum EntrySource {
  Bytes(Vec<u8>),
  File(PathBuf),
}

#[derive(Clone, Default)]
struct SharedBuffer(Arc<Mutex<Vec<u8>>>);

impl SharedBuffer {
  fn take(&self) -> Vec<u8> {
    mem::take(&mut *self.0.lock())
  }
}

impl Write for SharedBuffer {
  fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
    self.0.lock().extend_from_slice(bytes);
    Ok(bytes.len())
  }

  fn flush(&mut self) -> io::Result<()> {
    Ok(())
  }
}

struct PackagePaths {
  moved: BTreeMap<PathBuf, PathBuf>,
}

impl PackagePaths {
  fn of(command: AppCommand, settings: &Map<String, Value>, inputs: &[RunInput]) -> Result<Self, Report> {
    let mut taken: BTreeSet<PathBuf> = inputs
      .iter()
      .filter(|input| input.path.is_relative())
      .map(|input| input.path.clone())
      .collect();
    let mut moved = BTreeMap::new();
    for property in setting_properties(command.config_schema().as_value())? {
      let Some(value) = settings.get(&property.key) else {
        continue;
      };
      match property.path_role {
        Some(PathRole::Input) => {
          for path in path_strings(value).into_iter().map(Path::new) {
            if path.is_absolute() && !moved.contains_key(path) {
              let packaged = free_input_path(&taken, path)?;
              taken.insert(packaged.clone());
              moved.insert(path.to_path_buf(), packaged);
            }
          }
        },
        Some(PathRole::InputTemplate) => {
          if let Value::String(template) = value
            && Path::new(template).is_absolute()
          {
            moved.extend(template_paths(&property.key, template, settings)?);
          }
        },
        Some(PathRole::Output) | None => {},
      }
    }
    Ok(Self { moved })
  }

  fn package(&self, path: &Path) -> PathBuf {
    self.moved.get(path).cloned().unwrap_or_else(|| path.to_path_buf())
  }
}

fn free_input_path(taken: &BTreeSet<PathBuf>, path: &Path) -> Result<PathBuf, Report> {
  let name = path
    .file_name()
    .and_then(|name| name.to_str())
    .ok_or_else(|| invalid(format!("input '{}' has no file name", path.display())))?;
  validate_upload_name(name).wrap_err_with(|| format!("When packing the input '{}'", path.display()))?;
  let inputs = Path::new(INPUTS_DIR);
  let plain = inputs.join(name);
  if !taken.contains(&plain) {
    return Ok(plain);
  }
  (1..=taken.len() + 1)
    .map(|number| inputs.join(number.to_string()).join(name))
    .find(|candidate| !taken.contains(candidate))
    .ok_or_else(|| make_internal_report!("no free path for the input '{}' in the run package", path.display()))
}

fn template_paths(
  setting: &str,
  template: &str,
  settings: &Map<String, Value>,
) -> Result<Vec<(PathBuf, PathBuf)>, Report> {
  let cdses = cds_names(settings)?;
  let expanded = cdses.iter().map(|cds| translation_path(template, cds)).collect_vec();
  let base = Path::new(template)
    .ancestors()
    .skip(1)
    .find(|dir| expanded.iter().all(|path| path.starts_with(dir)))
    .ok_or_else(|| make_internal_report!("the path template '{template}' has no folder"))?;
  let packaged = Path::new(INPUTS_DIR)
    .join(setting)
    .join(Path::new(template).strip_prefix(base)?);
  let packaged_template = packaged.to_string_lossy().into_owned();
  Ok(
    iter::once((PathBuf::from(template), packaged))
      .chain(
        cdses
          .iter()
          .zip(expanded)
          .map(|(cds, path)| (path, translation_path(&packaged_template, cds))),
      )
      .collect(),
  )
}

fn package_inputs(
  manager: &RunManager,
  run_dir: &Path,
  inputs: &[RunInput],
  paths: &PackagePaths,
) -> Result<Vec<PackageEntry>, Report> {
  inputs
    .iter()
    .unique_by(|input| paths.package(&input.path))
    .map(|input| {
      let source = if input.path.is_absolute() {
        let (size, sha256) = file_sha256(&input.path)?;
        if (size, &sha256) != (input.size, &input.sha256) {
          return Err(conflict(format!(
            "input '{}' changed since the run read it: it has {size} bytes and SHA-256 {sha256}, the run read {} bytes \
             with SHA-256 {}",
            input.path.display(),
            input.size,
            input.sha256
          )));
        }
        input.path.clone()
      } else {
        input
          .stored()
          .map_or_else(|| run_dir.join(&input.path), |key| manager.input_store().path(&key))
      };
      if !source.is_file() {
        return Err(invalid(format!(
          "input `{}` of the run is missing: no file '{}'",
          input.path.display(),
          source.display()
        )));
      }
      Ok(PackageEntry::file(relative_name(&paths.package(&input.path)), source))
    })
    .collect()
}

fn package_script(config_command: &str) -> String {
  format!(
    "#!/bin/sh\n# Reproduces the run with the TreeTime command-line tool. Run it in this folder: the outputs replace the \
     files in out/.\n{config_command}\n"
  )
}

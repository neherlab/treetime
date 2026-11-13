use crate::runs::errors::{PackageTooLarge, invalid};
use crate::runs::input_store::format_size;
use crate::runs::record::RunRecord;
use crate::runs::store::RUN_FILE;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::fmt::Display;
use std::io::{Read, Seek};
use std::path::{Component, Path};
use treetime_utils::io::json::json_read_str;
use zip::ZipArchive;

pub const RECORD_SIZE_LIMIT: usize = 8 << 20;

const COPY_BUFFER_SIZE: usize = 1 << 16;

pub struct PackageReader<R> {
  archive: ZipArchive<R>,
  files: Vec<PackageFile>,
}

impl<R: Read + Seek> PackageReader<R> {
  pub fn open(reader: R) -> Result<Self, Report> {
    let mut archive = ZipArchive::new(reader)
      .map_err(|err| invalid(format!("the file is not a run package: it is not a zip archive: {err}")))?;
    let mut entries = vec![];
    for index in 0..archive.len() {
      let entry = archive.by_index_raw(index).map_err(unreadable)?;
      let name = entry.name().map_err(unreadable)?.into_owned();
      if !is_inside(Path::new(name.trim_end_matches('/'))) {
        return Err(outside_package(&name));
      }
      if entry.is_symlink() {
        return Err(invalid(format!(
          "the run package holds `{name}`, a symbolic link; a run package holds files and folders only"
        )));
      }
      if entry.is_file() {
        entries.push(PackageFile {
          name,
          index,
          size: entry.size(),
        });
      }
    }
    let root = package_root(&entries)?;
    let files = entries
      .into_iter()
      .filter_map(|file| match &root {
        Some(root) => file.name.strip_prefix(root.as_str()).map(|name| PackageFile {
          name: name.to_owned(),
          ..file
        }),
        None => Some(file),
      })
      .collect();
    Ok(Self { archive, files })
  }

  pub fn files(&self) -> &[PackageFile] {
    &self.files
  }

  pub fn ensure_within(&self, limit: usize) -> Result<(), Report> {
    let total = self
      .files
      .iter()
      .try_fold(0_u64, |total, file| total.checked_add(file.size));
    let limit_bytes = u64::try_from(limit)?;
    if total.is_none_or(|total| total > limit_bytes) {
      return Err(Report::new(PackageTooLarge {
        limit: format_size(limit),
      }));
    }
    Ok(())
  }

  pub fn record(&mut self) -> Result<RunRecord, Report> {
    let file = self
      .files
      .iter()
      .find(|file| file.name == RUN_FILE)
      .cloned()
      .ok_or_else(|| invalid(format!("the file is not a run package: it holds no `{RUN_FILE}`")))?;
    if file.size > u64::try_from(RECORD_SIZE_LIMIT)? {
      return Err(invalid(format!(
        "`{RUN_FILE}` of the run package is larger than {}, the limit of a run record",
        format_size(RECORD_SIZE_LIMIT)
      )));
    }
    let mut bytes = vec![];
    self.copy(&file, |chunk| {
      bytes.extend_from_slice(chunk);
      Ok(())
    })?;
    let text = String::from_utf8(bytes)
      .map_err(|err| invalid(format!("`{RUN_FILE}` of the run package is not UTF-8 text: {err}")))?;
    json_read_str::<RunRecord>(&text).map_err(|report| {
      invalid(format!(
        "`{RUN_FILE}` of the run package is not a run record: {report:#}"
      ))
    })
  }

  pub fn copy(&mut self, file: &PackageFile, mut write: impl FnMut(&[u8]) -> Result<(), Report>) -> Result<(), Report> {
    let mut entry = self.archive.by_index(file.index).map_err(unreadable)?;
    let mut buffer = vec![0_u8; COPY_BUFFER_SIZE];
    loop {
      let read = entry.read(&mut buffer).map_err(|err| {
        invalid(format!(
          "the file `{}` of the run package cannot be read: {err}",
          file.name
        ))
      })?;
      if read == 0 {
        return Ok(());
      }
      write(&buffer[..read]).wrap_err_with(|| format!("When unpacking `{}` of the run package", file.name))?;
    }
  }
}

#[derive(Clone, Debug)]
pub struct PackageFile {
  name: String,
  index: usize,
  size: u64,
}

impl PackageFile {
  pub fn name(&self) -> &str {
    &self.name
  }
}

pub fn is_inside(path: &Path) -> bool {
  let mut components = path.components().peekable();
  components.peek().is_some() && components.all(|component| matches!(component, Component::Normal(_)))
}

pub fn outside_package(path: &str) -> Report {
  invalid(format!(
    "the run package names `{path}`, a path outside the package; every path of a run package is relative and inside it"
  ))
}

fn package_root(entries: &[PackageFile]) -> Result<Option<String>, Report> {
  if entries.iter().any(|file| file.name == RUN_FILE) {
    return Ok(None);
  }
  let folders = entries
    .iter()
    .filter_map(|file| file.name.split_once('/').map(|(folder, _)| folder))
    .unique()
    .collect_vec();
  match folders.as_slice() {
    [folder] if entries.iter().any(|file| file.name == format!("{folder}/{RUN_FILE}")) => {
      Ok(Some(format!("{folder}/")))
    },
    _ => Err(invalid(format!(
      "the file is not a run package: it holds no `{RUN_FILE}`"
    ))),
  }
}

fn unreadable(err: impl Display) -> Report {
  invalid(format!("the run package cannot be read: {err}"))
}

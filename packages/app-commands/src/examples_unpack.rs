use eyre::{Report, WrapErr};
use smart_default::SmartDefault;
use std::fs::{self, OpenOptions};
use std::io::{self, Read, Seek};
use std::path::{Path, PathBuf};
use treetime_utils::make_error;
use zip::ZipArchive;
use zip::read::ZipFile;

const WINDOWS_DEVICE_NAMES: [&str; 32] = [
  "CON", "PRN", "AUX", "NUL", "COM0", "COM1", "COM2", "COM3", "COM4", "COM5", "COM6", "COM7", "COM8", "COM9", "COM¹",
  "COM²", "COM³", "LPT0", "LPT1", "LPT2", "LPT3", "LPT4", "LPT5", "LPT6", "LPT7", "LPT8", "LPT9", "LPT¹", "LPT²",
  "LPT³", "CONIN$", "CONOUT$",
];

const WINDOWS_RESERVED_CHARACTERS: [char; 8] = ['<', '>', ':', '"', '\\', '|', '?', '*'];

const FILE_TYPE_MASK: u32 = 0o170_000;
const FILE_TYPE_DIRECTORY: u32 = 0o040_000;
const FILE_TYPE_REGULAR: u32 = 0o100_000;

pub(crate) fn unpack_examples<R: Read + Seek>(
  archive: &mut ZipArchive<R>,
  folder: &Path,
  limits: UnpackLimits,
) -> Result<(), Report> {
  archive
    .validate_header_consistency()
    .wrap_err("the example archive is invalid")?;
  let count = archive.len();
  if count > limits.entries {
    return make_error!(
      "the example archive has {count} entries, more than the limit of {}",
      limits.entries
    );
  }
  let mut remaining = limits.bytes;
  for index in 0..count {
    let mut entry = archive
      .by_index(index)
      .wrap_err_with(|| format!("When reading entry {index} of the example archive"))?;
    let shown = shown_name(&entry)?;
    match entry_kind(&entry, &shown)? {
      EntryKind::Folder(path) => {
        let path = folder.join(path);
        fs::create_dir_all(&path).wrap_err_with(|| format!("When unpacking entry '{shown}'"))?;
      },
      EntryKind::File(path) => {
        let written = write_file(&mut entry, &folder.join(path), remaining)
          .wrap_err_with(|| format!("When unpacking entry '{shown}'"))?;
        if written > remaining {
          return make_error!(
            "the example archive unpacks to more than the limit of {} bytes (at entry '{shown}')",
            limits.bytes
          );
        }
        remaining -= written;
      },
    }
  }
  Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, SmartDefault)]
pub(crate) struct UnpackLimits {
  #[default(10_000)]
  pub(crate) entries: usize,
  #[default(1 << 30)]
  pub(crate) bytes: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum EntryKind {
  Folder(PathBuf),
  File(PathBuf),
}

fn shown_name<R: Read>(entry: &ZipFile<'_, R>) -> Result<String, Report> {
  let raw = entry
    .name()
    .wrap_err("When reading the name of an entry of the example archive")?;
  Ok(
    raw
      .chars()
      .map(|c| {
        if matches!(c, '\\' | '\'' | '"') {
          c.to_string()
        } else {
          c.escape_debug().to_string()
        }
      })
      .collect(),
  )
}

fn entry_kind<R: Read>(entry: &ZipFile<'_, R>, shown: &str) -> Result<EntryKind, Report> {
  let raw = entry
    .name()
    .wrap_err("When reading the name of an entry of the example archive")?;
  let path = match entry_path(&raw) {
    Ok(path) => path,
    Err(problem) => return make_error!("the example archive is invalid: entry '{shown}' {problem}"),
  };
  if entry.enclosed_name().is_none() {
    return make_error!("the example archive is invalid: entry '{shown}' has a path outside the folder");
  }
  if entry.is_symlink() {
    return make_error!("the example archive is invalid: entry '{shown}' is a symbolic link");
  }
  let file_type = entry.unix_mode().map(|mode| mode & FILE_TYPE_MASK);
  if entry.is_dir() || file_type == Some(FILE_TYPE_DIRECTORY) {
    if entry.size() != 0 {
      return make_error!("the example archive is invalid: entry '{shown}' is a folder with content");
    }
    return Ok(EntryKind::Folder(path));
  }
  if !matches!(file_type, None | Some(0 | FILE_TYPE_REGULAR)) {
    return make_error!("the example archive is invalid: entry '{shown}' is not a regular file or a folder");
  }
  if path.as_os_str().is_empty() {
    return make_error!("the example archive is invalid: entry '{shown}' is a file without a name");
  }
  Ok(EntryKind::File(path))
}

fn entry_path(raw: &str) -> Result<PathBuf, String> {
  let escapes = raw
    .split('/')
    .try_fold(0_usize, |depth, component| match component {
      "" | "." => Some(depth),
      ".." => depth.checked_sub(1),
      _ => Some(depth + 1),
    })
    .is_none();
  if raw.starts_with(['/', '\\']) || escapes {
    return Err("has a path outside the folder".to_owned());
  }
  let names = raw
    .split('/')
    .filter(|component| !matches!(*component, "" | "." | ".."));
  for name in names {
    if let Some(c) = name
      .chars()
      .find(|c| c.is_control() || WINDOWS_RESERVED_CHARACTERS.contains(c))
    {
      return Err(format!(
        "has a path component with the character '{}'",
        c.escape_debug()
      ));
    }
    if name.ends_with(['.', ' ']) {
      return Err("has a path component that ends with '.' or a space".to_owned());
    }
    if is_windows_device_name(name) {
      return Err("has a path component that is a device name on Windows".to_owned());
    }
  }
  Ok(raw.split('/').fold(PathBuf::new(), |mut path, component| {
    match component {
      "" | "." => {},
      ".." => {
        path.pop();
      },
      _ => path.push(component),
    }
    path
  }))
}

fn is_windows_device_name(name: &str) -> bool {
  let base = name.split('.').next().unwrap_or(name).trim_end_matches(' ');
  WINDOWS_DEVICE_NAMES
    .iter()
    .any(|device| base.eq_ignore_ascii_case(device))
}

fn write_file(entry: &mut impl Read, path: &Path, remaining: u64) -> Result<u64, Report> {
  if let Some(parent) = path.parent() {
    fs::create_dir_all(parent).wrap_err_with(|| format!("When creating the folder '{}'", parent.display()))?;
  }
  let mut file = OpenOptions::new()
    .write(true)
    .create_new(true)
    .open(path)
    .wrap_err_with(|| format!("When creating the file '{}'", path.display()))?;
  io::copy(&mut entry.take(remaining.saturating_add(1)), &mut file)
    .wrap_err_with(|| format!("When writing the file '{}'", path.display()))
}

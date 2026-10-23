use crate::examples_unpack::{UnpackLimits, unpack_examples};
use crate::version::{BUILD_MODE, LONG_VERSION, NIGHTLY_BUILD_MODE};
use eyre::{Report, WrapErr};
use std::fs;
#[cfg(unix)]
use std::fs::Permissions;
use std::io::{ErrorKind, Seek, SeekFrom};
#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;
use std::path::Path;
use tempfile::Builder;
use treetime_utils::make_error;
use zip::ZipArchive;

#[cfg(not(target_family = "wasm"))]
use crate::examples_fetch::fetch;
#[cfg(target_family = "wasm")]
use std::fs::File;

pub const EXAMPLES_ASSET: &str = "examples.zip";

const RELEASES_URL: &str = "https://github.com/neherlab/treetime-nightly/releases";

const TEMPORARY_PREFIX: &str = ".treetime-examples-";

const FINDER_METADATA_FILE: &str = ".DS_Store";

pub fn examples_url() -> String {
  if BUILD_MODE == NIGHTLY_BUILD_MODE {
    format!(
      "{RELEASES_URL}/download/{}/{EXAMPLES_ASSET}",
      LONG_VERSION.replace('+', "%2B")
    )
  } else {
    format!("{RELEASES_URL}/latest/download/{EXAMPLES_ASSET}")
  }
}

pub fn download_examples(url: &str, target: &Path, mut progress: impl FnMut(DownloadProgress)) -> Result<(), Report> {
  ensure_empty_target(target)?;
  let parent = target
    .parent()
    .filter(|parent| !parent.as_os_str().is_empty())
    .unwrap_or_else(|| Path::new("."));
  fs::create_dir_all(parent).wrap_err_with(|| format!("When creating the folder '{}'", parent.display()))?;
  remove_stale_temporaries(parent)?;

  let mut archive = Builder::new()
    .prefix(TEMPORARY_PREFIX)
    .suffix(".zip")
    .tempfile_in(parent)
    .wrap_err_with(|| format!("When creating a temporary file in '{}'", parent.display()))?;
  fetch(url, archive.as_file_mut(), &mut progress).wrap_err_with(|| format!("When downloading '{url}'"))?;
  archive
    .as_file_mut()
    .seek(SeekFrom::Start(0))
    .wrap_err_with(|| format!("When rewinding the downloaded archive '{}'", archive.path().display()))?;

  let mut folder_builder = Builder::new();
  folder_builder.prefix(TEMPORARY_PREFIX);
  #[cfg(unix)]
  folder_builder.permissions(Permissions::from_mode(0o700));
  let unpacked = folder_builder
    .tempdir_in(parent)
    .wrap_err_with(|| format!("When creating a temporary folder in '{}'", parent.display()))?;
  ZipArchive::new(archive.as_file_mut())
    .map_err(Report::new)
    .and_then(|mut zip| unpack_examples(&mut zip, unpacked.path(), UnpackLimits::default()))
    .wrap_err_with(|| format!("When unpacking the example datasets from '{url}'"))?;

  remove_empty_target(target)?;
  fs::rename(unpacked.path(), target).wrap_err_with(|| {
    format!(
      "When moving the example datasets from '{}' to '{}'",
      unpacked.path().display(),
      target.display()
    )
  })?;
  let _moved = unpacked.keep();
  Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DownloadProgress {
  pub received: u64,
  pub total: Option<u64>,
}

pub fn ensure_empty_target(target: &Path) -> Result<(), Report> {
  let entries = match fs::read_dir(target) {
    Ok(entries) => entries,
    Err(error) if error.kind() == ErrorKind::NotFound => return Ok(()),
    Err(error) => {
      return Err(Report::new(error)).wrap_err_with(|| format!("When reading the folder '{}'", target.display()));
    },
  };
  for entry in entries {
    let entry = entry.wrap_err_with(|| format!("When reading the folder '{}'", target.display()))?;
    if entry.file_name() != FINDER_METADATA_FILE {
      return make_error!(
        "the folder '{}' is not empty; the example datasets go into an empty or missing folder",
        target.display()
      );
    }
  }
  Ok(())
}

fn remove_stale_temporaries(parent: &Path) -> Result<(), Report> {
  let entries = fs::read_dir(parent).wrap_err_with(|| format!("When reading the folder '{}'", parent.display()))?;
  for entry in entries {
    let entry = entry.wrap_err_with(|| format!("When reading the folder '{}'", parent.display()))?;
    if !entry.file_name().to_string_lossy().starts_with(TEMPORARY_PREFIX) {
      continue;
    }
    let path = entry.path();
    let file_type = entry
      .file_type()
      .wrap_err_with(|| format!("When reading the type of '{}'", path.display()))?;
    let removed = if file_type.is_dir() {
      fs::remove_dir_all(&path)
    } else {
      fs::remove_file(&path)
    };
    removed.wrap_err_with(|| format!("When removing the stale temporary '{}'", path.display()))?;
  }
  Ok(())
}

#[cfg(target_family = "wasm")]
fn fetch(_url: &str, _file: &mut File, _progress: &mut impl FnMut(DownloadProgress)) -> Result<(), Report> {
  make_error!("the example datasets download only in the CLI and in the desktop app")
}

fn remove_empty_target(target: &Path) -> Result<(), Report> {
  match fs::remove_file(target.join(FINDER_METADATA_FILE)) {
    Ok(()) => {},
    Err(error) if error.kind() == ErrorKind::NotFound => {},
    Err(error) => {
      return Err(Report::new(error)).wrap_err_with(|| format!("When cleaning the folder '{}'", target.display()));
    },
  }
  match fs::remove_dir(target) {
    Ok(()) => Ok(()),
    Err(error) if error.kind() == ErrorKind::NotFound => Ok(()),
    Err(error) => Err(Report::new(error)).wrap_err_with(|| {
      format!(
        "When replacing the folder '{}', which is no longer empty",
        target.display()
      )
    }),
  }
}

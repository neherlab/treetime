use crate::io::file::read_file_with;
use crate::make_report;
use eyre::{Report, WrapErr};
use std::ffi::OsStr;
use std::io::{ErrorKind, Read};
use std::path::{Path, PathBuf};
use std::{env, fs};

pub fn ensure_dir(filepath: impl AsRef<Path>) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  {
    let parent_dir = filepath
      .parent()
      .ok_or_else(|| make_report!("Unable to get parent path for '{}'", filepath.display()))?;

    let parent_path = absolute_path(parent_dir)?;

    fs::create_dir_all(&parent_path).wrap_err_with(|| format!("When creating directory '{}'", parent_path.display()))
  }
  .wrap_err_with(|| format!("When ensuring parent directory for '{}'", filepath.display()))
}

pub fn absolute_path(path: impl AsRef<Path>) -> Result<PathBuf, Report> {
  let path = path.as_ref();

  let absolute_path = if path.is_absolute() {
    path.to_path_buf()
  } else {
    env::current_dir()
      .wrap_err("When resolving the current directory")?
      .join(path)
  };

  Ok(absolute_path)
}

pub fn extension(filepath: impl AsRef<Path>) -> Option<String> {
  let filepath = filepath.as_ref();
  filepath.extension().and_then(OsStr::to_str).map(str::to_owned)
}

pub fn read_file_to_string(filepath: impl AsRef<Path>) -> Result<String, Report> {
  read_file_with(filepath, |mut reader| {
    let mut data = String::new();
    reader.read_to_string(&mut data).wrap_err("When reading UTF-8 text")?;
    Ok(data)
  })
}

pub fn read_file_to_string_if_exists(filepath: impl AsRef<Path>) -> Result<Option<String>, Report> {
  let filepath = filepath.as_ref();
  match fs::read_to_string(filepath) {
    Ok(data) => Ok(Some(data)),
    Err(error) if error.kind() == ErrorKind::NotFound => Ok(None),
    Err(error) => Err(Report::new(error)).wrap_err_with(|| format!("When reading file '{}'", filepath.display())),
  }
}

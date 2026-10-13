use eyre::{Report, WrapErr};
use std::fs::File;
use std::path::Path;
use tempfile::NamedTempFile;

pub fn write_atomically(destination: &Path, write: impl FnOnce(&mut File) -> Result<(), Report>) -> Result<(), Report> {
  let dir = destination
    .parent()
    .filter(|dir| !dir.as_os_str().is_empty())
    .unwrap_or_else(|| Path::new("."));
  let mut file = NamedTempFile::new_in(dir).wrap_err_with(|| format!("When creating a file in '{}'", dir.display()))?;
  write(file.as_file_mut())?;
  file
    .persist(destination)
    .map_err(|err| Report::new(err.error))
    .wrap_err_with(|| format!("When saving '{}'", destination.display()))?;
  Ok(())
}

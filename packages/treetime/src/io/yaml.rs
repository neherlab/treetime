use crate::io::file::{create_file_or_stdout, open_file_or_stdin};
use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};
use std::io::Write;
use std::path::Path;

pub fn yaml_read_file<T: for<'de> Deserialize<'de>, P: AsRef<Path>>(filepath: &Option<P>) -> Result<T, Report> {
  yaml_read(open_file_or_stdin(filepath)?)
    .wrap_err_with(|| format!("When reading YAML file: {:#?}", filepath.as_ref().map(AsRef::as_ref)))
}

pub fn yaml_read_str<T: for<'de> Deserialize<'de>>(s: impl AsRef<str>) -> Result<T, Report> {
  let result = serde_yaml::from_str(s.as_ref()).wrap_err("When reading YAML string")?;
  Ok(result)
}

pub fn yaml_read<T: for<'de> Deserialize<'de>>(reader: impl std::io::Read) -> Result<T, Report> {
  let result = serde_yaml::from_reader(reader).wrap_err("When reading YAML")?;
  Ok(result)
}

pub fn yaml_write_file<T: Serialize>(filepath: impl AsRef<Path>, obj: &T) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  yaml_write(create_file_or_stdout(filepath)?, &obj).wrap_err_with(|| format!("When writing YAML file: {filepath:#?}"))
}

pub fn yaml_write_str<T: Serialize>(obj: &T) -> Result<String, Report> {
  serde_yaml::to_string(obj).wrap_err("When writing YAML string")
}

pub fn yaml_write<W: Write, T: Serialize>(writer: W, obj: &T) -> Result<(), Report> {
  serde_yaml::to_writer(writer, &obj).wrap_err("When writing YAML")
}

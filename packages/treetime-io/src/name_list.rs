use eyre::{Report, WrapErr};
use std::io::BufRead;
use std::path::Path;
use treetime_utils::io::file::read_file_with;

pub fn name_list_read_file(filepath: impl AsRef<Path>) -> Result<Vec<String>, Report> {
  read_file_with(filepath, name_list_read)
}

pub fn name_list_read(reader: impl BufRead) -> Result<Vec<String>, Report> {
  let mut names = Vec::new();
  for (index, line) in reader.split(b'\n').enumerate() {
    let line = line.wrap_err("When reading the name list")?;
    let line = String::from_utf8(line).wrap_err_with(|| format!("When reading name {} of the list", index + 1))?;
    let name = line.trim();
    if !name.is_empty() {
      names.push(name.to_owned());
    }
  }
  Ok(names)
}

use eyre::{Report, WrapErr};
use std::io::BufRead;
use std::path::Path;
use treetime_utils::io::file::read_file_with;

pub fn name_list_read_file(filepath: impl AsRef<Path>, delimiter: u8) -> Result<Vec<String>, Report> {
  read_file_with(filepath, |reader| name_list_read(reader, delimiter))
}

pub fn name_list_read_str(s: &str, delimiter: u8) -> Result<Vec<String>, Report> {
  name_list_read(s.as_bytes(), delimiter)
}

pub fn name_list_read(reader: impl BufRead, delimiter: u8) -> Result<Vec<String>, Report> {
  let mut names = Vec::new();
  for (index, field) in reader.split(delimiter).enumerate() {
    let field = field.wrap_err("When reading the name list")?;
    let field = String::from_utf8(field).wrap_err_with(|| format!("When reading name {} of the list", index + 1))?;
    let name = field.trim();
    if !name.is_empty() {
      names.push(name.to_owned());
    }
  }
  Ok(names)
}

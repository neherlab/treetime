use crate::csv::{delimiter_from_path, detect_csv_delimiter, table_records};
use eyre::{Report, WrapErr};
use std::io::BufRead;
use std::path::Path;
use treetime_utils::io::file::read_file_with;

pub fn state_weights_read_file(filepath: impl AsRef<Path>, delimiters: &[char]) -> Result<Vec<(String, f64)>, Report> {
  let filepath = filepath.as_ref();
  read_file_with(filepath, |reader| {
    state_weights_read(reader, delimiter_from_path(filepath), delimiters)
  })
}

pub fn state_weights_read(
  mut reader: impl BufRead,
  path_delimiter: Option<u8>,
  delimiters: &[char],
) -> Result<Vec<(String, f64)>, Report> {
  let delimiter = detect_csv_delimiter(&mut reader, path_delimiter, delimiters, |headers| headers.len() >= 2)
    .wrap_err("When detecting the table delimiter")?;
  let mut reader = table_records(reader, delimiter);
  reader.read::<Vec<String>>()?;

  reader
    .iter::<Vec<String>>()
    .enumerate()
    .map(|(index, record)| state_weight(index, &record?))
    .filter_map(Result::transpose)
    .collect()
}

fn state_weight(index: usize, record: &[String]) -> Result<Option<(String, f64)>, Report> {
  let (Some(state), Some(weight)) = (record.first(), record.get(1)) else {
    return Ok(None);
  };
  if weight.is_empty() {
    return Ok(None);
  }
  let row = index + 2;
  let weight = weight
    .parse::<f64>()
    .map_err(Report::from)
    .wrap_err_with(|| format!("Row {row}: weight of state '{state}' is '{weight}', expected a number"))?;
  Ok((!weight.is_nan()).then(|| (state.clone(), weight)))
}

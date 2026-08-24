#[cfg(test)]
mod __tests__;

use eyre::Report;
use std::{
  io::{BufRead, Cursor},
  path::Path,
};
use treetime_utils::error::make_report;
use treetime_utils::io::file::open_file_or_stdin;

pub fn parse_delimited<R: BufRead>(reader: R, delimiter: u8) -> impl Iterator<Item = Result<String, Report>> {
  reader.split(delimiter).map(|chunk| {
    let bytes = chunk.map_err(|e| make_report!("Failed to read chunk: {e}"))?;
    String::from_utf8(bytes).map_err(|e| make_report!("Invalid UTF-8 in input: {e}"))
  })
}

pub fn parse_delimited_str(s: &str, delimiter: u8) -> impl Iterator<Item = Result<String, Report>> {
  parse_delimited(Cursor::new(s.as_bytes()), delimiter)
}

pub fn parse_delimited_file(
  file_path: impl AsRef<Path>,
  delimiter: u8,
) -> Result<impl Iterator<Item = Result<String, Report>>, Report> {
  Ok(parse_delimited(open_file_or_stdin(&Some(file_path))?, delimiter))
}

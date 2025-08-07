use crate::io::file::open_file_or_stdin;
use eyre::Report;
use std::{
  io::{BufRead, Cursor},
  path::Path,
};

pub fn parse_delimited<R: BufRead>(reader: R, delimiter: u8) -> impl Iterator<Item = String> {
  reader
    .split(delimiter)
    .map(|chunk| String::from_utf8(chunk.expect("Failed to read chunk")).expect("Invalid UTF-8 in input"))
}

pub fn parse_delimited_str(s: &str, delimiter: u8) -> impl Iterator<Item = String> {
  parse_delimited(Cursor::new(s.as_bytes()), delimiter)
}

pub fn parse_delimited_file(
  file_path: impl AsRef<Path>,
  delimiter: u8,
) -> Result<impl Iterator<Item = String>, Report> {
  Ok(parse_delimited(open_file_or_stdin(&Some(file_path))?, delimiter))
}

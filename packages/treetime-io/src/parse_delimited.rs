use eyre::Report;
use std::{
  io::{BufRead, Cursor},
  path::Path,
};
use treetime_utils::io::file::open_file_or_stdin;

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

#[cfg(test)]
mod tests {
  use super::*;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_parse_delimited_str_comma_delimiter() {
    let result = parse_delimited_str("A,B,C", b',').collect_vec();
    assert_eq!(result, vec!["A", "B", "C"]);
  }

  #[test]
  fn test_parse_delimited_str_semicolon_delimiter() {
    let result = parse_delimited_str("node1;node2;node3", b';').collect_vec();
    assert_eq!(result, vec!["node1", "node2", "node3"]);
  }

  #[test]
  fn test_parse_delimited_str_pipe_delimiter() {
    let result = parse_delimited_str("A|B|C", b'|').collect_vec();
    assert_eq!(result, vec!["A", "B", "C"]);
  }

  #[test]
  fn test_parse_delimited_str_newline_delimiter() {
    let result = parse_delimited_str("line1\nline2\nline3", b'\n').collect_vec();
    assert_eq!(result, vec!["line1", "line2", "line3"]);
  }

  #[test]
  fn test_parse_delimited_str_whitespace_in_names() {
    let result = parse_delimited_str("node with space,another node", b',').collect_vec();
    assert_eq!(result, vec!["node with space", "another node"]);
  }

  #[test]
  fn test_parse_delimited_str_empty_string_produces_empty_vec() {
    let result = parse_delimited_str("", b',').collect_vec();
    let expected: Vec<String> = vec![];
    assert_eq!(result, expected);
  }

  #[test]
  fn test_parse_delimited_str_single_element_no_delimiter() {
    let result = parse_delimited_str("single", b',').collect_vec();
    assert_eq!(result, vec!["single"]);
  }

  #[test]
  fn test_parse_delimited_str_empty_elements_between_delimiters() {
    let result = parse_delimited_str("A,,B", b',').collect_vec();
    assert_eq!(result, vec!["A", "", "B"]);
  }

  #[test]
  fn test_parse_delimited_str_trailing_delimiter() {
    // BufRead::split does not produce trailing empty element
    let result = parse_delimited_str("A,B,", b',').collect_vec();
    assert_eq!(result, vec!["A", "B"]);
  }

  #[test]
  fn test_parse_delimited_str_leading_delimiter() {
    let result = parse_delimited_str(",A,B", b',').collect_vec();
    assert_eq!(result, vec!["", "A", "B"]);
  }

  #[test]
  fn test_parse_delimited_str_only_delimiters() {
    // BufRead::split does not include trailing empty element
    let result = parse_delimited_str(",,", b',').collect_vec();
    assert_eq!(result, vec!["", ""]);
  }

  #[test]
  fn test_parse_delimited_str_tab_delimiter() {
    let result = parse_delimited_str("col1\tcol2\tcol3", b'\t').collect_vec();
    assert_eq!(result, vec!["col1", "col2", "col3"]);
  }

  #[test]
  fn test_parse_delimited_str_space_delimiter() {
    let result = parse_delimited_str("word1 word2 word3", b' ').collect_vec();
    assert_eq!(result, vec!["word1", "word2", "word3"]);
  }

  #[test]
  fn test_parse_delimited_str_unicode_content() {
    let result = parse_delimited_str("α,β,γ", b',').collect_vec();
    assert_eq!(result, vec!["α", "β", "γ"]);
  }

  #[test]
  fn test_parse_delimited_str_mixed_whitespace_content() {
    let result = parse_delimited_str("has\ttab,has space,normal", b',').collect_vec();
    assert_eq!(result, vec!["has\ttab", "has space", "normal"]);
  }
}

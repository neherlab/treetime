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

#[cfg(test)]
mod tests {
  use super::*;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_parse_delimited_str_comma_delimiter() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str("A,B,C", b',').try_collect()?;
    assert_eq!(result, vec!["A", "B", "C"]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_semicolon_delimiter() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str("node1;node2;node3", b';').try_collect()?;
    assert_eq!(result, vec!["node1", "node2", "node3"]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_pipe_delimiter() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str("A|B|C", b'|').try_collect()?;
    assert_eq!(result, vec!["A", "B", "C"]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_newline_delimiter() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str("line1\nline2\nline3", b'\n').try_collect()?;
    assert_eq!(result, vec!["line1", "line2", "line3"]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_whitespace_in_names() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str("node with space,another node", b',').try_collect()?;
    assert_eq!(result, vec!["node with space", "another node"]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_empty_string_produces_empty_vec() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str("", b',').try_collect()?;
    let expected: Vec<String> = vec![];
    assert_eq!(result, expected);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_single_element_no_delimiter() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str("single", b',').try_collect()?;
    assert_eq!(result, vec!["single"]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_empty_elements_between_delimiters() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str("A,,B", b',').try_collect()?;
    assert_eq!(result, vec!["A", "", "B"]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_trailing_delimiter() -> Result<(), Report> {
    // BufRead::split does not produce trailing empty element
    let result: Vec<String> = parse_delimited_str("A,B,", b',').try_collect()?;
    assert_eq!(result, vec!["A", "B"]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_leading_delimiter() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str(",A,B", b',').try_collect()?;
    assert_eq!(result, vec!["", "A", "B"]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_only_delimiters() -> Result<(), Report> {
    // BufRead::split does not include trailing empty element
    let result: Vec<String> = parse_delimited_str(",,", b',').try_collect()?;
    assert_eq!(result, vec!["", ""]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_tab_delimiter() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str("col1\tcol2\tcol3", b'\t').try_collect()?;
    assert_eq!(result, vec!["col1", "col2", "col3"]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_space_delimiter() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str("word1 word2 word3", b' ').try_collect()?;
    assert_eq!(result, vec!["word1", "word2", "word3"]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_unicode_content() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str("α,β,γ", b',').try_collect()?;
    assert_eq!(result, vec!["α", "β", "γ"]);
    Ok(())
  }

  #[test]
  fn test_parse_delimited_str_mixed_whitespace_content() -> Result<(), Report> {
    let result: Vec<String> = parse_delimited_str("has\ttab,has space,normal", b',').try_collect()?;
    assert_eq!(result, vec!["has\ttab", "has space", "normal"]);
    Ok(())
  }
}

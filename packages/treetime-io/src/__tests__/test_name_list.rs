#[cfg(test)]
mod tests {
  use crate::name_list::{name_list_read, name_list_read_str};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::comma(            "A,B,C",                    b',',  vec!["A", "B", "C"])]
  #[case::semicolon(        "node1;node2;node3",        b';',  vec!["node1", "node2", "node3"])]
  #[case::newline(          "line1\nline2\nline3",      b'\n', vec!["line1", "line2", "line3"])]
  #[case::crlf(             "line1\r\nline2\r\n",       b'\n', vec!["line1", "line2"])]
  #[case::inner_whitespace( "node with space,another",  b',',  vec!["node with space", "another"])]
  #[case::trimmed(          " A , B ",                  b',',  vec!["A", "B"])]
  #[case::trailing(         "A,B,",                     b',',  vec!["A", "B"])]
  #[case::leading(          ",A,B",                     b',',  vec!["A", "B"])]
  #[case::interior_empty(   "A,,B",                     b',',  vec!["A", "B"])]
  #[case::blank_field(      "A,  ,B",                   b',',  vec!["A", "B"])]
  #[case::only_delimiters(  ",,",                       b',',  vec![])]
  #[case::empty(            "",                         b',',  vec![])]
  #[case::other_delimiter(  "A,B",                      b';',  vec!["A,B"])]
  #[trace]
  fn test_name_list_read_str(
    #[case] input: &str,
    #[case] delimiter: u8,
    #[case] expected: Vec<&str>,
  ) -> Result<(), Report> {
    let actual = name_list_read_str(input, delimiter)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_name_list_read_rejects_invalid_utf8() {
    let result = name_list_read(&b"A,\xff,B"[..], b',');
    assert_error!(
      result,
      "When reading name 2 of the list: invalid utf-8 sequence of 1 bytes from index 0"
    );
  }
}

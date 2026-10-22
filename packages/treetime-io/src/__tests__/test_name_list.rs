#[cfg(test)]
mod tests {
  use crate::name_list::name_list_read;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::one_per_line(     "A\nB\nC",                  vec!["A", "B", "C"])]
  #[case::crlf(             "line1\r\nline2\r\n",       vec!["line1", "line2"])]
  #[case::inner_whitespace( "node with space\nanother", vec!["node with space", "another"])]
  #[case::trimmed(          " A \n\tB ",                vec!["A", "B"])]
  #[case::trailing_newline( "A\nB\n",                   vec!["A", "B"])]
  #[case::leading_newline(  "\nA\nB",                   vec!["A", "B"])]
  #[case::blank_lines(      "A\n\n  \nB",               vec!["A", "B"])]
  #[case::only_newlines(    "\n\n",                     vec![])]
  #[case::empty(            "",                         vec![])]
  #[case::comma_in_name(    "A,B\nC;D",                 vec!["A,B", "C;D"])]
  #[trace]
  fn test_name_list_read(#[case] input: &str, #[case] expected: Vec<&str>) -> Result<(), Report> {
    let actual = name_list_read(input.as_bytes())?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_name_list_read_rejects_invalid_utf8() {
    let result = name_list_read(&b"A\n\xff\nB"[..]);
    assert_error!(
      result,
      "When reading name 2 of the list: invalid utf-8 sequence of 1 bytes from index 0"
    );
  }
}

#[cfg(test)]
mod tests {
  use crate::csv::{default_name_candidates, detect_csv_delimiter, get_col_name};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::io::{BufReader, Cursor, Read};
  use treetime_utils::assert_error;
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::strain_only(
    vec![o!("strain"), o!("date")],
    0,
  )]
  #[case::accession_leftmost(
    vec![o!("accession"), o!("date"), o!("strain")],
    0,
  )]
  #[case::name_leftmost(
    vec![o!("name"), o!("accession"), o!("strain")],
    0,
  )]
  #[case::strain_leftmost(
    vec![o!("strain"), o!("accession"), o!("name")],
    0,
  )]
  #[case::accession_after_non_matching(
    vec![o!("region"), o!("country"), o!("accession"), o!("strain")],
    2,
  )]
  #[case::case_insensitive(
    vec![o!("region"), o!("ACCESSION"), o!("Strain")],
    1,
  )]
  #[trace]
  fn test_csv_get_col_name_header_position_priority(
    #[case] headers: Vec<String>,
    #[case] expected_idx: usize,
  ) -> Result<(), Report> {
    let candidates = default_name_candidates();
    let actual = get_col_name(&headers, &candidates, &None)?;
    assert_eq!(expected_idx, actual);
    Ok(())
  }

  #[test]
  fn test_csv_get_col_name_provided_name() -> Result<(), Report> {
    let headers = vec![o!("accession"), o!("date"), o!("strain")];
    let candidates = default_name_candidates();
    let actual = get_col_name(&headers, &candidates, &Some(o!("strain")))?;
    assert_eq!(2, actual);
    Ok(())
  }

  #[test]
  fn test_csv_get_col_name_no_match() {
    let headers = vec![o!("region"), o!("country"), o!("date")];
    let candidates = default_name_candidates();
    let result = get_col_name(&headers, &candidates, &None);
    assert_error!(
      result,
      "Unable to find column:\n  Looking for: strain, name, accession\n  Available columns are: region, country, date"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::csv(       "metadata.csv.xz", b"strain,date\nA,2020\n",       b',')]
  #[case::tsv(       "metadata.tsv.xz", b"strain\tdate\nA\t2020\n",     b'\t')]
  #[case::ssv(       "metadata.ssv.gz", b"strain;date\nA;2020\n",       b';')]
  #[case::content_wins("metadata.csv.xz", b"strain\tdate\nA\t2020\n", b'\t')]
  #[trace]
  fn test_csv_detect_delimiter_from_bounded_content(
    #[case] filepath: &str,
    #[case] content: &[u8],
    #[case] expected: u8,
  ) -> Result<(), Report> {
    let mut reader = BufReader::new(Cursor::new(content));
    let actual = detect_csv_delimiter(&mut reader, filepath, &[',', '\t', ';'], |headers| {
      headers == ["strain", "date"]
    })?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_csv_detect_delimiter_does_not_consume_reader() -> Result<(), Report> {
    let content = b"strain\tdate\nA\t2020\n";
    let mut reader = BufReader::new(Cursor::new(content));
    detect_csv_delimiter(&mut reader, "metadata.tsv.xz", &[',', '\t', ';'], |headers| {
      headers == ["strain", "date"]
    })?;

    let mut actual = Vec::new();
    reader.read_to_end(&mut actual)?;
    assert_eq!(content, actual.as_slice());
    Ok(())
  }

  #[test]
  fn test_csv_detect_delimiter_uses_only_explicit_candidate() -> Result<(), Report> {
    let mut reader = BufReader::new(Cursor::new(b"not,a,header"));
    let actual = detect_csv_delimiter(&mut reader, "metadata.tsv.xz", &['|'], |_| false)?;
    assert_eq!(b'|', actual);
    Ok(())
  }

  #[test]
  fn test_csv_detect_delimiter_uses_logical_extension_for_invalid_header() -> Result<(), Report> {
    let mut reader = BufReader::new(Cursor::new(b"unexpected header"));
    let actual = detect_csv_delimiter(&mut reader, "metadata.tsv.xz", &[',', '\t', ';'], |_| false)?;
    assert_eq!(b'\t', actual);
    Ok(())
  }

  #[test]
  fn test_csv_detect_delimiter_rejects_empty_candidates() {
    let mut reader = BufReader::new(Cursor::new(b"strain,date"));
    let result = detect_csv_delimiter(&mut reader, "metadata.csv", &[], |_| true);
    assert_error!(result, "At least one metadata delimiter is required");
  }

  #[test]
  fn test_csv_detect_delimiter_rejects_multibyte_candidate() {
    let mut reader = BufReader::new(Cursor::new(b"strain,date"));
    let result = detect_csv_delimiter(&mut reader, "metadata.csv", &['‣'], |_| true);
    assert_error!(
      result,
      "Metadata delimiter '‣' must fit in one byte: out of range integral type conversion attempted"
    );
  }
}

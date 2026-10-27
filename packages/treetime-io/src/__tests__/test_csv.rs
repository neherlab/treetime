#[cfg(test)]
mod tests {
  use crate::csv::{
    ColumnRule, CsvWriter, TableFormat, csv_read, default_name_candidates, delimiter_from_path, detect_csv_delimiter,
    find_column,
  };
  use eyre::{Report, WrapErr};
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
  fn test_csv_find_column_header_position_priority(
    #[case] headers: Vec<String>,
    #[case] expected_idx: usize,
  ) -> Result<(), Report> {
    let candidates = default_name_candidates();
    let actual = find_column(&headers, ColumnRule::OneOf(&candidates))?;
    assert_eq!(expected_idx, actual);
    Ok(())
  }

  #[test]
  fn test_csv_find_column_provided_name() -> Result<(), Report> {
    let headers = vec![o!("accession"), o!("date"), o!("strain")];
    let candidates = default_name_candidates();
    let actual = find_column(
      &headers,
      ColumnRule::requested_or(Some("strain"), ColumnRule::OneOf(&candidates)),
    )?;
    assert_eq!(2, actual);
    Ok(())
  }

  #[test]
  fn test_csv_find_column_no_match() {
    let headers = vec![o!("region"), o!("country"), o!("date")];
    let candidates = default_name_candidates();
    let result = find_column(&headers, ColumnRule::OneOf(&candidates));
    assert_error!(
      result,
      "Unable to find column:\n  Looking for: strain, name, accession\n  Available columns are: region, country, date"
    );
  }

  #[test]
  fn test_csv_find_column_hash_name_as_candidate() -> Result<(), Report> {
    let headers = vec![o!("name"), o!("date"), o!("country")];
    let actual = find_column(&headers, ColumnRule::OneOf(&[o!("#name")]))?;
    assert_eq!(0, actual);
    Ok(())
  }

  #[test]
  fn test_csv_find_column_hash_name_as_provided_name() -> Result<(), Report> {
    let headers = vec![o!("country"), o!("name"), o!("date")];
    let actual = find_column(&headers, ColumnRule::Exact("#name"))?;
    assert_eq!(1, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::written_match_wins(vec![o!("name"), o!("date"), o!("#date")], "#date",   2)]
  #[case::normalized_match(  vec![o!("name"), o!("#date")],              "date",    1)]
  #[case::padded_request(    vec![o!("#name"), o!("date")],              " #name ", 0)]
  #[case::hash_and_space(    vec![o!("# name"), o!("date")],             "name",    0)]
  #[trace]
  fn test_csv_find_column_provided_name_as_written_then_normalized(
    #[case] headers: Vec<String>,
    #[case] name: &str,
    #[case] expected_idx: usize,
  ) -> Result<(), Report> {
    let actual = find_column(&headers, ColumnRule::Exact(name))?;
    assert_eq!(expected_idx, actual);
    Ok(())
  }

  #[test]
  fn test_csv_find_column_padded_candidate() -> Result<(), Report> {
    let headers = vec![o!("country"), o!("#name")];
    let actual = find_column(&headers, ColumnRule::OneOf(&[o!(" #name ")]))?;
    assert_eq!(1, actual);
    Ok(())
  }

  #[test]
  fn test_csv_find_column_provided_name_ambiguous_after_normalization() {
    let headers = vec![o!("name"), o!("#date"), o!("date#")];
    let result = find_column(&headers, ColumnRule::Exact("date"));
    assert_error!(
      result,
      "Column 'date' is ambiguous: the columns '#date', 'date#' all read as 'date'. Request one of them as written"
    );
  }

  #[test]
  fn test_csv_find_column_provided_name_is_case_sensitive() {
    let headers = vec![o!("strain"), o!("date")];
    let result = find_column(&headers, ColumnRule::Exact("Date"));
    assert_error!(
      result,
      "Unable to find column 'Date'. Available columns are: strain, date"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::date(           vec![o!("strain"), o!("date")],                              1)]
  #[case::num_date(       vec![o!("strain"), o!("num_date")],                          1)]
  #[case::collection_date(vec![o!("strain"), o!("country"), o!("collection_date")],    2)]
  #[case::date_submitted( vec![o!("strain"), o!("date_submitted")],                    1)]
  #[case::capitalized(    vec![o!("strain"), o!("Date")],                              1)]
  #[case::leftmost_wins(  vec![o!("strain"), o!("date_submitted"), o!("date")],        1)]
  #[trace]
  fn test_csv_find_column_containing(#[case] headers: Vec<String>, #[case] expected_idx: usize) -> Result<(), Report> {
    let actual = find_column(&headers, ColumnRule::Containing("date"))?;
    assert_eq!(expected_idx, actual);
    Ok(())
  }

  #[test]
  fn test_csv_find_column_containing_no_match() {
    let headers = vec![o!("strain"), o!("country")];
    let result = find_column(&headers, ColumnRule::Containing("date"));
    assert_error!(
      result,
      "Unable to find column:\n  Looking for: a name that contains 'date'\n  Available columns are: strain, country"
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
    let actual = detect_csv_delimiter(&mut reader, delimiter_from_path(filepath), &[',', '\t', ';'], |headers| {
      headers == ["strain", "date"]
    })?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_csv_detect_delimiter_does_not_consume_reader() -> Result<(), Report> {
    let content = b"strain\tdate\nA\t2020\n";
    let mut reader = BufReader::new(Cursor::new(content));
    detect_csv_delimiter(
      &mut reader,
      delimiter_from_path("metadata.tsv.xz"),
      &[',', '\t', ';'],
      |headers| headers == ["strain", "date"],
    )?;

    let mut actual = Vec::new();
    reader
      .read_to_end(&mut actual)
      .wrap_err("When reading the buffered input")?;
    assert_eq!(content, actual.as_slice());
    Ok(())
  }

  #[test]
  fn test_csv_detect_delimiter_uses_only_explicit_candidate() -> Result<(), Report> {
    let mut reader = BufReader::new(Cursor::new(b"not,a,header"));
    let actual = detect_csv_delimiter(&mut reader, delimiter_from_path("metadata.tsv.xz"), &['|'], |_| false)?;
    assert_eq!(b'|', actual);
    Ok(())
  }

  #[test]
  fn test_csv_detect_delimiter_uses_logical_extension_for_invalid_header() -> Result<(), Report> {
    let mut reader = BufReader::new(Cursor::new(b"unexpected header"));
    let actual = detect_csv_delimiter(
      &mut reader,
      delimiter_from_path("metadata.tsv.xz"),
      &[',', '\t', ';'],
      |_| false,
    )?;
    assert_eq!(b'\t', actual);
    Ok(())
  }

  #[test]
  fn test_csv_detect_delimiter_rejects_empty_candidates() {
    let mut reader = BufReader::new(Cursor::new(b"strain,date"));
    let result = detect_csv_delimiter(&mut reader, delimiter_from_path("metadata.csv"), &[], |_| true);
    assert_error!(result, "At least one metadata delimiter is required");
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::beyond_one_byte('‣', "Metadata delimiter '‣' must be an ASCII character")]
  #[case::latin1(         'é', "Metadata delimiter 'é' must be an ASCII character")]
  #[trace]
  fn test_csv_detect_delimiter_rejects_non_ascii_candidate(#[case] delimiter: char, #[case] expected: &str) {
    let mut reader = BufReader::new(Cursor::new(b"strain,date"));
    let result = detect_csv_delimiter(&mut reader, delimiter_from_path("metadata.csv"), &[delimiter], |_| true);
    assert_error!(result, expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::csv(TableFormat::Csv, "name,value\n\"A,1\",0.5\nB,2.0\n")]
  #[case::tsv(TableFormat::Tsv, "name\tvalue\nA,1\t0.5\nB\t2.0\n")]
  #[trace]
  fn test_csv_writer_round_trips_rows(#[case] format: TableFormat, #[case] expected: &str) -> Result<(), Report> {
    let rows = vec![
      helpers::Row { name: o!("A,1"), value: 0.5 },
      helpers::Row { name: o!("B"), value: 2.0 },
    ];

    let mut buf = Vec::new();
    let mut csv = CsvWriter::new(&mut buf, format);
    rows.iter().try_for_each(|row| csv.write_row(row))?;
    csv.into_inner()?;

    assert_eq!(expected, String::from_utf8(buf.clone())?);
    assert_eq!(rows, csv_read::<helpers::Row>(buf.as_slice(), format)?);
    Ok(())
  }

  #[test]
  fn test_csv_writer_round_trips_missing_values_as_empty_fields() -> Result<(), Report> {
    let rows = vec![
      helpers::OptionalRow {
        name: o!("A"),
        label: None,
        score: None,
      },
      helpers::OptionalRow {
        name: o!("B"),
        label: Some(o!("x")),
        score: Some(1.5),
      },
      helpers::OptionalRow {
        name: o!(""),
        label: None,
        score: Some(2.0),
      },
    ];

    let mut buf = Vec::new();
    let mut csv = CsvWriter::new(&mut buf, TableFormat::Csv);
    rows.iter().try_for_each(|row| csv.write_row(row))?;
    csv.into_inner()?;

    assert_eq!(
      "name,label,score\nA,,\nB,x,1.5\n,,2.0\n",
      String::from_utf8(buf.clone())?
    );
    assert_eq!(
      rows,
      csv_read::<helpers::OptionalRow>(buf.as_slice(), TableFormat::Csv)?
    );
    Ok(())
  }

  #[test]
  fn test_csv_writer_writes_records_of_any_width() -> Result<(), Report> {
    let mut buf = Vec::new();
    let mut csv = CsvWriter::new(&mut buf, TableFormat::Csv);
    csv.write_record(["node", "Congo, DR", "usa"])?;
    csv.write_record(["A", "0.25", "0.75"])?;
    csv.into_inner()?;

    assert_eq!("node,\"Congo, DR\",usa\nA,0.25,0.75\n", String::from_utf8(buf)?);
    Ok(())
  }

  mod helpers {
    use deser::adapters::SkipBlank;
    use deser::{Deserialize, Serialize};

    #[derive(Debug, PartialEq, Serialize, Deserialize)]
    pub(super) struct Row {
      pub(super) name: String,
      pub(super) value: f64,
    }

    #[derive(Debug, PartialEq, Serialize, Deserialize)]
    pub(super) struct OptionalRow {
      pub(super) name: String,
      #[deser(as = SkipBlank<Option<_>>)]
      pub(super) label: Option<String>,
      pub(super) score: Option<f64>,
    }
  }
}

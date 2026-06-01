#[cfg(test)]
mod tests {
  use crate::csv::default_name_candidates;
  use crate::discrete_states_csv::*;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::tsv_delimiter(b'\t', "name\tlocation\nA\tusa\nB\teurope")]
  #[case::csv_delimiter(b',',  "name,location\nA,usa\nB,europe")]
  #[trace]
  fn test_discrete_states_csv_delimiter_handling(
    #[case] delimiter: u8,
    #[case] content: &str,
  ) -> Result<(), Report> {
    let (values, attr_name) =
      read_discrete_attrs_from_str(content, delimiter, &default_name_candidates(), &None, &Some(o!("location")), |s| Ok(s.to_owned()))?;

    let expected = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("europe"),
    };

    assert_eq!(expected, values);
    assert_eq!("location", attr_name);

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::default_name_column(None,                    "name\tlocation\nA\tusa\nB\teurope")]
  #[case::strain_column(      Some(o!("strain")),      "strain\tlocation\nA\tusa\nB\teurope")]
  #[case::accession_column(   Some(o!("accession")),   "accession\tlocation\nA\tusa\nB\teurope")]
  #[case::custom_column(      Some(o!("sample_id")),   "sample_id\tlocation\nA\tusa\nB\teurope")]
  #[trace]
  fn test_discrete_states_csv_name_column_detection(
    #[case] name_column: Option<String>,
    #[case] content: &str,
  ) -> Result<(), Report> {
    let (values, _) =
      read_discrete_attrs_from_str(content, b'\t', &default_name_candidates(), &name_column, &Some(o!("location")), |s| Ok(s.to_owned()))?;

    let expected = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("europe"),
    };

    assert_eq!(expected, values);

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::explicit_value_column(   "location", "name\tlocation\tregion\nA\tusa\tnorth_america")]
  #[case::column_among_multiple(   "region",   "name\tlocation\tregion\nA\tusa\tnorth_america")]
  #[trace]
  fn test_discrete_states_csv_value_column_selection(
    #[case] value_col_name: &str,
    #[case] content: &str,
  ) -> Result<(), Report> {
    let (values, attr_name) =
      read_discrete_attrs_from_str(content, b'\t', &default_name_candidates(), &None, &Some(o!(value_col_name)), |s| Ok(s.to_owned()))?;

    assert_eq!(attr_name, value_col_name);
    assert!(values.contains_key(&o!("A")));

    Ok(())
  }

  #[test]
  fn test_discrete_states_csv_header_normalization() -> Result<(), Report> {
    let content = "#name#\t#location#\nA\tusa\nB\teurope";

    let (values, attr_name) =
      read_discrete_attrs_from_str(content, b'\t', &[], &Some(o!("name")), &Some(o!("location")), |s| {
        Ok(s.to_owned())
      })?;

    let expected = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("europe"),
    };

    assert_eq!(expected, values);
    assert_eq!("location", attr_name);

    Ok(())
  }

  #[test]
  fn test_discrete_states_csv_custom_parser() -> Result<(), Report> {
    let content = "name\tweight\nA\t1.5\nB\t2.0\nC\t0.5";

    let (values, attr_name) = read_discrete_attrs_from_str(
      content,
      b'\t',
      &default_name_candidates(),
      &None,
      &Some(o!("weight")),
      |s| Ok(s.parse::<f64>()?),
    )?;

    let expected = btreemap! {
      o!("A") => 1.5,
      o!("B") => 2.0,
      o!("C") => 0.5,
    };

    assert_eq!(expected, values);
    assert_eq!("weight", attr_name);

    Ok(())
  }

  #[test]
  fn test_discrete_states_csv_whitespace_trimming() -> Result<(), Report> {
    let content = "name\tlocation\n  A  \t  usa  \n  B  \t  europe  ";

    let (values, _) = read_discrete_attrs_from_str(
      content,
      b'\t',
      &default_name_candidates(),
      &None,
      &Some(o!("location")),
      |s| Ok(s.to_owned()),
    )?;

    let expected = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("europe"),
    };

    assert_eq!(expected, values);

    Ok(())
  }
}

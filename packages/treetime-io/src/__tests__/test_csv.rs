#[cfg(test)]
mod tests {
  use crate::csv::{default_name_candidates, get_col_name};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
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
}

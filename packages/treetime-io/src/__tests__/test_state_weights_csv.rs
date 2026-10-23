#[cfg(test)]
mod tests {
  use crate::csv::default_metadata_delimiters;
  use crate::state_weights_csv::state_weights_read;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::{assert_error, o};

  #[rustfmt::skip]
  #[rstest]
  #[case::header_names_any_columns( "state,weight\nusa,0.5\nuk,2\n")]
  #[case::header_with_hash_marks(   "#country,weight\nusa,0.5\nuk,2\n")]
  #[case::tab_delimited(            "location\tpopulation\nusa\t0.5\nuk\t2\n")]
  #[case::extra_columns_ignored(    "state,weight,note\nusa,0.5,a\nuk,2,b\n")]
  #[trace]
  fn test_state_weights_csv_reads_state_and_weight_by_position(#[case] content: &str) -> Result<(), Report> {
    let expected = vec![(o!("usa"), 0.5), (o!("uk"), 2.0)];

    let actual = state_weights_read(content.as_bytes(), None, &default_metadata_delimiters())?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_state_weights_csv_drops_rows_without_a_weight() -> Result<(), Report> {
    let content = indoc! {"
      state,weight
      usa,0.5
      fiji,
      tonga,NaN
      uk,2
    "};
    let expected = vec![(o!("usa"), 0.5), (o!("uk"), 2.0)];

    let actual = state_weights_read(content.as_bytes(), None, &default_metadata_delimiters())?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_state_weights_csv_rejects_a_weight_that_is_not_a_number() {
    let content = "state,weight\nusa,0.5\nuk,many\n";

    assert_error!(
      state_weights_read(content.as_bytes(), None, &default_metadata_delimiters()),
      "Row 3: weight of state 'uk' is 'many', expected a number: invalid float literal"
    );
  }
}

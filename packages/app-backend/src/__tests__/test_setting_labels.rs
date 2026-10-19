#[cfg(test)]
mod tests {
  use crate::config::labels::setting_label;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::words_of_the_key(  "clock_rate",            "Clock rate")]
  #[case::fixed_label(       "clock_std_dev",         "Clock rate std. dev.")]
  #[case::abbreviation(      "output_tree_nwk",       "Output tree Newick")]
  #[case::leading_acronym(   "gtr_iterations",        "GTR iterations")]
  #[trace]
  fn test_setting_label(#[case] key: &str, #[case] expected: &str) {
    assert_eq!(expected, setting_label(key));
  }
}

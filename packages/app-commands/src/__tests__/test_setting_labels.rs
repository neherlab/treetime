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
  #[case::nested_key(        "branch_split.n_points", "Branch split: Number of points")]
  #[trace]
  fn test_setting_label(#[case] key: &str, #[case] expected: &str) {
    let key_path = key.split('.').map(str::to_owned).collect::<Vec<_>>();
    assert_eq!(expected, setting_label(&key_path));
  }
}

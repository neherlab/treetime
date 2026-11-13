#[cfg(test)]
mod tests {
  use crate::config::labels::{option_label, setting_label};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::words_of_the_key(  "clock_rate",            "Clock rate")]
  #[case::fixed_label(       "clock_std_dev",         "Clock rate std. dev.")]
  #[case::abbreviation(      "output_tree_nwk",       "Output tree Newick")]
  #[case::leading_acronym(   "gtr_model",             "GTR model")]
  #[trace]
  fn test_setting_label(#[case] key: &str, #[case] expected: &str) {
    assert_eq!(expected, setting_label(key));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::words_of_the_value(  "least-squares",            "Least squares")]
  #[case::fixed_label(         "min-dev",                  "Minimum deviation")]
  #[case::model_name(          "hky85",                    "HKY85")]
  #[case::abbreviation(        "mat-pb",                   "MAT protobuf")]
  #[case::reversed(            "descendant-count-reverse", "Descendant count, reversed")]
  #[case::reversed_fixed(      "label-reverse",            "Label, reversed")]
  #[case::parts(               "enewick,plain",            "Extended Newick, plain")]
  #[case::parts_keep_acronyms( "classic,beast",            "Classic Newick, BEAST")]
  #[case::parts_keep_mixed(    "rich,mrbayes",             "Rich Newick, MrBayes")]
  #[trace]
  fn test_option_label(#[case] value: &str, #[case] expected: &str) {
    assert_eq!(expected, option_label(value));
  }
}

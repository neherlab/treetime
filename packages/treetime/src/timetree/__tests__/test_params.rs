#[cfg(test)]
mod tests {
  use crate::timetree::params::{TimeMarginalMode, compute_effective_time_marginal};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::never_no_confidence(            (TimeMarginalMode::Never,     false, false, None),      TimeMarginalMode::Never)]
  #[case::never_confidence_covariation(   (TimeMarginalMode::Never,     true,  true,  None),      TimeMarginalMode::OnlyFinal)]
  #[case::never_confidence_clock_std(     (TimeMarginalMode::Never,     true,  false, Some(0.1)), TimeMarginalMode::OnlyFinal)]
  #[case::never_confidence_no_prereqs(    (TimeMarginalMode::Never,     true,  false, None),      TimeMarginalMode::Never)]
  #[case::always_no_confidence(           (TimeMarginalMode::Always,    false, false, None),      TimeMarginalMode::Always)]
  #[case::always_confidence_covariation(  (TimeMarginalMode::Always,    true,  true,  None),      TimeMarginalMode::Always)]
  #[case::always_confidence_clock_std(    (TimeMarginalMode::Always,    true,  false, Some(0.1)), TimeMarginalMode::Always)]
  #[case::only_final_no_confidence(       (TimeMarginalMode::OnlyFinal, false, false, None),      TimeMarginalMode::OnlyFinal)]
  #[case::only_final_confidence(          (TimeMarginalMode::OnlyFinal, true,  true,  None),      TimeMarginalMode::OnlyFinal)]
  #[trace]
  fn test_timetree_compute_effective_time_marginal(
    #[case] (mode, confidence, covariation, clock_std_dev): (TimeMarginalMode, bool, bool, Option<f64>),
    #[case] expected: TimeMarginalMode,
  ) {
    assert_eq!(expected, compute_effective_time_marginal(mode, confidence, clock_std_dev, covariation));
  }
}

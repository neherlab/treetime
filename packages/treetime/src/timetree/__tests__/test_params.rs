#[cfg(test)]
mod tests {
  use crate::progress::NoopProgress;
  use crate::timetree::params::{TimeMarginalMode, compute_effective_time_marginal};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::never_no_confidence(       (TimeMarginalMode::Never,     false), TimeMarginalMode::Never)]
  #[case::never_confidence(          (TimeMarginalMode::Never,     true),  TimeMarginalMode::OnlyFinal)]
  #[case::always_no_confidence(      (TimeMarginalMode::Always,    false), TimeMarginalMode::Always)]
  #[case::always_confidence(         (TimeMarginalMode::Always,    true),  TimeMarginalMode::Always)]
  #[case::only_final_no_confidence(  (TimeMarginalMode::OnlyFinal, false), TimeMarginalMode::OnlyFinal)]
  #[case::only_final_confidence(     (TimeMarginalMode::OnlyFinal, true),  TimeMarginalMode::OnlyFinal)]
  #[trace]
  fn test_timetree_compute_effective_time_marginal(
    #[case] (mode, confidence): (TimeMarginalMode, bool),
    #[case] expected: TimeMarginalMode,
  ) {
    assert_eq!(expected, compute_effective_time_marginal(mode, confidence, &NoopProgress));
  }
}

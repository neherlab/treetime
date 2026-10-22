#[cfg(test)]
mod tests {
  use helpers::PropAssertion;
  use proptest::test_runner::TestCaseError;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::ulps(          PropAssertion::Ulps)]
  #[case::relative(      PropAssertion::Relative)]
  #[case::array_ulps(    PropAssertion::ArrayUlps)]
  #[case::array_relative(PropAssertion::ArrayRelative)]
  #[case::map_ulps(      PropAssertion::MapUlps)]
  #[case::map_relative(  PropAssertion::MapRelative)]
  #[trace]
  fn test_assert_prop_rejects_tiny_values_far_apart(#[case] assertion: PropAssertion) {
    assert!(matches!(assertion.check(1e-20, 2e-20), Err(TestCaseError::Fail(_))));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::ulps(          PropAssertion::Ulps)]
  #[case::relative(      PropAssertion::Relative)]
  #[case::array_ulps(    PropAssertion::ArrayUlps)]
  #[case::array_relative(PropAssertion::ArrayRelative)]
  #[case::map_ulps(      PropAssertion::MapUlps)]
  #[case::map_relative(  PropAssertion::MapRelative)]
  #[trace]
  fn test_assert_prop_accepts_adjacent_tiny_values(#[case] assertion: PropAssertion) {
    assertion.check(1e-20, 1e-20_f64.next_up()).unwrap();
  }

  mod helpers {
    use crate::{
      prop_assert_array_relative_eq, prop_assert_array_ulps_eq, prop_assert_map_relative_eq, prop_assert_map_ulps_eq,
      prop_assert_relative_eq, prop_assert_ulps_eq,
    };
    use maplit::btreemap;
    use ndarray::array;
    use proptest::test_runner::TestCaseError;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum PropAssertion {
      Ulps,
      Relative,
      ArrayUlps,
      ArrayRelative,
      MapUlps,
      MapRelative,
    }

    impl PropAssertion {
      pub(super) fn check(self, lhs: f64, rhs: f64) -> Result<(), TestCaseError> {
        match self {
          Self::Ulps => prop_assert_ulps_eq!(lhs, rhs, max_ulps = 4),
          Self::Relative => prop_assert_relative_eq!(lhs, rhs, max_relative = 1e-10),
          Self::ArrayUlps => prop_assert_array_ulps_eq!(array![1.0, lhs], array![1.0, rhs], max_ulps = 4),
          Self::ArrayRelative => {
            prop_assert_array_relative_eq!(array![1.0, lhs], array![1.0, rhs], max_relative = 1e-10);
          },
          Self::MapUlps => prop_assert_map_ulps_eq!(btreemap! {1 => lhs}, btreemap! {1 => rhs}, max_ulps = 4),
          Self::MapRelative => {
            prop_assert_map_relative_eq!(btreemap! {1 => lhs}, btreemap! {1 => rhs}, max_relative = 1e-10);
          },
        }
        Ok(())
      }
    }
  }
}

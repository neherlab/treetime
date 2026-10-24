#[cfg(test)]
mod tests {
  use crate::split_fraction::SplitFraction;
  use proptest::prelude::*;
  use treetime_utils::prop_assert_ulps_eq;

  proptest! {
    #[test]
    fn test_prop_split_fraction_lengths_preserve_the_length(
      fraction in 0.0..=1.0_f64,
      length in 0.0..1e3_f64,
    ) {
      let lengths = SplitFraction::new(fraction).unwrap().lengths(length);

      prop_assert!(lengths.parent_side.is_finite() && lengths.parent_side >= 0.0);
      prop_assert!(lengths.child_side.is_finite() && lengths.child_side >= 0.0);
      prop_assert_ulps_eq!(length, lengths.parent_side + lengths.child_side, max_ulps = 4);
    }
  }
}

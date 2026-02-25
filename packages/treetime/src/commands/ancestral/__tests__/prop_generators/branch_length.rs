//! Proptest generators for branch lengths.

use proptest::prelude::*;

/// Minimum branch length for numerical stability.
pub const MIN_BRANCH_LENGTH: f64 = 0.001;

/// Generate branch lengths with diverse coverage.
///
/// Avoids near-zero values (< 0.001) that cause numerical instability.
/// Uses `prop_filter` to prevent proptest shrinking below the minimum.
pub fn arb_branch_length() -> impl Strategy<Value = f64> {
  prop_oneof![
    1 => Just(MIN_BRANCH_LENGTH),  // minimum stable (won't shrink)
    4 => 0.001..0.1_f64,           // short (most common)
    3 => 0.1..0.5_f64,             // moderate
    1 => 0.5..2.0_f64,             // long
  ]
  .prop_filter("branch length >= MIN_BRANCH_LENGTH", |&bl| bl >= MIN_BRANCH_LENGTH)
}

#[cfg(test)]
mod tests {
  use super::*;
  use proptest::proptest;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn test_prop_branch_length_arb_branch_length_positive(bl in arb_branch_length()) {
      prop_assert!(bl >= 0.001, "Branch length must be >= 0.001: {bl}");
      prop_assert!(bl.is_finite(), "Branch length must be finite: {bl}");
    }

    #[test]
    fn test_prop_branch_length_arb_branch_length_bounded(bl in arb_branch_length()) {
      prop_assert!(bl <= 2.0, "Branch length should be bounded: {bl}");
    }
  }
}

#[cfg(test)]
mod tests {
  use crate::cast::{F64_EXACT_INTEGER_LIMIT, ToF64Exact, f64_to_usize_saturating, u32_to_usize, usize_to_u64};
  use proptest::prelude::*;
  use rstest::rstest;

  #[rstest]
  #[case::zero(0, 0.0)]
  #[case::count(42, 42.0)]
  #[case::limit(1 << 53, 9_007_199_254_740_992.0)]
  fn test_cast_usize_converts_exactly(#[case] value: usize, #[case] expected: f64) {
    assert_eq!(expected.to_bits(), value.to_f64_exact().to_bits());
  }

  #[rstest]
  #[case::negative_limit(-(1 << 53), -9_007_199_254_740_992.0)]
  #[case::negative(-3, -3.0)]
  fn test_cast_i64_converts_negative_values_exactly(#[case] value: i64, #[case] expected: f64) {
    assert_eq!(expected.to_bits(), value.to_f64_exact().to_bits());
  }

  #[test]
  fn test_cast_u32_converts_its_full_range() {
    assert_eq!(4_294_967_295.0_f64.to_bits(), u32::MAX.to_f64_exact().to_bits());
  }

  #[test]
  #[should_panic(expected = "integer 9007199254740993 has no exact f64 representation")]
  fn test_cast_usize_above_limit_panics_in_debug_builds() {
    let above_limit = usize::try_from(F64_EXACT_INTEGER_LIMIT + 1).unwrap();
    let _ = above_limit.to_f64_exact();
  }

  #[test]
  #[should_panic(expected = "integer -9007199254740993 has no exact f64 representation")]
  fn test_cast_i64_below_negative_limit_panics_in_debug_builds() {
    let below_limit = -i64::try_from(F64_EXACT_INTEGER_LIMIT + 1).unwrap();
    let _ = below_limit.to_f64_exact();
  }

  #[test]
  fn test_cast_usize_to_u64_keeps_the_largest_value() {
    assert_eq!(u64::MAX, usize_to_u64(usize::MAX));
  }

  #[test]
  fn test_cast_u32_to_usize_keeps_the_largest_value() {
    let expected = usize::try_from(u32::MAX).unwrap();
    assert_eq!(expected, u32_to_usize(u32::MAX));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::fraction_truncates(2.7,           2)]
  #[case::negative_fraction( -0.5,          0)]
  #[case::negative(          -3.0,          0)]
  #[case::nan(               f64::NAN,      0)]
  #[case::too_large(         1e30,          usize::MAX)]
  #[case::infinity(          f64::INFINITY, usize::MAX)]
  #[trace]
  fn test_cast_f64_to_usize_saturating_matches_as_semantics(#[case] value: f64, #[case] expected: usize) {
    assert_eq!(expected, f64_to_usize_saturating(value));
  }

  proptest! {
    #[test]
    fn test_prop_cast_u64_roundtrips_up_to_limit(value in 0..=F64_EXACT_INTEGER_LIMIT) {
      let converted = value.to_f64_exact();
      prop_assert_eq!(converted.to_string(), value.to_string());
    }
  }
}

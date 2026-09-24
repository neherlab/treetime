use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use num_traits::clamp_min;

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn fix_branch_length(seq_length: usize, branch_length: f64) -> f64 {
  let one_mutation = 1.0 / seq_length as f64;
  let min_branch_len = MIN_BRANCH_LENGTH_FRACTION * one_mutation;
  clamp_min(branch_length, min_branch_len)
}

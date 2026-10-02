use crate::branch_lengths::one_mutation;
use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use num_traits::clamp_min;

pub(crate) fn fix_branch_length(seq_length: usize, branch_length: f64) -> f64 {
  let min_branch_len = MIN_BRANCH_LENGTH_FRACTION * one_mutation(seq_length);
  clamp_min(branch_length, min_branch_len)
}

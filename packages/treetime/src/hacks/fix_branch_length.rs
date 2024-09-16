use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use num_traits::clamp_min;

/// Clamp minimum branch length to a fraction of "1 mutation" which depends on sequence length
///
/// HACK: this is required to avoid NaNs when branch length is zero. Decide whether it is worth to fix
/// the incorrect input trees with branch lengths of 0, or to simply report them as errors instead. Perhaps users would
/// prefer to know that their inputs are broken and to act on that, rather than hide it. (And if not, the errors can
/// potentially be turned off with a CLI flag, then the "fixup" will be enabled)
pub fn fix_branch_length(seq_length: usize, branch_length: f64) -> f64 {
  let one_mutation = 1.0 / seq_length as f64;
  let min_branch_len = MIN_BRANCH_LENGTH_FRACTION * one_mutation;
  clamp_min(branch_length, min_branch_len)
}

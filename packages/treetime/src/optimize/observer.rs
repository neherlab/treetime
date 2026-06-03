use argmin::core::observers::Observe;
use argmin::core::{Error, KV, State};
use log::debug;

/// Logging observer for argmin optimization progress.
///
/// Shared across optimization sites (Tc time scale, Brent root search, golden
/// section root search). Logs at `debug` level on every 10th iteration plus the
/// first `early_threshold` iterations, surfacing early convergence behavior
/// while keeping later output sparse.
pub struct OptimizationObserver {
  /// Human-readable name of the optimization, used as the log prefix.
  pub label: &'static str,
  /// Number of leading iterations to log unconditionally.
  pub early_threshold: u64,
}

impl<I> Observe<I> for OptimizationObserver
where
  I: State,
  <I as State>::Param: std::fmt::Debug,
  <I as State>::Float: std::fmt::LowerExp,
{
  fn observe_iter(&mut self, state: &I, _kv: &KV) -> Result<(), Error> {
    let iter = state.get_iter();
    if should_log(iter, self.early_threshold) {
      debug!(
        "{} iteration {}: best_param = {:?}, best_cost = {:.6e}",
        self.label,
        iter,
        state.get_best_param(),
        state.get_best_cost()
      );
    }
    Ok(())
  }
}

/// Decides whether iteration `iter` should be logged.
///
/// True on every 10th iteration, or while still within the first
/// `early_threshold` iterations.
fn should_log(iter: u64, early_threshold: u64) -> bool {
  iter.is_multiple_of(10) || iter <= early_threshold
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  // iter 0 always logs (multiple of 10 and within any threshold)
  #[case::zero_threshold3(          (0,  3),  true)]
  #[case::early_within_threshold3(  (3,  3),  true)]
  #[case::first_after_threshold3(   (4,  3),  false)]
  #[case::just_before_ten_thr3(     (9,  3),  false)]
  #[case::multiple_of_ten_thr3(     (10, 3),  true)]
  #[case::between_multiples_thr3(   (15, 3),  false)]
  #[case::early_within_threshold5(  (5,  5),  true)]
  #[case::first_after_threshold5(   (6,  5),  false)]
  #[case::just_before_ten_thr5(     (9,  5),  false)]
  #[case::multiple_of_ten_thr5(     (20, 5),  true)]
  #[case::between_multiples_thr5(   (23, 5),  false)]
  #[trace]
  fn test_observer_should_log(#[case] (iter, early_threshold): (u64, u64), #[case] expected: bool) {
    assert_eq!(expected, should_log(iter, early_threshold));
  }
}

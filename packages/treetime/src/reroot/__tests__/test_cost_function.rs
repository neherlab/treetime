#[cfg(test)]
mod tests {
  use crate::reroot::cost_function::EdgeCostFn;
  use crate::reroot::div_stats::DivStats;
  use crate::reroot::traits::RootStats;
  use approx::assert_ulps_eq;
  use argmin::core::CostFunction;

  fn two_tip_cost_fn() -> EdgeCostFn<DivStats> {
    EdgeCostFn {
      to_parent: DivStats::new(1.0, 0.1, 0.01),
      to_child: DivStats::new(1.0, 0.3, 0.09),
      branch_length: 0.4,
      branch_variance: 0.0,
      is_leaf: false,
      leaf_time: None,
      variance_offset_leaf: 1.0,
    }
  }

  #[test]
  fn test_cost_function_evaluate_at_endpoints() {
    let cost_fn = two_tip_cost_fn();

    // x=0: root at source. Parent message propagates 0, child propagates full bl.
    let at_source = cost_fn.evaluate(0.0);
    // x=0: to_child propagates by 0 (stays at parent), to_parent propagates by bl.
    assert_ulps_eq!(at_source.count(), 2.0, max_ulps = 4);

    // x=1: root at target. Parent propagates full bl, child propagates 0.
    let at_target = cost_fn.evaluate(1.0);
    assert_ulps_eq!(at_target.count(), 2.0, max_ulps = 4);

    // The two endpoints should give different scores (asymmetric tips).
    let diff = (at_source.score() - at_target.score()).abs();
    assert!(diff > 1e-10, "expected different scores at endpoints, diff={diff}");
  }

  #[test]
  fn test_cost_function_evaluate_at_midpoint() {
    let cost_fn = two_tip_cost_fn();
    let at_mid = cost_fn.evaluate(0.5);
    assert_ulps_eq!(at_mid.count(), 2.0, max_ulps = 4);
    assert!(at_mid.score().is_finite());
  }

  #[test]
  fn test_cost_function_out_of_range_returns_infinity() {
    let cost_fn = two_tip_cost_fn();
    let neg = (&cost_fn).cost(&-0.1).unwrap();
    let over = (&cost_fn).cost(&1.1).unwrap();
    assert!(neg.is_infinite() && neg > 0.0, "expected +INFINITY for x<0, got {neg}");
    assert!(
      over.is_infinite() && over > 0.0,
      "expected +INFINITY for x>1, got {over}"
    );
  }

  #[test]
  fn test_cost_function_in_range_returns_finite() {
    let cost_fn = two_tip_cost_fn();
    let at_zero = (&cost_fn).cost(&0.0).unwrap();
    let at_half = (&cost_fn).cost(&0.5).unwrap();
    let at_one = (&cost_fn).cost(&1.0).unwrap();
    assert!(at_zero.is_finite());
    assert!(at_half.is_finite());
    assert!(at_one.is_finite());
  }
}

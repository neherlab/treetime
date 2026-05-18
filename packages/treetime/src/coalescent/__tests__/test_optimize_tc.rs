#[cfg(test)]
mod tests {
  use super::super::helpers::setup_graph;
  use crate::coalescent::optimize_tc::optimize_tc;
  use eyre::Report;
  use rstest::rstest;

  #[test]
  fn test_optimize_tc_converges() -> Result<(), Report> {
    let graph = setup_graph()?;
    let result = optimize_tc(&graph, 1.0)?;

    assert!(result.success, "Tc optimization should succeed");
    assert!(result.tc > 0.0, "Optimized Tc should be positive");
    assert!(result.tc.is_finite(), "Optimized Tc should be finite");
    assert!(result.likelihood.is_finite(), "Likelihood should be finite");

    Ok(())
  }

  #[rstest]
  #[case(0.1, 1.0)]
  #[case(0.1, 10.0)]
  #[case(1.0, 10.0)]
  fn test_optimize_tc_convergence_from_different_starts(
    #[case] initial_tc_a: f64,
    #[case] initial_tc_b: f64,
  ) -> Result<(), Report> {
    let graph = setup_graph()?;

    let result_a = optimize_tc(&graph, initial_tc_a)?;
    let result_b = optimize_tc(&graph, initial_tc_b)?;

    assert!(result_a.success, "optimization from {initial_tc_a} failed");
    assert!(result_b.success, "optimization from {initial_tc_b} failed");

    let ratio = result_b.tc / result_a.tc;
    assert!(
      ratio > 0.5 && ratio < 2.0,
      "Optimized Tc values should be similar regardless of initial value: {:.4e} vs {:.4e}",
      result_a.tc,
      result_b.tc
    );

    Ok(())
  }
}

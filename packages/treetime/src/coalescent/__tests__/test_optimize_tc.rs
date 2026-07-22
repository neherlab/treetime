#[cfg(test)]
mod tests {
  use super::super::helpers::setup_graph;
  use crate::coalescent::optimize_tc::optimize_tc;
  use eyre::Report;

  #[test]
  fn test_optimize_tc_returns_positive_finite_optimum() -> Result<(), Report> {
    let graph = setup_graph()?;
    let result = optimize_tc(&graph)?;

    assert!(result.tc > 0.0, "Optimized Tc should be positive");
    assert!(result.tc.is_finite(), "Optimized Tc should be finite");
    assert!(result.likelihood.is_finite(), "Likelihood should be finite");

    Ok(())
  }

  #[test]
  fn test_optimize_tc_is_deterministic() -> Result<(), Report> {
    // The optimum is a closed form (Tc = I/M), so repeated calls are identical.
    let graph = setup_graph()?;

    let a = optimize_tc(&graph)?;
    let b = optimize_tc(&graph)?;

    assert_eq!(a.tc, b.tc);
    assert_eq!(a.likelihood, b.likelihood);

    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use super::super::helpers::setup_graph;
  use crate::coalescent::optimize_tc::optimize_tc;
  use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
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

  #[test]
  fn test_optimize_tc_matches_single_segment_skyline() -> Result<(), Report> {
    // A constant Tc is defined as the one-segment skyline, and `optimize_tc` is a
    // facade over it. Lock that contract so the two cannot silently diverge.
    let graph = setup_graph()?;

    let constant = optimize_tc(&graph)?;
    let skyline = optimize_skyline(&graph, &SkylineParams { n_points: 1, ..SkylineParams::default() })?;

    assert_eq!(constant.tc, skyline.tc_values[0]);
    assert_eq!(constant.likelihood, skyline.log_likelihood);

    Ok(())
  }
}

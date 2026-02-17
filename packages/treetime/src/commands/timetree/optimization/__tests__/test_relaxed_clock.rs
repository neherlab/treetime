#[cfg(test)]
mod tests {
  use crate::commands::timetree::optimization::relaxed_clock::apply_relaxed_clock;
  use crate::representation::partition::timetree::GraphTimetree;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use treetime_io::nwk::nwk_read_str;

  /// Build a simple tree with time_length set on edges
  fn build_simple_tree() -> Result<GraphTimetree, Report> {
    let graph: GraphTimetree = nwk_read_str("(A:0.1,B:0.2)root:0.0;")?;

    // Set time_length on edges (simulating clock inference)
    for edge in graph.get_edges() {
      let edge = edge.write_arc();
      let branch_length = edge.payload().read_arc().base.branch_length.unwrap_or(0.0);
      // time_length = branch_length * 100 (arbitrary scaling for testing)
      edge.payload().write_arc().time_length = Some(branch_length * 100.0);
    }

    Ok(graph)
  }

  /// Build a deeper tree with time_length set
  fn build_deep_tree() -> Result<GraphTimetree, Report> {
    let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2)AB:0.15,(C:0.05,D:0.1)CD:0.08)root:0.0;")?;

    for edge in graph.get_edges() {
      let edge = edge.write_arc();
      let branch_length = edge.payload().read_arc().base.branch_length.unwrap_or(0.0);
      edge.payload().write_arc().time_length = Some(branch_length * 100.0);
    }

    Ok(graph)
  }

  #[test]
  fn test_relaxed_clock_default_params_produce_reasonable_gamma() -> Result<(), Report> {
    let graph = build_simple_tree()?;
    let one_mutation = 0.01;
    let params = [1.0, 1.0];

    apply_relaxed_clock(&graph, &params, one_mutation);

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let gamma = edge.payload().read_arc().gamma;

      assert!(gamma >= 0.1, "gamma={gamma} should be >= 0.1 (minimum bound)");
      assert!(gamma < 10.0, "gamma={gamma} should be reasonable (< 10.0)");
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_all_gamma_above_minimum() -> Result<(), Report> {
    let graph = build_deep_tree()?;
    let one_mutation = 0.001;
    let params = [1.0, 1.0];

    apply_relaxed_clock(&graph, &params, one_mutation);

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let gamma = edge.payload().read_arc().gamma;
      assert!(gamma >= 0.1, "gamma={gamma} must be >= 0.1 (algorithm minimum bound)");
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_uniform_branches_produce_similar_gamma() -> Result<(), Report> {
    // Tree with uniform branch lengths
    let graph: GraphTimetree = nwk_read_str("(A:0.1,B:0.1,C:0.1)root:0.0;")?;

    for edge in graph.get_edges() {
      let edge = edge.write_arc();
      edge.payload().write_arc().time_length = Some(10.0);
    }

    let one_mutation = 0.01;
    let params = [1.0, 1.0];
    apply_relaxed_clock(&graph, &params, one_mutation);

    let gammas: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().gamma)
      .collect();

    let mean_gamma: f64 = gammas.iter().sum::<f64>() / gammas.len() as f64;
    for gamma in &gammas {
      assert_ulps_eq!(*gamma, mean_gamma, max_ulps = 1000);
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_empty_params_uses_defaults() -> Result<(), Report> {
    let graph = build_simple_tree()?;
    let one_mutation = 0.01;

    apply_relaxed_clock(&graph, &[], one_mutation);

    for edge in graph.get_edges() {
      let gamma = edge.read_arc().payload().read_arc().gamma;
      assert!(gamma >= 0.1, "gamma should be >= minimum bound");
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_gamma_stored_in_edges() -> Result<(), Report> {
    let graph = build_simple_tree()?;
    let one_mutation = 0.01;
    let params = [1.0, 1.0];

    // Before: gamma should be default 1.0
    for edge in graph.get_edges() {
      let gamma = edge.read_arc().payload().read_arc().gamma;
      assert_ulps_eq!(gamma, 1.0, max_ulps = 4);
    }

    apply_relaxed_clock(&graph, &params, one_mutation);

    // After: gamma values should be computed (may differ from 1.0)
    let mut any_changed = false;
    for edge in graph.get_edges() {
      let gamma = edge.read_arc().payload().read_arc().gamma;
      if (gamma - 1.0).abs() > 1e-6 {
        any_changed = true;
      }
    }

    assert!(any_changed, "At least one gamma should differ from default 1.0");

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_high_slack_pulls_toward_one() -> Result<(), Report> {
    let graph = build_deep_tree()?;
    let one_mutation = 0.01;

    // Compare low slack vs high slack - high slack should have gammas closer to 1.0
    let params_low = [1.0, 1.0];
    apply_relaxed_clock(&graph, &params_low, one_mutation);

    let gammas_low: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().gamma)
      .collect();
    let deviation_low: f64 = gammas_low.iter().map(|g| (g - 1.0).abs()).sum();

    let params_high = [100.0, 1.0];
    apply_relaxed_clock(&graph, &params_high, one_mutation);

    let gammas_high: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().gamma)
      .collect();
    let deviation_high: f64 = gammas_high.iter().map(|g| (g - 1.0).abs()).sum();

    // High slack should reduce total deviation from 1.0
    assert!(
      deviation_high <= deviation_low,
      "High slack should pull gammas toward 1.0: high_dev={deviation_high}, low_dev={deviation_low}"
    );

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_high_coupling_reduces_variation() -> Result<(), Report> {
    let graph = build_deep_tree()?;
    let one_mutation = 0.01;

    // Run with low coupling
    let params_low = [1.0, 0.1];
    apply_relaxed_clock(&graph, &params_low, one_mutation);

    let gammas_low: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().gamma)
      .collect();
    let variance_low = compute_variance(&gammas_low);

    // Run with high coupling
    let params_high = [1.0, 10.0];
    apply_relaxed_clock(&graph, &params_high, one_mutation);

    let gammas_high: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().gamma)
      .collect();
    let variance_high = compute_variance(&gammas_high);

    assert!(
      variance_high <= variance_low,
      "High coupling should reduce gamma variance: high={variance_high}, low={variance_low}"
    );

    Ok(())
  }

  fn compute_variance(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n
  }
}

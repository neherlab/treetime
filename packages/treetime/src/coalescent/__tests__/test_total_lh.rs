#[cfg(test)]
mod tests {
  use super::super::helpers::{coalescent_node_times, setup_graph};
  use crate::coalescent::edge_data::collect_coalescent_edges;
  use crate::coalescent::optimize_tc::optimize_tc;
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use rstest::rstest;
  use treetime_distribution::{Distribution, NegLog};

  #[test]
  fn test_total_lh_returns_finite_value() -> Result<(), Report> {
    let (graph, names, constraints) = setup_graph()?;
    let tc = Distribution::constant(1.0);

    let lh = compute_coalescent_total_lh(&graph, &tc, &coalescent_node_times(&graph, &constraints))?.value();

    assert!(lh.is_finite(), "Total coalescent LH should be finite, got {lh}");
    Ok(())
  }

  #[test]
  fn test_total_lh_negative_for_reasonable_tc() -> Result<(), Report> {
    let (graph, names, constraints) = setup_graph()?;
    let tc = Distribution::constant(1.0);

    let lh = compute_coalescent_total_lh(&graph, &tc, &coalescent_node_times(&graph, &constraints))?.value();

    assert!(lh < 0.0, "Coalescent log-likelihood should be negative, got {lh}");
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::small( 0.1)]
  #[case::medium(1.0)]
  #[case::large( 10.0)]
  #[case::very_large(100.0)]
  #[trace]
  fn test_total_lh_finite_for_tc(#[case] tc_value: f64) -> Result<(), Report> {
    let (graph, names, constraints) = setup_graph()?;
    let tc = Distribution::constant(tc_value);

    let lh = compute_coalescent_total_lh(&graph, &tc, &coalescent_node_times(&graph, &constraints))?.value();

    assert!(lh.is_finite(), "LH should be finite for Tc={tc_value}");
    Ok(())
  }

  #[test]
  fn test_total_lh_differs_across_tc_values() -> Result<(), Report> {
    let (graph, names, constraints) = setup_graph()?;
    let node_times = coalescent_node_times(&graph, &constraints);

    let lh_small = compute_coalescent_total_lh(&graph, &Distribution::constant(0.1), &node_times)?.value();
    let lh_large = compute_coalescent_total_lh(&graph, &Distribution::constant(100.0), &node_times)?.value();

    assert!(
      (lh_small - lh_large).abs() > 1e-6,
      "LH should differ across Tc values: lh(0.1)={lh_small:.6}, lh(100.0)={lh_large:.6}"
    );
    Ok(())
  }

  #[test]
  fn test_total_lh_monotonic_near_optimum() -> Result<(), Report> {
    let (graph, names, constraints) = setup_graph()?;
    let node_times = coalescent_node_times(&graph, &constraints);
    let opt = optimize_tc(&graph, &node_times)?;

    let tc_opt = opt.tc;
    let lh_opt = compute_coalescent_total_lh(&graph, &Distribution::constant(tc_opt), &node_times)?.value();
    let lh_low = compute_coalescent_total_lh(&graph, &Distribution::constant(tc_opt * 0.01), &node_times)?.value();
    let lh_high = compute_coalescent_total_lh(&graph, &Distribution::constant(tc_opt * 100.0), &node_times)?.value();

    assert!(
      lh_opt >= lh_low,
      "LH at optimum ({lh_opt:.4}) should be >= LH below ({lh_low:.4})"
    );
    assert!(
      lh_opt >= lh_high,
      "LH at optimum ({lh_opt:.4}) should be >= LH above ({lh_high:.4})"
    );

    Ok(())
  }

  #[test]
  fn test_total_lh_matches_optimize_tc_likelihood() -> Result<(), Report> {
    let (graph, names, constraints) = setup_graph()?;
    let node_times = coalescent_node_times(&graph, &constraints);
    let opt = optimize_tc(&graph, &node_times)?;

    let lh = compute_coalescent_total_lh(&graph, &Distribution::constant(opt.tc), &node_times)?.value();

    pretty_assert_ulps_eq!(opt.likelihood.value(), lh, max_ulps = 10);

    Ok(())
  }

  #[test]
  fn test_total_lh_with_formula_distribution() -> Result<(), Report> {
    let (graph, names, constraints) = setup_graph()?;
    let node_times = coalescent_node_times(&graph, &constraints);
    let tc_value = 5.0;
    let tc_const = Distribution::constant(tc_value);
    let tc_formula = Distribution::Formula(treetime_distribution::DistributionFormula::new(
      move |_t| Ok(tc_value),
      -1000.0,
      1000.0,
    ));

    let lh_const = compute_coalescent_total_lh(&graph, &tc_const, &node_times)?.value();
    let lh_formula = compute_coalescent_total_lh(&graph, &tc_formula, &node_times)?.value();

    pretty_assert_ulps_eq!(lh_const, lh_formula, max_ulps = 10);

    Ok(())
  }

  #[test]
  fn test_total_lh_rejects_child_older_than_parent() -> Result<(), Report> {
    let (graph, names, constraints) = setup_graph()?;
    let leaf_key = find_node_key_by_name(&graph, &names, "leaf1").expect("leaf1 not found");
    let mut node_times = coalescent_node_times(&graph, &constraints);
    node_times
      .get_mut(&leaf_key)
      .expect("leaf1 in node times")
      .time_dist_likely = Distribution::<NegLog>::point(1990.0, 1.0).likely_time();

    let error = collect_coalescent_edges(&graph, &node_times).unwrap_err();

    assert!(error.to_string().contains("child older than parent"));
    Ok(())
  }
}

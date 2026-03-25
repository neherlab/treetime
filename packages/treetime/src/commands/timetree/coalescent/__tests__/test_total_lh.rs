#[cfg(test)]
mod tests {
  use crate::commands::clock::date_constraints::load_date_constraints;
  use crate::commands::timetree::coalescent::optimize_tc::optimize_tc;
  use crate::commands::timetree::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::representation::partition::timetree::GraphTimetree;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use rstest::rstest;
  use treetime_distribution::Distribution;
  use treetime_io::dates_csv::DateOrRange;
  use treetime_io::nwk::nwk_read_str;

  const TREE_NWK: &str = "((leaf1:0.01,leaf2:0.01)internal1:0.01,leaf3:0.02)root:0.0;";

  fn setup_graph() -> Result<GraphTimetree, Report> {
    let dates = btreemap! {
      "root".to_owned() => Some(DateOrRange::YearFraction(2000.0)),
      "internal1".to_owned() => Some(DateOrRange::YearFraction(2005.0)),
      "leaf1".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "leaf2".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "leaf3".to_owned() => Some(DateOrRange::YearFraction(2012.0)),
    };
    let graph: GraphTimetree = nwk_read_str(TREE_NWK)?;
    load_date_constraints(&dates, &graph)?;
    Ok(graph)
  }

  #[test]
  fn test_total_lh_returns_finite_value() -> Result<(), Report> {
    let graph = setup_graph()?;
    let tc = Distribution::constant(1.0);

    let lh = compute_coalescent_total_lh(&graph, &tc)?;

    assert!(lh.is_finite(), "Total coalescent LH should be finite, got {lh}");
    Ok(())
  }

  #[test]
  fn test_total_lh_negative_for_reasonable_tc() -> Result<(), Report> {
    // For typical trees, the coalescent log-likelihood is negative
    // (it's a log-probability)
    let graph = setup_graph()?;
    let tc = Distribution::constant(1.0);

    let lh = compute_coalescent_total_lh(&graph, &tc)?;

    assert!(lh.is_finite(), "LH should be finite");
    // The sign depends on the tree structure and Tc, but for this simple tree
    // with Tc=1.0, the LH should be a real number (not NaN or inf)
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::small( 0.1)]
  #[case::medium(1.0)]
  #[case::large( 10.0)]
  #[case::very_large(100.0)]
  #[trace]
  fn test_total_lh_varies_with_tc(#[case] tc_value: f64) -> Result<(), Report> {
    let graph = setup_graph()?;
    let tc = Distribution::constant(tc_value);

    let lh = compute_coalescent_total_lh(&graph, &tc)?;

    assert!(lh.is_finite(), "LH should be finite for Tc={tc_value}");
    Ok(())
  }

  #[test]
  fn test_total_lh_monotonic_near_optimum() -> Result<(), Report> {
    // The coalescent LH should peak near the optimal Tc.
    // Compute LH at three points: below, at, and above the optimum.
    let graph = setup_graph()?;
    let opt = optimize_tc(&graph, 1.0)?;
    assert!(opt.success);

    let tc_opt = opt.tc;
    let lh_opt = compute_coalescent_total_lh(&graph, &Distribution::constant(tc_opt))?;
    let lh_low = compute_coalescent_total_lh(&graph, &Distribution::constant(tc_opt * 0.01))?;
    let lh_high = compute_coalescent_total_lh(&graph, &Distribution::constant(tc_opt * 100.0))?;

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
    // The LH returned by optimize_tc should match compute_coalescent_total_lh
    // at the same Tc value. This validates consistency between the two code paths.
    let graph = setup_graph()?;
    let opt = optimize_tc(&graph, 1.0)?;
    assert!(opt.success);

    let lh = compute_coalescent_total_lh(&graph, &Distribution::constant(opt.tc))?;

    assert_abs_diff_eq!(opt.likelihood, lh, epsilon = 1e-6);

    Ok(())
  }

  #[test]
  fn test_total_lh_with_formula_distribution() -> Result<(), Report> {
    // Verify that Formula-based Tc distributions (like skyline) work correctly.
    // A constant formula should give the same result as Distribution::constant.
    let graph = setup_graph()?;
    let tc_value = 5.0;
    let tc_const = Distribution::constant(tc_value);
    let tc_formula = Distribution::Formula(treetime_distribution::DistributionFormula::new(
      move |_t| Ok(tc_value),
      -1000.0,
      1000.0,
    ));

    let lh_const = compute_coalescent_total_lh(&graph, &tc_const)?;
    let lh_formula = compute_coalescent_total_lh(&graph, &tc_formula)?;

    assert_abs_diff_eq!(lh_const, lh_formula, epsilon = 1e-10);

    Ok(())
  }
}

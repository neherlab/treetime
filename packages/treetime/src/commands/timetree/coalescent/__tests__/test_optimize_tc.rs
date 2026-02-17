#[cfg(test)]
mod tests {
  use crate::commands::clock::date_constraints::load_date_constraints;
  use crate::commands::timetree::coalescent::optimize_tc::optimize_tc;
  use crate::representation::partition::timetree::GraphTimetree;
  use eyre::Report;
  use maplit::btreemap;
  use treetime_io::dates_csv::DateOrRange;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_optimize_tc_converges() -> Result<(), Report> {
    const TREE_NWK: &str = "((leaf1:0.01,leaf2:0.01)internal1:0.01,leaf3:0.02)root:0.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateOrRange::YearFraction(2000.0)),
      "internal1".to_owned() => Some(DateOrRange::YearFraction(2005.0)),
      "leaf1".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "leaf2".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "leaf3".to_owned() => Some(DateOrRange::YearFraction(2012.0)),
    };

    let graph: GraphTimetree = nwk_read_str(TREE_NWK)?;
    load_date_constraints(&dates, &graph)?;

    // Run optimization with a reasonable initial Tc
    let initial_tc = 1.0;
    let result = optimize_tc(&graph, initial_tc)?;

    // Optimization should succeed
    assert!(result.success, "Tc optimization should succeed");

    // Optimized Tc should be positive and finite
    assert!(result.tc > 0.0, "Optimized Tc should be positive");
    assert!(result.tc.is_finite(), "Optimized Tc should be finite");

    // Likelihood should be finite
    assert!(result.likelihood.is_finite(), "Likelihood should be finite");

    Ok(())
  }

  #[test]
  fn test_optimize_tc_different_initial_values() -> Result<(), Report> {
    const TREE_NWK: &str = "((leaf1:0.01,leaf2:0.01)internal1:0.01,leaf3:0.02)root:0.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateOrRange::YearFraction(2000.0)),
      "internal1".to_owned() => Some(DateOrRange::YearFraction(2005.0)),
      "leaf1".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "leaf2".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "leaf3".to_owned() => Some(DateOrRange::YearFraction(2012.0)),
    };

    let graph: GraphTimetree = nwk_read_str(TREE_NWK)?;
    load_date_constraints(&dates, &graph)?;

    // Test with different initial Tc values
    let initial_values = [0.1, 1.0, 10.0];
    let mut optimized_values = Vec::new();

    for initial_tc in initial_values {
      let result = optimize_tc(&graph, initial_tc)?;
      if result.success {
        optimized_values.push(result.tc);
      }
    }

    // All successful optimizations should converge to similar values
    if optimized_values.len() >= 2 {
      let first = optimized_values[0];
      for tc in &optimized_values[1..] {
        let ratio = tc / first;
        assert!(
          ratio > 0.5 && ratio < 2.0,
          "Optimized Tc values should be similar regardless of initial value: {first:.4e} vs {tc:.4e}"
        );
      }
    }

    Ok(())
  }
}

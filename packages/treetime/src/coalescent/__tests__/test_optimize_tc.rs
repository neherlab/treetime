#[cfg(test)]
mod tests {
  use super::super::helpers::setup_graph;
  use crate::clock::date_constraints::load_date_constraints;
  use crate::coalescent::optimize_tc::optimize_tc;
  use crate::partition::timetree::GraphTimetree;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::o;

  fn graph_with_dates(nwk: &str, dates: &DatesMap) -> Result<GraphTimetree, Report> {
    let graph = nwk_read_str(nwk)?;
    load_date_constraints(dates, &graph)?;
    Ok(graph)
  }

  /// Isochronous 3-tip tree `((a,b)x, c)root` with all tips sampled at `t_tip`, the
  /// inner coalescence `x` at `t_x`, and the root at `t_root`. The lineage count is
  /// `k = 2` on `(t_root, t_x)` and `k = 3` on `(t_x, t_tip)`, so the pairwise-rate
  /// integral is `I = (t_x - t_root)·1 + (t_tip - t_x)·3` and there are two binary
  /// mergers (`M = 2`). The analytic constant optimum is `Tc* = I/M`.
  fn tree3(t_root: f64, t_x: f64, t_tip: f64) -> Result<GraphTimetree, Report> {
    let dates = btreemap! {
      o!("root") => Some(DateConstraint::exact(t_root)),
      o!("x") => Some(DateConstraint::exact(t_x)),
      o!("a") => Some(DateConstraint::exact(t_tip)),
      o!("b") => Some(DateConstraint::exact(t_tip)),
      o!("c") => Some(DateConstraint::exact(t_tip)),
    };
    graph_with_dates("((a:1,b:1)x:1,c:1)root:0;", &dates)
  }

  /// Hand-derived analytic optimum `Tc* = I/M` for [`tree3`].
  fn tree3_analytic_tc(t_root: f64, t_x: f64, t_tip: f64) -> f64 {
    let i = (t_x - t_root) * 1.0 + (t_tip - t_x) * 3.0;
    i / 2.0
  }

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
  fn test_optimize_tc_equals_analytic_i_over_m() -> Result<(), Report> {
    // Oracle: compare the optimizer against the hand-derived closed form Tc* = I/M
    // on trees of varying depth and inner-node placement (see `tree3`).
    for &(t_x, t_tip) in &[(2005.0, 2010.0), (2003.0, 2015.0), (2010.0, 2020.0)] {
      let graph = tree3(2000.0, t_x, t_tip)?;
      let result = optimize_tc(&graph)?;
      pretty_assert_ulps_eq!(tree3_analytic_tc(2000.0, t_x, t_tip), result.tc, max_ulps = 8);
    }

    Ok(())
  }

  #[test]
  fn test_optimize_tc_scale_equivariant() -> Result<(), Report> {
    // Tc carries units of time: I scales with the time span, M is a pure count, so
    // Tc = I/M scales linearly. Scaling the tree's timescale by s (about the root)
    // scales Tc by s.
    let base = optimize_tc(&tree3(2000.0, 2005.0, 2010.0)?)?.tc;
    let scaled = optimize_tc(&tree3(2000.0, 2010.0, 2020.0)?)?.tc; // s = 2

    pretty_assert_ulps_eq!(2.0 * base, scaled, max_ulps = 8);

    Ok(())
  }

  #[test]
  fn test_optimize_tc_degenerate_tree_errors() -> Result<(), Report> {
    // A tree with no time span (all nodes coincident) is degenerate for the
    // coalescent: the pairwise-rate integral is zero, so no finite positive Tc
    // exists. Optimization must error rather than invent a timescale.
    let graph = tree3(2000.0, 2000.0, 2000.0)?;
    assert!(
      optimize_tc(&graph).is_err(),
      "a zero-span tree must not yield a coalescent Tc"
    );

    Ok(())
  }
}

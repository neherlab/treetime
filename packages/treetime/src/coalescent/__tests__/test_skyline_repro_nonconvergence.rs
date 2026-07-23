#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::load_date_constraints;
  use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
  use crate::partition::timetree::GraphTimetree;
  use eyre::Report;
  use maplit::btreemap;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::nwk::nwk_read_str;

  const SMALL_TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0)root:1.0;";

  fn small_tree_dates() -> DatesMap {
    btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "internal1".to_owned() => Some(DateConstraint::exact(2005.0)),
      "leaf1".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf3".to_owned() => Some(DateConstraint::exact(2012.0)),
    }
  }

  fn small_tree() -> Result<GraphTimetree, Report> {
    let graph = nwk_read_str(SMALL_TREE_NWK)?;
    load_date_constraints(&small_tree_dates(), &graph)?;
    Ok(graph)
  }

  // Counterexample for "the solve is exact and cannot fail numerically".
  //
  // `solve_log_tc` returns `Ok(z)` unconditionally after its iteration loop. With
  // `max_iter = 0` it performs no Newton step and returns the decoupled per-segment
  // warm start, which for the default (positive) stiffness is not the regularized
  // optimum. `optimize_skyline` still returns a successful `SkylineResult`, so a
  // caller cannot distinguish an optimized skyline from an un-optimized one, and
  // the pipeline's fail-loud contract cannot fire because no error is raised.
  //
  // Expected: a solve that ran no iterations (or that exhausts `max_iter` above the
  // gradient tolerance) returns an error rather than reporting success.
  #[test]
  fn test_skyline_zero_iterations_reported_as_success() -> Result<(), Report> {
    let params = SkylineParams {
      n_points: 3,
      max_iter: 0,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&small_tree()?, &params);

    assert!(
      result.is_err(),
      "optimize_skyline with max_iter=0 must not report success on a multi-segment regularized solve, but returned Ok"
    );
    Ok(())
  }
}

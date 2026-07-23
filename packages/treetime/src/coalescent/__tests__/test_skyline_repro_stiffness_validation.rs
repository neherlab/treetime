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

  // Counterexample for the unenforced `stiffness > 0` precondition.
  //
  // `SkylineParams::stiffness` is documented "Must be positive for more than one
  // segment", but `optimize_skyline` validates only `n_points` and the CLI accepts
  // any `f64`. This branch also removed the previous `stiffness <= 0` short-circuit
  // in `solve_log_tc`. A negative stiffness turns the smoothing penalty into
  // `-(|s|/2) * sum (ln(Tc_{i+1}/Tc_i))^2`, which is unbounded below, so the
  // objective is no longer convex and the Thomas solve loses the positive-definite
  // Hessian it assumes (no pivoting).
  //
  // Expected: multi-segment optimization rejects a non-positive stiffness with a
  // typed error. Once the guard is restored, tighten this to
  // `assert_err!(result, "<message>")`.
  #[test]
  fn test_skyline_negative_stiffness_is_rejected() -> Result<(), Report> {
    let params = SkylineParams {
      n_points: 5,
      stiffness: -1.0,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&small_tree()?, &params);

    assert!(
      result.is_err(),
      "multi-segment optimize_skyline must reject negative stiffness, but returned Ok"
    );
    Ok(())
  }
}

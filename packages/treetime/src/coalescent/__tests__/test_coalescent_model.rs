#[cfg(test)]
mod tests {
  use crate::coalescent::coalescent::CoalescentModel;
  use crate::coalescent::edge_data::CoalescentEdgeData;
  use crate::coalescent::time_coordinate::CalendarTime;
  use crate::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};
  use eyre::Report;
  use ndarray::array;
  use proptest::prelude::*;
  use treetime_distribution::Distribution;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_utils::{assert_error, prop_assert_abs_diff_eq};

  #[test]
  fn test_coalescent_model_node_costs_follow_telescoped_objective() -> Result<(), Report> {
    let model = model(array![0.0, 5.0, 10.0], array![1.0, 2.0, 3.0, 0.0], 2.0)?;

    // Analytical oracle: κ=1/4 on [0,5], κ=1/2 on [5,10].
    pretty_assert_ulps_eq!(-3.75, model.leaf_contribution(0.0), max_ulps = 4);
    pretty_assert_abs_diff_eq!(
      2.5 - 1.5_f64.ln(),
      model.internal_contribution(5.0, 2)?,
      epsilon = 1e-10
    );
    pretty_assert_abs_diff_eq!(7.5 - 0.5_f64.ln(), model.root_contribution(0.0, 2)?, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_coalescent_model_branch_merger_rate_schedule_uses_all_breakpoints() -> Result<(), Report> {
    let model = model(array![0.0, 5.0], array![1.0, 3.0, 0.0], 2.0)?;
    let tc_schedule = PiecewiseConstantFn::new(array![2.0], array![2.0, 4.0]);

    let actual = model.branch_merger_rate_schedule(&tc_schedule)?;

    pretty_assert_ulps_eq!(actual.breakpoints(), &array![0.0, 2.0, 5.0], max_ulps = 4);
    pretty_assert_ulps_eq!(actual.values(), &array![0.125, 0.5, 0.25, 0.0625], max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_coalescent_model_branch_merger_rate_schedule_rejects_invalid_tc() -> Result<(), Report> {
    let model = model(array![0.0, 5.0], array![1.0, 3.0, 0.0], 2.0)?;
    let tc_schedule = PiecewiseConstantFn::new(array![], array![f64::NAN]);

    assert_error!(
      model.branch_merger_rate_schedule(&tc_schedule),
      "Coalescent Tc region 0 must be finite and positive, got NaN"
    );
    Ok(())
  }

  proptest! {
    #[test]
    fn test_prop_coalescent_model_node_and_edge_objectives_match(
      n_children in 2_usize..8,
      root_time in 1900.0_f64..2000.0,
      duration in 0.1_f64..100.0,
      tc in 0.001_f64..100.0,
    ) {
      let child_time = root_time + duration;
      let model = model(
        array![root_time, child_time],
        array![1.0, n_children as f64, 0.0],
        tc,
      ).unwrap();
      let edge = CoalescentEdgeData::new(
        CalendarTime::new(child_time),
        CalendarTime::new(root_time),
        n_children as f64,
      );

      let node_cost = model.root_contribution(root_time, n_children).unwrap()
        + n_children as f64 * model.leaf_contribution(child_time);
      let edge_cost = n_children as f64 * model.edge_contribution(&edge).unwrap();

      // Algebraic oracle: telescoping branch survival leaves the grouped
      // leaf/internal/root objective exactly (Kingman 1982).
      prop_assert_abs_diff_eq!(node_cost, edge_cost, epsilon = 1e-10);
    }
  }

  fn model(
    breakpoints: ndarray::Array1<f64>,
    values: ndarray::Array1<f64>,
    tc: f64,
  ) -> Result<CoalescentModel, Report> {
    CoalescentModel::new(
      &PiecewiseConstantFn::new(breakpoints, values),
      &Distribution::constant(tc),
    )
  }
}

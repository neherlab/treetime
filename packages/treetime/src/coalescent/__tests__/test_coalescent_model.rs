#[cfg(test)]
mod tests {
  use crate::coalescent::coalescent::CoalescentModel;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use ndarray::array;
  use treetime_distribution::Distribution;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_utils::assert_error;

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

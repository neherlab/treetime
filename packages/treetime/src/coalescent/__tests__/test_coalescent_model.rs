#[cfg(test)]
mod tests {
  use crate::coalescent::coalescent::CoalescentModel;
  use crate::coalescent::edge_data::CoalescentEdgeData;
  use crate::coalescent::precomputed::CoalescentPrecomputed;
  use crate::coalescent::time_coordinate::CalendarTime;
  use crate::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};
  use eyre::Report;
  use ndarray::array;
  use proptest::prelude::*;
  use treetime_distribution::{Distribution, DistributionFormula};
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_utils::prop_assert_abs_diff_eq;

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
  fn test_coalescent_model_leaf_application_preserves_range_support() -> Result<(), Report> {
    let model = model(array![0.0, 10.0], array![1.0, 2.0, 0.0], 2.0)?;
    let distribution = Distribution::range((0.0, 10.0), 1.0);

    let actual = model.apply_leaf_cost(&distribution)?;

    let expected_times = array![0.0, 10.0];
    assert_eq!(expected_times, actual.t());
    assert!(matches!(actual, Distribution::Function(_)));
    assert!(actual.eval(0.0)? > actual.eval(10.0)?);
    Ok(())
  }

  #[test]
  fn test_coalescent_model_leaf_application_preserves_point() -> Result<(), Report> {
    let model = model(array![0.0, 10.0], array![1.0, 2.0, 0.0], 2.0)?;
    let distribution = Distribution::point(5.0, 0.25);

    let actual = model.apply_leaf_cost(&distribution)?;

    // A one-point likelihood normalizes to unit peak without changing support.
    let expected = Distribution::point(5.0, 1.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_coalescent_model_leaf_application_preserves_empty() -> Result<(), Report> {
    let model = model(array![0.0, 10.0], array![1.0, 2.0, 0.0], 2.0)?;

    let actual = model.apply_leaf_cost(&Distribution::Empty)?;

    assert_eq!(Distribution::Empty, actual);
    Ok(())
  }

  #[test]
  fn test_coalescent_model_application_preserves_function_grid() -> Result<(), Report> {
    let model = model(array![0.0, 10.0], array![1.0, 2.0, 0.0], 2.0)?;
    let expected_times = array![0.0, 2.5, 5.0, 7.5, 10.0];
    let distribution = Distribution::function(expected_times.clone(), array![0.2, 0.5, 1.0, 0.5, 0.2])?;

    let actual = model.apply_leaf_cost(&distribution)?;

    assert_eq!(expected_times, actual.t());
    assert!(!matches!(actual, Distribution::Formula(_)));
    Ok(())
  }

  #[test]
  fn test_coalescent_model_rejects_formula_contribution_input() -> Result<(), Report> {
    let model = model(array![0.0, 10.0], array![1.0, 2.0, 0.0], 2.0)?;
    let distribution = Distribution::Formula(DistributionFormula::new(|_| Ok(1.0), 0.0, 10.0));

    let error = model.apply_leaf_cost(&distribution).unwrap_err();

    assert!(error.to_string().contains("concrete Point, Range, or Function"));
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
    let precomputed = CoalescentPrecomputed::from_lineage_counts(PiecewiseConstantFn::new(breakpoints, values));
    CoalescentModel::new(&precomputed, &Distribution::constant(tc))
  }
}

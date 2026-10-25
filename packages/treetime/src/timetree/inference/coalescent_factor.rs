use crate::coalescent::mean_field::CoalescentPrior;
use crate::timetree::inference::runner::GRID_POINTS;
use eyre::Report;
use ndarray::Array1;
use treetime_distribution::NegLog;
use treetime_distribution::distribution_multiply_by_fn;
use treetime_distribution::{BoundaryBehavior, Distribution, DistributionFunction};
use treetime_graph::node::GraphNodeKey;
use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};

pub(crate) fn apply_coalescent_prior(
  coalescent_model: Option<&CoalescentPrior>,
  key: GraphNodeKey,
  n_children: usize,
  distribution: Distribution<NegLog>,
) -> Result<Distribution<NegLog>, Report> {
  let Some(model) = coalescent_model else {
    return Ok(distribution);
  };
  if n_children == 0 || matches!(distribution, Distribution::Empty) {
    return Ok(distribution);
  }
  distribution_multiply_by_fn(
    &distribution,
    |time| model.node_contribution(key, n_children, time),
    |side, edge| model.node_contribution_outward_slope(key, n_children, side, edge),
  )
}

pub(crate) fn apply_leaf_coalescent_factor(
  coalescent_model: Option<&CoalescentPrior>,
  key: GraphNodeKey,
  distribution: Distribution<NegLog>,
) -> Result<Distribution<NegLog>, Report> {
  let Some(model) = coalescent_model else {
    return Ok(distribution);
  };
  let distribution = match distribution.time_bounds() {
    Some((start, end)) if matches!(distribution, Distribution::Range(_)) => {
      let n_points =
        (f64_to_usize_saturating(((end - start) / model.time_resolution()).ceil()) + 1).clamp(2, GRID_POINTS);
      let dx = (end - start) / (n_points - 1).to_f64_exact();
      let values = Array1::from_elem(n_points, distribution.eval(start)?);
      Distribution::Function(
        DistributionFunction::from_start_dx_values(start, dx, values)?
          .with_left_extrap(BoundaryBehavior::Hard)?
          .with_right_extrap(BoundaryBehavior::Hard)?,
      )
    },
    _ => distribution,
  };
  distribution_multiply_by_fn(
    &distribution,
    |time| model.node_contribution(key, 0, time),
    |side, edge| model.node_contribution_outward_slope(key, 0, side, edge),
  )
}

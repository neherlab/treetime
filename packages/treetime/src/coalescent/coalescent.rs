use crate::coalescent::edge_data::CoalescentEdgeData;
use crate::coalescent::integration::{compute_integral_merger_rate, compute_merger_rate_total_scalar};
use crate::coalescent::precomputed::CoalescentPrecomputed;
use crate::payload::traits::TimetreeNode;
use eyre::{Context, Report};
use ndarray::Array1;
use treetime_distribution::Distribution;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_grid::piecewise_linear_fn::PiecewiseLinearFn;
use treetime_utils::make_error;

/// Builds one calendar-coordinate Kingman coalescent model for a tree and $T_c$.
pub fn compute_coalescent_model<N, E, D>(graph: &Graph<N, E, D>, tc: &Distribution) -> Result<CoalescentModel, Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let precomputed = CoalescentPrecomputed::from_graph(graph)?;
  CoalescentModel::new(&precomputed, tc)
}

/// Calendar-coordinate Kingman coalescent rates and cumulative hazard.
#[derive(Clone, Debug)]
pub struct CoalescentModel {
  lineage_counts: PiecewiseConstantFn,
  tc: Distribution,
  cumulative_branch_rate: PiecewiseLinearFn,
}

impl CoalescentModel {
  pub(crate) fn new(precomputed: &CoalescentPrecomputed, tc: &Distribution) -> Result<Self, Report> {
    let cumulative_branch_rate = compute_integral_merger_rate(tc, precomputed.lineage_counts())?;
    Ok(Self {
      lineage_counts: precomputed.lineage_counts().clone(),
      tc: tc.clone(),
      cumulative_branch_rate,
    })
  }

  pub fn leaf_cost(&self, time: f64) -> f64 {
    -self.cumulative_branch_rate.eval(time)
  }

  pub fn internal_cost(&self, time: f64, n_children: usize) -> Result<f64, Report> {
    let multiplicity = n_children.saturating_sub(1) as f64;
    let total_rate = self.total_rate(time)?;
    Ok(multiplicity * (self.cumulative_branch_rate.eval(time) - total_rate.ln()))
  }

  pub fn root_cost(&self, time: f64, n_children: usize) -> Result<f64, Report> {
    Ok(self.internal_cost(time, n_children)? + self.cumulative_branch_rate.eval(time))
  }

  pub(crate) fn edge_cost(&self, edge: &CoalescentEdgeData) -> Result<f64, Report> {
    let parent_time = edge.parent_time().value();
    let child_time = edge.child_time().value();
    let survival_cost = self.cumulative_branch_rate.eval(parent_time) - self.cumulative_branch_rate.eval(child_time);
    let merger_credit = self.total_rate(parent_time)?.ln() * (edge.multiplicity() - 1.0) / edge.multiplicity();
    Ok(survival_cost - merger_credit)
  }

  pub(crate) fn apply_leaf_cost(&self, distribution: &Distribution) -> Result<Distribution, Report> {
    self.apply_cost(distribution, |time| Ok(self.leaf_cost(time)))
  }

  pub(crate) fn apply_internal_cost(
    &self,
    distribution: &Distribution,
    n_children: usize,
  ) -> Result<Distribution, Report> {
    self.apply_cost(distribution, |time| self.internal_cost(time, n_children))
  }

  pub(crate) fn apply_root_cost(&self, distribution: &Distribution, n_children: usize) -> Result<Distribution, Report> {
    self.apply_cost(distribution, |time| self.root_cost(time, n_children))
  }

  fn total_rate(&self, time: f64) -> Result<f64, Report> {
    // Calendar right-continuity gives the number of lineages immediately on
    // the sampled-tree side of a merger, equivalent to TBP eval_left().
    let k = self.lineage_counts.eval(time);
    let tc = self
      .tc
      .eval(time)
      .wrap_err_with(|| format!("When evaluating coalescent Tc at calendar time {time:.6e}"))?;
    if !k.is_finite() {
      return make_error!("Coalescent lineage count must be finite at calendar time {time:.6e}, got {k:.6e}");
    }
    if !tc.is_finite() || tc <= 0.0 {
      return make_error!("Coalescent Tc must be finite and positive at calendar time {time:.6e}, got {tc:.6e}");
    }
    Ok(compute_merger_rate_total_scalar(k, tc))
  }

  fn apply_cost<F>(&self, distribution: &Distribution, cost: F) -> Result<Distribution, Report>
  where
    F: Fn(f64) -> Result<f64, Report>,
  {
    match distribution {
      Distribution::Empty => Ok(Distribution::Empty),
      Distribution::Formula(_) => {
        make_error!("Coalescent contributions require a concrete Point, Range, or Function distribution")
      },
      Distribution::Point(_) | Distribution::Range(_) | Distribution::Function(_) => {
        let times = distribution.t();
        let costs = times
          .iter()
          .map(|&time| cost(time))
          .collect::<Result<Vec<_>, Report>>()?;
        let minimum = costs
          .iter()
          .copied()
          .reduce(f64::min)
          .filter(|minimum| minimum.is_finite())
          .ok_or_else(|| eyre::eyre!("Coalescent contribution has no finite cost"))?;
        let amplitudes = Array1::from_iter(
          distribution
            .y()
            .iter()
            .zip(costs)
            .map(|(&amplitude, cost)| amplitude * (minimum - cost).exp()),
        );
        Ok(Distribution::function(times, amplitudes)?.normalize())
      },
    }
  }
}

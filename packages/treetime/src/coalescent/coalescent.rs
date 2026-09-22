use crate::coalescent::edge_data::CoalescentEdgeData;
use crate::coalescent::integration::{
  compute_integral_merger_rate, compute_merger_rate_per_lineage_scalar, compute_merger_rate_total_scalar,
};
use eyre::{Context, Report};
use treetime_distribution::Distribution;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_grid::piecewise_linear_fn::PiecewiseLinearFn;
use treetime_utils::make_error;

#[derive(Clone, Debug)]
pub struct CoalescentModel {
  lineage_counts: PiecewiseConstantFn,
  tc: Distribution,
  expected_mergers: PiecewiseLinearFn,
}

impl CoalescentModel {
  pub(crate) fn new(lineage_counts: &PiecewiseConstantFn, tc: &Distribution) -> Result<Self, Report> {
    let expected_mergers = compute_integral_merger_rate(tc, lineage_counts)?;
    Ok(Self {
      lineage_counts: lineage_counts.clone(),
      tc: tc.clone(),
      expected_mergers,
    })
  }

  pub(crate) fn leaf_contribution(&self, time: f64) -> f64 {
    -self.expected_mergers.eval(time)
  }

  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  pub(crate) fn internal_contribution(&self, time: f64, n_children: usize) -> Result<f64, Report> {
    let n_mergers = n_children.saturating_sub(1) as f64;
    let total_merger_rate = self.total_merger_rate(time)?;
    Ok(n_mergers * (self.expected_mergers.eval(time) - total_merger_rate.ln()))
  }

  pub(crate) fn root_contribution(&self, time: f64, n_children: usize) -> Result<f64, Report> {
    Ok(self.internal_contribution(time, n_children)? + self.expected_mergers.eval(time))
  }

  pub(crate) fn edge_contribution(&self, edge: &CoalescentEdgeData) -> Result<f64, Report> {
    let parent_time = edge.parent_time().value();
    let child_time = edge.child_time().value();
    let survival_term = self.expected_mergers.eval(parent_time) - self.expected_mergers.eval(child_time);
    let n_siblings = edge.n_siblings();
    let merger_credit = self.total_merger_rate(parent_time)?.ln() * (n_siblings - 1.0) / n_siblings;
    Ok(survival_term - merger_credit)
  }

  fn total_merger_rate(&self, time: f64) -> Result<f64, Report> {
    let (k, tc) = self.lineage_count_and_tc(time)?;
    Ok(compute_merger_rate_total_scalar(k, tc))
  }

  pub fn branch_merger_rate(&self, time: f64) -> Result<f64, Report> {
    let (k, tc) = self.lineage_count_and_tc(time)?;
    Ok(compute_merger_rate_per_lineage_scalar(k, tc))
  }

  pub(crate) fn branch_merger_rate_schedule(
    &self,
    tc_schedule: &PiecewiseConstantFn,
  ) -> Result<PiecewiseConstantFn, Report> {
    for (index, &k) in self.lineage_counts.values().iter().enumerate() {
      if !k.is_finite() {
        return make_error!("Coalescent lineage count region {index} must be finite, got {k:.6e}");
      }
    }
    for (index, &tc) in tc_schedule.values().iter().enumerate() {
      if !tc.is_finite() || tc <= 0.0 {
        return make_error!("Coalescent Tc region {index} must be finite and positive, got {tc:.6e}");
      }
    }

    Ok(
      self
        .lineage_counts
        .zip_map(tc_schedule, compute_merger_rate_per_lineage_scalar),
    )
  }

  fn lineage_count_and_tc(&self, time: f64) -> Result<(f64, f64), Report> {
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
    Ok((k, tc))
  }
}

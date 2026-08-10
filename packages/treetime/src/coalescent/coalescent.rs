use crate::coalescent::edge_data::CoalescentEdgeData;
use crate::coalescent::integration::{
  compute_integral_merger_rate, compute_merger_rate_per_lineage_scalar, compute_merger_rate_total_scalar,
};
use eyre::{Context, Report};
use treetime_distribution::Distribution;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_grid::piecewise_linear_fn::PiecewiseLinearFn;
use treetime_utils::make_error;

/// Calendar-coordinate Kingman coalescent rates and expected merger counts.
///
/// The three quantities enter differently, which is why the model is assembled from them rather
/// than read off a tree: $k(t)$ is a fixed input, $T_c$ is estimated, and $H(t)$ is a compound
/// of the two, materialised once here because it is an integral.
///
/// $k(t)$ has two roles, and only one of them may track the times being inferred. As the
/// **prior** it sets the merger rate imposed on node times, and must be held fixed: deriving it
/// from the current times makes the prior self-referential, so each pass moves the times, which
/// moves $k(t)$, which moves the prior the next pass is inferred under, and the loop chases a
/// receding target. As the **statistic** that $T_c$ is estimated from and the likelihood is
/// evaluated against, it is read from the live tree, in
/// [`optimize_skyline`](crate::coalescent::skyline::optimize_skyline) and
/// [`compute_coalescent_total_lh`](crate::coalescent::total_lh::compute_coalescent_total_lh).
///
/// Correspondence to v0 (`packages/legacy/treetime/treetime/merger_models.py`):
///
/// | v1 (this struct)          | v0                     | Quantity                                                    |
/// | ------------------------- | ---------------------- | ----------------------------------------------------------- |
/// | `lineage_counts`          | `nlineages`            | number of extant lineages $k(t)$                            |
/// | `tc`                      | `Tc`                   | coalescent time scale $T_c(t)$                              |
/// | `expected_mergers`        | `integral_merger_rate` | $H(t)=\int_0^t \kappa(s)\,ds$, expected mergers on a branch |
/// | `total_merger_rate(t)`    | `total_merger_rate`    | total pairwise merger rate $\lambda(t)$                     |
/// | (`integration.rs`: $\kappa$) | `branch_merger_rate` | per-branch merger rate $\kappa(t)$                          |
#[derive(Clone, Debug)]
pub struct CoalescentModel {
  lineage_counts: PiecewiseConstantFn,
  tc: Distribution,
  /// $H(t)=\int_0^t \kappa(s)\,ds$: the expected number of coalescent merger
  /// events a branch experiences from the present to calendar time $t$. v0's
  /// `integral_merger_rate`.
  expected_mergers: PiecewiseLinearFn,
}

impl CoalescentModel {
  /// Assemble a model from lineage counts and a timescale.
  ///
  /// Cheap enough to rebuild whenever $T_c$ is re-estimated: it is one integral over the
  /// breakpoints of `lineage_counts`, next to the solve that produced `tc`.
  pub fn new(lineage_counts: &PiecewiseConstantFn, tc: &Distribution) -> Result<Self, Report> {
    let expected_mergers = compute_integral_merger_rate(tc, lineage_counts)?;
    Ok(Self {
      lineage_counts: lineage_counts.clone(),
      tc: tc.clone(),
      expected_mergers,
    })
  }

  // Per-node and per-edge additive terms of the Kingman coalescent objective,
  // matching v0's signed `node_contribution`. Each is the term's contribution to
  // the coalescent cost (negative log-likelihood); a value can be negative, as a
  // leaf's branch-survival credit is. `coalescent_log_likelihood` sums the edge
  // contributions and negates the total.
  pub fn leaf_contribution(&self, time: f64) -> f64 {
    -self.expected_mergers.eval(time)
  }

  pub fn internal_contribution(&self, time: f64, n_children: usize) -> Result<f64, Report> {
    let n_mergers = n_children.saturating_sub(1) as f64;
    let total_merger_rate = self.total_merger_rate(time)?;
    Ok(n_mergers * (self.expected_mergers.eval(time) - total_merger_rate.ln()))
  }

  pub fn root_contribution(&self, time: f64, n_children: usize) -> Result<f64, Report> {
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

  /// Per-branch merger rate $\kappa(t)$: the rate at which one lineage merges with any
  /// other. v0's `branch_merger_rate`.
  ///
  /// Distinct from [`Self::total_merger_rate`], which is the rate over all pairs. Callers
  /// simulating mergers within one local group scale this by their own lineage count rather
  /// than the tree-wide $k(t)$ baked into the total.
  pub fn branch_merger_rate(&self, time: f64) -> Result<f64, Report> {
    let (k, tc) = self.lineage_count_and_tc(time)?;
    Ok(compute_merger_rate_per_lineage_scalar(k, tc))
  }

  /// Piecewise-constant per-branch merger rate over all lineage-count and Tc changes.
  pub fn branch_merger_rate_schedule(&self, tc_schedule: &PiecewiseConstantFn) -> Result<PiecewiseConstantFn, Report> {
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

    // The breakpoint union keeps every rate discontinuity visible to event sampling.
    Ok(
      self
        .lineage_counts
        .zip_map(tc_schedule, compute_merger_rate_per_lineage_scalar),
    )
  }

  fn lineage_count_and_tc(&self, time: f64) -> Result<(f64, f64), Report> {
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
    Ok((k, tc))
  }
}

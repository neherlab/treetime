use crate::coalescent::edge_data::CoalescentEdgeData;
use crate::coalescent::integration::{compute_integral_merger_rate, compute_merger_rate_total_scalar};
use crate::coalescent::precomputed::CoalescentPrecomputed;
use crate::payload::traits::TimetreeNode;
use eyre::{Context, Report};
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

/// Calendar-coordinate Kingman coalescent rates and expected merger counts.
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
  pub(crate) fn new(precomputed: &CoalescentPrecomputed, tc: &Distribution) -> Result<Self, Report> {
    let expected_mergers = compute_integral_merger_rate(tc, precomputed.lineage_counts())?;
    Ok(Self {
      lineage_counts: precomputed.lineage_counts().clone(),
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
}

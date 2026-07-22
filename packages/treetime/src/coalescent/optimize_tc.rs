use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::edge_data::{CoalescentEdgeData, coalescent_log_likelihood, collect_coalescent_edges};
use crate::coalescent::integration::compute_integral_merger_rate;
use crate::coalescent::precomputed::CoalescentPrecomputed;
use crate::payload::traits::TimetreeNode;
use eyre::Report;
use log::{info, warn};
use treetime_distribution::Distribution;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

/// Result of Tc optimization.
pub struct OptimizeTcResult {
  /// Optimized coalescence time scale.
  pub tc: f64,
  /// Total coalescent likelihood at optimized Tc.
  pub likelihood: f64,
  /// Whether optimization succeeded.
  pub success: bool,
}

/// Computes the coalescence time scale Tc that maximizes the coalescent likelihood.
///
/// # Algorithm
///
/// The Kingman coalescent likelihood factorizes over intervals of constant lineage
/// count `k`. Each merger event contributes a factor `(k(k-1)/2 / Tc)` and each
/// interval of duration `Δt` a survival factor `exp(-Δt · k(k-1)/2 / Tc)`. Taking
/// the product over the whole tree collapses to
///
/// ```text
///   L(Tc) ∝ Tc^(-M) · exp(-I / Tc)
/// ```
///
/// where `I = ∫ k(k-1)/2 dt` is the pairwise-merger-rate integral over the tree
/// (with `Tc` factored out) and `M` is the total number of merger events.
/// Maximizing `log L = -M·ln(Tc) - I/Tc` yields the closed form
///
/// ```text
///   Tc* = I / M.
/// ```
///
/// Both `I` and `M` are accumulated per edge rather than as global quantities: this
/// keeps them consistent with the bad branches that [`collect_coalescent_edges`]
/// excludes, and generalizes cleanly to a piecewise-constant `Tc` (where per-interval
/// mergers matter). Summed over the `k` edges spanning an interval, the per-edge
/// per-lineage integral `∫ (k-1)/2 dt` reproduces `∫ k(k-1)/2 dt`; summed over a
/// node's `n_siblings` edges, `(n_siblings - 1)/n_siblings` reproduces the node's
/// merger count `n_siblings - 1`.
///
/// # Returns
///
/// Returns the optimized Tc value, or falls back to the initial Tc on failure.
pub fn optimize_tc<N, E, D>(graph: &Graph<N, E, D>, initial_tc: f64) -> Result<OptimizeTcResult, Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge,
  D: Sync + Send,
{
  info!("Optimizing coalescent time scale Tc (initial Tc = {initial_tc:.6e})");

  let precomputed = CoalescentPrecomputed::from_graph(graph)?;
  let edges = collect_coalescent_edges(graph)?;

  // Per-lineage merger integral H₀(t) = ∫ₜᴾ (k-1)/2 dt with Tc factored out, i.e.
  // the survival integrand at Tc = 1. Within any edge span k ≥ 2, so the
  // `max(0.5, k-1)` clamp inside `compute_integral_merger_rate` is inert here and
  // H₀ equals the textbook ∫ (k-1)/2 dt.
  let bare_integral = compute_integral_merger_rate(&Distribution::constant(1.0), precomputed.lineage_counts())?;

  let mut integral = 0.0;
  let mut n_mergers = 0.0;
  for edge in &edges {
    integral += bare_integral.eval(edge.parent_time().value()) - bare_integral.eval(edge.child_time().value());
    let n_siblings = edge.n_siblings();
    // n_siblings is the number of children of the edge's parent node, so this term
    // is added n_siblings times per parent and sums to (n_siblings - 1) there.
    n_mergers += (n_siblings - 1.0) / n_siblings;
  }
  info!("Tc optimization: integral = {integral:.6e}, mergers = {n_mergers:.6e}");
  let tc = integral / n_mergers;
  let success = n_mergers > 0.0 && integral.is_finite() && tc.is_finite() && tc > 0.0;

  let tc = if success {
    tc
  } else {
    warn!(
      "Analytic Tc optimum unavailable (integral = {integral:.6e}, mergers = {n_mergers:.6e}); \
       falling back to initial Tc {initial_tc:.6e}"
    );
    initial_tc
  };

  let likelihood = compute_total_lh(&precomputed, &edges, tc)?;

  info!("Tc optimization completed: Tc = {tc:.6e}, LH = {likelihood:.4}");

  Ok(OptimizeTcResult {
    tc,
    likelihood,
    success,
  })
}

/// Evaluates the total coalescent log-likelihood for a constant Tc, reusing the
/// shared per-edge model so the reported likelihood matches [`compute_coalescent_total_lh`].
///
/// [`compute_coalescent_total_lh`]: crate::coalescent::total_lh::compute_coalescent_total_lh
fn compute_total_lh(precomputed: &CoalescentPrecomputed, edges: &[CoalescentEdgeData], tc: f64) -> Result<f64, Report> {
  let tc_dist = Distribution::constant(tc);
  let model = CoalescentModel::new(precomputed, &tc_dist)?;
  coalescent_log_likelihood(edges, &model)
}

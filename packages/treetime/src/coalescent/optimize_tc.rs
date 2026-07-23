use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
use crate::payload::traits::TimetreeNode;
use eyre::Report;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};

/// Result of Tc optimization.
pub struct OptimizeTcResult {
  /// Optimized coalescence time scale.
  pub tc: f64,
  /// Total coalescent likelihood at optimized Tc.
  pub likelihood: f64,
}

/// Computes the constant coalescence time scale Tc that maximizes the coalescent
/// likelihood.
///
/// A constant Tc is exactly a one-segment skyline, so this is a thin wrapper over
/// [`optimize_skyline`] with `n_points = 1`. For a single segment the smoothing
/// penalty is inert and the skyline's decoupled per-segment optimum reduces to the
/// closed form `Tc* = I / M` — the maximizer of the Kingman likelihood, where
/// `I = ∫ k(k-1)/2 dt` is the pairwise-merger-rate integral over the tree and `M`
/// is the total number of merger events.
///
/// See [`optimize_skyline`] for the shared machinery (per-edge `I`/`M` accumulation,
/// self-consistent likelihood reporting) and the degenerate-tree error returned when
/// the tree has no time span or no mergers.
pub fn optimize_tc<N, E, D>(graph: &Graph<N, E, D>) -> Result<OptimizeTcResult, Report>
where
  N: GraphNode + TimetreeNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  let result = optimize_skyline(graph, &SkylineParams { n_points: 1, ..SkylineParams::default() })?;
  Ok(OptimizeTcResult {
    tc: result.tc_values[0],
    likelihood: result.log_likelihood,
  })
}

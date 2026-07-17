use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::edge_data::{CoalescentEdgeData, coalescent_log_likelihood, collect_coalescent_edges};
use crate::coalescent::precomputed::CoalescentPrecomputed;
use crate::optimize::observer::OptimizationObserver;
use crate::payload::traits::TimetreeNode;
use crate::{make_error, make_report};
use argmin::core::observers::ObserverMode;
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::neldermead::NelderMead;
use eyre::Report;
use log::info;
use ndarray::Array1;
use std::sync::Arc;
use treetime_distribution::{Distribution, DistributionFunction};
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_grid::grid::Grid;

/// Parameters for skyline optimization.
#[derive(Debug, Clone)]
pub struct SkylineParams {
  /// Number of grid points for the piecewise linear Tc(t).
  pub n_points: usize,
  /// Penalty for rapid changes in log(Tc).
  pub stiffness: f64,
  /// Penalty for log(Tc) outside [-100, 0].
  pub regularization: f64,
  /// Optimization tolerance.
  pub tolerance: f64,
  /// Maximum iterations.
  pub max_iter: u64,
}

impl Default for SkylineParams {
  fn default() -> Self {
    Self {
      n_points: 10,
      stiffness: 2.0,
      regularization: 10.0,
      tolerance: 0.03,
      max_iter: 1000,
    }
  }
}

/// Result of skyline optimization.
#[derive(Debug, Clone)]
pub struct SkylineResult {
  /// Optimized Tc distribution (piecewise linear).
  pub tc_distribution: Distribution,
  /// Time grid points.
  pub time_grid: Array1<f64>,
  /// Optimized log(Tc) values at grid points.
  pub log_tc_values: Array1<f64>,
  /// Final coalescent log-likelihood.
  pub log_likelihood: f64,
}

/// Optimizes the skyline coalescent model to find the best Tc(t) trajectory.
///
/// This estimates a piecewise linear Tc(t) history that maximizes coalescent
/// likelihood with smoothness regularization.
///
/// The cost function minimizes:
/// - `-total_coalescent_LH` (negative log likelihood)
/// - `+ stiffness * sum(diff(logTc)^2)` (smoothness penalty)
/// - `+ regularization * boundary_penalty` (keep logTc in [-100, 0])
pub fn optimize_skyline<N, E, D>(graph: &Graph<N, E, D>, params: &SkylineParams) -> Result<SkylineResult, Report>
where
  N: GraphNode + TimetreeNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  debug_assert!(
    params.n_points >= 2,
    "skyline requires at least 2 grid points, got {}",
    params.n_points
  );

  info!(
    "Starting skyline optimization with {} grid points, stiffness={}, regularization={}",
    params.n_points, params.stiffness, params.regularization
  );

  let pre = CoalescentPrecomputed::from_graph(graph)?;

  // Create time grid spanning the tree event range
  let breakpoints = pre.lineage_counts().breakpoints();
  if breakpoints.len() < 2 {
    return make_error!(
      "Skyline optimization requires at least 2 breakpoints, got {}",
      breakpoints.len()
    );
  }
  let t_min = breakpoints[0];
  let t_max = breakpoints[breakpoints.len() - 1];

  let time_grid = Grid::from_range_n_points(t_min, t_max, params.n_points)?.to_array();

  // Initial guess: constant log(Tc) = 0 (Tc = 1)
  let initial_log_tc = Array1::zeros(params.n_points);

  // Create cost function
  let cost_fn = SkylineCostFunction {
    edges: Arc::from(collect_coalescent_edges(graph)?),
    precomputed: Arc::new(pre),
    time_grid: time_grid.clone(),
    stiffness: params.stiffness,
    regularization: params.regularization,
  };

  // Set up Nelder-Mead solver with initial simplex
  let simplex = create_initial_simplex(&initial_log_tc);
  let solver = NelderMead::new(simplex);

  // Run optimization
  let result = Executor::new(&cost_fn, solver)
    .configure(|cfg| cfg.max_iters(params.max_iter).target_cost(params.tolerance))
    .add_observer(
      OptimizationObserver {
        label: "Skyline",
        early_threshold: 5,
      },
      ObserverMode::Always,
    )
    .run()
    .map_err(|e| make_report!("Skyline optimization failed: {e}"))?;

  let opt_log_tc = result
    .state
    .best_param
    .ok_or_else(|| make_report!("Skyline optimization returned no result"))?;
  let opt_log_tc = Array1::from_vec(opt_log_tc);

  let opt_log_tc_clamped = Array1::from_iter(opt_log_tc.iter().map(|&x| x.clamp(-200.0, 100.0)));
  let tc_distribution = build_tc_distribution(&time_grid, &opt_log_tc_clamped)?;
  let model = CoalescentModel::new(&cost_fn.precomputed, &tc_distribution)?;
  let log_likelihood = coalescent_log_likelihood(&cost_fn.edges, &model)?;

  info!(
    "Skyline optimization completed: final_cost={:.4}, log_likelihood={log_likelihood:.4}",
    result.state.best_cost
  );

  Ok(SkylineResult {
    tc_distribution,
    time_grid,
    log_tc_values: opt_log_tc_clamped,
    log_likelihood,
  })
}

/// Cost function for skyline optimization.
struct SkylineCostFunction {
  edges: Arc<[CoalescentEdgeData]>,
  precomputed: Arc<CoalescentPrecomputed>,
  time_grid: Array1<f64>,
  stiffness: f64,
  regularization: f64,
}

impl CostFunction for &SkylineCostFunction {
  type Param = Vec<f64>;
  type Output = f64;

  fn cost(&self, log_tc: &Self::Param) -> Result<Self::Output, Error> {
    // Clamp log_tc to avoid overflow/underflow
    let log_tc_clamped = Array1::from_iter(log_tc.iter().map(|&x| x.clamp(-200.0, 100.0)));

    // Build Tc distribution from log values
    let tc_dist = build_tc_distribution(&self.time_grid, &log_tc_clamped)
      .map_err(|e| Error::msg(format!("Failed to build Tc distribution: {e}")))?;

    let model = CoalescentModel::new(&self.precomputed, &tc_dist)
      .map_err(|e| Error::msg(format!("Failed to construct coalescent model: {e}")))?;
    let neg_log_lh = -coalescent_log_likelihood(&self.edges, &model)
      .map_err(|e| Error::msg(format!("Failed to compute coalescent likelihood: {e}")))?;

    // Add smoothness penalty: stiffness * sum(diff(logTc)^2)
    // Use clamped values so optimizer gradients match the cost surface
    let smoothness_penalty: f64 = log_tc_clamped
      .windows(2)
      .into_iter()
      .map(|w| (w[1] - w[0]).powi(2))
      .sum::<f64>()
      * self.stiffness;

    // Add boundary penalty for log_tc outside [-100, 0]
    let boundary_penalty = compute_boundary_penalty(log_tc, self.regularization);

    Ok(neg_log_lh + smoothness_penalty + boundary_penalty)
  }
}

/// Computes boundary penalty for log_tc values outside [-100, 0].
fn compute_boundary_penalty(log_tc: &[f64], regularization: f64) -> f64 {
  log_tc
    .iter()
    .map(|&x| {
      let upper_penalty = if x > 0.0 { x } else { 0.0 };
      let lower_penalty = if x < -100.0 { -x - 100.0 } else { 0.0 };
      (upper_penalty + lower_penalty) * regularization
    })
    .sum()
}

/// Builds a piecewise linear Tc distribution from time grid and log values.
///
/// `Tc(t)` is `exp(log_tc)` interpolated linearly between grid points, with
/// constant extrapolation outside the grid (clamped to the first/last value).
pub(crate) fn build_tc_distribution(time_grid: &Array1<f64>, log_tc: &Array1<f64>) -> Result<Distribution, Report> {
  let tc_values = log_tc.mapv(f64::exp);
  let function = DistributionFunction::from_range_values((time_grid[0], time_grid[time_grid.len() - 1]), tc_values)?;
  Ok(Distribution::Function(function))
}

/// Creates initial simplex for Nelder-Mead optimization.
fn create_initial_simplex(initial: &Array1<f64>) -> Vec<Vec<f64>> {
  let n = initial.len();
  let mut simplex = Vec::with_capacity(n + 1);

  // First vertex is the initial point
  simplex.push(initial.to_vec());

  // Other vertices are perturbations along each dimension.
  // Perturbation of 0.5 in log-space corresponds to ~1.65x multiplicative factor in Tc,
  // providing reasonable exploration without extreme values.
  for i in 0..n {
    let mut vertex = initial.to_vec();
    vertex[i] += 0.5;
    simplex.push(vertex);
  }

  simplex
}

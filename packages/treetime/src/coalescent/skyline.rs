use crate::coalescent::events::collect_tree_events;
use crate::coalescent::integration::compute_integral_merger_rate;
use crate::coalescent::lineage_dynamics::compute_lineage_count_distribution;
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
use treetime_distribution::{Distribution, DistributionFormula};
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_grid::piecewise_linear_fn::PiecewiseLinearFn;

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

  // Collect tree events and compute lineage counts
  let (present_time, events_calendar) = collect_tree_events(graph)?;
  let events_tbp: Vec<_> = events_calendar
    .iter()
    .map(|(t, delta)| (t.to_tbp(present_time), *delta))
    .collect();

  let lineage_counts = compute_lineage_count_distribution(&events_tbp)?;

  // Create time grid spanning the tree event range
  let breakpoints = lineage_counts.breakpoints();
  if breakpoints.len() < 2 {
    return make_error!(
      "Skyline optimization requires at least 2 breakpoints, got {}",
      breakpoints.len()
    );
  }
  let t_min = breakpoints[0];
  let t_max = breakpoints[breakpoints.len() - 1];

  let time_grid = Array1::linspace(t_min, t_max, params.n_points);

  // Initial guess: constant log(Tc) = 0 (Tc = 1)
  let initial_log_tc = Array1::zeros(params.n_points);

  // Create cost function
  let cost_fn = SkylineCostFunction {
    lineage_counts: Arc::new(lineage_counts),
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

  // Compute final likelihood using clamped values (consistent with cost function)
  let final_cost = result.state.best_cost;
  let opt_log_tc_clamped = Array1::from_iter(opt_log_tc.iter().map(|&x| x.clamp(-200.0, 100.0)));
  let smoothness_penalty: f64 = opt_log_tc_clamped
    .windows(2)
    .into_iter()
    .map(|w| (w[1] - w[0]).powi(2))
    .sum::<f64>()
    * params.stiffness;
  let boundary_penalty = compute_boundary_penalty(opt_log_tc.as_slice().unwrap(), params.regularization);
  let log_likelihood = -(final_cost - smoothness_penalty - boundary_penalty);

  info!("Skyline optimization completed: final_cost={final_cost:.4}, log_likelihood={log_likelihood:.4}");

  // Build piecewise constant Tc distribution
  let tc_distribution = build_tc_distribution(&time_grid, &opt_log_tc);

  Ok(SkylineResult {
    tc_distribution,
    time_grid,
    log_tc_values: opt_log_tc,
    log_likelihood,
  })
}

/// Cost function for skyline optimization.
struct SkylineCostFunction {
  lineage_counts: Arc<PiecewiseConstantFn>,
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
    let tc_dist = build_tc_distribution(&self.time_grid, &log_tc_clamped);

    // Compute integral merger rate and total coalescent likelihood
    let integral_merger_rate = compute_integral_merger_rate(&tc_dist, &self.lineage_counts)
      .map_err(|e| Error::msg(format!("Failed to compute integral merger rate: {e}")))?;

    // Compute negative log-likelihood (sum over all breakpoints)
    let neg_log_lh = compute_total_neg_log_lh(&integral_merger_rate, &self.lineage_counts, &tc_dist);

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
pub(crate) fn build_tc_distribution(time_grid: &Array1<f64>, log_tc: &Array1<f64>) -> Distribution {
  let tc_values = log_tc.mapv(f64::exp);

  let t_min = time_grid[0];
  let t_max = time_grid[time_grid.len() - 1];

  let pwl = PiecewiseLinearFn::new(time_grid.clone(), tc_values);

  Distribution::Formula(DistributionFormula::new(move |t| Ok(pwl.eval(t)), t_min, t_max))
}

/// Computes total negative log-likelihood from the coalescent model.
///
/// This sums contributions from all tree segments based on the integral merger rate.
fn compute_total_neg_log_lh(
  integral_merger_rate: &PiecewiseLinearFn,
  lineage_counts: &PiecewiseConstantFn,
  tc_dist: &Distribution,
) -> f64 {
  let breakpoints = lineage_counts.breakpoints();
  let n = breakpoints.len();

  // Sum contributions from each segment
  let mut neg_log_lh = 0.0;

  for i in 1..n {
    let t0 = breakpoints[i - 1];
    let t1 = breakpoints[i];

    // Integral contribution: I(t1) - I(t0)
    let i0 = integral_merger_rate.eval(t0);
    let i1 = integral_merger_rate.eval(t1);
    neg_log_lh += i1 - i0;

    // Merger event contribution at t1 (if it's a merger, not a sample)
    let k = lineage_counts.eval(t1);
    if k < lineage_counts.eval(t0) {
      // This is a merger event - compute total merger rate (Kingman coalescent)
      // Rate proportional to k*(k-1)/2 pairs of lineages
      if let Ok(tc) = tc_dist.eval(t1) {
        let nlineages = f64::max(0.5, k - 1.0);
        let lambda = 0.5 * nlineages * (nlineages + 1.0) / tc;
        neg_log_lh -= lambda.ln();
      }
    }
  }

  neg_log_lh
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

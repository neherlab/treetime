use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::edge_data::{CoalescentEdgeData, coalescent_log_likelihood, collect_coalescent_edges};
use crate::coalescent::precomputed::CoalescentPrecomputed;
use crate::make_error;
use crate::payload::traits::TimetreeNode;
use eyre::Report;
use log::{info, warn};
use ndarray::Array1;
use ordered_float::OrderedFloat;
use treetime_distribution::{Distribution, DistributionFormula};
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;

/// Parameters for skyline (piecewise-constant Tc) optimization.
#[derive(Debug, Clone)]
pub struct SkylineParams {
  /// Number of piecewise-constant Tc segments over the tree's time span.
  pub n_points: usize,
  /// Smoothing penalty on adjacent inverse time scales: `stiffness * Σ (xᵢ₊₁ - xᵢ)²`
  /// where `xᵢ = 1/Tc_i`. Must be positive for more than one segment.
  pub stiffness: f64,
  /// Newton convergence tolerance on the gradient infinity-norm.
  pub tolerance: f64,
  /// Maximum Newton iterations.
  pub max_iter: u64,
}

impl Default for SkylineParams {
  fn default() -> Self {
    Self {
      n_points: 10,
      stiffness: 2.0,
      tolerance: 1e-8,
      max_iter: 100,
    }
  }
}

/// Result of skyline optimization.
#[derive(Debug, Clone)]
pub struct SkylineResult {
  /// Optimized piecewise-constant Tc(t).
  pub tc_distribution: Distribution,
  /// Segment boundaries in calendar time (length `n_points + 1`, ascending).
  pub segment_boundaries: Array1<f64>,
  /// Optimized Tc value per segment (length `n_points`).
  pub tc_values: Array1<f64>,
  /// Coalescent log-likelihood at the optimized Tc(t).
  pub log_likelihood: f64,
}

/// Optimizes a piecewise-constant Tc(t) trajectory that maximizes the coalescent
/// likelihood with a smoothness penalty.
///
/// # Algorithm
///
/// The tree's time span is split into `n_points` segments whose boundaries are the
/// quantiles of the merger (coalescence) times, so every segment — including the
/// two boundary segments spanning to the root and the tips — contains mergers.
/// Within segment `i`, Tc is constant; writing the coalescent rate `xᵢ = 1/Tc_i`,
/// the negative log-likelihood plus penalty is
///
/// ```text
///   C(x) = Σ_i (Iᵢ xᵢ - Mᵢ ln xᵢ) + (stiffness/2) Σ_i (xᵢ₊₁ - xᵢ)²
/// ```
///
/// where `Iᵢ = ∫_seg_i k(k-1)/2 dt` (Tc-independent pairwise-rate integral) and
/// `Mᵢ` is the merger count in segment `i`. In the `x` parametrization every term
/// is convex, so `C` has a unique positive minimizer. It is found with Newton's
/// method on the symmetric tridiagonal Hessian, warm-started from the decoupled
/// per-segment optimum `xᵢ = Mᵢ / Iᵢ` and damped to keep every `xᵢ > 0`. Because
/// each segment owns at least one merger, `Mᵢ > 0` keeps the `-Mᵢ ln xᵢ` barrier
/// active, so no segment collapses to `Tc → ∞`.
///
/// `Iᵢ` and `Mᵢ` are attributed to segments using the same interval-midpoint and
/// node-time conventions as [`CoalescentModel`], so the analytic optimum coincides
/// with the maximizer of the model-evaluated likelihood.
pub fn optimize_skyline<N, E, D>(graph: &Graph<N, E, D>, params: &SkylineParams) -> Result<SkylineResult, Report>
where
  N: GraphNode + TimetreeNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  if params.n_points < 1 {
    return make_error!(
      "Skyline optimization requires at least 1 segment, got {}",
      params.n_points
    );
  }

  info!(
    "Starting skyline optimization with {} segments, stiffness={}",
    params.n_points, params.stiffness
  );

  let precomputed = CoalescentPrecomputed::from_graph(graph)?;
  let edges = collect_coalescent_edges(graph)?;

  let breakpoints = precomputed.lineage_counts().breakpoints();
  if breakpoints.len() < 2 {
    return make_error!(
      "Skyline optimization requires at least 2 breakpoints, got {}",
      breakpoints.len()
    );
  }
  let t_min = breakpoints[0];
  let t_max = breakpoints[breakpoints.len() - 1];

  // Merger times drive the boundary quantiles; a non-finite time would make the
  // grid meaningless and is rejected here rather than propagated into the sort.
  if let Some(bad) = edges
    .iter()
    .map(|edge| edge.parent_time().value())
    .find(|t| !t.is_finite())
  {
    return make_error!("Coalescent merger time must be finite, got {bad}");
  }

  let boundaries = merger_quantile_boundaries(&edges, t_min, t_max, params.n_points);

  let (i_seg, m_seg) = accumulate_segment_terms(precomputed.lineage_counts(), &edges, &boundaries);
  let x = solve_inverse_tc(&i_seg, &m_seg, params.stiffness, params.tolerance, params.max_iter)?;

  let tc_values = Array1::from_iter(x.iter().map(|&xi| 1.0 / xi));
  let tc_distribution = build_tc_distribution(&boundaries, &tc_values);

  // Report the likelihood via the shared model so it matches `compute_coalescent_total_lh`.
  let model = CoalescentModel::new(&precomputed, &tc_distribution)?;
  let log_likelihood = coalescent_log_likelihood(&edges, &model)?;

  info!("Skyline optimization completed: log_likelihood={log_likelihood:.4}");
  info!("Skyline Tc(t) trajectory ({} segments):", tc_values.len());
  for (i, &tc) in tc_values.iter().enumerate() {
    info!(
      "  segment {i}: [{:.4}, {:.4}]  Tc = {tc:.6e}",
      boundaries[i],
      boundaries[i + 1]
    );
  }

  Ok(SkylineResult {
    tc_distribution,
    segment_boundaries: Array1::from(boundaries),
    tc_values,
    log_likelihood,
  })
}

/// Segment index containing calendar time `t`, clamped to `[0, n_seg - 1]`.
///
/// `boundaries` is ascending with length `n_seg + 1`; segment `i` covers
/// `[boundaries[i], boundaries[i + 1])`.
fn segment_index(boundaries: &[f64], t: f64) -> usize {
  let n_seg = boundaries.len() - 1;
  let above = boundaries.partition_point(|&b| b <= t);
  above.saturating_sub(1).min(n_seg - 1)
}

/// Computes `n_seg + 1` ascending segment boundaries at quantiles of the merger
/// times, with the outer boundaries pinned to `t_min` and `t_max`.
///
/// Interior boundaries fall between merger events so that each segment owns roughly
/// the same number of mergers; the first and last segments therefore contain the
/// oldest and youngest mergers rather than an empty sampling tail.
fn merger_quantile_boundaries(edges: &[CoalescentEdgeData], t_min: f64, t_max: f64, n_seg: usize) -> Vec<f64> {
  if n_seg <= 1 {
    return vec![t_min, t_max];
  }

  let mut times: Vec<f64> = edges.iter().map(|e| e.parent_time().value()).collect();
  times.sort_by_key(|&t| OrderedFloat(t));

  let n = times.len();
  let mut boundaries = Vec::with_capacity(n_seg + 1);
  boundaries.push(t_min);
  for k in 1..n_seg {
    // Split between the (k·n/n_seg)-th sorted merger and its predecessor.
    let pos = (k * n) / n_seg;
    let candidate = if n == 0 {
      t_min + (t_max - t_min) * (k as f64 / n_seg as f64)
    } else {
      let lo = times[pos.saturating_sub(1).min(n - 1)];
      let hi = times[pos.min(n - 1)];
      if hi > lo { f64::midpoint(lo, hi) } else { lo }
    };
    // Keep boundaries strictly inside (t_min, t_max) and non-decreasing.
    let prev = *boundaries.last().unwrap();
    boundaries.push(candidate.clamp(prev, t_max));
  }
  boundaries.push(t_max);
  boundaries
}

/// Accumulates the per-segment pairwise-rate integral `Iᵢ` and merger count `Mᵢ`.
///
/// `Iᵢ` sums, over lineage-count intervals whose midpoint falls in segment `i`, the
/// interval's per-lineage merger integral (at Tc = 1) times the number of collected
/// edges covering it — matching the model's per-edge survival term. `Mᵢ` sums
/// `(n_siblings - 1)/n_siblings` over edges whose parent (merger) time lies in
/// segment `i`.
fn accumulate_segment_terms(
  lineage_counts: &PiecewiseConstantFn,
  edges: &[CoalescentEdgeData],
  boundaries: &[f64],
) -> (Vec<f64>, Vec<f64>) {
  let n_seg = boundaries.len() - 1;
  let breakpoints = lineage_counts.breakpoints();
  let n_int = breakpoints.len() - 1;

  // Per-interval midpoint and per-lineage integral (Tc = 1).
  let mids: Vec<f64> = (0..n_int)
    .map(|j| f64::midpoint(breakpoints[j], breakpoints[j + 1]))
    .collect();
  let rate: Vec<f64> = (0..n_int)
    .map(|j| {
      let dt = breakpoints[j + 1] - breakpoints[j];
      let k = lineage_counts.eval(mids[j]);
      dt * 0.5 * f64::max(0.5, k - 1.0)
    })
    .collect();

  // Number of collected edges covering each interval, via a difference array.
  let mut coverage = vec![0_i64; n_int + 1];
  for edge in edges {
    let parent_time = edge.parent_time().value();
    let child_time = edge.child_time().value();
    // Intervals covered are those whose midpoint lies within the edge's span.
    let lo = mids.partition_point(|&m| m < parent_time);
    let hi = mids.partition_point(|&m| m < child_time);
    coverage[lo] += 1;
    coverage[hi] -= 1;
  }
  let mut running = 0_i64;
  let mut i_seg = vec![0.0; n_seg];
  for j in 0..n_int {
    running += coverage[j];
    i_seg[segment_index(boundaries, mids[j])] += running as f64 * rate[j];
  }

  let mut m_seg = vec![0.0; n_seg];
  for edge in edges {
    let n_siblings = edge.n_siblings();
    m_seg[segment_index(boundaries, edge.parent_time().value())] += (n_siblings - 1.0) / n_siblings;
  }

  for i in 0..n_seg {
    info!("Skyline segment {i}: I = {:.6e}, M = {:.6e}", i_seg[i], m_seg[i]);
  }

  (i_seg, m_seg)
}

/// Minimizes `C(x) = Σ (Iᵢ xᵢ - Mᵢ ln xᵢ) + (γ/2) Σ (xᵢ₊₁ - xᵢ)²` over `xᵢ > 0`.
///
/// Convex in `x`, solved by damped Newton on the symmetric tridiagonal Hessian.
fn solve_inverse_tc(
  i_seg: &[f64],
  m_seg: &[f64],
  stiffness: f64,
  tolerance: f64,
  max_iter: u64,
) -> Result<Vec<f64>, Report> {
  let n = i_seg.len();

  // Decoupled per-segment optimum, with a pooled fallback for empty segments.
  let i_tot: f64 = i_seg.iter().sum();
  let m_tot: f64 = m_seg.iter().sum();
  let x_pooled = if i_tot > 0.0 && m_tot > 0.0 { m_tot / i_tot } else { 1.0 };
  let mut x: Vec<f64> = (0..n)
    .map(|k| {
      let xk = m_seg[k] / i_seg[k];
      if xk.is_finite() && xk > 0.0 { xk } else { x_pooled }
    })
    .collect();

  // Single segment: the decoupled solution is already optimal, no smoothing applies.
  if n == 1 {
    return Ok(x);
  }

  // Multi-segment smoothing requires well-formed control parameters. A non-positive
  // or non-finite stiffness was previously treated as "no smoothing" silently,
  // masking misconfiguration; the parameter is documented as positive.
  if !(stiffness.is_finite() && stiffness > 0.0) {
    return make_error!(
      "Skyline smoothing stiffness must be finite and positive for a multi-segment fit, got {stiffness}"
    );
  }
  if !(tolerance.is_finite() && tolerance > 0.0) {
    return make_error!("Skyline Newton tolerance must be finite and positive, got {tolerance}");
  }
  if max_iter == 0 {
    return make_error!("Skyline Newton max_iter must be at least 1");
  }

  let mut converged = false;
  let mut last_g_norm = f64::INFINITY;
  for _ in 0..max_iter {
    // Gradient and tridiagonal Hessian (diagonal `diag`, off-diagonal `off = -γ`).
    let mut g = vec![0.0; n];
    let mut diag = vec![0.0; n];
    for i in 0..n {
      g[i] = i_seg[i] - m_seg[i] / x[i];
      diag[i] = m_seg[i] / (x[i] * x[i]);
    }
    let mut off = vec![0.0; n - 1];
    for i in 0..n - 1 {
      let d = x[i] - x[i + 1];
      g[i] += stiffness * d;
      g[i + 1] -= stiffness * d;
      diag[i] += stiffness;
      diag[i + 1] += stiffness;
      off[i] = -stiffness;
    }

    let g_norm = g.iter().fold(0.0_f64, |acc, &v| acc.max(v.abs()));
    last_g_norm = g_norm;
    if g_norm < tolerance {
      converged = true;
      break;
    }

    let neg_g: Vec<f64> = g.iter().map(|&v| -v).collect();
    let dx = solve_symmetric_tridiagonal(&diag, &off, &neg_g);

    // Cap the step so every xᵢ stays positive, then backtrack for cost decrease.
    let mut alpha: f64 = 1.0;
    for i in 0..n {
      if dx[i] < 0.0 {
        alpha = alpha.min(-0.99 * x[i] / dx[i]);
      }
    }
    let c0 = skyline_cost(&x, i_seg, m_seg, stiffness);
    let slope: f64 = g.iter().zip(&dx).map(|(&gi, &di)| gi * di).sum();
    loop {
      let x_new: Vec<f64> = (0..n).map(|i| x[i] + alpha * dx[i]).collect();
      if x_new.iter().all(|&xi| xi > 0.0) && skyline_cost(&x_new, i_seg, m_seg, stiffness) <= c0 + 1e-4 * alpha * slope
      {
        x = x_new;
        break;
      }
      alpha *= 0.5;
      if alpha < 1e-12 {
        break;
      }
    }
  }

  // A non-converged fit was previously indistinguishable from a converged one.
  // Surface it so a returned iterate that failed the gradient tolerance is visible.
  if !converged {
    warn!(
      "Skyline Newton did not reach gradient tolerance {tolerance:.3e} within {max_iter} iterations \
       (final gradient infinity-norm {last_g_norm:.3e}); returning the last iterate"
    );
  }

  Ok(x)
}

/// Value of the skyline objective `C(x)` (constants dropped).
fn skyline_cost(x: &[f64], i_seg: &[f64], m_seg: &[f64], stiffness: f64) -> f64 {
  let data: f64 = (0..x.len()).map(|i| i_seg[i] * x[i] - m_seg[i] * x[i].ln()).sum();
  let penalty: f64 = x.windows(2).map(|w| (w[1] - w[0]).powi(2)).sum::<f64>() * 0.5 * stiffness;
  data + penalty
}

/// Solves a symmetric tridiagonal system `A x = rhs` via the Thomas algorithm.
///
/// `diag` is the main diagonal (length `n`); `off` holds the shared sub/super
/// diagonal (length `n - 1`). `A` is positive-definite here, so no pivoting.
fn solve_symmetric_tridiagonal(diag: &[f64], off: &[f64], rhs: &[f64]) -> Vec<f64> {
  let n = diag.len();
  if n == 1 {
    return vec![rhs[0] / diag[0]];
  }
  let mut c = vec![0.0; n];
  let mut d = vec![0.0; n];
  c[0] = off[0] / diag[0];
  d[0] = rhs[0] / diag[0];
  for i in 1..n {
    let denom = diag[i] - off[i - 1] * c[i - 1];
    if i < n - 1 {
      c[i] = off[i] / denom;
    }
    d[i] = (rhs[i] - off[i - 1] * d[i - 1]) / denom;
  }
  let mut x = vec![0.0; n];
  x[n - 1] = d[n - 1];
  for i in (0..n - 1).rev() {
    x[i] = d[i] - c[i] * x[i + 1];
  }
  x
}

/// Builds a piecewise-constant Tc(t) distribution from per-segment Tc values.
///
/// Uses the same segment lookup as the optimizer, so model evaluation reproduces
/// the optimized per-segment rates. Times outside the grid clamp to the first/last
/// segment.
fn build_tc_distribution(boundaries: &[f64], tc_values: &Array1<f64>) -> Distribution {
  let t_min = boundaries[0];
  let t_max = boundaries[boundaries.len() - 1];
  let boundaries = boundaries.to_vec();
  let values = tc_values.to_vec();
  Distribution::Formula(DistributionFormula::new(
    move |t| Ok(values[segment_index(&boundaries, t)]),
    t_min,
    t_max,
  ))
}

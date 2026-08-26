use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::edge_data::{CoalescentEdgeData, coalescent_log_likelihood, collect_coalescent_edges};
use crate::coalescent::lineage_counts::compute_lineage_counts;
use crate::make_error;
use crate::payload::traits::TimetreeNode;
use eyre::{Report, WrapErr};
use log::{info, warn};
use ndarray::{Array1, Array2, array};
use ndarray_linalg::layout::MatrixLayout;
use ndarray_linalg::{SolveTridiagonal, Tridiagonal};
use treetime_distribution::{Distribution, DistributionFormula};
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_primitives::LogLh;
use treetime_utils::array::ndarray::exp;

/// Parameters for skyline (piecewise-constant Tc) optimization.
#[derive(Debug, Clone)]
pub struct SkylineParams {
  /// Number of piecewise-constant Tc segments over the tree's time span.
  pub n_points: usize,
  /// Smoothing penalty on adjacent log time scales: `(stiffness/2) * Σ (zᵢ₊₁ - zᵢ)²`
  /// where `zᵢ = ln Tc_i`, i.e. it penalizes squared log-fold-changes of Tc and is
  /// therefore scale-independent. Must be positive for more than one segment.
  pub stiffness: f64,
  /// Newton convergence tolerance on the gradient infinity-norm.
  pub tolerance: f64,
  /// Maximum Newton iterations.
  pub max_iter: u64,
  /// Confidence level, in standard deviations, for the Tc(t) bands reported alongside
  /// the optimized time scale. The band spans `Tc_i * exp(±n_std * σ_i)`, where `σ_i`
  /// is the standard deviation of `ln Tc_i` recovered from the optimization curvature.
  pub n_std: f64,
}

impl Default for SkylineParams {
  fn default() -> Self {
    Self {
      n_points: 20,
      stiffness: 2.0,
      tolerance: 1e-8,
      max_iter: 100,
      n_std: 2.0,
    }
  }
}

/// Result of skyline optimization.
#[derive(Debug, Clone)]
pub struct SkylineResult {
  /// Optimized piecewise-constant Tc(t).
  pub tc_distribution: Distribution,
  /// Exact step schedule for consumers that integrate across every Tc change.
  pub tc_schedule: PiecewiseConstantFn,
  /// Segment boundaries in calendar time (length `n_points + 1`, ascending).
  pub segment_boundaries: Array1<f64>,
  /// Optimized Tc value per segment (length `n_points`).
  pub tc_values: Array1<f64>,
  /// Diagonal of the inverse Hessian in `ln Tc` coordinates.
  pub log_tc_variances: Array1<f64>,
  /// Lower confidence bound per segment.
  pub tc_lower_bounds: Array1<f64>,
  /// Upper confidence bound per segment.
  pub tc_upper_bounds: Array1<f64>,
  /// Coalescent log-likelihood at the optimized Tc(t).
  pub log_likelihood: LogLh,
}

/// Optimizes a piecewise-constant Tc(t) trajectory that maximizes the coalescent
/// likelihood with a smoothness penalty.
///
/// # Algorithm
///
/// The tree's time span is split into `n_points` equal-width segments. Within
/// segment `i`, Tc is constant; writing `zᵢ = ln Tc_i` (so the coalescent
/// rate is `1/Tc_i = e^{-zᵢ}`), the negative log-likelihood plus penalty is
///
/// ```text
///   C(z) = Σ_i (Iᵢ e^{-zᵢ} + Mᵢ zᵢ) + (stiffness/2) Σ_i (zᵢ₊₁ - zᵢ)²
/// ```
///
/// where `Iᵢ = ∫_seg_i k(k-1)/2 dt` (Tc-independent pairwise-rate integral) and
/// `Mᵢ` is the merger count in segment `i`. Modeling `z = ln Tc` makes the penalty
/// scale-independent — it charges squared *log-fold-changes* `zᵢ₊₁ - zᵢ =
/// ln(Tc_{i+1}/Tc_i)` — and guarantees `Tc = e^z > 0` with no constraint. Every
/// term is convex in `z`, so `C` has a unique minimizer, found with Newton's method
/// on the symmetric tridiagonal Hessian, warm-started from the decoupled per-segment
/// optimum `zᵢ = ln(Iᵢ / Mᵢ)` and globalized with an Armijo line search. A segment in
/// a merger-sparse region may own no mergers (`Mᵢ = 0`); then the linear `Mᵢ zᵢ` term
/// vanishes and `zᵢ` is pinned by the smoothing prior alone, so `stiffness > 0` is
/// required for more than one segment to keep the Hessian positive-definite and every
/// `zᵢ` finite (a lone data term `Iᵢ e^{-zᵢ}` would otherwise drive `Tc → ∞`).
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
  if !(params.n_std.is_finite() && params.n_std >= 0.0) {
    return make_error!(
      "Skyline confidence must be finite and nonnegative, got {}",
      params.n_std
    );
  }

  info!(
    "Starting skyline optimization with {} segments, stiffness={}",
    params.n_points, params.stiffness
  );

  let lineage_counts = compute_lineage_counts(graph)?;
  let edges = collect_coalescent_edges(graph)?;

  let breakpoints = lineage_counts.breakpoints();
  if breakpoints.len() < 2 {
    return make_error!(
      "Skyline optimization requires at least 2 breakpoints, got {}",
      breakpoints.len()
    );
  }
  let t_min = breakpoints[0];
  let t_max = breakpoints[breakpoints.len() - 1];

  let boundaries = equal_width_boundaries(t_min, t_max, params.n_points);

  let (i_seg, m_seg) = accumulate_segment_terms(&lineage_counts, &edges, &boundaries);

  // The whole-tree pairwise-rate integral and merger count must be positive and
  // finite for a coalescent Tc to exist. A tree with no time span or no internal
  // mergers is degenerate for the coalescent; fail loudly rather than letting the
  // per-segment pooled fallback silently substitute Tc = 1. (Per-segment emptiness
  // is still handled inside `solve_log_tc`, and merger-quantile boundaries keep
  // every segment non-empty in practice.)
  let i_tot: f64 = i_seg.iter().sum();
  let m_tot: f64 = m_seg.iter().sum();
  if !(m_tot > 0.0 && i_tot > 0.0 && i_tot.is_finite()) {
    return make_error!(
      "Cannot estimate a coalescent Tc: the tree is degenerate for the coalescent \
       (pairwise-rate integral = {i_tot:.6e}, mergers = {m_tot:.6e}). This means the \
       tree has effectively no time span or no internal mergers. Provide a fixed Tc via \
       --coalescent, or run without a coalescent prior."
    );
  }

  let (z, hessian) = solve_log_tc(&i_seg, &m_seg, params.stiffness, params.tolerance, params.max_iter)?;

  let tc_values = Array1::from_iter(z.iter().map(|&zi| zi.exp()));
  let confidence = skyline_confidence_band(&hessian, params.n_std, &tc_values)?;
  let tc_distribution = build_tc_distribution(&boundaries, &tc_values);
  // Internal segment boundaries are the exact discontinuities of the clamped Tc schedule.
  let tc_schedule = PiecewiseConstantFn::new(
    Array1::from(boundaries[1..boundaries.len() - 1].to_vec()),
    tc_values.clone(),
  );

  // Report the likelihood via the shared model so it matches `compute_coalescent_total_lh`.
  let model = CoalescentModel::new(&lineage_counts, &tc_distribution)?;
  let log_likelihood = coalescent_log_likelihood(&edges, &model)?;

  info!(
    "Skyline optimization completed: log_likelihood={:.4}",
    log_likelihood.value()
  );
  info!("Skyline Tc(t) trajectory ({} segments):", tc_values.len());
  for (i, &tc) in tc_values.iter().enumerate() {
    info!(
      "  segment {i}: [{:.4}, {:.4}]  Tc = {tc:.6e} [{:.6e}, {:.6e}]",
      boundaries[i],
      boundaries[i + 1],
      confidence.tc_lower_bounds[i],
      confidence.tc_upper_bounds[i]
    );
  }

  Ok(SkylineResult {
    tc_distribution,
    tc_schedule,
    segment_boundaries: Array1::from(boundaries),
    tc_values,
    log_tc_variances: confidence.log_tc_variances,
    tc_lower_bounds: confidence.tc_lower_bounds,
    tc_upper_bounds: confidence.tc_upper_bounds,
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

/// Computes `n_seg + 1` equally spaced ascending segment boundaries spanning
/// `[t_min, t_max]`.
///
/// Uniform widths give the smoothing penalty `Σ (zᵢ₊₁ - zᵢ)²` a clean, grid-
/// independent meaning — a consistent discretization of the squared log-Tc gradient
/// in time — so the stiffness has a well-defined scale. Unlike merger-quantile
/// boundaries this does not guarantee every segment owns a merger: segments in
/// merger-sparse regions can be empty (`Mᵢ = 0`) and are then pinned by the
/// smoothing prior, which is why `stiffness > 0` is required for `n_seg > 1`.
fn equal_width_boundaries(t_min: f64, t_max: f64, n_seg: usize) -> Vec<f64> {
  let n_seg = n_seg.max(1);
  let mut boundaries: Vec<f64> = (0..=n_seg)
    .map(|k| t_min + (t_max - t_min) * (k as f64 / n_seg as f64))
    .collect();
  // Pin the endpoints exactly, guarding against floating-point drift at the edges.
  boundaries[0] = t_min;
  boundaries[n_seg] = t_max;
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

/// Minimizes `C(z) = Σ (Iᵢ e^{-zᵢ} + Mᵢ zᵢ) + (γ/2) Σ (zᵢ₊₁ - zᵢ)²` over `zᵢ = ln Tc_i`.
///
/// Convex in `z`, solved by Newton on the symmetric tridiagonal Hessian with an
/// Armijo line search. `Tc = e^z` is positive by construction, so no step capping.
///
/// Returns the optimum `z` and the Hessian evaluated at that `z`. The confidence band
/// is the inverse of this same operator, so returning it lets the band reuse the final
/// Newton Hessian instead of rebuilding it.
fn solve_log_tc(
  i_seg: &[f64],
  m_seg: &[f64],
  stiffness: f64,
  tolerance: f64,
  max_iter: u64,
) -> Result<(Vec<f64>, Tridiagonal<f64>), Report> {
  let n = i_seg.len();

  // Decoupled per-segment optimum zᵢ = ln(Iᵢ / Mᵢ), with a pooled fallback for
  // empty/degenerate segments (where Iᵢ or Mᵢ is zero and the ratio's log is
  // non-finite).
  let i_tot: f64 = i_seg.iter().sum();
  let m_tot: f64 = m_seg.iter().sum();
  let z_pooled = if i_tot > 0.0 && m_tot > 0.0 {
    (i_tot / m_tot).ln()
  } else {
    0.0
  };
  let mut z: Vec<f64> = (0..n)
    .map(|k| {
      let zk = (i_seg[k] / m_seg[k]).ln();
      if zk.is_finite() { zk } else { z_pooled }
    })
    .collect();

  // Single segment or no smoothing: the decoupled solution is already optimal.
  if n == 1 {
    let hessian = skyline_hessian(&z, i_seg, stiffness)?;
    return Ok((z, hessian));
  }
  // error if stiffness is non-positive, which would make the Hessian indefinite.
  if stiffness <= 0.0 {
    return make_error!(
      "Skyline optimization requires positive stiffness for more than one segment, got {}",
      stiffness
    );
  }

  for _ in 0..max_iter {
    let g = skyline_gradient(&z, i_seg, m_seg, stiffness);
    let hessian = skyline_hessian(&z, i_seg, stiffness)?;

    let g_norm = g.iter().fold(0.0_f64, |acc, &v| acc.max(v.abs()));
    if g_norm < tolerance {
      // Converged: this Hessian is evaluated at the returned `z`, so hand it back
      // for the confidence band rather than rebuilding the same operator.
      return Ok((z, hessian));
    }

    // Solve the Hessian system with the gradient, then negate in place: the
    // Newton step is dz = -H⁻¹g.
    let mut dz = hessian
      .solve_tridiagonal(&g)
      .wrap_err("Failed to solve the skyline Hessian system")?;
    dz.mapv_inplace(|d| -d);

    // Armijo backtracking line search; Tc = e^z stays positive for any step.
    let mut alpha: f64 = 1.0;
    let c0 = skyline_cost(&z, i_seg, m_seg, stiffness);
    let slope: f64 = g.iter().zip(&dz).map(|(&gi, &di)| gi * di).sum();
    loop {
      let z_new: Vec<f64> = (0..n).map(|i| z[i] + alpha * dz[i]).collect();
      if skyline_cost(&z_new, i_seg, m_seg, stiffness) <= c0 + 1e-4 * alpha * slope {
        z = z_new;
        break;
      }
      alpha *= 0.5;
      if alpha < 1e-12 {
        break;
      }
    }
  }

  // Reached the iteration cap without meeting the gradient tolerance. The last step
  // advanced `z` past the loop's Hessian, so evaluate a fresh one at the final iterate.
  warn!("Skyline optimization did not converge within {max_iter} iterations");
  let hessian = skyline_hessian(&z, i_seg, stiffness)?;
  Ok((z, hessian))
}

/// Returns the objective gradient at `z`.
fn skyline_gradient(z: &[f64], i_seg: &[f64], m_seg: &[f64], stiffness: f64) -> Array1<f64> {
  let n = z.len();
  let mut gradient = Array1::from_iter((0..n).map(|i| -i_seg[i] * (-z[i]).exp() + m_seg[i]));

  for i in 0..n - 1 {
    let difference = z[i] - z[i + 1];
    gradient[i] += stiffness * difference;
    gradient[i + 1] -= stiffness * difference;
  }
  gradient
}

/// Returns the symmetric tridiagonal objective Hessian at `z`.
///
/// Exposed to the crate so tests can build a known Hessian and check the marginal
/// variance recovery against a hand-inverted oracle.
pub(crate) fn skyline_hessian(z: &[f64], i_seg: &[f64], stiffness: f64) -> Result<Tridiagonal<f64>, Report> {
  let n = z.len();
  let matrix_size = i32::try_from(n).wrap_err("Skyline segment count exceeds the linear algebra limit")?;
  let mut diagonal: Vec<f64> = (0..n).map(|i| i_seg[i] * (-z[i]).exp()).collect();
  for i in 0..n - 1 {
    diagonal[i] += stiffness;
    diagonal[i + 1] += stiffness;
  }
  let off_diagonal = vec![-stiffness; n - 1];

  Ok(Tridiagonal {
    l: MatrixLayout::C {
      row: matrix_size,
      lda: matrix_size,
    },
    dl: off_diagonal.clone(),
    d: diagonal,
    du: off_diagonal,
  })
}

/// Computes the local Gaussian confidence band from the Hessian at the optimum.
///
/// The inverse Hessian is the covariance of the Laplace approximation in `ln Tc`
/// coordinates. See https://doi.org/10.1080/01621459.1986.10478240.
fn skyline_confidence_band(
  hessian: &Tridiagonal<f64>,
  n_std: f64,
  tc_values: &Array1<f64>,
) -> Result<SkylineConfidenceBand, Report> {
  let log_tc_variances = marginal_log_tc_variances(hessian)?;
  let band_factors = exp(&(log_tc_variances.mapv(f64::sqrt) * n_std));
  Ok(SkylineConfidenceBand {
    log_tc_variances,
    tc_lower_bounds: tc_values / &band_factors,
    tc_upper_bounds: tc_values * &band_factors,
  })
}

/// Per-segment marginal variances of `ln Tc`: the diagonal of the inverse Hessian.
///
/// Inverting the full tridiagonal Hessian against the identity keeps the off-diagonal
/// stiffness coupling between adjacent segments, so each marginal variance is strictly
/// larger than the diagonal-only `1/H_ii` it would carry in isolation whenever the
/// smoothing couples the segments. Exposed to the crate for the analytic oracle test.
pub(crate) fn marginal_log_tc_variances(hessian: &Tridiagonal<f64>) -> Result<Array1<f64>, Report> {
  let n = hessian.d.len();
  let log_tc_variances = if n == 1 {
    // A 1x1 Hessian inverts to the scalar reciprocal.
    array![1.0 / hessian.d[0]]
  } else {
    hessian
      .solve_tridiagonal(&Array2::eye(n))
      .wrap_err("Failed to invert the skyline Hessian")?
      .diag()
      .to_owned()
  };
  if log_tc_variances
    .iter()
    .any(|variance| !variance.is_finite() || *variance <= 0.0)
  {
    return make_error!("Skyline Hessian inverse contains a nonpositive or non-finite variance");
  }
  Ok(log_tc_variances)
}

struct SkylineConfidenceBand {
  log_tc_variances: Array1<f64>,
  tc_lower_bounds: Array1<f64>,
  tc_upper_bounds: Array1<f64>,
}

/// Value of the skyline objective `C(z)` (constants dropped).
fn skyline_cost(z: &[f64], i_seg: &[f64], m_seg: &[f64], stiffness: f64) -> f64 {
  let data: f64 = (0..z.len()).map(|i| i_seg[i] * (-z[i]).exp() + m_seg[i] * z[i]).sum();
  let penalty: f64 = z.windows(2).map(|w| (w[1] - w[0]).powi(2)).sum::<f64>() * 0.5 * stiffness;
  data + penalty
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

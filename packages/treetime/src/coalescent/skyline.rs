use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::edge_data::{CoalescentEdgeData, coalescent_log_likelihood, collect_coalescent_edges};
use crate::coalescent::lineage_counts::compute_lineage_counts;
use crate::coalescent::node_time::CoalescentNodeTimes;
use crate::make_error;
use eyre::{Report, WrapErr};
use log::{info, warn};
use ndarray::{Array1, Array2, array};
use ndarray_linalg::layout::MatrixLayout;
use ndarray_linalg::{SolveTridiagonal, Tridiagonal};
use treetime_distribution::{Distribution, DistributionFormula};
use treetime_graph::graph::Graph;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_primitives::LogLh;
use treetime_utils::array::ndarray::exp;

pub(crate) fn optimize_skyline(
  graph: &Graph,
  params: &SkylineParams,
  node_times: &CoalescentNodeTimes,
) -> Result<SkylineResult, Report> {
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

  let lineage_counts = compute_lineage_counts(graph, node_times)?;
  let edges = collect_coalescent_edges(graph, node_times)?;

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
  let tc_schedule = PiecewiseConstantFn::new(
    Array1::from(boundaries[1..boundaries.len() - 1].to_vec()),
    tc_values.clone(),
  );

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

#[derive(Debug, Clone)]
pub struct SkylineParams {
  pub(crate) n_points: usize,
  pub(crate) stiffness: f64,
  pub(crate) tolerance: f64,
  pub(crate) max_iter: u64,
  pub(crate) n_std: f64,
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

#[derive(Debug, Clone)]
pub struct SkylineResult {
  pub(crate) tc_distribution: Distribution,
  pub(crate) tc_schedule: PiecewiseConstantFn,
  pub(crate) segment_boundaries: Array1<f64>,
  pub(crate) tc_values: Array1<f64>,
  pub(crate) log_tc_variances: Array1<f64>,
  pub(crate) tc_lower_bounds: Array1<f64>,
  pub(crate) tc_upper_bounds: Array1<f64>,
  pub(crate) log_likelihood: LogLh,
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn equal_width_boundaries(t_min: f64, t_max: f64, n_seg: usize) -> Vec<f64> {
  let n_seg = n_seg.max(1);
  let mut boundaries: Vec<f64> = (0..=n_seg)
    .map(|k| t_min + (t_max - t_min) * (k as f64 / n_seg as f64))
    .collect();
  boundaries[0] = t_min;
  boundaries[n_seg] = t_max;
  boundaries
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn accumulate_segment_terms(
  lineage_counts: &PiecewiseConstantFn,
  edges: &[CoalescentEdgeData],
  boundaries: &[f64],
) -> (Vec<f64>, Vec<f64>) {
  let n_seg = boundaries.len() - 1;
  let breakpoints = lineage_counts.breakpoints();
  let n_int = breakpoints.len() - 1;

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

  let mut coverage = vec![0_i64; n_int + 1];
  for edge in edges {
    let parent_time = edge.parent_time().value();
    let child_time = edge.child_time().value();
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

fn solve_log_tc(
  i_seg: &[f64],
  m_seg: &[f64],
  stiffness: f64,
  tolerance: f64,
  max_iter: u64,
) -> Result<(Vec<f64>, Tridiagonal<f64>), Report> {
  let n = i_seg.len();

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

  if n == 1 {
    let hessian = skyline_hessian(&z, i_seg, stiffness)?;
    return Ok((z, hessian));
  }
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
      return Ok((z, hessian));
    }

    let mut dz = hessian
      .solve_tridiagonal(&g)
      .wrap_err("Failed to solve the skyline Hessian system")?;
    dz.mapv_inplace(|d| -d);

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

  warn!("Skyline optimization did not converge within {max_iter} iterations");
  let hessian = skyline_hessian(&z, i_seg, stiffness)?;
  Ok((z, hessian))
}

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

pub(crate) fn marginal_log_tc_variances(hessian: &Tridiagonal<f64>) -> Result<Array1<f64>, Report> {
  let n = hessian.d.len();
  let log_tc_variances = if n == 1 {
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

fn skyline_cost(z: &[f64], i_seg: &[f64], m_seg: &[f64], stiffness: f64) -> f64 {
  let data: f64 = (0..z.len()).map(|i| i_seg[i] * (-z[i]).exp() + m_seg[i] * z[i]).sum();
  let penalty: f64 = z.windows(2).map(|w| (w[1] - w[0]).powi(2)).sum::<f64>() * 0.5 * stiffness;
  data + penalty
}

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

fn segment_index(boundaries: &[f64], t: f64) -> usize {
  let n_seg = boundaries.len() - 1;
  let above = boundaries.partition_point(|&b| b <= t);
  above.saturating_sub(1).min(n_seg - 1)
}

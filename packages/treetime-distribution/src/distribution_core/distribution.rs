use crate::DistributionFunction;
use crate::distribution_core::formula::DistributionFormula;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use crate::distribution_ops::negate::{distribution_negation, distribution_negation_inplace};
use crate::policy::{NegLog, Plain, YAxisPolicy};
use approx::ulps_eq;
use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use treetime_utils::make_error;

pub const TIME_LIMIT: f64 = 1e10;
pub const TIME_EPSILON: f64 = 1e-10;
const FORMULA_GRID_SIZE: usize = 200;

#[must_use]
#[allow(variant_size_differences)]
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum Distribution<Y: YAxisPolicy = Plain> {
  #[default]
  Empty,
  Point(DistributionPoint<f64, Y>),
  Range(DistributionRange<f64, Y>),
  Function(DistributionFunction<f64, Y>),
  Formula(DistributionFormula<Y>),
}

impl<Y: YAxisPolicy> Distribution<Y> {
  pub fn empty() -> Self {
    Self::Empty
  }

  pub fn point(x: f64, y: f64) -> Self {
    Self::Point(DistributionPoint::new(x, y))
  }

  pub fn range((x1, x2): (f64, f64), y: f64) -> Self {
    Self::Range(DistributionRange::new((x1, x2), y))
  }

  #[allow(clippy::needless_pass_by_value)]
  pub fn function(x: Array1<f64>, y: Array1<f64>) -> Result<Self, Report> {
    assert_eq!(x.shape(), y.shape());

    if x.is_empty() {
      return Ok(Self::empty());
    }

    if x.len() == 1 {
      return Ok(Self::point(x[0], y[0]));
    }

    if x.len() == 2 && ulps_eq!(y[0], y[1], max_ulps = 10) {
      return Ok(Self::range((x[0], x[1]), y[1]));
    }

    Ok(Self::Function(DistributionFunction::from_arrays(&x, y)?))
  }

  pub fn constant(amplitude: f64) -> Self {
    Distribution::range((-TIME_LIMIT, TIME_LIMIT), amplitude)
  }

  pub fn likely_time(&self) -> Option<f64> {
    match self {
      Self::Empty => {
        None //
      },
      Self::Point(p) => {
        Some(p.t()) //
      },
      Self::Range(r) => {
        Some(f64::midpoint(r.start(), r.end())) //
      },
      Self::Function(f) => {
        f.likely_time() //
      },
      Self::Formula(f) => {
        Some(f.likely_time()) //
      },
    }
  }

  pub fn t(&self) -> Array1<f64> {
    match self {
      Self::Empty => {
        ndarray::array![] //
      },
      Self::Point(p) => {
        ndarray::array![p.t()] //
      },
      Self::Range(r) => {
        ndarray::array![r.start(), r.end()] //
      },
      Self::Function(f) => {
        f.t().to_owned() //
      },
      Self::Formula(f) => {
        ndarray::array![f.t_min(), f.t_max()] //
      },
    }
  }

  pub fn y(&self) -> Array1<f64> {
    match self {
      Self::Point(p) => ndarray::array![p.amplitude()],
      Self::Range(r) => ndarray::array![r.amplitude(), r.amplitude()],
      Self::Function(f) => f.y().clone(),
      Self::Formula(f) => {
        let t = ndarray::array![f.t_min(), f.t_max()];
        f.eval_many(&t).unwrap_or_else(|_| ndarray::array![0.0, 0.0])
      },
      Self::Empty => ndarray::array![],
    }
  }

  pub fn negate(&self) -> Result<Self, Report> {
    distribution_negation(self)
  }

  pub fn negate_inplace(&mut self) -> Result<(), Report> {
    distribution_negation_inplace(self)
  }

  pub fn time_bounds(&self) -> (f64, f64) {
    let t = self.t();
    debug_assert!(!t.is_empty(), "Cannot extract time bounds from empty distribution");
    (t[0], t[t.len() - 1])
  }

  pub fn eval(&self, t: f64) -> Result<f64, Report> {
    match self {
      Self::Function(f) => Ok(f.interp(t)),
      Self::Formula(f) => f.eval_single(t),
      Self::Point(p) => {
        if ulps_eq!(t, p.t(), max_ulps = 10) {
          Ok(p.amplitude())
        } else {
          make_error!("Cannot evaluate point distribution outside its support")
        }
      },
      Self::Range(r) => {
        if t >= r.start() && t <= r.end() {
          Ok(r.amplitude())
        } else {
          make_error!("Cannot evaluate range distribution outside its support")
        }
      },
      Self::Empty => make_error!("Cannot evaluate empty distribution"),
    }
  }

  pub fn eval_many(&self, t: &Array1<f64>) -> Result<Array1<f64>, Report> {
    match self {
      Self::Function(f) => Ok(f.interp_many(t)),
      Self::Formula(f) => f.eval_many(t),
      Self::Point(p) => {
        let results = t
          .iter()
          .map(|&ti| {
            if ulps_eq!(ti, p.t(), max_ulps = 10) {
              Ok(p.amplitude())
            } else {
              make_error!("Cannot evaluate point distribution outside its support")
            }
          })
          .collect::<Result<Vec<f64>, Report>>()?;
        Ok(Array1::from(results))
      },
      Self::Range(r) => {
        let results = t
          .iter()
          .map(|&ti| {
            if ti >= r.start() && ti <= r.end() {
              Ok(r.amplitude())
            } else {
              make_error!("Cannot evaluate range distribution outside its support")
            }
          })
          .collect::<Result<Vec<f64>, Report>>()?;
        Ok(Array1::from(results))
      },
      Self::Empty => make_error!("Cannot evaluate empty distribution"),
    }
  }
}

impl Distribution<Plain> {
  /// Returns the maximum value of the distribution.
  pub fn max_value(&self) -> f64 {
    match self {
      Distribution::Empty => 0.0,
      Distribution::Point(p) => p.amplitude(),
      Distribution::Range(r) => r.amplitude(),
      Distribution::Function(f) => f.y().max().ok().copied().unwrap_or(0.0),
      Distribution::Formula(f) => discretize_formula(f).map_or(0.0, |df| df.y().max().ok().copied().unwrap_or(0.0)),
    }
  }

  /// Returns a new distribution with all values multiplied by factor.
  pub fn scale_by(&self, factor: f64) -> Self {
    match self {
      Distribution::Empty => Distribution::Empty,
      Distribution::Point(p) => Distribution::point(p.t(), p.amplitude() * factor),
      Distribution::Range(r) => Distribution::range((r.start(), r.end()), r.amplitude() * factor),
      Distribution::Function(f) => f.scale_y(factor).map_or(Distribution::Empty, Distribution::Function),
      Distribution::Formula(f) => discretize_formula(f)
        .and_then(|df| df.scale_y(factor))
        .map_or(Distribution::Empty, Distribution::Function),
    }
  }

  /// Returns a new distribution with max value = 1.0.
  /// Divides all values by max_value().
  /// Returns Empty if max_value() <= 0.
  pub fn normalize(&self) -> Self {
    let max_val = self.max_value();
    if max_val <= 0.0 || !max_val.is_finite() {
      return Distribution::Empty;
    }
    self.scale_by(1.0 / max_val)
  }

  /// Compute quantile (inverse CDF) for the distribution.
  ///
  /// Given probability p in [0, 1], returns the x value where CDF(x) = p.
  /// Uses trapezoidal integration to compute CDF and linear interpolation to find quantile.
  ///
  /// For Point distributions, returns the point location for any p.
  /// For Range distributions, returns linear interpolation between bounds.
  /// For Function distributions, computes from discrete CDF.
  #[allow(clippy::many_single_char_names)]
  pub fn quantile(&self, p: f64) -> Option<f64> {
    if !(0.0..=1.0).contains(&p) {
      return None;
    }

    match self {
      Distribution::Empty => None,
      Distribution::Point(point) => Some(point.t()),
      Distribution::Range(range) => {
        let start = range.start();
        let end = range.end();
        Some(start + p * (end - start))
      },
      Distribution::Function(f) => {
        let t = f.t();
        let y = f.y();
        let n = t.len();
        if n == 0 {
          return None;
        }
        if n == 1 {
          return Some(t[0]);
        }

        let Some(cdf) = compute_normalized_cdf(y, f.dx()) else {
          return self.likely_time();
        };

        // Find where CDF crosses p using linear interpolation
        if p <= 0.0 {
          return Some(t[0]);
        }
        if p >= 1.0 {
          return Some(t[n - 1]);
        }

        for i in 1..n {
          if cdf[i] >= p {
            let t0 = t[i - 1];
            let t1 = t[i];
            let c0 = cdf[i - 1];
            let c1 = cdf[i];
            if ulps_eq!(c0, c1, max_ulps = 10) {
              return Some(t0);
            }
            let frac = (p - c0) / (c1 - c0);
            return Some(t0 + frac * (t1 - t0));
          }
        }

        Some(t[n - 1])
      },
      Distribution::Formula(_) => {
        // Formula distributions don't support quantile computation directly
        self.likely_time()
      },
    }
  }

  /// Compute confidence interval bounds at given probabilities.
  ///
  /// Returns (lower, upper) where lower = quantile(p_lower) and upper = quantile(p_upper).
  pub fn confidence_interval(&self, p_lower: f64, p_upper: f64) -> Option<(f64, f64)> {
    let lower = self.quantile(p_lower)?;
    let upper = self.quantile(p_upper)?;
    Some((lower, upper))
  }

  /// Compute the highest posterior density (HPD) region containing `fraction`
  /// of the probability mass.
  ///
  /// The HPD region is the shortest interval containing the specified fraction
  /// of probability mass. For unimodal distributions, this is the interval
  /// where the PDF exceeds a threshold, found by bisection over thresholds.
  ///
  /// v0: `get_max_posterior_region(node, fraction=0.9)` in
  /// `clock_tree.py:1146-1230`.
  ///
  /// Boundary handling (v0 `clock_tree.py:1175-1178`):
  /// - Peak at left boundary: `[t_min, quantile(fraction)]`
  /// - Peak at right boundary: `[quantile(1-fraction), t_max]`
  /// - Interior peak: bisection to find probability threshold where the
  ///   level-set interval contains exactly `fraction` of the mass.
  ///
  /// For symmetric distributions, HPD equals equal-tailed CI. For skewed
  /// distributions (nodes near tree boundaries), HPD is narrower and
  /// centered on the peak.
  pub fn hpd_region(&self, fraction: f64) -> Option<(f64, f64)> {
    if !(0.0..=1.0).contains(&fraction) {
      return None;
    }

    match self {
      Distribution::Empty => None,
      Distribution::Point(p) => Some((p.t(), p.t())),
      Distribution::Range(r) => {
        // Uniform: any sub-interval of correct length has equal density.
        // Use centered interval (same as equal-tailed).
        let width = r.end() - r.start();
        let margin = (1.0 - fraction) * 0.5 * width;
        Some((r.start() + margin, r.end() - margin))
      },
      Distribution::Function(f) => hpd_region_function(f, fraction),
      Distribution::Formula(_) => {
        // Fall back to equal-tailed for Formula distributions
        let p_lo = (1.0 - fraction) * 0.5;
        self.confidence_interval(p_lo, 1.0 - p_lo)
      },
    }
  }

  pub fn to_neglog(&self) -> Distribution<NegLog> {
    match self {
      Self::Empty => Distribution::Empty,
      Self::Point(p) => Distribution::point(p.t(), NegLog::from_plain(p.amplitude())),
      Self::Range(r) => Distribution::range((r.start(), r.end()), NegLog::from_plain(r.amplitude())),
      Self::Function(f) => {
        let y_neglog = f.y().mapv(NegLog::from_plain);
        Distribution::Function(DistributionFunction::from_grid_fn(
          treetime_grid::GridFn::from_start_dx_values(f.x_min(), f.dx(), y_neglog)
            .expect("Grid construction should not fail for valid input"),
        ))
      },
      Self::Formula(f) => match discretize_formula(f) {
        Ok(df) => {
          let y_neglog = df.y().mapv(NegLog::from_plain);
          Distribution::Function(DistributionFunction::from_grid_fn(
            treetime_grid::GridFn::from_start_dx_values(df.x_min(), df.dx(), y_neglog)
              .expect("Grid construction should not fail for valid input"),
          ))
        },
        Err(_) => Distribution::Empty,
      },
    }
  }
}

impl Distribution<NegLog> {
  pub fn to_plain(&self) -> Distribution<Plain> {
    match self {
      Self::Empty => Distribution::Empty,
      Self::Point(p) => Distribution::point(p.t(), NegLog::to_plain(p.amplitude())),
      Self::Range(r) => Distribution::range((r.start(), r.end()), NegLog::to_plain(r.amplitude())),
      Self::Function(f) => {
        let y_plain = f.y().mapv(NegLog::to_plain);
        Distribution::Function(DistributionFunction::from_grid_fn(
          treetime_grid::GridFn::from_start_dx_values(f.x_min(), f.dx(), y_plain)
            .expect("Grid construction should not fail for valid input"),
        ))
      },
      Self::Formula(f) => {
        let t_min = f.t_min();
        let t_max = f.t_max();
        let f = f.clone();
        let eval_fn = move |t: f64| -> Result<f64, Report> {
          let y = f.eval_single(t)?;
          Ok(NegLog::to_plain(y))
        };
        Distribution::Formula(DistributionFormula::new(eval_fn, t_min, t_max))
      },
    }
  }
}

fn discretize_formula(f: &DistributionFormula<Plain>) -> Result<DistributionFunction<f64, Plain>, Report> {
  let n_points = FORMULA_GRID_SIZE;
  let t = Array1::from_shape_fn(n_points, |i| {
    f.t_min() + (f.t_max() - f.t_min()) * (i as f64 / (n_points - 1) as f64)
  });
  let values = f.eval_many(&t)?;
  DistributionFunction::from_range_values((f.t_min(), f.t_max()), values)
}

/// HPD region for a discretized distribution on a uniform grid.
///
/// Finds the shortest interval containing `fraction` of the probability mass
/// by bisecting on a probability threshold. All points where the PDF exceeds
/// the threshold form the HPD region.
///
/// Assumes unimodality: the crossing search radiates outward from a single
/// argmax peak, so the result is always one contiguous interval. For a
/// multimodal distribution the interval spans the valley between modes,
/// producing a conservative (wider than true HPD) estimate. Typical
/// phylogenetic posteriors are approximately Gaussian, so this is not
/// a practical concern. Matches v0 (`clock_tree.py:1185`).
///
/// v0 equivalent: the optimization loop in `get_max_posterior_region`
/// (`clock_tree.py:1184-1226`), adapted from neg-log space to plain
/// probability space.
fn hpd_region_function(f: &DistributionFunction<f64, Plain>, fraction: f64) -> Option<(f64, f64)> {
  let t = f.t();
  let y = f.y();
  let n = t.len();
  let dx = f.dx();
  let x_min = f.x_min();

  if n == 0 {
    return None;
  }
  if n == 1 {
    return Some((t[0], t[0]));
  }

  let cdf = compute_normalized_cdf(y, dx)?;

  // Peak = argmax of probability density
  let pidx = y.argmax().ok()?;

  // Boundary cases (v0 clock_tree.py:1175-1178, 1194):
  // Peak on boundary or array too short for interior interpolation.
  // HPD for a boundary peak is a one-sided interval from the peak boundary
  // to the quantile capturing `fraction` of the mass.
  if n < 3 || pidx == 0 {
    let upper = interp_cdf_inverse(&t, &cdf, fraction);
    return Some((t[0], upper));
  }
  if pidx == n - 1 {
    let lower = interp_cdf_inverse(&t, &cdf, 1.0 - fraction);
    return Some((lower, t[n - 1]));
  }

  // Interior peak: bisect on probability threshold p_thresh in [0, peak_val].
  //
  // For each threshold, the HPD region is [left(p_thresh), right(p_thresh)]
  // where left/right are where the PDF crosses the threshold on each side
  // of the peak. As p_thresh decreases from peak_val toward 0, the interval
  // widens and captures more mass.
  //
  // Find p_thresh where CDF(right) - CDF(left) = fraction.
  let peak_val = y[pidx];
  let mut lo = 0.0_f64; // wide interval, mass ~ 1
  let mut hi = peak_val; // narrow interval, mass ~ 0

  // ~50 iterations gives ~15 decimal digits; 64 is comfortable
  for _ in 0..64 {
    let mid = 0.5 * (lo + hi);
    let left_pos = interp_crossing_left(&t, y, pidx, mid);
    let right_pos = interp_crossing_right(&t, y, pidx, mid);
    let mass = interp_cdf_at_uniform(&cdf, x_min, dx, right_pos) - interp_cdf_at_uniform(&cdf, x_min, dx, left_pos);

    if mass > fraction {
      lo = mid; // raise threshold to narrow interval
    } else {
      hi = mid; // lower threshold to widen interval
    }

    if (hi - lo) < 1e-14 * peak_val.max(1e-30) {
      break;
    }
  }

  let p_thresh = 0.5 * (lo + hi);
  let left_pos = interp_crossing_left(&t, y, pidx, p_thresh);
  let right_pos = interp_crossing_right(&t, y, pidx, p_thresh);
  Some((left_pos, right_pos))
}

/// Compute normalized CDF from probability density values on a uniform grid.
///
/// Uses the trapezoidal rule for integration and normalizes to [0, 1].
/// Returns `None` if the total integral is non-positive or non-finite.
fn compute_normalized_cdf(y: &Array1<f64>, dx: f64) -> Option<Array1<f64>> {
  let n = y.len();
  debug_assert!(n >= 2, "CDF requires at least 2 points, got {n}");
  let mut cdf = Array1::<f64>::zeros(n);
  for i in 1..n {
    cdf[i] = cdf[i - 1] + 0.5 * (y[i - 1] + y[i]) * dx;
  }
  let total = cdf[n - 1];
  if total <= 0.0 || !total.is_finite() {
    return None;
  }
  cdf.mapv_inplace(|v| v / total);
  Some(cdf)
}

/// Find position left of peak where y crosses `threshold` (linear interpolation).
///
/// Searches from peak toward the left boundary. Returns `t[0]` if the
/// threshold is never reached (PDF stays above threshold across the full
/// left side).
fn interp_crossing_left(t: &Array1<f64>, y: &Array1<f64>, pidx: usize, threshold: f64) -> f64 {
  // Search from peak leftward: find first interval where y drops below threshold
  for i in (1..=pidx).rev() {
    if y[i - 1] <= threshold && y[i] > threshold {
      let frac = (threshold - y[i - 1]) / (y[i] - y[i - 1]);
      return t[i - 1] + frac * (t[i] - t[i - 1]);
    }
  }
  t[0]
}

/// Find position right of peak where y crosses `threshold` (linear interpolation).
///
/// Searches from peak toward the right boundary. Returns `t[n-1]` if the
/// threshold is never reached.
fn interp_crossing_right(t: &Array1<f64>, y: &Array1<f64>, pidx: usize, threshold: f64) -> f64 {
  let n = t.len();
  // Search from peak rightward: find first interval where y drops below threshold
  for i in pidx..n - 1 {
    if y[i] > threshold && y[i + 1] <= threshold {
      // y[i] > threshold > y[i+1], both diffs negative, frac in [0,1]
      let frac = (threshold - y[i]) / (y[i + 1] - y[i]);
      return t[i] + frac * (t[i + 1] - t[i]);
    }
  }
  t[n - 1]
}

/// Interpolate CDF value at an arbitrary position on a uniform grid.
///
/// O(1) via direct index computation from uniform grid spacing.
fn interp_cdf_at_uniform(cdf: &Array1<f64>, x_min: f64, dx: f64, pos: f64) -> f64 {
  let n = cdf.len();
  if n == 0 {
    return 0.0;
  }
  let idx_f = (pos - x_min) / dx;
  if idx_f <= 0.0 {
    return cdf[0];
  }
  let last = (n - 1) as f64;
  if idx_f >= last {
    return cdf[n - 1];
  }
  let idx = (idx_f.floor() as usize).min(n - 2);
  let frac = idx_f - idx as f64;
  cdf[idx] + frac * (cdf[idx + 1] - cdf[idx])
}

/// Inverse CDF: find position where CDF = p (linear interpolation).
fn interp_cdf_inverse(t: &Array1<f64>, cdf: &Array1<f64>, p: f64) -> f64 {
  let n = t.len();
  if p <= 0.0 {
    return t[0];
  }
  if p >= 1.0 {
    return t[n - 1];
  }

  for i in 1..n {
    if cdf[i] >= p {
      let c0 = cdf[i - 1];
      let c1 = cdf[i];
      if (c1 - c0).abs() < 1e-30 {
        return t[i - 1];
      }
      let frac = (p - c0) / (c1 - c0);
      return t[i - 1] + frac * (t[i] - t[i - 1]);
    }
  }
  t[n - 1]
}

pub type DistributionPlain = Distribution<Plain>;
pub type DistributionNegLog = Distribution<NegLog>;

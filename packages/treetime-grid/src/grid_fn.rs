use std::fmt::Debug;

use crate::InterpElem;
use crate::grid::Grid;
use crate::interp_nonuniform::interp_nonuniform;
use approx::{UlpsEq, ulps_eq};
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::{Array1, s};
use ndarray_stats::QuantileExt;
use num::Float;
use serde::{Deserialize, Serialize};
use treetime_utils::array::ndarray::has_uniform_spacing;
use treetime_utils::array::serde::{array1_as_vec, array1_from_vec};
use treetime_utils::make_error;

/// Power-law interpolation between a hard boundary and the nearest grid point.
///
/// Near a hard boundary the density follows `p(t) = C * |t - t_hard|^b` in
/// plain-probability space, equivalently `y(t) = a - b * log|t - t_hard|` in
/// neg-log space where `C = exp(-a)`.
///
/// - `b = 0`: constant at the boundary (zero-mutation branch, density is maximal
///   at t=0). The grid endpoint stores the finite boundary value directly.
/// - `b > 0`: power-law decay toward zero (n-mutation branch, `p(t) ~ t^n`). The
///   grid stays strictly inside the hard boundary and the approach law covers the
///   gap `[t_hard, t_first)`.
///
/// Under multiplication the exponents add (`b_result = b_a + b_b`) because the
/// product of two power laws is a power law with summed exponents.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApproachLaw {
  /// Location of the hard boundary (e.g. t=0 for branch lengths).
  pub t_hard: f64,
  /// Amplitude coefficient C in `p(t) = C * |t - t_hard|^b`. Always positive.
  pub coeff: f64,
  /// Power-law exponent b >= 0.
  pub exponent: f64,
}

impl ApproachLaw {
  /// Fit an approach law from the innermost grid points adjacent to a hard boundary.
  ///
  /// Uses least-squares regression on `(log|t - t_hard|, log(p(t)))` from the `n_fit`
  /// innermost points. The slope gives `b`, the intercept gives `log(C)`.
  ///
  /// `side` selects which end of the grid is near the hard boundary:
  /// - `Side::Left`: the hard boundary is to the left of `x_min`, fit from the leftmost points
  /// - `Side::Right`: the hard boundary is to the right of `x_max`, fit from the rightmost points
  ///
  /// Returns `None` when the fit is not meaningful (all y values zero, or fewer than 2
  /// positive points near the boundary).
  pub fn fit(grid_fn: &GridFn<f64>, t_hard: f64, side: Side, n_fit: usize) -> Option<Self> {
    let n = grid_fn.n_points();
    let n_fit = n_fit.min(n);

    let (xs, ys): (Vec<f64>, Vec<f64>) = match side {
      Side::Left => (0..n_fit)
        .filter_map(|i| {
          let x = grid_fn.grid().x_at(i);
          let y = grid_fn.y()[i];
          (y > 0.0).then(|| ((x - t_hard).abs().ln(), y.ln()))
        })
        .collect(),
      Side::Right => (0..n_fit)
        .filter_map(|i| {
          let idx = n - 1 - i;
          let x = grid_fn.grid().x_at(idx);
          let y = grid_fn.y()[idx];
          (y > 0.0).then(|| ((t_hard - x).abs().ln(), y.ln()))
        })
        .collect(),
    };

    if xs.len() < 2 {
      return None;
    }

    let (slope, intercept) = least_squares_fit(&xs, &ys);

    // b >= 0: a negative fit means the density increases into the boundary faster
    // than any power law, which is unphysical; clamp to b = 0.
    let exponent = slope.max(0.0);
    let coeff = intercept.exp();

    if !coeff.is_finite() || coeff <= 0.0 {
      return None;
    }

    Some(ApproachLaw {
      t_hard,
      coeff,
      exponent,
    })
  }

  /// Evaluate the approach law at a point between `t_hard` and the nearest grid point.
  ///
  /// Returns `C * |t - t_hard|^b` in plain-probability space.
  pub fn eval(&self, t: f64) -> f64 {
    let dt = (t - self.t_hard).abs();
    if dt == 0.0 {
      if self.exponent == 0.0 {
        return self.coeff;
      }
      return 0.0;
    }
    self.coeff * dt.powf(self.exponent)
  }

  /// Compose two approach laws under multiplication (exponents add).
  ///
  /// The product `p_a(t) * p_b(t) = C_a * dt^b_a * C_b * dt^b_b = (C_a * C_b) * dt^(b_a + b_b)`.
  /// Both laws must share the same `t_hard`.
  #[must_use]
  pub fn compose_multiply(&self, other: &ApproachLaw) -> ApproachLaw {
    ApproachLaw {
      t_hard: self.t_hard,
      coeff: self.coeff * other.coeff,
      exponent: self.exponent + other.exponent,
    }
  }

  /// Mass in the approach region between `t_hard` and `t_grid` (the nearest grid point).
  ///
  /// Closed form: `C * |t_grid - t_hard|^(b+1) / (b + 1)`.
  pub fn mass(&self, t_grid: f64) -> f64 {
    let dt = (t_grid - self.t_hard).abs();
    self.coeff * dt.powf(self.exponent + 1.0) / (self.exponent + 1.0)
  }

  /// Transform the approach law when the argument is negated: `f(x) -> f(-x)`.
  /// The hard boundary moves to `-t_hard`.
  #[must_use]
  pub fn negate_arg(&self) -> Self {
    ApproachLaw {
      t_hard: -self.t_hard,
      coeff: self.coeff,
      exponent: self.exponent,
    }
  }
}

/// Which side of the grid a hard boundary is on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
  Left,
  Right,
}

/// Simple least-squares linear regression: y = slope * x + intercept.
fn least_squares_fit(xs: &[f64], ys: &[f64]) -> (f64, f64) {
  let n = xs.len() as f64;
  let sum_x: f64 = xs.iter().sum();
  let sum_y: f64 = ys.iter().sum();
  let sum_xx: f64 = xs.iter().map(|x| x * x).sum();
  let sum_xy: f64 = xs.iter().zip(ys).map(|(x, y)| x * y).sum();

  let sum_x_sq = sum_x * sum_x;
  let denom = n * sum_xx - sum_x_sq;
  if denom.abs() < 1e-30 {
    return (0.0, sum_y / n);
  }

  let slope = (n * sum_xy - sum_x * sum_y) / denom;
  let intercept = (sum_y - slope * sum_x) / n;
  (slope, intercept)
}

/// Behavior of a [`GridFn`] when evaluated outside its grid support.
///
/// A bare grid function is a generic interpolant with no probabilistic meaning, so the
/// default is [`BoundaryBehavior::Error`]: a query outside the grid is a programming error
/// unless the caller has declared how the tail should behave. The two non-default variants
/// are explicit opt-ins used by finite-support distributions and by the timetree message
/// passes, which assign a per-side tail policy to each message.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub enum BoundaryBehavior {
  /// Out-of-support evaluation is an error.
  #[default]
  Error,
  /// Hard boundary: the domain terminates and probability is zero beyond the grid edge.
  /// An optional [`ApproachLaw`] provides power-law interpolation between the hard boundary
  /// and the nearest grid point. In plain-probability space a bare `Hard(None)` returns
  /// `0.0`; a negative-log representation cannot express zero probability as `0.0` (it
  /// would be `+inf`), so that combination is rejected at the distribution layer.
  Hard(Option<ApproachLaw>),
  /// Return the nearest boundary value (`y[0]` to the left, `y[n-1]` to the right),
  /// i.e. a flat tail. Use when the function is genuinely uninformative beyond the edge.
  Constant,
}

impl BoundaryBehavior {
  /// Whether this boundary is *soft*: the distribution continues past the grid edge, so the
  /// grid boundary is a representation choice (where interpolation stopped), not a fact about
  /// the distribution. A soft boundary is evaluable beyond the grid via its declared tail law
  /// and therefore *extends* the evaluable domain under multiplication.
  ///
  /// The complement is a *hard* boundary: the grid edge terminates the domain (`Hard`: zero
  /// probability beyond; `Error`: out-of-support evaluation is undefined). A hard boundary
  /// *restricts* the result domain and is never silently extended.
  ///
  /// Currently `Constant` is the only soft tail. This predicate is what the multiplication
  /// rule keys off, so a future soft tail law extends the domain without touching that rule.
  pub fn is_soft(self) -> bool {
    matches!(self, BoundaryBehavior::Constant)
  }

  pub fn approach_law(self) -> Option<ApproachLaw> {
    match self {
      BoundaryBehavior::Hard(law) => law,
      _ => None,
    }
  }
}

/// Function represented on a uniform grid for piecewise linear interpolation
///
/// Represents a function as a set of (x, y) points on a uniformly-spaced grid, providing
/// linear interpolation between points. Behavior outside the grid is governed by the
/// per-side [`BoundaryBehavior`] policies `left_extrap` and `right_extrap`, which default
/// to [`BoundaryBehavior::Error`].
///
/// # Invariants
///
/// - Grid must be uniformly spaced with spacing `dx`
/// - Grid starts at `x_min`
/// - `y` array must contain at least 2 points
/// - `dx` must be positive
/// - These invariants are enforced by the type system and cannot be violated
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(bound(serialize = "T: Serialize", deserialize = "T: Deserialize<'de>"))]
pub struct GridFn<T: InterpElem> {
  grid: Grid<T>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  y: Array1<T>,
  // Out-of-support tail policy is runtime behavior, not persisted data: skip it so serialized
  // output (auspice node data, snapshots) is unchanged. Deserialization restores the default.
  #[serde(skip)]
  left_extrap: BoundaryBehavior,
  #[serde(skip)]
  right_extrap: BoundaryBehavior,
}

impl<T: InterpElem> GridFn<T> {
  pub fn from_grid_array(grid: Grid<T>, y: Array1<T>) -> Result<Self, Report> {
    if grid.n_points() != y.len() {
      return make_error!(
        "Grid has {} points but y array has {} elements",
        grid.n_points(),
        y.len()
      );
    }
    Ok(Self {
      grid,
      y,
      left_extrap: BoundaryBehavior::default(),
      right_extrap: BoundaryBehavior::default(),
    })
  }

  /// Set the out-of-support behavior for the left (below `x_min`) tail.
  #[must_use]
  pub fn with_left_extrap(mut self, behavior: BoundaryBehavior) -> Self {
    self.left_extrap = behavior;
    self
  }

  /// Set the out-of-support behavior for the right (above `x_max`) tail.
  #[must_use]
  pub fn with_right_extrap(mut self, behavior: BoundaryBehavior) -> Self {
    self.right_extrap = behavior;
    self
  }

  /// Set the same out-of-support behavior for both tails.
  #[must_use]
  pub fn with_extrap(self, behavior: BoundaryBehavior) -> Self {
    self.with_left_extrap(behavior).with_right_extrap(behavior)
  }

  pub fn left_extrap(&self) -> BoundaryBehavior {
    self.left_extrap
  }

  pub fn right_extrap(&self) -> BoundaryBehavior {
    self.right_extrap
  }

  pub fn from_grid_fn<F>(grid: Grid<T>, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    let n_points = grid.n_points();
    let y = Array1::from_shape_fn(n_points, |i| y_fn(grid.x_at(i)));
    Self::from_grid_array(grid, y)
  }

  pub fn from_arrays(x: &Array1<T>, y: Array1<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    if x.len() != y.len() {
      return make_error!(
        "x and y arrays must have the same length, got {} and {}",
        x.len(),
        y.len()
      );
    }
    let grid = Grid::from_array(x)?;
    Self::from_grid_array(grid, y)
  }

  /// Constructs GridFn from non-uniformly spaced arrays by resampling to uniform grid
  ///
  /// Takes non-uniform (x, y) arrays and resamples them to a uniform grid using linear
  /// interpolation. The grid spacing is determined by the smallest spacing in the input
  /// to preserve detail.
  ///
  /// # Arguments
  ///
  /// * `x` - Non-uniform x coordinates (must be sorted ascending)
  /// * `y` - Corresponding y values
  ///
  /// # Returns
  ///
  /// GridFn with uniformly spaced grid covering the same range as input
  pub fn from_arrays_nonuniform(x: &Array1<T>, y: &Array1<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    if x.len() < 2 {
      return make_error!("Grid must have at least 2 points, got {}", x.len());
    }
    if x.len() != y.len() {
      return make_error!(
        "x and y arrays must have the same length, got {} and {}",
        x.len(),
        y.len()
      );
    }

    if has_uniform_spacing(x) {
      let grid = Grid::from_array(x)?;
      return Self::from_grid_array(grid, y.clone());
    }

    let x_min = x[0];
    let x_max = x[x.len() - 1];
    let dx = find_min_spacing(x)?;
    let n_points = ((x_max - x_min) / dx).ceil().to_usize().unwrap() + 1;
    let dx = (x_max - x_min) / T::from(n_points - 1).unwrap();

    if n_points > 1_000_000 {
      return make_error!("Resampling would require {n_points} points, which exceeds safety limit");
    }

    let grid = Grid::from_start_dx(x_min, dx, n_points)?;
    let y_uniform = interp_nonuniform(x, y, n_points, |i| grid.x_at(i))?;
    Self::from_grid_array(grid, y_uniform)
  }

  pub fn from_n_points<F>((x_min, x_max): (T, T), n_points: usize, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    let grid = Grid::from_range_n_points(x_min, x_max, n_points)?;
    Self::from_grid_fn(grid, y_fn)
  }

  pub fn from_grid<F>((x_min, x_max): (T, T), dx: T, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    let grid = Grid::from_range_dx(x_min, x_max, dx)?;
    Self::from_grid_fn(grid, y_fn)
  }

  pub fn from_start_dx_values(x_min: T, dx: T, y: Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid = Grid::from_start_dx(x_min, dx, y.len())?;
    Self::from_grid_array(grid, y)
  }

  pub fn from_range_values((x_min, x_max): (T, T), y: Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid = Grid::from_range_n_points(x_min, x_max, y.len())?;
    Self::from_grid_array(grid, y)
  }

  pub fn constant((x_min, x_max): (T, T), n_points: usize, value: T) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid = Grid::from_range_n_points(x_min, x_max, n_points)?;
    let y = Array1::from_elem(n_points, value);
    Self::from_grid_array(grid, y)
  }

  pub fn zeros((x_min, x_max): (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    Self::constant((x_min, x_max), n_points, T::zero())
  }

  pub fn ones((x_min, x_max): (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    Self::constant((x_min, x_max), n_points, T::one())
  }

  // TODO: inefficient. Try to remove this method
  pub fn x(&self) -> Array1<T>
  where
    T: Float,
  {
    self.grid.to_array()
  }

  pub fn y(&self) -> &Array1<T> {
    &self.y
  }

  pub fn grid(&self) -> &Grid<T> {
    &self.grid
  }

  pub fn x_min(&self) -> T {
    self.grid.x_min()
  }

  pub fn x_max(&self) -> T
  where
    T: Float,
  {
    self.grid.x_max()
  }

  pub fn x_range(&self) -> (T, T)
  where
    T: Float,
  {
    self.grid.x_range()
  }

  pub fn n_points(&self) -> usize {
    self.grid.n_points()
  }

  pub fn len(&self) -> usize {
    self.grid.len()
  }

  pub fn is_empty(&self) -> bool {
    self.grid.is_empty()
  }

  pub fn dx(&self) -> T {
    self.grid.dx()
  }

  /// Returns the minimum y value.
  ///
  /// Uses `fold` starting from first element since `GridFn` invariants guarantee at least 2 points.
  pub fn y_min(&self) -> T {
    self
      .y
      .iter()
      .copied()
      .skip(1)
      .fold(self.y[0], |a, b| if a < b { a } else { b })
  }

  /// Returns the maximum y value.
  ///
  /// Uses `fold` starting from first element since `GridFn` invariants guarantee at least 2 points.
  pub fn y_max(&self) -> T {
    self
      .y
      .iter()
      .copied()
      .skip(1)
      .fold(self.y[0], |a, b| if a > b { a } else { b })
  }

  pub fn y_range(&self) -> (T, T) {
    (self.y_min(), self.y_max())
  }

  pub fn to_pairs(&self) -> Vec<(T, T)>
  where
    T: Float,
  {
    izip!(self.grid.iter(), self.y.iter())
      .map(|(x, &y)| (x, y))
      .collect_vec()
  }

  /// Interpolate function value at a single point
  ///
  /// Uses piecewise linear interpolation within the grid bounds. Outside the bounds the
  /// per-side [`BoundaryBehavior`] applies: `Error` (default) rejects the query, `Hard`
  /// returns `0.0`, `Constant` returns the nearest boundary value.
  ///
  /// # Arguments
  ///
  /// * `xi` - Point at which to evaluate the function
  ///
  /// # Returns
  ///
  /// Interpolated value at `xi`, or an error when `xi` is outside the support and the
  /// relevant tail policy is [`BoundaryBehavior::Error`].
  pub fn interp(&self, xi: T) -> Result<T, Report>
  where
    T: Float + UlpsEq,
  {
    let x_min = self.grid.x_min();
    let x_max = self.grid.x_max();

    // A query at the nominal support boundary can land a few ulps outside the grid, because
    // x_max is reconstructed as x_min + (n-1)*dx and need not reproduce the value originally
    // passed in. Treat such a query as the boundary grid point, not as extrapolation, so that
    // evaluating a distribution exactly at its own endpoint is in-support.
    if xi < x_min {
      if ulps_eq!(xi, x_min, max_ulps = 4) {
        return Ok(self.y[0]);
      }
      return self.extrapolate(self.left_extrap, self.y[0], xi, x_min, "below");
    }

    if xi > x_max {
      let n = self.grid.n_points();
      if ulps_eq!(xi, x_max, max_ulps = 4) {
        return Ok(self.y[n - 1]);
      }
      return self.extrapolate(self.right_extrap, self.y[n - 1], xi, x_max, "above");
    }

    let idx = self.grid.find_interval_index(xi);
    Ok(self.interpolate_at(xi, idx))
  }

  /// Interpolate function values at multiple points.
  ///
  /// Applies [`GridFn::interp`] to each query; a single out-of-support query under an
  /// `Error` tail fails the whole call.
  pub fn interp_many(&self, queries: &Array1<T>) -> Result<Array1<T>, Report>
  where
    T: Float + UlpsEq,
  {
    let values = queries
      .iter()
      .map(|&q| self.interp(q))
      .collect::<Result<Vec<T>, Report>>()?;
    Ok(Array1::from_vec(values))
  }

  fn extrapolate(&self, behavior: BoundaryBehavior, boundary_value: T, xi: T, bound: T, side: &str) -> Result<T, Report>
  where
    T: Float,
  {
    match behavior {
      BoundaryBehavior::Constant => Ok(boundary_value),
      BoundaryBehavior::Hard(Some(law)) => {
        let xi_f64 = xi.to_f64().unwrap();
        // Beyond the hard boundary: zero probability
        if (side == "below" && xi_f64 < law.t_hard) || (side != "below" && xi_f64 > law.t_hard) {
          return Ok(T::zero());
        }
        // Between hard boundary and grid edge: use approach law
        Ok(T::from(law.eval(xi_f64)).unwrap())
      },
      BoundaryBehavior::Hard(None) => Ok(T::zero()),
      BoundaryBehavior::Error => make_error!(
        "GridFn evaluated at {xi:?}, {side} the support boundary {bound:?}, but no extrapolation policy is set for that side"
      ),
    }
  }

  fn interpolate_at(&self, q: T, idx: usize) -> T
  where
    T: Float,
  {
    let n = self.grid.n_points();
    if idx >= n - 1 {
      return self.y[n - 1];
    }
    let x0 = self.grid.x_at(idx);
    let y0 = self.y[idx];
    let y1 = self.y[idx + 1];
    let t = (q - x0) / self.grid.dx();
    y0 + t * (y1 - y0)
  }

  #[must_use]
  pub fn mapv<F>(&self, f: F) -> Self
  where
    F: Fn(T) -> T,
  {
    Self {
      grid: self.grid,
      y: self.y.mapv(f),
      // Arbitrary y-transforms invalidate approach laws because the coefficient
      // cannot be updated without knowing the transform's structure.
      left_extrap: strip_approach_law(self.left_extrap),
      right_extrap: strip_approach_law(self.right_extrap),
    }
  }

  /// Scale all y-values by a multiplicative factor, preserving approach laws.
  ///
  /// Unlike [`Self::mapv`] with an arbitrary function, multiplication by a constant
  /// preserves the power-law shape: `(s * C) * dt^b` is still a power law with the
  /// same exponent. The approach law coefficient is scaled by the same factor.
  #[must_use]
  pub fn scale_y(&self, factor: f64) -> Self
  where
    T: Float,
  {
    Self {
      grid: self.grid,
      y: self.y.mapv(|v| v * T::from(factor).unwrap()),
      left_extrap: scale_approach_law(self.left_extrap, factor),
      right_extrap: scale_approach_law(self.right_extrap, factor),
    }
  }

  pub fn mapv_inplace<F>(&mut self, f: F)
  where
    F: Fn(T) -> T,
  {
    self.y.mapv_inplace(f);
  }

  /// Negates the argument of the function: f(x) -> f(-x).
  /// This reflects the function across the y-axis.
  /// The domain [x_min, x_max] becomes [-x_max, -x_min].
  /// The y-values are reversed.
  pub fn negate_arg(&self) -> Result<Self, Report>
  where
    T: Float,
  {
    let mut result = self.clone();
    result.negate_arg_inplace()?;
    Ok(result)
  }

  /// Negates the argument of the function in-place: f(x) -> f(-x).
  /// This reflects the function across the y-axis.
  /// The domain [x_min, x_max] becomes [-x_max, -x_min].
  /// The y-values are reversed.
  pub fn negate_arg_inplace(&mut self) -> Result<(), Report>
  where
    T: Float,
  {
    let x_max = self.grid.x_max();
    let dx = self.grid.dx();
    let n_points = self.grid.n_points();
    self.grid = Grid::from_start_dx(-x_max, dx, n_points)?;

    // Reverse y in-place
    let n = self.y.len();
    for i in 0..n / 2 {
      self.y.swap(i, n - 1 - i);
    }

    // Negating the argument reflects the domain, so the left and right tails swap sides.
    // Approach laws inside Hard variants also need their t_hard negated.
    self.left_extrap = negate_approach_law(self.left_extrap);
    self.right_extrap = negate_approach_law(self.right_extrap);
    std::mem::swap(&mut self.left_extrap, &mut self.right_extrap);
    Ok(())
  }

  /// Resamples function to a new grid
  ///
  /// Creates a new GridFn on the specified grid, using linear interpolation
  /// (and extrapolation if needed) from the current function.
  ///
  /// # Arguments
  ///
  /// * `grid` - Target grid for resampling
  ///
  /// # Returns
  ///
  /// New GridFn with values interpolated onto the target grid
  pub fn resample(&self, grid: &Grid<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let n_points = grid.n_points();
    let y_new = (0..n_points)
      .map(|i| self.interp(grid.x_at(i)))
      .collect::<Result<Vec<T>, Report>>()?;
    let resampled = Self::from_grid_array(*grid, Array1::from_vec(y_new))?;
    Ok(
      resampled
        .with_left_extrap(self.left_extrap)
        .with_right_extrap(self.right_extrap),
    )
  }

  /// Resamples function to a new uniform grid with specified start, spacing, and length
  ///
  /// Creates a new GridFn on a uniform grid defined by starting point, spacing,
  /// and number of points, using linear interpolation (and extrapolation if needed)
  /// from the current function.
  ///
  /// # Arguments
  ///
  /// * `x_min` - Starting x coordinate
  /// * `dx` - Grid spacing
  /// * `n_points` - Number of points in the new grid
  ///
  /// # Returns
  ///
  /// New GridFn with uniformly spaced grid
  pub fn resample_start_dx(&self, x_min: T, dx: T, n_points: usize) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid = Grid::from_start_dx(x_min, dx, n_points)?;
    self.resample(&grid)
  }

  /// Resamples function to a new uniform grid with specified range and number of points
  ///
  /// Creates a new GridFn on a uniform grid with the given number of points, using
  /// linear interpolation (and extrapolation if needed) from the current function.
  ///
  /// # Arguments
  ///
  /// * `x_range` - New grid range (x_min, x_max)
  /// * `n_points` - Number of points in the new grid
  ///
  /// # Returns
  ///
  /// New GridFn with uniformly spaced grid covering the specified range
  pub fn resample_range_n_points(&self, x_range: (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let (x_min, x_max) = x_range;
    let grid = Grid::from_range_n_points(x_min, x_max, n_points)?;
    self.resample(&grid)
  }

  /// Resamples function to a new uniform grid with specified range and spacing
  ///
  /// Creates a new GridFn on a uniform grid with the given spacing, using linear
  /// interpolation (and extrapolation if needed) from the current function.
  ///
  /// # Arguments
  ///
  /// * `x_range` - New grid range (x_min, x_max)
  /// * `dx` - Grid spacing for the new grid
  ///
  /// # Returns
  ///
  /// New GridFn with uniformly spaced grid covering the specified range
  pub fn resample_range_dx(&self, x_range: (T, T), dx: T) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let (x_min, x_max) = x_range;
    let grid = Grid::from_range_dx(x_min, x_max, dx)?;
    self.resample(&grid)
  }
}

fn strip_approach_law(behavior: BoundaryBehavior) -> BoundaryBehavior {
  match behavior {
    BoundaryBehavior::Hard(_) => BoundaryBehavior::Hard(None),
    other => other,
  }
}

fn scale_approach_law(behavior: BoundaryBehavior, factor: f64) -> BoundaryBehavior {
  match behavior {
    BoundaryBehavior::Hard(Some(law)) => BoundaryBehavior::Hard(Some(ApproachLaw {
      coeff: law.coeff * factor,
      ..law
    })),
    other => other,
  }
}

fn negate_approach_law(behavior: BoundaryBehavior) -> BoundaryBehavior {
  match behavior {
    BoundaryBehavior::Hard(Some(law)) => BoundaryBehavior::Hard(Some(law.negate_arg())),
    other => other,
  }
}

/// Finds the smallest spacing between consecutive points in a sorted array
///
/// Used to determine optimal grid spacing when resampling non-uniform data.
/// Preserves maximum detail by using the finest resolution present in the input.
///
/// # Arguments
///
/// * `x` - Sorted array of x coordinates (must be ascending)
///
/// # Returns
///
/// Minimum spacing between consecutive points
fn find_min_spacing<T>(x: &Array1<T>) -> Result<T, Report>
where
  T: Float + Debug,
{
  if x.len() < 2 {
    return make_error!("Array must have at least 2 points");
  }

  let diffs = &x.slice(s![1..]) - &x.slice(s![..-1]);

  if diffs.iter().any(|&dx| dx <= T::zero()) {
    return make_error!("x array must be sorted in ascending order");
  }

  let min_dx = *diffs.min().unwrap();

  if !min_dx.is_finite() || min_dx <= T::zero() {
    return make_error!("Invalid spacing in input array: {min_dx:?}");
  }

  Ok(min_dx)
}

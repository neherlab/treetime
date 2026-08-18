use std::fmt::Debug;

use crate::InterpElem;
use crate::boundary_behavior::BoundaryBehavior;
use crate::grid::Grid;
use crate::interp_nonuniform::interp_nonuniform;
use crate::soft_tail_law::SoftTailLaw;
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
  /// Construct a fresh grid function from a grid and matching y-array.
  ///
  /// Both tails default to [`BoundaryBehavior::Error`]: a raw grid and y-array carry no declared
  /// out-of-support policy, so evaluating beyond the grid is a programming error until a caller
  /// opts in with [`Self::with_left_extrap`] / [`Self::with_right_extrap`]. Defaulting to a soft
  /// tail here would fabricate a boundary law the data never declared and silently corrupt the
  /// quantile and HPD integrals that depend on the declared tail.
  ///
  /// This is the *fresh-construction* entry point. Regridding an existing function is different:
  /// the per-side policy and any fitted boundary law are properties of the distribution, not of
  /// the grid, so they must survive every regrid. That carry lives in [`Self::regridded`], which
  /// every regridding method routes through, so a regrid can never silently drop the policy by
  /// forgetting to re-apply it.
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

  /// Rebuild this function on a new grid and y-array, carrying the per-side boundary policy and
  /// any fitted boundary law across the regrid.
  ///
  /// [`Self::from_grid_array`] deliberately resets both tails to [`BoundaryBehavior::Error`]
  /// because it constructs a fresh function. Regridding is the opposite case: the declared domain
  /// (hard versus soft) and the fitted [`SoftTailLaw`] / [`HardApproachLaw`](crate::HardApproachLaw) describe the
  /// distribution, not the grid, so they must be preserved. Centralizing the carry here makes
  /// preservation structural: every regridding method (`resample`, and any future addition)
  /// builds through this helper and therefore cannot lose the policy.
  ///
  /// Both the hard approach law and the soft-tail law are edge-relative: they store only shape (the
  /// exponent `b`, the slope) and read the resampled edge ordinate on evaluation, so both stay valid
  /// across a regrid without refitting. Refitting is required only after operations that change the
  /// tail's shape (convolution), never after a pure regrid.
  fn regridded(&self, grid: Grid<T>, y: Array1<T>) -> Result<Self, Report> {
    Ok(
      Self::from_grid_array(grid, y)?
        .with_left_extrap(self.left_extrap)
        .with_right_extrap(self.right_extrap),
    )
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
  /// returns `0.0`, `HardApproach` follows the fitted approach law up to its hard boundary, and
  /// `Linear` continues the density along its fitted log-linear tail.
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
      // Soft log-linear tail: a straight neg-log line anchored on the live grid edge, evaluated on
      // the stored ordinate axis so it meets the grid continuously at `boundary_value`.
      BoundaryBehavior::Linear(law) => {
        let value = law.eval(
          boundary_value.to_f64().unwrap(),
          bound.to_f64().unwrap(),
          xi.to_f64().unwrap(),
        );
        Ok(T::from(value).unwrap())
      },
      BoundaryBehavior::HardApproach(law) => {
        let xi_f64 = xi.to_f64().unwrap();
        // Beyond the hard boundary: zero probability
        if (side == "below" && xi_f64 < law.t_hard) || (side != "below" && xi_f64 > law.t_hard) {
          return Ok(T::zero());
        }
        // Between hard boundary and grid edge: use the edge-relative approach law, anchored on the
        // live grid edge (`boundary_value` at `bound`), exactly like the soft-tail arm above.
        let value = law.eval(boundary_value.to_f64().unwrap(), bound.to_f64().unwrap(), xi_f64);
        Ok(T::from(value).unwrap())
      },
      BoundaryBehavior::Hard => Ok(T::zero()),
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
      // An arbitrary y-transform invalidates any fitted boundary law and can move the tail off the
      // log-linear family the laws describe, so the result carries no declared tail: both sides
      // reset to `Error`, exactly as a fresh `from_grid_array` does. A caller that still wants a tail
      // re-declares it (and refits) with a known transform (`scale_y`, `shift_y`) or an explicit
      // `with_*_extrap`. Structure-preserving transforms keep their law through their own methods.
      left_extrap: BoundaryBehavior::default(),
      right_extrap: BoundaryBehavior::default(),
    }
  }

  /// Scale all y-values by a multiplicative factor, preserving boundary laws.
  ///
  /// Under NegLog storage, multiplying stored ordinates by `factor` is a pointwise scale of the
  /// neg-log values (`p -> p^factor`). The hard approach exponent `b` scales by `factor`, and a
  /// soft-tail slope scales by `factor`: the stored ordinates steepen by `factor`, so both neg-log
  /// laws do as well, while the edge-relative tails read the already-scaled edge value.
  #[must_use]
  pub fn scale_y(&self, factor: f64) -> Self
  where
    T: Float,
  {
    Self {
      grid: self.grid,
      y: self.y.mapv(|v| v * T::from(factor).unwrap()),
      left_extrap: scale_tail_law(self.left_extrap, factor),
      right_extrap: scale_tail_law(self.right_extrap, factor),
    }
  }

  /// Add a constant `delta` to every y-value, preserving boundary policy and fitted laws.
  ///
  /// This is the additive counterpart of [`Self::scale_y`]. Under `NegLog` storage the ordinate is
  /// `-ln(probability)`, so adding a constant is normalization by a pure shift (subtracting the
  /// peak ordinate moves the mode to zero): exact, no scaling, and it preserves likelihood ratios.
  ///
  /// Both fitted laws carry through unchanged. They are edge-relative: a hard approach law
  /// `y = y_edge - b·ln|Δt/Δt_edge|` and a soft-tail law both store only shape (the exponent `b`,
  /// the slope `d(-ln p)/dt`), which is invariant under a vertical shift, while evaluation reads the
  /// shifted edge ordinate.
  ///
  /// Unlike the arbitrary [`Self::mapv`], a shift is a known transform whose effect on both laws is
  /// closed form, so the law is updated rather than dropped.
  #[must_use]
  pub fn shift_y(&self, delta: T) -> Self
  where
    T: Float,
  {
    let delta_f64 = delta.to_f64().unwrap();
    Self {
      grid: self.grid,
      y: self.y.mapv(|v| v + delta),
      left_extrap: shift_tail_law(self.left_extrap, delta_f64),
      right_extrap: shift_tail_law(self.right_extrap, delta_f64),
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

    // The reflection swaps the left and right tails and reflects each fitted law's argument: a
    // hard approach law negates its `t_hard`, a soft-tail slope flips sign.
    self.left_extrap = negate_tail_law(self.left_extrap);
    self.right_extrap = negate_tail_law(self.right_extrap);
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
    self.regridded(*grid, Array1::from_vec(y_new))
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

  /// Resamples onto a uniform grid over `x_range`, clamping any target point that grid-construction
  /// rounding pushes marginally outside this function's own support back to the nearest boundary.
  ///
  /// [`Grid::from_range_dx`] rounds the point count to the nearest integer, so the final target
  /// point can land up to `dx / 2` beyond `x_range.1`. When a caller regrids a function onto its own
  /// support (or a sub-window of it), that overshoot is a gridding artifact, not a genuine
  /// out-of-support query: clamping the query into `[x_min, x_max]` reads the boundary value there,
  /// exactly as evaluating the function at its own endpoint would. Unlike [`Self::resample`] this
  /// never consults the per-side tail policy, so it neither errors on an `Error` tail nor invents
  /// density from a soft one; the caller restores the real tails on the result.
  pub fn resample_range_dx_clamped(&self, x_range: (T, T), dx: T) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid = Grid::from_range_dx(x_range.0, x_range.1, dx)?;
    let x_min = self.grid.x_min();
    let x_max = self.grid.x_max();
    let y_new = (0..grid.n_points())
      .map(|i| self.interp(grid.x_at(i).max(x_min).min(x_max)))
      .collect::<Result<Vec<T>, Report>>()?;
    self.regridded(grid, Array1::from_vec(y_new))
  }
}

/// Scale a fitted boundary law when every ordinate is multiplied by `factor`. Scaling every neg-log
/// ordinate by `factor` raises the probability to a power (`p -> p^factor`), so both shape terms of
/// the hard approach law scale by `factor`: `y = y_edge - b*ln|dt/dt_edge| + slope*(t - t_edge)`
/// becomes `factor*y = factor*y_edge - (factor*b)*ln|dt/dt_edge| + (factor*slope)*(t - t_edge)`, with
/// the edge read live. A soft-tail slope scales with the ordinates the same way: the neg-log tail
/// line `y_edge + slope*(t - t_edge)` steepens by `factor` (see [`GridFn::scale_y`]).
fn scale_tail_law(behavior: BoundaryBehavior, factor: f64) -> BoundaryBehavior {
  match behavior {
    BoundaryBehavior::HardApproach(law) => BoundaryBehavior::HardApproach(law.scale(factor)),
    BoundaryBehavior::Linear(law) => BoundaryBehavior::Linear(SoftTailLaw {
      slope: law.slope * factor,
    }),
    other => other,
  }
}

/// Shift a fitted boundary law when a constant `delta` is added to every ordinate. Both the hard
/// approach law and the soft-tail law are edge-relative: they read the shifted edge ordinate on
/// evaluation and store only shape (the exponent `b`, the slope), which is invariant under a
/// vertical shift. So both carry through unchanged (see [`GridFn::shift_y`]); the `delta` argument
/// is retained for symmetry with [`scale_tail_law`].
fn shift_tail_law(behavior: BoundaryBehavior, _delta: f64) -> BoundaryBehavior {
  behavior
}

/// Reflect a fitted boundary law for `f(x) -> f(-x)`.
fn negate_tail_law(behavior: BoundaryBehavior) -> BoundaryBehavior {
  match behavior {
    BoundaryBehavior::HardApproach(law) => BoundaryBehavior::HardApproach(law.negate_arg()),
    BoundaryBehavior::Linear(law) => BoundaryBehavior::Linear(law.negate_arg()),
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

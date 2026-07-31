use crate::policy::{Plain, PolicyMarker, YAxisPolicy};
use approx::UlpsEq;
use eyre::Report;
use ndarray::{Array1, Axis, concatenate};
use ndarray_stats::QuantileExt;
use num::Float;
use serde::{Deserialize, Serialize};
use treetime_grid::grid::Grid;
use treetime_grid::{BoundaryBehavior, GridFn, InterpElem};
use treetime_utils::array::ndarray::reverse;
use treetime_utils::{make_error, make_report};

/// Side of a distribution's grid on which a tail can be extended outward.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TailSide {
  /// Below `x_min` (toward smaller x).
  Left,
  /// Above `x_max` (toward larger x).
  Right,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DistributionFunction<T: InterpElem, Y: YAxisPolicy = Plain> {
  grid_fn: GridFn<T>,
  #[serde(skip)]
  _policy: PolicyMarker<Y>,
}

impl<T: InterpElem, Y: YAxisPolicy> DistributionFunction<T, Y> {
  pub fn from_arrays(x: &Array1<T>, y: Array1<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = GridFn::from_arrays(x, y)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn from_arrays_nonuniform(x: &Array1<T>, y: &Array1<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = GridFn::from_arrays_nonuniform(x, y)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn from_range_values(x_range: (T, T), y: Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn from_start_dx_values(x_min: T, dx: T, y: Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::from_start_dx_values(x_min, dx, y)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn from_n_points<F>((x_min, x_max): (T, T), n_points: usize, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    let grid_fn = GridFn::from_n_points((x_min, x_max), n_points, y_fn)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn from_grid<F>((x_min, x_max): (T, T), dx: T, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    let grid_fn = GridFn::from_grid((x_min, x_max), dx, y_fn)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn constant((x_min, x_max): (T, T), n_points: usize, value: T) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::constant((x_min, x_max), n_points, value)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn zeros((x_min, x_max): (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::zeros((x_min, x_max), n_points)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn ones((x_min, x_max): (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::ones((x_min, x_max), n_points)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn from_grid_fn(grid_fn: GridFn<T>) -> Self {
    Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    }
  }

  pub fn t(&self) -> Array1<T>
  where
    T: Float,
  {
    self.grid_fn.x()
  }

  pub fn x_min(&self) -> T {
    self.grid_fn.x_min()
  }

  pub fn x_max(&self) -> T
  where
    T: Float,
  {
    self.grid_fn.x_max()
  }

  pub fn dx(&self) -> T {
    self.grid_fn.dx()
  }

  pub fn y(&self) -> &Array1<T> {
    self.grid_fn.y()
  }

  pub fn grid(&self) -> &Grid<T> {
    self.grid_fn.grid()
  }

  pub fn interp(&self, x: T) -> Result<T, Report>
  where
    T: Float + UlpsEq,
  {
    self.grid_fn.interp(x)
  }

  pub fn interp_many(&self, xs: &Array1<T>) -> Result<Array1<T>, Report>
  where
    T: Float + UlpsEq,
  {
    self.grid_fn.interp_many(xs)
  }

  pub fn left_extrap(&self) -> BoundaryBehavior {
    self.grid_fn.left_extrap()
  }

  pub fn right_extrap(&self) -> BoundaryBehavior {
    self.grid_fn.right_extrap()
  }

  /// Set the left (below `x_min`) out-of-support tail policy.
  ///
  /// Rejects a [`BoundaryBehavior::Zero`] tail when the representation cannot express zero
  /// probability as `0.0` (negative-log), where zero probability is `+inf`.
  pub fn with_left_extrap(self, behavior: BoundaryBehavior) -> Result<Self, Report> {
    Self::check_zero_boundary(behavior)?;
    Ok(Self::from_grid_fn(self.grid_fn.with_left_extrap(behavior)))
  }

  /// Set the right (above `x_max`) out-of-support tail policy. See [`Self::with_left_extrap`].
  pub fn with_right_extrap(self, behavior: BoundaryBehavior) -> Result<Self, Report> {
    Self::check_zero_boundary(behavior)?;
    Ok(Self::from_grid_fn(self.grid_fn.with_right_extrap(behavior)))
  }

  /// Set the same out-of-support tail policy on both sides.
  pub fn with_extrap(self, behavior: BoundaryBehavior) -> Result<Self, Report> {
    self.with_left_extrap(behavior)?.with_right_extrap(behavior)
  }

  fn check_zero_boundary(behavior: BoundaryBehavior) -> Result<(), Report> {
    if behavior == BoundaryBehavior::Zero && !Y::supports_zero_boundary() {
      return make_error!(
        "Refusing a Zero boundary tail: it writes 0.0 outside support, which is the multiplicative identity (probability one), not zero probability, under this distribution's negative-log representation"
      );
    }
    Ok(())
  }

  pub fn resample(&self, grid: &Grid<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = self.grid_fn.resample(grid)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn resample_start_dx(&self, x_min: T, dx: T, n_points: usize) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = self.grid_fn.resample_start_dx(x_min, dx, n_points)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn resample_range_n_points(&self, x_range: (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = self.grid_fn.resample_range_n_points(x_range, n_points)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn resample_range_dx(&self, x_range: (T, T), dx: T) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = self.grid_fn.resample_range_dx(x_range, dx)?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn resample_dx(&self, dx: T) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    // Re-grid onto the function's own support. When dx does not divide the range evenly,
    // the final uniform grid point can fall marginally beyond x_max; hold the boundary
    // value there instead of erroring, since this is a gridding artifact, not a genuine
    // out-of-support query. The resampled result keeps this function's own tail policy.
    let regridded = self.grid_fn.clone().with_extrap(BoundaryBehavior::Constant);
    let resampled = regridded.resample_range_dx((self.x_min(), self.x_max()), dx)?;
    let resampled = resampled
      .with_left_extrap(self.grid_fn.left_extrap())
      .with_right_extrap(self.grid_fn.right_extrap());
    Ok(Self::from_grid_fn(resampled))
  }

  pub fn len(&self) -> usize {
    self.grid_fn.len()
  }

  pub fn is_empty(&self) -> bool {
    self.grid_fn.len() == 0
  }

  pub fn negate_arg(&self) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = self.grid_fn.negate_arg()?;
    Ok(Self {
      grid_fn,
      _policy: PolicyMarker::new(),
    })
  }

  pub fn negate_arg_inplace(&mut self) -> Result<(), Report>
  where
    T: Float,
  {
    self.grid_fn.negate_arg_inplace()
  }

  /// Find the most likely time point (x-value corresponding to maximum y-value)
  /// Returns the x-coordinate where the function reaches its maximum y-value
  pub fn likely_time(&self) -> Option<T>
  where
    T: Float,
  {
    let y_values = self.y();
    if let Ok(max_idx) = y_values.argmax() {
      Some(self.t()[max_idx])
    } else {
      None
    }
  }

  /// Create a new distribution function with y values scaled by factor.
  ///
  /// Preserves the grid parameters exactly (via `fn GridFn.mapv()`), avoiding floating point
  /// issues that can occur when regenerating the x array, and preserves the per-side tail
  /// policies. Scaling by a positive factor does not change out-of-support behavior, so the
  /// tails carry through. This makes `fn Distribution.normalize()` tail-preserving.
  pub fn scale_y(&self, factor: T) -> Result<Self, Report>
  where
    T: Float,
  {
    Ok(Self::from_grid_fn(self.grid_fn.mapv(|v| v * factor)))
  }
}

impl<T: InterpElem> DistributionFunction<T, Plain> {
  /// Extend the grid outward on `side` with a log-linear (exponential-in-probability) tail.
  ///
  /// Branch-length likelihoods decay exponentially for long branches: under a Poisson
  /// substitution model the branch contribution is $\mathcal{L}(t) \propto e^{-\mu t}(\mu t)^k$,
  /// whose logarithm is asymptotically linear in $t$, so beyond the computed grid the density
  /// should continue that decay rather than hold flat. See Felsenstein 1981
  /// (<https://doi.org/10.1007/bf01734359>).
  ///
  /// The outward decay rate is estimated in log space from the outermost `min(3, n/3)` points and
  /// the density is continued as `f = f_boundary * exp(-rate * distance)`. Points are appended at
  /// the existing spacing until the extrapolated value falls below `rel_floor * peak`, the number
  /// of appended points reaches `max_added_points`, or the total grid reaches the 1e6-point cap.
  /// The tail is extended only when the density strictly decays outward on that side (boundary
  /// value below the interior anchor); a flat or rising boundary is returned unchanged. Both
  /// out-of-support tail policies are preserved.
  ///
  /// Restricted to the [`Plain`] axis: the extrapolation takes logarithms of the ordinate, which
  /// is a probability density here, not a negative-log value.
  pub fn extend_log_linear_tail(&self, side: TailSide, rel_floor: T, max_added_points: usize) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let y = self.y();
    let n = y.len();
    let dx = self.dx();
    if n < 2 || max_added_points == 0 || !dx.is_finite() || dx <= T::zero() {
      return Ok(self.clone());
    }

    let margin = (n / 3).clamp(1, 3);
    let (f_boundary, f_anchor, boundary_x) = match side {
      TailSide::Right => (y[n - 1], y[n - 1 - margin], self.x_max()),
      TailSide::Left => (y[0], y[margin], self.x_min()),
    };

    let f_peak = *y
      .max()
      .map_err(|e| make_report!("Tail extension requires a non-empty grid: {e}"))?;
    let floor = rel_floor * f_peak;
    // Log-linear extrapolation needs positive, finite endpoints, and extension only makes sense
    // when the density is strictly lower at the boundary than at the interior anchor (decaying
    // outward) and the boundary is still above the negligibility floor.
    let can_extend = f_boundary.is_finite()
      && f_anchor.is_finite()
      && f_boundary > T::zero()
      && f_anchor > f_boundary
      && f_boundary > floor;
    if !can_extend {
      return Ok(self.clone());
    }

    // Outward decay rate per unit x (positive by the check above).
    let rate = (f_anchor.ln() - f_boundary.ln()) / (T::from(margin).unwrap() * dx);
    if !rate.is_finite() || rate <= T::zero() {
      return Ok(self.clone());
    }

    let cap = max_added_points.min(1_000_000_usize.saturating_sub(n));
    let mut tail = Vec::new();
    for k in 1..=cap {
      let distance = T::from(k).unwrap() * dx;
      let value = f_boundary * (-rate * distance).exp();
      if value < floor {
        break;
      }
      tail.push(value);
    }
    if tail.is_empty() {
      return Ok(self.clone());
    }
    let tail = Array1::from_vec(tail);

    let (new_x_min, new_y) = match side {
      TailSide::Right => (
        self.x_min(),
        concatenate(Axis(0), &[y.view(), tail.view()]).map_err(|e| make_report!("Tail concatenation failed: {e}"))?,
      ),
      TailSide::Left => {
        let tail = reverse(&tail);
        let new_x_min = boundary_x - T::from(tail.len()).unwrap() * dx;
        (
          new_x_min,
          concatenate(Axis(0), &[tail.view(), y.view()]).map_err(|e| make_report!("Tail concatenation failed: {e}"))?,
        )
      },
    };

    let grid_fn = GridFn::from_start_dx_values(new_x_min, dx, new_y)?
      .with_left_extrap(self.left_extrap())
      .with_right_extrap(self.right_extrap());
    Ok(Self::from_grid_fn(grid_fn))
  }
}

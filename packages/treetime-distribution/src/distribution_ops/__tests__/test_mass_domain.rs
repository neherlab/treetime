#[cfg(test)]
mod tests {
  use crate::distribution_ops::mass_domain::{mass_bounded_domain, rewindow_to_mass, total_mass};
  use crate::policy::NegLog;
  use crate::{Distribution, DistributionFunction};
  use approx::assert_abs_diff_eq;
  use ndarray::Array1;
  use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, Side, SoftTailLaw};

  /// Per-side tail mass fraction used across these tests (design D4, proposal Axis 2).
  const EPS: f64 = 5e-4;
  const GRID_POINTS: usize = 300;

  /// Total mass of an `Exponential(lambda)` distribution is `1/lambda`.
  ///
  /// Oracle: analytic `integral_0^inf lambda*exp(-lambda t) dt / lambda = 1/lambda` for the
  /// peak-relative density `exp(-lambda t)` (the neg-log peak sits at `t = 0`). The grid trapezoid
  /// plus the closed-form right-tail mass must recover it.
  #[test]
  fn test_mass_domain_total_mass_exponential_matches_analytic() {
    let lambda = 1.0;
    let f = helpers::exponential_neglog(lambda, 5.0, 8001);
    // 1/lambda = 1.0; trapezoid + tail closed form recover it to grid-refinement order.
    assert_abs_diff_eq!(total_mass(&f).unwrap(), 1.0 / lambda, epsilon = 1e-6);
  }

  /// The right mass edge of an `Exponential(lambda)` truncated below its analytic `eps` quantile
  /// extends into the fitted tail to `-ln(eps)/lambda`; the hard left edge stays at `0`.
  ///
  /// Oracle: `integral_hi^inf exp(-lambda t) dt = eps * Z` with `Z = 1/lambda` gives
  /// `hi = -ln(eps)/lambda`. Here `t_max = 5 < hi`, so the edge is found in the closed-form tail.
  #[test]
  fn test_mass_domain_edge_extends_into_tail_to_analytic_quantile() {
    let lambda = 1.0;
    let f = helpers::exponential_neglog(lambda, 5.0, 8001);
    let (lo, hi) = mass_bounded_domain(&f, EPS).unwrap();
    assert_abs_diff_eq!(lo, 0.0, epsilon = 1e-12);
    assert_abs_diff_eq!(hi, -EPS.ln() / lambda, epsilon = 1e-6);
  }

  /// The same analytic edge is recovered by trapezoid-CDF inversion when the grid already extends
  /// past the `eps` quantile (`t_max = 10 > hi`), so the edge is trimmed inward rather than extended.
  ///
  /// Oracle: same `hi = -ln(eps)/lambda`.
  #[test]
  fn test_mass_domain_edge_trims_inward_to_analytic_quantile() {
    let lambda = 1.0;
    let f = helpers::exponential_neglog(lambda, 10.0, 16001);
    let (lo, hi) = mass_bounded_domain(&f, EPS).unwrap();
    assert_abs_diff_eq!(lo, 0.0, epsilon = 1e-12);
    assert_abs_diff_eq!(hi, -EPS.ln() / lambda, epsilon = 1e-6);
  }

  /// A symmetric two-sided exponential (Laplace, slope `s`) trims both soft edges symmetrically.
  ///
  /// Oracle: `Z = 2/s`, and `integral_hi^inf exp(-s t) dt = eps * Z` gives `hi = -ln(2 eps)/s`; by
  /// symmetry `lo = -hi`.
  #[test]
  fn test_mass_domain_symmetric_edges_hold_equal_tail_mass() {
    let slope = 1.0;
    let f = helpers::laplace_neglog(slope, 12.0, 24001);
    let (lo, hi) = mass_bounded_domain(&f, EPS).unwrap();
    let expected = -(2.0 * EPS).ln() / slope;
    assert_abs_diff_eq!(hi, expected, epsilon = 1e-6);
    assert_abs_diff_eq!(lo, -expected, epsilon = 1e-6);
  }

  /// Re-windowing a fixed distribution repeatedly conserves total mass and holds the mode.
  ///
  /// Oracle: mass conservation is an invariant of the re-window (it re-fits the soft tail rather than
  /// discarding it), so drift after the fixed point is bounded by the tail-fit residual, not by
  /// arithmetic. Compares the once-windowed distribution against the 100-times-windowed one.
  #[test]
  fn test_mass_domain_rewindow_conserves_mass_and_mode() {
    let f = helpers::exponential_neglog(1.0, 5.0, 4001);
    let once = helpers::as_function(rewindow_to_mass(&Distribution::Function(f), EPS, GRID_POINTS).unwrap());

    let mass_once = total_mass(&once).unwrap();
    let mode_once = once.likely_time().unwrap();
    let dx = once.dx();

    let mut current = once;
    for _ in 0..100 {
      current = helpers::as_function(rewindow_to_mass(&Distribution::Function(current), EPS, GRID_POINTS).unwrap());
    }

    let mass_100 = total_mass(&current).unwrap();
    let mode_100 = current.likely_time().unwrap();
    assert_abs_diff_eq!(mass_100, mass_once, epsilon = 1e-6 * mass_once);
    assert!(
      (mode_100 - mode_once).abs() <= dx,
      "mode moved {} beyond one spacing {dx}",
      (mode_100 - mode_once).abs()
    );
  }

  /// A distribution whose mode sits on a finite hard boundary keeps that mode on the lower edge after
  /// re-windowing, with no `+inf` ordinate stored.
  ///
  /// Oracle: a nullary `Hard` lower boundary is the exact bound, so re-windowing leaves the lower edge
  /// on it (probability is zero beyond, no `eps` trim) and the mode stays on the first grid point.
  #[test]
  fn test_mass_domain_rewindow_keeps_mode_on_hard_bound() {
    let f = helpers::hard_bound_mode_neglog(1.0, 5.0, 500);
    let rewindowed = helpers::as_function(rewindow_to_mass(&Distribution::Function(f), EPS, GRID_POINTS).unwrap());

    assert_abs_diff_eq!(rewindowed.x_min(), 0.0, epsilon = 1e-12);
    let mode_index =
      rewindowed
        .y()
        .iter()
        .copied()
        .enumerate()
        .fold(
          (0_usize, f64::INFINITY),
          |(best_i, best), (i, y)| if y < best { (i, y) } else { (best_i, best) },
        );
    assert_eq!(mode_index.0, 0, "mode must sit on the lower hard edge");
    assert!(
      rewindowed.y().iter().all(|y| y.is_finite()),
      "no +inf ordinate may be stored"
    );
  }

  /// The stored grid holds at least `grid_points` points (design D3 resolution floor).
  #[test]
  fn test_mass_domain_rewindow_holds_at_least_grid_points() {
    let f = helpers::exponential_neglog(1.0, 5.0, 50);
    let rewindowed = helpers::as_function(rewindow_to_mass(&Distribution::Function(f), EPS, GRID_POINTS).unwrap());
    assert!(rewindowed.len() >= GRID_POINTS, "got {} points", rewindowed.len());
  }

  /// A non-`Function` distribution has no grid to re-window, so re-windowing is the shift-only
  /// normalize: a `Point` keeps its location and its amplitude drops to the peak identity `0`.
  #[test]
  fn test_mass_domain_rewindow_point_is_shift_only_normalize() {
    let point: Distribution<NegLog> = Distribution::point(3.0, 7.5);
    let rewindowed = rewindow_to_mass(&point, EPS, GRID_POINTS).unwrap();
    assert_eq!(rewindowed, point.normalize());
  }

  mod helpers {
    use super::*;

    pub fn as_function(dist: Distribution<NegLog>) -> DistributionFunction<f64, NegLog> {
      match dist {
        Distribution::Function(f) => f,
        other => panic!("expected a Function distribution, got {other:?}"),
      }
    }

    /// `Exponential(lambda)` on `[0, t_max]` in neg-log storage: `y(t) = -ln(lambda) + lambda*t`,
    /// a hard left boundary at `t = 0`, and a fitted soft right tail (slope `lambda`).
    pub fn exponential_neglog(lambda: f64, t_max: f64, n: usize) -> DistributionFunction<f64, NegLog> {
      let t = Array1::linspace(0.0, t_max, n);
      let y = t.mapv(|ti| -lambda.ln() + lambda * ti);
      let f = DistributionFunction::from_arrays(&t, y).unwrap();
      let right = SoftTailLaw::fit(f.grid_fn(), Side::Right, DEFAULT_TAIL_FIT_POINTS).unwrap();
      f.with_left_extrap(BoundaryBehavior::Hard)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Linear(right))
        .unwrap()
    }

    /// Symmetric two-sided exponential (Laplace, slope `s`) on `[-t_max, t_max]`: `y(t) = s*|t|`,
    /// fitted soft tails on both sides.
    pub fn laplace_neglog(slope: f64, t_max: f64, n: usize) -> DistributionFunction<f64, NegLog> {
      let t = Array1::linspace(-t_max, t_max, n);
      let y = t.mapv(|ti| slope * ti.abs());
      let f = DistributionFunction::from_arrays(&t, y).unwrap();
      let left = SoftTailLaw::fit(f.grid_fn(), Side::Left, DEFAULT_TAIL_FIT_POINTS).unwrap();
      let right = SoftTailLaw::fit(f.grid_fn(), Side::Right, DEFAULT_TAIL_FIT_POINTS).unwrap();
      f.with_left_extrap(BoundaryBehavior::Linear(left))
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Linear(right))
        .unwrap()
    }

    /// A density whose mode sits on a finite hard boundary at `t = 0`: `y(t) = slope*t` on
    /// `[0, t_max]` with a nullary `Hard` left boundary (the grid edge is the exact hard bound,
    /// carrying the mode) and a fitted soft right tail.
    pub fn hard_bound_mode_neglog(slope: f64, t_max: f64, n: usize) -> DistributionFunction<f64, NegLog> {
      let t = Array1::linspace(0.0, t_max, n);
      let y = t.mapv(|ti| slope * ti);
      let f = DistributionFunction::from_arrays(&t, y).unwrap();
      let right = SoftTailLaw::fit(f.grid_fn(), Side::Right, DEFAULT_TAIL_FIT_POINTS).unwrap();
      f.with_left_extrap(BoundaryBehavior::Hard)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Linear(right))
        .unwrap()
    }
  }
}

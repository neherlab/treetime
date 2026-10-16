#[cfg(test)]
mod tests {
  use crate::optimize::method_newton::{SqrtInterval, SqrtSolution, newton_sqrt_maximize};
  use crate::optimize::solver_outcome::Convergence;
  use crate::pretty_assert_abs_diff_eq;
  use helpers::{flat, poisson, shifted_poisson};
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use treetime_utils::prop_assert_relative_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::start_above_maximum(   (3.0,  1000.0, 0.3 ))]
  #[case::start_below_maximum(   (3.0,  1000.0, 1e-5))]
  #[case::many_events(           (50.0, 2e4,    0.5 ))]
  #[case::long_branch(           (1.0,  10.0,   2.0 ))]
  #[trace]
  fn test_optimize_method_newton_safeguard_reaches_poisson_maximum(#[case] (k, sites, t0): (f64, f64, f64)) {
    let interval = SqrtInterval { lower: 1e-4, upper: 3.0 };

    let solution = newton_sqrt_maximize(t0.sqrt(), interval, 50, poisson(k, sites)).unwrap();

    assert_eq!(Convergence::Converged, solution.convergence);
    pretty_assert_abs_diff_eq!(
      (k / sites).sqrt(),
      solution.s,
      epsilon = 1e-8,
      "maximum of 2k ln(s) - L s^2 at s = sqrt(k/L) (derivative 2k/s - 2Ls = 0)"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::maximum_above_interval((0.0001, 0.01), 0.01)]
  #[case::maximum_below_interval((0.1,    1.0 ), 0.1 )]
  #[trace]
  fn test_optimize_method_newton_safeguard_returns_interval_end_nearest_an_outside_maximum(
    #[case] (lower, upper): (f64, f64),
    #[case] expected: f64,
  ) {
    let interval = SqrtInterval { lower, upper };

    let solution = newton_sqrt_maximize(f64::midpoint(lower, upper), interval, 50, poisson(3.0, 1000.0)).unwrap();

    let expected = SqrtSolution { s: expected, convergence: Convergence::Converged };
    assert_eq!(expected, solution, "maximum s = sqrt(3/1000) = 0.0548 lies outside [{lower}, {upper}]");
  }

  #[test]
  fn test_optimize_method_newton_safeguard_leaves_zero_when_zero_is_a_minimum() {
    let interval = SqrtInterval { lower: 0.0, upper: 1.0 };

    let solution = newton_sqrt_maximize(0.0, interval, 50, shifted_poisson(3.0, 1000.0, 1e-3)).unwrap();

    assert_eq!(Convergence::Converged, solution.convergence);
    pretty_assert_abs_diff_eq!(
      (3.0_f64 / 1000.0 - 1e-3).sqrt(),
      solution.s,
      epsilon = 1e-8,
      "maximum of k ln(t + a) - L t at t = k/L - a; at s = 0 the slope is 0 and the curvature 2k/a - 2L > 0"
    );
  }

  #[test]
  fn test_optimize_method_newton_safeguard_returns_start_on_flat_surface() {
    let interval = SqrtInterval { lower: 0.0, upper: 1.0 };

    let solution = newton_sqrt_maximize(0.2, interval, 50, flat).unwrap();

    let expected = SqrtSolution {
      s: 0.2,
      convergence: Convergence::Converged,
    };
    assert_eq!(expected, solution, "no evaluated point is better than the start");
  }

  #[test]
  fn test_optimize_method_newton_safeguard_reports_iteration_limit() {
    let interval = SqrtInterval {
      lower: 1e-4,
      upper: 3.0,
    };

    let solution = newton_sqrt_maximize(1.0, interval, 1, poisson(3.0, 1000.0)).unwrap();

    assert_eq!(Convergence::IterationLimit, solution.convergence);
  }

  mod generators {
    use proptest::prelude::*;

    pub(super) fn gen_events() -> impl Strategy<Value = f64> {
      1.0_f64..1e3_f64
    }

    pub(super) fn gen_sites() -> impl Strategy<Value = f64> {
      1e2_f64..1e6_f64
    }

    pub(super) fn gen_start() -> impl Strategy<Value = f64> {
      1e-4_f64..2.0_f64
    }

    pub(super) fn gen_max_iter() -> impl Strategy<Value = usize> {
      0_usize..6_usize
    }
  }

  proptest! {
    #[test]
    fn test_prop_optimize_method_newton_safeguard_never_worse_than_start(
      k in generators::gen_events(),
      sites in generators::gen_sites(),
      s0 in generators::gen_start(),
      max_iter in generators::gen_max_iter(),
    ) {
      let interval = SqrtInterval { lower: 1e-4, upper: 3.0 };
      let evaluate = poisson(k, sites);

      let solution = newton_sqrt_maximize(s0, interval, max_iter, &evaluate).unwrap();

      let log_lh_start = evaluate(s0).unwrap().log_lh;
      let log_lh_result = evaluate(solution.s).unwrap().log_lh;
      prop_assert!(
        log_lh_result >= log_lh_start,
        "result s={} has log_lh {log_lh_result} below the start s={s0} with {log_lh_start}", solution.s
      );
    }

    #[test]
    fn test_prop_optimize_method_newton_safeguard_matches_poisson_maximum(
      k in generators::gen_events(),
      sites in generators::gen_sites(),
      s0 in generators::gen_start(),
    ) {
      let interval = SqrtInterval { lower: 1e-4, upper: 3.0 };

      let solution = newton_sqrt_maximize(s0, interval, 50, poisson(k, sites)).unwrap();

      prop_assert_eq!(Convergence::Converged, solution.convergence);
      let expected = (k / sites).sqrt().clamp(interval.lower, interval.upper);
      prop_assert_relative_eq!(expected, solution.s, max_relative = 1e-6);
    }
  }

  mod helpers {
    use crate::optimize::method_newton::SqrtPoint;
    use eyre::Report;

    pub(super) fn poisson(k: f64, sites: f64) -> impl Fn(f64) -> Result<SqrtPoint, Report> {
      move |s| {
        Ok(SqrtPoint {
          log_lh: 2.0 * k * s.ln() - sites * s * s,
          slope: 2.0 * k / s - 2.0 * sites * s,
          curvature: -2.0 * k / (s * s) - 2.0 * sites,
        })
      }
    }

    pub(super) fn shifted_poisson(k: f64, sites: f64, shift: f64) -> impl Fn(f64) -> Result<SqrtPoint, Report> {
      move |s| {
        let u = s * s + shift;
        Ok(SqrtPoint {
          log_lh: k * u.ln() - sites * s * s,
          slope: 2.0 * k * s / u - 2.0 * sites * s,
          curvature: 2.0 * k * (shift - s * s) / (u * u) - 2.0 * sites,
        })
      }
    }

    pub(super) fn flat(_s: f64) -> Result<SqrtPoint, Report> {
      Ok(SqrtPoint {
        log_lh: 0.0,
        slope: 0.0,
        curvature: 0.0,
      })
    }
  }
}

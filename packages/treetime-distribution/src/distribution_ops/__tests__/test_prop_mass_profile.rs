#[cfg(test)]
mod tests {
  use crate::distribution_ops::mass_profile::{MassProfile, TailMass};
  use crate::policy::NegLog;
  use crate::{BoundaryBehavior, Distribution, DistributionFunction};
  use itertools::Itertools;
  use ndarray::Array1;
  use proptest::prelude::*;
  use treetime_grid::{GridEdge, HardApproachLaw, SoftTailLaw};
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{prop_assert_abs_diff_eq, prop_assert_ulps_eq};

  proptest! {
    #[test]
    fn test_prop_mass_profile_cdf_of_quantile_is_identity(
      f in generators::function(true),
      q in 0.0..=1.0_f64,
    ) {
      let profile = MassProfile::new(&f).unwrap();
      let fractions = profile.fractions(TailMass::Required).unwrap();
      let (below, above) = helpers::ulps_around(&f, fractions.quantile(q), 8.0);
      prop_assert!(
        fractions.cdf(below) - 1e-15 <= q && q <= fractions.cdf(above) + 1e-15,
        "cdf over [{below}, {above}] around quantile({q}) does not bracket {q}"
      );
      let (below, above) = helpers::ulps_around(&f, fractions.survival_inverse(q), 8.0);
      prop_assert!(
        fractions.survival(above) - 1e-15 <= q && q <= fractions.survival(below) + 1e-15,
        "survival over [{below}, {above}] around survival_inverse({q}) does not bracket {q}"
      );
    }

    #[test]
    fn test_prop_mass_profile_quantile_of_cdf_is_identity(
      f in generators::function(true),
      at in 0.0..=1.0_f64,
    ) {
      let profile = MassProfile::new(&f).unwrap();
      let fractions = profile.fractions(TailMass::Required).unwrap();
      let x = f.x_min() + at * (f.x_max() - f.x_min());
      let from_cdf = fractions.quantile(fractions.cdf(x));
      let from_survival = fractions.survival_inverse(fractions.survival(x));
      let tolerance = 1e-15 + 4.0 * (x.next_up() - x) / fractions.total();
      for (name, inverse) in [("quantile(cdf", from_cdf), ("survival_inverse(survival", from_survival)] {
        let between = fractions.mass_between(x.min(inverse), x.max(inverse)).unwrap() / fractions.total();
        prop_assert!(between <= tolerance, "{name}({x})) = {inverse} holds mass {between} between them");
      }
    }

    #[test]
    fn test_prop_mass_profile_cdf_is_monotone_within_unit_interval(f in generators::function(true)) {
      let profile = MassProfile::new(&f).unwrap();
      let fractions = profile.fractions(TailMass::Required).unwrap();
      let xs = Array1::linspace(f.x_min() - 3.0, f.x_max() + 3.0, 997);
      let cdf = xs.mapv(|x| fractions.cdf(x));
      let survival = xs.mapv(|x| fractions.survival(x));
      prop_assert!(cdf.iter().all(|c| (0.0..=1.0).contains(c)), "cdf outside [0, 1]: {cdf}");
      prop_assert!(cdf.iter().tuple_windows().all(|(a, b)| a <= b), "cdf decreases: {cdf}");
      prop_assert!(survival.iter().tuple_windows().all(|(a, b)| a >= b), "survival increases: {survival}");
      for (c, s) in cdf.iter().zip(&survival) {
        prop_assert_abs_diff_eq!(1.0, c + s, epsilon = 1e-15);
      }
    }

    #[test]
    fn test_prop_mass_profile_fractions_do_not_change_under_ordinate_shift(
      f in generators::function(false),
      shift in -50_i32..50,
      q in 0.0..=1.0_f64,
      at in -0.2..=1.2_f64,
    ) {
      let f = helpers::dyadic_ordinates(&f);
      let shifted = f.shift_y(f64::from(shift));
      let profile = MassProfile::new(&f).unwrap();
      let profile_shifted = MassProfile::new(&shifted).unwrap();
      let fractions = profile.fractions(TailMass::Required).unwrap();
      let fractions_shifted = profile_shifted.fractions(TailMass::Required).unwrap();
      let x = f.x_min() + at * (f.x_max() - f.x_min());
      prop_assert_ulps_eq!(fractions.cdf(x), fractions_shifted.cdf(x), max_ulps = 0);
      prop_assert_ulps_eq!(fractions.quantile(q), fractions_shifted.quantile(q), max_ulps = 0);
      prop_assert_ulps_eq!(fractions.survival_inverse(q), fractions_shifted.survival_inverse(q), max_ulps = 0);
    }

    #[test]
    fn test_prop_mass_profile_mass_between_is_additive(
      f in generators::function(true),
      cuts in proptest::array::uniform3(-0.3..=1.3_f64),
    ) {
      let profile = MassProfile::new(&f).unwrap();
      let fractions = profile.fractions(TailMass::Required).unwrap();
      let [a, b, c] = helpers::sorted_positions(&f, cuts);
      let whole = profile.mass_between(a, c).unwrap() / fractions.total();
      let parts = (profile.mass_between(a, b).unwrap() + profile.mass_between(b, c).unwrap()) / fractions.total();
      prop_assert_abs_diff_eq!(whole, parts, epsilon = 1e-15);
      let grid_point_rounding = helpers::position_ulp(&f) / fractions.total();
      prop_assert!(
        (fractions.cdf(c) - fractions.cdf(a) - whole).abs() <= 1e-15 + 4.0 * grid_point_rounding,
        "cdf({c}) - cdf({a}) differs from the mass between them {whole}"
      );
    }

    #[test]
    fn test_prop_log_mass_without_hard_cells_is_trapezoid_sum_with_tail_laws(
      f in generators::function(false),
      amplitude in -20.0..20.0_f64,
    ) {
      let f = f.with_amplitude(amplitude);
      let y_peak = f.y().iter().copied().fold(f64::INFINITY, f64::min);
      let p = f.y().mapv(|y| (y_peak - y).exp());
      let n = p.len();
      let trapezoid = f.dx() * (p.sum() - f64::midpoint(p[0], p[n - 1]));
      let left = helpers::tail_mass(f.left_extrap(), GridEdge { t: f.x_min(), y: f.y()[0] - y_peak });
      let right = helpers::tail_mass(f.right_extrap(), GridEdge { t: f.x_max(), y: f.y()[n - 1] - y_peak });
      let expected = -amplitude - y_peak + (trapezoid + left + right).ln();
      let actual = Distribution::Function(f).log_mass().unwrap();
      prop_assert_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    }
  }

  mod helpers {
    use super::*;

    pub(super) fn tail_mass(behavior: BoundaryBehavior, edge: GridEdge) -> f64 {
      match behavior {
        BoundaryBehavior::Linear(law) => (-edge.y).exp() / law.slope.abs(),
        BoundaryBehavior::HardApproach(law) => (-edge.y).exp() * (edge.t - law.t_hard).abs() / (law.b + 1.0),
        BoundaryBehavior::Hard | BoundaryBehavior::HardEdgeCell(_) | BoundaryBehavior::Error => 0.0,
      }
    }

    pub(super) fn position_ulp(f: &DistributionFunction<f64, NegLog>) -> f64 {
      let scale = f.x_min().abs().max(f.x_max().abs());
      scale.next_up() - scale
    }

    pub(super) fn ulps_around(f: &DistributionFunction<f64, NegLog>, x: f64, ulps: f64) -> (f64, f64) {
      let scale = x.abs().max(f.dx()).max(f.x_min().abs()).max(f.x_max().abs());
      let ulp = scale.next_up() - scale;
      (x - ulps * ulp, x + ulps * ulp)
    }

    pub(super) fn dyadic_ordinates(f: &DistributionFunction<f64, NegLog>) -> DistributionFunction<f64, NegLog> {
      let y = f.y().mapv(|y| (y * 1_048_576.0).round() / 1_048_576.0);
      DistributionFunction::from_start_dx_values(f.x_min(), f.dx(), y)
        .unwrap()
        .with_left_extrap(f.left_extrap())
        .unwrap()
        .with_right_extrap(f.right_extrap())
        .unwrap()
    }

    pub(super) fn sorted_positions(f: &DistributionFunction<f64, NegLog>, cuts: [f64; 3]) -> [f64; 3] {
      let width = f.x_max() - f.x_min();
      let mut xs = cuts.map(|c| f.x_min() + c * width);
      xs.sort_by(f64::total_cmp);
      xs
    }
  }

  mod generators {
    use super::*;

    pub(super) fn function(hard_cells: bool) -> impl Strategy<Value = DistributionFunction<f64, NegLog>> {
      let ordinate = if hard_cells {
        prop_oneof![9 => 0.0..30.0_f64, 1 => Just(f64::INFINITY)].boxed()
      } else {
        (0.0..30.0_f64).boxed()
      };
      (
        -10.0..10.0_f64,
        0.01..5.0_f64,
        proptest::collection::vec(ordinate, 2..20).prop_filter("a finite peak", |y| y.iter().any(|y| y.is_finite())),
        tail(),
        tail(),
      )
        .prop_map(|(x_min, dx, y, left, right)| {
          let n = y.len();
          let f = DistributionFunction::from_start_dx_values(x_min, dx, Array1::from_iter(y)).unwrap();
          let left = match left {
            BoundaryBehavior::Linear(law) => BoundaryBehavior::Linear(SoftTailLaw { slope: -law.slope }),
            BoundaryBehavior::HardApproach(law) => BoundaryBehavior::HardApproach(HardApproachLaw {
              t_hard: x_min - law.t_hard,
              b: law.b,
            }),
            other => other,
          };
          let right = match right {
            BoundaryBehavior::HardApproach(law) => BoundaryBehavior::HardApproach(HardApproachLaw {
              t_hard: x_min + dx * (n - 1).to_f64_exact() + law.t_hard,
              b: law.b,
            }),
            other => other,
          };
          f.with_left_extrap(left).unwrap().with_right_extrap(right).unwrap()
        })
    }

    fn tail() -> impl Strategy<Value = BoundaryBehavior> {
      prop_oneof![
        Just(BoundaryBehavior::Hard),
        (0.05..5.0_f64).prop_map(|slope| BoundaryBehavior::Linear(SoftTailLaw { slope })),
        (0.01..3.0_f64, 0.0..4.0_f64)
          .prop_map(|(t_hard, b)| BoundaryBehavior::HardApproach(HardApproachLaw { t_hard, b })),
      ]
    }
  }
}

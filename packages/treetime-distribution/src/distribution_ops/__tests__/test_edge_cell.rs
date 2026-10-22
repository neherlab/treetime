#[cfg(test)]
mod tests {
  use crate::distribution_ops::convolve::distribution_convolution;
  use crate::distribution_ops::edge_convolution::convolve_across_edge;
  use crate::distribution_ops::mass_domain::{mass_bounded_domain, trim_to_mass};
  use crate::distribution_ops::mass_profile::{MassProfile, TailMass};
  use crate::policy::NegLog;
  use crate::{Distribution, DistributionFunction, distribution_division, distribution_multiplication};
  use eyre::Report;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_grid::{BoundaryBehavior, HardApproachLaw, MaxGridPoints, Side};
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq, pretty_assert_relative_eq, pretty_assert_ulps_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::whole_cell(        (0.0,   1.0 ))]
  #[case::near_the_edge(     (0.0,   0.3 ))]
  #[case::near_the_anchor(   (0.3,   1.0 ))]
  #[case::at_the_offset_scale((1e-5, 2e-3))]
  #[case::across_two_cells(  (0.5,   1.5 ))]
  #[trace]
  fn test_mass_profile_edge_cell_mass_matches_quadrature_of_the_poisson_cost(#[case] (a, b): (f64, f64)) -> Result<(), Report> {
    let poisson = helpers::Poisson::one_mutation();
    let profile = MassProfile::new(&poisson.function()?)?;
    let (a, b) = (a * poisson.h, b * poisson.h);

    let expected = poisson.mass_by_quadrature(a, b.min(poisson.h)) + poisson.second_cell_mass_to(b);
    pretty_assert_relative_eq!(expected, profile.mass_between(a, b)?, max_relative = 1e-10);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::at_the_offset_scale(1e-5)]
  #[case::near_the_edge(      0.3)]
  #[case::near_the_anchor(    0.9)]
  #[trace]
  fn test_mass_fractions_in_the_edge_cell_invert_and_give_the_poisson_density(#[case] at: f64) -> Result<(), Report> {
    let poisson = helpers::Poisson::one_mutation();
    let profile = MassProfile::new(&poisson.function()?)?;
    let fractions = profile.fractions(TailMass::Required).expect("the function has a finite mass");
    let x = at * poisson.h;

    let Distribution::Function(negated) = Distribution::Function(poisson.function()?).negate()? else {
      panic!("the negation of a function is a function");
    };
    let negated_profile = MassProfile::new(&negated)?;
    let mirrored = negated_profile.fractions(TailMass::Required).expect("the function has a finite mass");

    pretty_assert_relative_eq!(x, fractions.quantile(fractions.cdf(x)), max_relative = 1e-9);
    pretty_assert_relative_eq!(-x, mirrored.survival_inverse(mirrored.survival(-x)), max_relative = 1e-9);
    pretty_assert_relative_eq!(
      (-poisson.cost(x)).exp(),
      fractions.density(x, Side::Left) * fractions.total(),
      max_relative = 1e-12
    );
    Ok(())
  }

  #[test]
  fn test_negated_edge_cell_follows_the_mirrored_law() -> Result<(), Report> {
    let poisson = helpers::Poisson::one_mutation();
    let f = poisson.function()?;
    let Distribution::Function(negated) = Distribution::Function(f.clone()).negate()? else {
      panic!("the negation of a function is a function");
    };
    let x = 0.3 * poisson.h;

    assert_eq!(
      BoundaryBehavior::HardEdgeCell(poisson.law().scale_arg(-1.0)),
      negated.right_extrap()
    );
    pretty_assert_relative_eq!(poisson.cost(x), negated.interp(-x)?, max_relative = 1e-12);
    pretty_assert_relative_eq!(
      MassProfile::new(&f)?.mass_between(0.0, x)?,
      MassProfile::new(&negated)?.mass_between(-x, 0.0)?,
      max_relative = 1e-12
    );
    Ok(())
  }

  #[test]
  fn test_scaled_edge_cell_follows_the_scaled_law() -> Result<(), Report> {
    let poisson = helpers::Poisson::one_mutation();
    let scaled = poisson.function()?.scale_time(2.0)?;
    pretty_assert_relative_eq!(
      poisson.cost(0.3 * poisson.h),
      scaled.interp(0.6 * poisson.h)?,
      max_relative = 1e-12
    );
    Ok(())
  }

  #[test]
  fn test_edge_cell_support_is_hard_at_the_grid_edge() -> Result<(), Report> {
    let poisson = helpers::Poisson::one_mutation();
    let f = poisson.function()?;
    let profile = MassProfile::new(&f)?;
    let distribution = Distribution::Function(f.clone());

    assert!(distribution.eval(-0.5 * poisson.offset)?.is_infinite());
    assert_eq!(Some((0.0, f.x_max())), distribution.support_bounds());
    pretty_assert_ulps_eq!(0.0, mass_bounded_domain(&profile, 1e-6)?.0, max_ulps = 0);
    pretty_assert_ulps_eq!(0.0, profile.mass_between(-poisson.h, 0.0)?, max_ulps = 0);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::quantile_inside_the_edge_cell(0.0001, true )]
  #[case::quantile_beyond_the_edge_cell(0.01,   false)]
  #[trace]
  fn test_mass_trim_keeps_the_edge_cell_only_at_its_edge(#[case] k: f64, #[case] keeps_edge: bool) -> Result<(), Report> {
    let poisson = helpers::Poisson::new(k, 1000);
    let Distribution::Function(trimmed) = trim_to_mass(Distribution::Function(poisson.function()?), 1e-6)? else {
      panic!("trimming a function gives a function");
    };
    let expected = if keeps_edge {
      BoundaryBehavior::HardEdgeCell(poisson.law())
    } else {
      BoundaryBehavior::Hard
    };
    assert_eq!(expected, trimmed.left_extrap());
    assert_eq!(keeps_edge, trimmed.x_min() <= 0.0);
    Ok(())
  }

  #[test]
  fn test_range_convolution_integrates_the_edge_cell() -> Result<(), Report> {
    let poisson = helpers::Poisson::one_mutation();
    let (r0, r1) = (1.0, 1.0 + 3.0 * poisson.h);
    let range = Distribution::range((r0, r1), 0.0);
    let Distribution::Function(g) = distribution_convolution(
      &range,
      &Distribution::Function(poisson.function()?),
      MaxGridPoints::default(),
    )?
    else {
      panic!("a range convolved with a function is a function");
    };

    let inside = g
      .t()
      .iter()
      .copied()
      .filter(|&t| t - r0 > 0.0 && t - r0 < poisson.h)
      .collect::<Vec<f64>>();
    assert!(inside.len() > 1, "the convolution samples the edge cell");
    for t in inside {
      let expected = -poisson.mass_by_quadrature(0.0, t - r0).ln();
      pretty_assert_relative_eq!(expected, g.interp(t)?, max_relative = 1e-10);
    }
    assert_eq!((r0, BoundaryBehavior::Hard), (g.x_min(), g.left_extrap()));
    Ok(())
  }

  #[test]
  fn test_point_convolution_shifts_the_edge_cell() -> Result<(), Report> {
    let poisson = helpers::Poisson::one_mutation();
    let point = Distribution::point(2.0, 0.0);
    let Distribution::Function(g) = distribution_convolution(
      &point,
      &Distribution::Function(poisson.function()?),
      MaxGridPoints::default(),
    )?
    else {
      panic!("a point convolved with a function is a function");
    };
    pretty_assert_relative_eq!(
      poisson.cost(0.3 * poisson.h),
      g.interp(2.0 + 0.3 * poisson.h)?,
      max_relative = 1e-10
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::edge_on_the_left( Side::Right, 1.0)]
  #[case::edge_on_the_right(Side::Left, -1.0)]
  #[trace]
  fn test_point_convolution_across_the_edge_keeps_the_shifted_poisson_cost(
    #[case] soft: Side,
    #[case] direction: f64,
  ) -> Result<(), Report> {
    let poisson = helpers::Poisson::one_mutation();
    let f = Distribution::Function(poisson.function()?);
    let f = if direction > 0.0 { f } else { f.negate()? };
    let at = 3.0 * poisson.h;
    let Distribution::Function(g) = convolve_across_edge(
      &Distribution::point(at, 0.0),
      &f,
      soft,
      1e-6,
      100,
      MaxGridPoints::default(),
    )?
    else {
      panic!("a point convolved with a function is a function");
    };

    let offsets = helpers::cell_offsets(poisson.h).mapv(|x| (at + direction * x - at) * direction);
    let reference = poisson.cost(poisson.h) - g.interp(at + direction * poisson.h)?;
    let expected = offsets.mapv(|x| poisson.cost(x));
    let actual = offsets.mapv(|x| g.interp(at + direction * x).expect("the offset lies on the grid") + reference);
    pretty_assert_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_product_keeps_the_edge_cell_law_beside_a_hard_edge() -> Result<(), Report> {
    let poisson = helpers::Poisson::one_mutation();
    let product = distribution_multiplication(
      &Distribution::Function(poisson.function()?),
      &Distribution::Function(helpers::line(&poisson, BoundaryBehavior::Hard)?),
      MaxGridPoints::default(),
    )?;
    let Distribution::Function(product) = product else {
      panic!("a product of functions is a function");
    };

    assert_eq!(BoundaryBehavior::HardEdgeCell(poisson.law()), product.left_extrap());
    let x = 0.3 * poisson.h;
    pretty_assert_relative_eq!(
      poisson.cost(x) + 2.0 + 5.0 * x,
      product.interp(x)?,
      max_relative = 1e-12
    );
    Ok(())
  }

  #[test]
  fn test_product_of_an_edge_cell_and_a_hard_approach_at_one_point_is_an_internal_error() -> Result<(), Report> {
    let poisson = helpers::Poisson::one_mutation();
    let approach = BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 1.0 });
    assert_error!(
      distribution_multiplication(
        &Distribution::Function(poisson.function()?),
        &Distribution::Function(helpers::line(&poisson, approach)?),
        MaxGridPoints::default(),
      ),
      "Cannot multiply an edge cell law at 0 by a hard-approach tail at the same point: the product follows neither law, and no inference path multiplies them. This is an internal error. Please report it to developers."
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::same_offset( 1e-3)]
  #[case::other_offset(3e-3)]
  #[trace]
  fn test_quotient_of_a_product_follows_the_dividend_law_between_grid_points(#[case] offset_share: f64) -> Result<(), Report> {
    let f = helpers::Poisson::one_mutation();
    let g = helpers::Poisson::with_offset(0.001, 1000, offset_share);
    let (f_fn, g_fn) = (Distribution::Function(f.function()?), Distribution::Function(g.function()?));
    let product = distribution_multiplication(&f_fn, &g_fn, MaxGridPoints::default())?;
    let Distribution::Function(quotient) = distribution_division(&product, &g_fn, MaxGridPoints::default())? else {
      panic!("a quotient of functions is a function");
    };

    let offsets = helpers::cell_offsets(f.h);
    let reference = f.cost(f.h) - quotient.interp(f.h)?;
    let expected = offsets.mapv(|x| f.cost(x));
    let actual = offsets.mapv(|x| quotient.interp(x).expect("the offset lies on the grid") + reference);
    pretty_assert_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_product_at_one_edge_with_other_offsets_stays_within_the_offset_bound() -> Result<(), Report> {
    let f = helpers::Poisson::one_mutation();
    let g = helpers::Poisson::with_offset(0.001, 1000, 3e-3);
    let Distribution::Function(product) = distribution_multiplication(
      &Distribution::Function(f.function()?),
      &Distribution::Function(g.function()?),
      MaxGridPoints::default(),
    )?
    else {
      panic!("a product of functions is a function");
    };

    let b = f.law().b;
    let ln_mean_offset = f64::midpoint(f.offset.ln(), g.offset.ln());
    let bound = b * ((f.offset.ln() - ln_mean_offset).powi(2) + (g.offset.ln() - ln_mean_offset).powi(2)) / 8.0;
    let offsets = helpers::cell_offsets(f.h);
    let error = offsets.mapv(|x| product.interp(x).expect("the offset lies on the grid") - f.cost(x) - g.cost(x));
    let BoundaryBehavior::HardEdgeCell(law) = product.left_extrap() else {
      panic!("the product keeps an edge cell law, got {:?}", product.left_extrap());
    };
    pretty_assert_relative_eq!(2.0 * b, law.b, max_relative = 1e-15);
    pretty_assert_relative_eq!(ln_mean_offset.exp(), law.offset, max_relative = 1e-14);
    assert!(
      error.iter().all(|e| e.abs() <= bound),
      "errors {error} exceed the bound {bound}"
    );
    Ok(())
  }

  #[test]
  fn test_product_of_edge_cells_at_different_edges_keeps_the_inner_law() -> Result<(), Report> {
    let f = helpers::Poisson::one_mutation();
    let shift = 2.0 * f.h;
    let f_fn = Distribution::Function(f.function()?);
    let shifted = distribution_convolution(&Distribution::point(shift, 0.0), &f_fn, MaxGridPoints::default())?;
    let Distribution::Function(product) = distribution_multiplication(&f_fn, &shifted, MaxGridPoints::default())?
    else {
      panic!("a product of functions is a function");
    };

    assert_eq!(
      BoundaryBehavior::HardEdgeCell(f.law().shift_arg(shift)),
      product.left_extrap()
    );
    let bound = f.law().b * f.h.powi(2) / (8.0 * (shift + f.offset).powi(2));
    let offsets = helpers::cell_offsets(f.h);
    let error =
      offsets.mapv(|x| product.interp(shift + x).expect("the offset lies on the grid") - f.cost(x) - f.cost(shift + x));
    assert!(
      error.iter().all(|e| e.abs() <= bound),
      "errors {error} exceed the bound {bound}"
    );
    Ok(())
  }

  mod helpers {
    use super::*;
    use treetime_grid::EdgeCellLaw;

    pub(super) struct Poisson {
      k: f64,
      one_mutation: f64,
      pub(super) offset: f64,
      pub(super) h: f64,
      peak: f64,
    }

    impl Poisson {
      pub(super) fn one_mutation() -> Self {
        Self::new(0.001, 1000)
      }

      pub(super) fn new(k: f64, sequence_length: usize) -> Self {
        Self::with_offset(k, sequence_length, 1e-3)
      }

      pub(super) fn with_offset(k: f64, sequence_length: usize, offset_share: f64) -> Self {
        let one_mutation = sequence_length.to_f64_exact().recip();
        let h = (k.max(one_mutation) * one_mutation).sqrt() / 10.0;
        let raw = |x: f64| (x - k * (x + offset_share * one_mutation).ln()) / one_mutation;
        let peak = (0..=CELLS)
          .map(|i| raw(h * i.to_f64_exact()))
          .fold(f64::INFINITY, f64::min);
        Self {
          k,
          one_mutation,
          offset: offset_share * one_mutation,
          h,
          peak,
        }
      }

      pub(super) fn cost(&self, x: f64) -> f64 {
        (x - self.k * (x + self.offset).ln()) / self.one_mutation - self.peak
      }

      pub(super) fn law(&self) -> EdgeCellLaw {
        EdgeCellLaw {
          t_edge: 0.0,
          offset: self.offset,
          b: self.k / self.one_mutation,
        }
      }

      pub(super) fn function(&self) -> Result<DistributionFunction<f64, NegLog>, Report> {
        let x_max = self.h * CELLS.to_f64_exact();
        let y = Array1::linspace(0.0, x_max, CELLS + 1).mapv(|x| self.cost(x));
        DistributionFunction::from_range_values((0.0, x_max), y)?
          .with_left_extrap(BoundaryBehavior::HardEdgeCell(self.law()))?
          .with_right_extrap(BoundaryBehavior::Hard)
      }

      pub(super) fn mass_by_quadrature(&self, a: f64, b: f64) -> f64 {
        if b <= a {
          return 0.0;
        }
        let intervals = 1 << 16;
        let (s_a, s_b) = ((a + self.offset).ln(), (b + self.offset).ln());
        let step = (s_b - s_a) / f64::from(intervals);
        let sum: f64 = (0..=intervals)
          .map(|i| {
            let s = s_a + step * f64::from(i);
            let weight = if i == 0 || i == intervals {
              1.0
            } else if i % 2 == 1 {
              4.0
            } else {
              2.0
            };
            weight * (-self.cost(s.exp() - self.offset) + s).exp()
          })
          .sum();
        sum * step / 3.0
      }

      pub(super) fn second_cell_mass_to(&self, b: f64) -> f64 {
        let (start, end) = (self.h, 2.0 * self.h);
        if b <= start {
          return 0.0;
        }
        let (p_start, p_end) = ((-self.cost(start)).exp(), (-self.cost(end)).exp());
        let line = |x: f64| p_start + (p_end - p_start) * (x - start) / (end - start);
        0.5 * (b - start) * (line(start) + line(b))
      }
    }

    pub(super) fn line(poisson: &Poisson, left: BoundaryBehavior) -> Result<DistributionFunction<f64, NegLog>, Report> {
      let f = poisson.function()?;
      let line = f.t().mapv(|t| 2.0 + 5.0 * t);
      DistributionFunction::from_range_values((f.x_min(), f.x_max()), line)?
        .with_left_extrap(left)?
        .with_right_extrap(BoundaryBehavior::Hard)
    }

    pub(super) fn cell_offsets(h: f64) -> Array1<f64> {
      array![1e-5, 1e-3, 0.1, 0.3, 0.7, 0.99].mapv(|u| u * h)
    }

    pub(super) const CELLS: usize = 40;
  }
}

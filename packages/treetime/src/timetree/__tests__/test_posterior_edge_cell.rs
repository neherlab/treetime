#[cfg(test)]
mod tests {
  use crate::coalescent::mean_field::NodeSurvival;
  use crate::timetree::max_posterior_region::{
    MaxPosteriorRegion, PeakLocation, max_posterior_region, posterior_survival,
  };
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_distribution::{BoundaryBehavior, Distribution, DistributionFunction, MassProfile, NegLog, TailMass};
  use treetime_grid::EdgeCellLaw;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_relative_eq};

  const ONE_MUTATION: f64 = 1e-3;
  const OFFSET: f64 = 1e-3 * ONE_MUTATION;
  const MUTATIONS: f64 = 0.3;
  const FRACTIONAL_MUTATIONS: f64 = 0.01;
  const CELL: f64 = 0.1 * ONE_MUTATION;
  const CELLS: usize = 40;

  #[test]
  fn test_max_posterior_region_crosses_the_edge_cell_at_the_level_of_its_law() -> Result<(), Report> {
    let distribution = helpers::poisson()?;
    let Some(MaxPosteriorRegion::Grid(region)) = max_posterior_region(&distribution, 0.7)? else {
      panic!("a function has a grid region");
    };
    assert!(
      region.lower > 0.0 && region.lower < CELL,
      "the region {region:?} starts inside the edge cell"
    );
    let level_lower = distribution.eval(region.lower)?;
    let level_upper = helpers::density_linear_level(region.upper);
    pretty_assert_relative_eq!(level_upper, level_lower, max_relative = 1e-12);
    pretty_assert_relative_eq!(helpers::cost(region.lower), level_lower, max_relative = 1e-12);
    Ok(())
  }

  #[test]
  fn test_posterior_survival_follows_the_mass_model_inside_the_edge_cell() -> Result<(), Report> {
    let distribution = helpers::poisson()?;
    let survival = posterior_survival(&distribution, 1e-9, 1e-6)?.expect("a function has a survival function");
    let Distribution::Function(function) = distribution.normalize()? else {
      panic!("the posterior is a function");
    };
    let profile = MassProfile::new(&function)?;
    let fractions = profile
      .fractions(TailMass::Required)
      .expect("the posterior has a finite mass");

    let interpolant = NodeSurvival::from(survival.clone());
    let knots_in_edge_cell = survival.t.iter().filter(|&&t| t > 0.0 && t < CELL).count();
    assert!(knots_in_edge_cell > 0, "the survival function refines the edge cell");
    let x: Array1<f64> = survival
      .t
      .iter()
      .tuple_windows()
      .filter(|&(&t0, _)| t0 < CELL)
      .flat_map(|(&t0, &t1)| (0..16).map(move |k| t0 + (t1 - t0) * (k.to_f64_exact() / 16.0)))
      .collect();
    let expected = x.mapv(|x| fractions.survival(x));
    let actual = x.mapv(|x| interpolant.at(x));
    pretty_assert_abs_diff_eq!(expected, actual, epsilon = 1e-6);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::edge_on_the_left(  1.0)]
  #[case::edge_on_the_right(-1.0)]
  #[trace]
  fn test_max_posterior_region_around_a_mode_inside_the_edge_cell(#[case] direction: f64) -> Result<(), Report> {
    let distribution = helpers::poisson_with(FRACTIONAL_MUTATIONS)?;
    let distribution = if direction > 0.0 { distribution } else { distribution.negate()? };
    let Some(MaxPosteriorRegion::Grid(region)) = max_posterior_region(&distribution, 0.01)? else {
      panic!("a function has a grid region");
    };

    let mode = FRACTIONAL_MUTATIONS * ONE_MUTATION - OFFSET;
    let (near, far) = if direction > 0.0 {
      (region.lower, region.upper)
    } else {
      (-region.upper, -region.lower)
    };
    assert_eq!(PeakLocation::Interior, region.peak);
    assert!(near > 0.0 && near < mode && far > mode && far < CELL, "the region {region:?} surrounds the mode in the edge cell");
    pretty_assert_relative_eq!(helpers::cost_with(FRACTIONAL_MUTATIONS, near), helpers::cost_with(FRACTIONAL_MUTATIONS, far), max_relative = 1e-12);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn cost(x: f64) -> f64 {
      cost_with(MUTATIONS, x)
    }

    pub(super) fn cost_with(mutations: f64, x: f64) -> f64 {
      (x - mutations * ONE_MUTATION * (x + OFFSET).ln()) / ONE_MUTATION
    }

    pub(super) fn poisson() -> Result<Distribution<NegLog>, Report> {
      poisson_with(MUTATIONS)
    }

    pub(super) fn poisson_with(mutations: f64) -> Result<Distribution<NegLog>, Report> {
      let x_max = CELL * CELLS.to_f64_exact();
      let y = Array1::linspace(0.0, x_max, CELLS + 1).mapv(|x| cost_with(mutations, x));
      let law = EdgeCellLaw {
        t_edge: 0.0,
        offset: OFFSET,
        b: mutations,
      };
      let function = DistributionFunction::from_range_values((0.0, x_max), y)?
        .with_left_extrap(BoundaryBehavior::HardEdgeCell(law))?
        .with_right_extrap(BoundaryBehavior::Hard)?;
      Ok(Distribution::Function(function))
    }

    pub(super) fn density_linear_level(t: f64) -> f64 {
      let cell = (t / CELL).floor();
      let (x0, x1) = (cell * CELL, (cell + 1.0) * CELL);
      let (p0, p1) = ((-cost(x0)).exp(), (-cost(x1)).exp());
      let share = (t - x0) / (x1 - x0);
      -(p0 + (p1 - p0) * share).ln()
    }
  }
}

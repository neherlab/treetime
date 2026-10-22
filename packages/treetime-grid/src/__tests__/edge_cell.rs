#[cfg(test)]
mod tests {
  use crate::*;
  use eyre::Report;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{assert_error, pretty_assert_relative_eq, pretty_assert_ulps_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::zero_length(            (0.0,    1000))]
  #[case::tenth_mutation(         (0.0001, 1000))]
  #[case::one_mutation(           (0.001,  1000))]
  #[case::three_mutations(        (0.003,  1000))]
  #[case::ten_mutations(          (0.01,   1000))]
  #[case::fifty_mutations(        (0.05,   1000))]
  #[case::long_genome(            (0.001, 30000))]
  #[case::below_gaussian_threshold((0.0499, 10000))]
  #[trace]
  fn test_edge_cell_follows_the_poisson_cost_inside_the_first_cell(#[case] (k, sequence_length): (f64, usize)) -> Result<(), Report> {
    let poisson = helpers::Poisson::new(k, sequence_length);
    let cell = poisson.edge_cell()?;
    let h = poisson.cell_width();

    pretty_assert_ulps_eq!(poisson.cost(0.0), cell.eval(0.0), max_ulps = 0);
    pretty_assert_ulps_eq!(poisson.cost(h), cell.eval(h), max_ulps = 0);
    let inside = helpers::points_inside(poisson.offset, h);
    let expected = inside.mapv(|x| poisson.cost(x));
    let actual = inside.mapv(|x| cell.eval(x));
    pretty_assert_relative_eq!(expected, actual, epsilon = 1e-9, max_relative = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero_length(            (0.0,    1000))]
  #[case::tenth_mutation(         (0.0001, 1000))]
  #[case::one_mutation(           (0.001,  1000))]
  #[case::three_mutations(        (0.003,  1000))]
  #[case::ten_mutations(          (0.01,   1000))]
  #[case::fifty_mutations(        (0.05,   1000))]
  #[case::long_genome(            (0.001, 30000))]
  #[case::below_gaussian_threshold((0.0499, 10000))]
  #[trace]
  fn test_edge_cell_mass_matches_quadrature_of_the_poisson_cost(#[case] (k, sequence_length): (f64, usize)) -> Result<(), Report> {
    let poisson = helpers::Poisson::new(k, sequence_length);
    let cell = poisson.edge_cell()?;
    let (d_edge, d_anchor) = (cell.edge_distance(), cell.anchor_distance());
    let bounds = [d_edge, 2.0 * d_edge, 0.01 * d_anchor, 0.5 * d_anchor, d_anchor];
    let intervals = [(0, 4), (0, 1), (1, 3), (2, 3), (3, 4)].map(|(near, far)| (bounds[near], bounds[far]));

    let expected = intervals.map(|(near, far)| poisson.ln_mass_by_quadrature(near, far));
    let actual = intervals.map(|(near, far)| cell.ln_relative_mass(near, far));
    pretty_assert_relative_eq!(Array1::from_iter(expected), Array1::from_iter(actual), epsilon = 1e-10, max_relative = 1e-10);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::edge_far_below(  ( 2.0, -40.0))]
  #[case::edge_far_above(  ( 2.0,  40.0))]
  #[case::flat_power(      ( 0.0,  -3.0))]
  #[case::large_exponent(  (60.0,   5.0))]
  #[trace]
  fn test_edge_cell_mass_matches_quadrature_for_any_slope(#[case] (b, edge_y): (f64, f64)) -> Result<(), Report> {
    let law = EdgeCellLaw { t_edge: 2.0, offset: 0.01, b };
    let cell = EdgeCell::new(law, GridEdge { t: 2.0, y: edge_y }, GridEdge { t: 2.5, y: 0.0 })?;
    let (d_edge, d_anchor) = (cell.edge_distance(), cell.anchor_distance());

    let expected = helpers::ln_quadrature(|d| cell.ln_relative_density(d), d_edge, d_anchor);
    pretty_assert_relative_eq!(expected, cell.ln_relative_mass(d_edge, d_anchor), epsilon = 1e-10, max_relative = 1e-10);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_mutation(  (0.001,  1000))]
  #[case::fifty_mutations((0.05,  1000))]
  #[case::long_genome(   (0.001, 30000))]
  #[trace]
  fn test_edge_cell_narrow_masses_match_quadrature(#[case] (k, sequence_length): (f64, usize)) -> Result<(), Report> {
    let poisson = helpers::Poisson::new(k, sequence_length);
    let cell = poisson.edge_cell()?;
    let (near, far) = (0.3 * cell.anchor_distance(), 0.30003 * cell.anchor_distance());
    let parts = Array1::linspace(near, far, 11);

    let sum: f64 = parts
      .windows(2)
      .into_iter()
      .map(|pair| cell.ln_relative_mass(pair[0], pair[1]).exp())
      .sum();
    let expected = poisson.ln_mass_by_quadrature(near, far);
    pretty_assert_relative_eq!(expected, sum.ln(), max_relative = 1e-12);
    pretty_assert_relative_eq!(expected, cell.ln_relative_mass(near, far), max_relative = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::from_edge_one_mutation(    (0.001,  1000), CellEnd::Edge,   0.3)]
  #[case::from_anchor_one_mutation(  (0.001,  1000), CellEnd::Anchor, 0.3)]
  #[case::from_edge_tenth_mutation(  (0.0001, 1000), CellEnd::Edge,   1e-3)]
  #[case::from_anchor_ten_mutations( (0.01,   1000), CellEnd::Anchor, 0.999)]
  #[trace]
  fn test_edge_cell_distance_at_mass_inverts_the_mass(
    #[case] (k, sequence_length): (f64, usize),
    #[case] from: CellEnd,
    #[case] share: f64,
  ) -> Result<(), Report> {
    let cell = helpers::Poisson::new(k, sequence_length).edge_cell()?;
    let (d_edge, d_anchor) = (cell.edge_distance(), cell.anchor_distance());
    let ln_mass = cell.ln_relative_mass(d_edge, d_anchor) + share.ln();

    let distance = cell.distance_at_mass(from, ln_mass);
    let actual = match from {
      CellEnd::Edge => cell.ln_relative_mass(d_edge, distance),
      CellEnd::Anchor => cell.ln_relative_mass(distance, d_anchor),
    };
    pretty_assert_relative_eq!(ln_mass, actual, epsilon = 1e-12, max_relative = 1e-12);
    Ok(())
  }

  #[test]
  fn test_edge_cell_distance_at_mass_beyond_the_cell_is_its_far_end() -> Result<(), Report> {
    let cell = helpers::Poisson::new(0.001, 1000).edge_cell()?;
    let total = cell.ln_relative_mass(cell.edge_distance(), cell.anchor_distance());
    assert_eq!(
      (cell.anchor_distance(), cell.edge_distance(), cell.edge_distance()),
      (
        cell.distance_at_mass(CellEnd::Edge, total + 1.0),
        cell.distance_at_mass(CellEnd::Anchor, total + 1.0),
        cell.distance_at_mass(CellEnd::Edge, f64::NEG_INFINITY),
      )
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::near_anchor(0.1)]
  #[case::middle(     0.5)]
  #[case::near_edge(  0.99)]
  #[trace]
  fn test_edge_cell_position_between_inverts_the_value(#[case] rise: f64) -> Result<(), Report> {
    let cell = helpers::Poisson::new(0.001, 1000).edge_cell()?;
    let level = cell.anchor.y + rise * (cell.edge.y - cell.anchor.y);

    let t = cell.position_between(level, cell.anchor.t, cell.edge.t);
    assert!(t > cell.edge.t && t < cell.anchor.t, "position {t} must lie inside the cell");
    pretty_assert_relative_eq!(level, cell.eval(t), max_relative = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::steep_fall_whole_cell(  ( 64.0,       (1.0, 2.0)), 64.0 + (-(-64.0_f64).exp()).ln_1p() - 64.0_f64.ln())]
  #[case::steep_fall_near_half(   ( 64.0,       (1.0, 1.5)), 64.0 + (-(-32.0_f64).exp()).ln_1p() - 64.0_f64.ln())]
  #[case::steep_fall_far_half(    ( 64.0,       (1.5, 2.0)), 32.0 + (-(-32.0_f64).exp()).ln_1p() - 64.0_f64.ln())]
  #[case::steep_rise_whole_cell(  (-64.0,       (1.0, 2.0)), (-(-64.0_f64).exp()).ln_1p() - 64.0_f64.ln())]
  #[case::extreme_rise_whole_cell((-2_000_000.0, (1.0, 2.0)), -2_000_000.0_f64.ln())]
  #[case::extreme_rise_far_end(   (-2_000_000.0, (2.0 - 1.0 / 1_048_576.0, 2.0)), (-(-2_000_000.0_f64 / 1_048_576.0).exp_m1()).ln() - 2_000_000.0_f64.ln())]
  #[case::extreme_fall_whole_cell(( 2_000_000.0, (1.0, 2.0)), 2_000_000.0 - 2_000_000.0_f64.ln())]
  #[trace]
  fn test_edge_cell_mass_of_a_steep_exponential_matches_its_closed_form(
    #[case] (anchor_y, (near, far)): (f64, (f64, f64)),
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let law = EdgeCellLaw { t_edge: 0.0, offset: 1.0, b: 0.0 };
    let cell = EdgeCell::new(law, GridEdge { t: 0.0, y: 0.0 }, GridEdge { t: 1.0, y: anchor_y })?;
    pretty_assert_relative_eq!(expected, cell.ln_relative_mass(near, far), epsilon = 1e-12, max_relative = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tenth_of_a_cell((1e-5,  1000))]
  #[case::half_a_cell(    (5e-5,  1000))]
  #[case::long_genome(    (1e-6, 30000))]
  #[trace]
  fn test_edge_cell_mode_lies_at_the_minimum_of_the_poisson_cost(#[case] (k, sequence_length): (f64, usize)) -> Result<(), Report> {
    let poisson = helpers::Poisson::new(k, sequence_length);
    let cell = poisson.edge_cell()?;
    let mode = cell.mode().expect("the cost has its minimum inside the edge cell");

    let expected = k - poisson.offset;
    pretty_assert_relative_eq!(expected, mode.t, max_relative = 1e-12);
    pretty_assert_relative_eq!(poisson.cost(expected), mode.y, max_relative = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::toward_the_edge(  true )]
  #[case::toward_the_anchor(false)]
  #[trace]
  fn test_edge_cell_position_between_finds_each_branch_around_the_mode(#[case] toward_edge: bool) -> Result<(), Report> {
    let poisson = helpers::Poisson::new(1e-5, 1000);
    let cell = poisson.edge_cell()?;
    let mode = cell.mode().expect("the cost has its minimum inside the edge cell");
    let outer = if toward_edge { cell.edge } else { cell.anchor };
    let level = f64::midpoint(mode.y, cell.edge.y.min(cell.anchor.y));

    let t = cell.position_between(level, mode.t, outer.t);
    assert!((t - mode.t) * (outer.t - mode.t) > 0.0, "position {t} lies between the mode and {}", outer.t);
    pretty_assert_relative_eq!(level, poisson.cost(t), max_relative = 1e-12);
    Ok(())
  }

  #[test]
  fn test_edge_cell_without_interior_minimum_has_no_mode() -> Result<(), Report> {
    let cell = helpers::Poisson::new(0.001, 1000).edge_cell()?;
    assert_eq!(None, cell.mode());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::same_offset( (1e-6, 1.5), 1e-6)]
  #[case::other_offset((4e-6, 1.5), 2e-6)]
  #[case::empty_factor((5e-6, 0.0), 1e-6)]
  #[trace]
  fn test_edge_cell_law_product_at_one_edge_takes_the_weighted_geometric_mean_offset(
    #[case] (other_offset, other_b): (f64, f64),
    #[case] expected_offset: f64,
  ) {
    let law = EdgeCellLaw { t_edge: 2.0, offset: 1e-6, b: 1.5 };
    let other = EdgeCellLaw { t_edge: 2.0, offset: other_offset, b: other_b };

    let product = law.compose_multiply(&other).expect("both laws share the edge");
    pretty_assert_relative_eq!(1.5 + other_b, product.b, max_relative = 1e-15);
    pretty_assert_relative_eq!(expected_offset, product.offset, max_relative = 1e-14);
  }

  #[test]
  fn test_edge_cell_law_quotient_inverts_the_product() -> Result<(), Report> {
    let law = EdgeCellLaw {
      t_edge: 0.0,
      offset: 1e-6,
      b: 0.7,
    };
    let other = EdgeCellLaw {
      t_edge: 0.0,
      offset: 3e-6,
      b: 2.0,
    };
    let product = law.compose_multiply(&other).expect("both laws share the edge");

    let quotient = product.compose_divide(&other)?.expect("both laws share the edge");
    pretty_assert_relative_eq!(law.b, quotient.b, max_relative = 1e-14);
    pretty_assert_relative_eq!(law.offset, quotient.offset, max_relative = 1e-12);
    assert_eq!(None, law.compose_multiply(&EdgeCellLaw { t_edge: 1.0, ..other }));
    assert_error!(
      law.compose_divide(&other),
      "Cannot divide an edge cell law of exponent 0.7 by one of exponent 2: the quotient has a negative exponent. This is an internal error. Please report it to developers."
    );
    Ok(())
  }

  #[test]
  fn test_edge_cell_without_exponent_is_the_line_between_its_ends() -> Result<(), Report> {
    let law = EdgeCellLaw {
      t_edge: 0.0,
      offset: 1e-6,
      b: 0.0,
    };
    let cell = EdgeCell::new(law, GridEdge { t: 0.0, y: 3.0 }, GridEdge { t: 0.5, y: 1.0 })?;
    pretty_assert_relative_eq!(2.0, cell.eval(0.25), max_relative = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::origin_inside_the_cell(EdgeCellLaw { t_edge: 0.0, offset: -0.1, b: 1.0 }, "The law of an edge cell needs a finite positive offset of its origin beyond the cell edge 0 and a cell of positive width, got the offset -0.1 and the other cell end 1")]
  #[case::origin_at_the_edge(    EdgeCellLaw { t_edge: 0.0, offset:  0.0, b: 1.0 }, "The law of an edge cell needs a finite positive offset of its origin beyond the cell edge 0 and a cell of positive width, got the offset 0 and the other cell end 1")]
  #[case::negative_exponent(     EdgeCellLaw { t_edge: 0.0, offset:  0.1, b: -1.0 }, "The law of an edge cell needs a finite non-negative exponent, got -1")]
  #[trace]
  fn test_edge_cell_rejects_an_invalid_law(#[case] law: EdgeCellLaw, #[case] message: &str) {
    assert_error!(
      EdgeCell::new(law, GridEdge { t: 0.0, y: 1.0 }, GridEdge { t: 1.0, y: 0.0 }),
      message
    );
  }

  #[test]
  fn test_gridfn_follows_the_edge_cell_law_in_the_edge_cells_only() -> Result<(), Report> {
    let poisson = helpers::Poisson::new(0.001, 1000);
    let h = poisson.cell_width();
    let x = Array1::linspace(0.0, 4.0 * h, 5);
    let y = x.mapv(|x| poisson.cost(x));
    let grid_fn = GridFn::from_range_values((0.0, 4.0 * h), y.clone())?;
    let left = BoundaryBehavior::HardEdgeCell(poisson.law());

    let in_edge_cell = grid_fn.interp_with_extrap(0.3 * h, left, BoundaryBehavior::Hard)?;
    let in_second_cell = grid_fn.interp_with_extrap(1.5 * h, left, BoundaryBehavior::Hard)?;
    pretty_assert_relative_eq!(poisson.cost(0.3 * h), in_edge_cell, max_relative = 1e-12);
    pretty_assert_relative_eq!(f64::midpoint(y[1], y[2]), in_second_cell, max_relative = 1e-12);
    Ok(())
  }

  #[test]
  fn test_gridfn_follows_the_mirrored_edge_cell_law_at_the_right_end() -> Result<(), Report> {
    let poisson = helpers::Poisson::new(0.001, 1000);
    let h = poisson.cell_width();
    let x = Array1::linspace(-4.0 * h, 0.0, 5);
    let y = x.mapv(|x| poisson.cost(-x));
    let grid_fn = GridFn::from_range_values((-4.0 * h, 0.0), y)?;
    let law = poisson.law().scale_arg(-1.0);

    let actual = grid_fn.interp_with_extrap(-0.3 * h, BoundaryBehavior::Hard, BoundaryBehavior::HardEdgeCell(law))?;
    pretty_assert_relative_eq!(poisson.cost(0.3 * h), actual, max_relative = 1e-12);
    Ok(())
  }

  #[test]
  fn test_gridfn_edge_cell_has_hard_support_beyond_its_edge() -> Result<(), Report> {
    let poisson = helpers::Poisson::new(0.001, 1000);
    let h = poisson.cell_width();
    let grid_fn = GridFn::from_range_values((0.0, h), array![poisson.cost(0.0), poisson.cost(h)])?;
    let left = BoundaryBehavior::HardEdgeCell(poisson.law());

    let actual = grid_fn.interp_with_extrap(-0.5 * poisson.offset, left, BoundaryBehavior::Hard)?;
    pretty_assert_ulps_eq!(0.0, actual, max_ulps = 0);
    assert!(!left.is_soft());
    Ok(())
  }

  #[test]
  fn test_gridfn_edge_cell_law_away_from_the_grid_edge_leaves_the_cell_linear() -> Result<(), Report> {
    let poisson = helpers::Poisson::new(0.001, 1000);
    let h = poisson.cell_width();
    let y = array![poisson.cost(h), poisson.cost(2.0 * h)];
    let grid_fn = GridFn::from_range_values((h, 2.0 * h), y.clone())?;
    let left = BoundaryBehavior::HardEdgeCell(poisson.law());

    assert_eq!(None, grid_fn.edge_cell(0, left, BoundaryBehavior::Hard)?);
    let actual = grid_fn.interp_with_extrap(1.5 * h, left, BoundaryBehavior::Hard)?;
    pretty_assert_relative_eq!(f64::midpoint(y[0], y[1]), actual, max_relative = 1e-15);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) struct Poisson {
      k: f64,
      one_mutation: f64,
      pub(super) offset: f64,
    }

    impl Poisson {
      pub(super) fn new(k: f64, sequence_length: usize) -> Self {
        let one_mutation = sequence_length.to_f64_exact().recip();
        Self {
          k,
          one_mutation,
          offset: 1e-3 * one_mutation,
        }
      }

      pub(super) fn cost(&self, x: f64) -> f64 {
        (x - self.k * (x + self.offset).ln()) / self.one_mutation
      }

      pub(super) fn cell_width(&self) -> f64 {
        (self.k.max(self.one_mutation) * self.one_mutation).sqrt() / 10.0
      }

      pub(super) fn law(&self) -> EdgeCellLaw {
        EdgeCellLaw {
          t_edge: 0.0,
          offset: self.offset,
          b: self.k / self.one_mutation,
        }
      }

      pub(super) fn edge_cell(&self) -> Result<EdgeCell, Report> {
        let h = self.cell_width();
        EdgeCell::new(
          self.law(),
          GridEdge {
            t: 0.0,
            y: self.cost(0.0),
          },
          GridEdge { t: h, y: self.cost(h) },
        )
      }

      pub(super) fn ln_mass_by_quadrature(&self, near: f64, far: f64) -> f64 {
        let y_anchor = self.cost(self.cell_width());
        ln_quadrature(|d| y_anchor - self.cost(d - self.offset), near, far)
      }
    }

    pub(super) fn points_inside(offset: f64, h: f64) -> Array1<f64> {
      array![
        1e-3 * offset,
        0.5 * offset,
        offset,
        10.0 * offset,
        0.1 * h,
        0.5 * h,
        0.9 * h
      ]
    }

    pub(super) fn ln_quadrature(ln_density: impl Fn(f64) -> f64, near: f64, far: f64) -> f64 {
      let intervals = 1 << 20;
      let (s_near, s_far) = (near.ln(), far.ln());
      let step = (s_far - s_near) / intervals.to_f64_exact();
      let ln_integrand = |i: u32| {
        let s = s_near + step * i.to_f64_exact();
        ln_density(s.exp()) + s
      };
      let reference = (0..=intervals).map(ln_integrand).fold(f64::NEG_INFINITY, f64::max);
      let sum: f64 = (0..=intervals)
        .map(|i| {
          let weight = if i == 0 || i == intervals {
            1.0
          } else if i % 2 == 1 {
            4.0
          } else {
            2.0
          };
          weight * (ln_integrand(i) - reference).exp()
        })
        .sum();
      reference + (sum * step / 3.0).ln()
    }
  }
}

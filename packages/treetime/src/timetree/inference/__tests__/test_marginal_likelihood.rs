#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::coalescent::mean_field::CoalescentPrior;
  use crate::test_utils::{find_node_key_by_name, parent_edge_key, point_coalescent_prior, unknown_branches};
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::result::{BranchLikelihood, log_marginal_likelihood};
  use crate::timetree::optimization::marginal_clock_rate::branches_at_rate;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use itertools::{Itertools, izip};
  use ndarray::Array1;
  use ordered_float::OrderedFloat;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::{BoundaryBehavior, Distribution, DistributionFunction, NegLog};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::MaxGridPoints;
  use treetime_io::nwk::nwk_read;

  #[rstest]
  #[case::no_prior_reference_rate(1e-3, None)]
  #[case::no_prior_slower(0.5e-3, None)]
  #[case::coalescent_reference_rate(1e-3, Some(5.0))]
  #[case::coalescent_slower(0.5e-3, Some(5.0))]
  #[case::coalescent_faster_long_tc(2e-3, Some(20.0))]
  #[trace]
  fn test_marginal_likelihood_matches_piecewise_integration(
    #[case] rate: f64,
    #[case] tc: Option<f64>,
  ) -> Result<(), Report> {
    let setup = helpers::setup()?;
    let branches = &branches_at_rate(&setup.branches, helpers::REFERENCE_RATE, rate)?;
    let model = tc.map(|tc| helpers::coalescent(&setup, tc)).transpose()?;
    let backward = propagate_distributions_backward(
      &setup.graph,
      &setup.constraints,
      model.as_ref(),
      &setup.bad_branches,
      branches,
      MaxGridPoints::default(),
    )?;
    let actual = log_marginal_likelihood(&backward, setup.graph.root_key()?, None)?.expect("root has a distribution");
    let expected = helpers::piecewise_log_z(&setup, branches, model.as_ref())?;
    assert_abs_diff_eq!(expected, actual, epsilon = 1e-6);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) const REFERENCE_RATE: f64 = 1e-3;
    const SITES: f64 = 1000.0;
    const TIP_DATE_A: f64 = 2015.0;
    const TIP_DATE_B: f64 = 2014.0;
    const TIP_DATE_C: f64 = 2013.0;
    const TIPS: [(&str, f64); 3] = [("A", TIP_DATE_A), ("B", TIP_DATE_B), ("C", TIP_DATE_C)];
    const PRIOR_TIMES: [(&str, f64); 5] = [
      ("A", TIP_DATE_A),
      ("B", TIP_DATE_B),
      ("C", TIP_DATE_C),
      ("I", 2011.0),
      ("root", 2005.0),
    ];
    const COUNTS: [(&str, f64); 4] = [("A", 3.0), ("B", 2.0), ("I", 2.0), ("C", 5.0)];

    pub(super) struct Setup {
      pub graph: Graph,
      pub constraints: DateConstraints,
      pub bad_branches: BTreeMap<GraphNodeKey, bool>,
      pub branches: BTreeMap<GraphEdgeKey, BranchLikelihood>,
      pub keys: BTreeMap<&'static str, GraphNodeKey>,
      pub names: BTreeMap<GraphNodeKey, Option<String>>,
    }

    pub(super) fn setup() -> Result<Setup, Report> {
      let parsed = nwk_read(b"((A:3,B:2)I:2,C:5)root;".as_slice())?;
      let names = parsed.names();
      let graph = parsed.graph;
      let keys: BTreeMap<&'static str, GraphNodeKey> = ["A", "B", "C", "I", "root"]
        .into_iter()
        .map(|name| {
          (
            name,
            find_node_key_by_name(&graph, &names, name).expect("node not found"),
          )
        })
        .collect();
      let mut constraints = DateConstraints::default();
      for (name, date) in TIPS {
        constraints
          .by_node
          .insert(keys[name], Some(Arc::new(Distribution::<NegLog>::point(date, 0.0))));
      }
      let mut branches = unknown_branches(&graph);
      for (name, k) in COUNTS {
        branches.insert(parent_edge_key(&graph, keys[name]), poisson_branch(k)?);
      }
      let bad_branches = graph.get_nodes().map(|node| (node.key(), false)).collect();
      Ok(Setup {
        graph,
        constraints,
        bad_branches,
        branches,
        keys,
        names,
      })
    }

    fn poisson_branch(k: f64) -> Result<BranchLikelihood, Report> {
      let mean_time = k / (SITES * REFERENCE_RATE);
      let t = Array1::linspace(1e-4, mean_time * 6.0 + 4.0, 6001);
      let y = t.mapv(|t| {
        let m = t * SITES * REFERENCE_RATE;
        m - k * m.ln() - (k - k * k.ln())
      });
      let f = DistributionFunction::from_start_dx_values(t[0], t[1] - t[0], y)?
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Hard)?;
      Ok(BranchLikelihood {
        distribution: Some(Arc::new(Distribution::Function(f))),
        time_length: Some(mean_time),
      })
    }

    pub(super) fn coalescent(setup: &Setup, tc: f64) -> Result<CoalescentPrior, Report> {
      point_coalescent_prior(&setup.graph, &setup.names, &PRIOR_TIMES, tc)
    }

    pub(super) fn piecewise_log_z(
      setup: &Setup,
      branches: &BTreeMap<GraphEdgeKey, BranchLikelihood>,
      model: Option<&CoalescentPrior>,
    ) -> Result<f64, Report> {
      let branch = |name: &str| -> Result<NegLogPieces, Report> {
        let distribution = branches[&parent_edge_key(&setup.graph, setup.keys[name])]
          .distribution
          .as_deref()
          .expect("branch distribution");
        let x = distribution.t().to_vec();
        let y = x.iter().map(|&t| distribution.eval(t)).collect::<Result<_, Report>>()?;
        Ok(NegLogPieces { x, y })
      };
      let contribution = |name: &str, n_children: usize, t: f64| -> Result<f64, Report> {
        model.map_or(Ok(0.0), |m| m.node_contribution(setup.keys[name], n_children, t))
      };
      let kinks: Vec<f64> = PRIOR_TIMES.iter().map(|&(_, t)| t).collect();
      let leaf_cost = TIPS
        .iter()
        .map(|&(name, date)| contribution(name, 0, date))
        .sum::<Result<f64, Report>>()?;

      let (branch_a, branch_b) = (branch("A")?, branch("B")?);
      let (branch_i, branch_c) = (branch("I")?, branch("C")?);
      let root_side = branch_c
        .reflected(TIP_DATE_C)
        .with_kinks(&kinks)
        .plus(|t| contribution("root", 2, t))?;
      let tip_side_b = branch_b.reflected(TIP_DATE_B);
      let tip_side = branch_a
        .reflected(TIP_DATE_A)
        .intersect_kinks(&tip_side_b)
        .with_kinks(&kinks)
        .map_values(|time, neglog| Ok(neglog + tip_side_b.eval(time) + contribution("I", 2, time)? + leaf_cost))?;

      let i_reflected = branch_i.reflected(0.0);
      let shifted = ShiftedProduct::new(&root_side, &i_reflected);
      let inner = |t_i: f64| shifted.exp_integral(t_i);
      let reference = tip_side.y.iter().copied().fold(f64::INFINITY, f64::min);
      let mut total = 0.0;
      for (k, (&x0, &x1)) in tip_side.x.iter().tuple_windows().enumerate() {
        let (y0, y1) = (tip_side.y[k], tip_side.y[k + 1]);
        if y0.min(y1) > reference + NEGLIGIBLE_NEGLOG {
          continue;
        }
        for (node, weight) in GAUSS_LEGENDRE {
          let time = x0 + 0.5 * (x1 - x0) * (1.0 + node);
          let neglog = y0 + (time - x0) / (x1 - x0) * (y1 - y0);
          total += 0.5 * (x1 - x0) * weight * (reference - neglog).exp() * inner(time);
        }
      }
      Ok(total.ln() - reference)
    }

    const NEGLIGIBLE_NEGLOG: f64 = 60.0;

    const GAUSS_LEGENDRE: [(f64, f64); 3] = [
      (-0.774_596_669_241_483_4, 5.0 / 9.0),
      (0.0, 8.0 / 9.0),
      (0.774_596_669_241_483_4, 5.0 / 9.0),
    ];

    struct NegLogPieces {
      x: Vec<f64>,
      y: Vec<f64>,
    }

    impl NegLogPieces {
      fn reflected(&self, origin: f64) -> Self {
        Self {
          x: self.x.iter().rev().map(|&t| origin - t).collect(),
          y: self.y.iter().rev().copied().collect(),
        }
      }

      fn lines(&self) -> Vec<(f64, f64)> {
        izip!(self.x.iter().tuple_windows(), self.y.iter().tuple_windows())
          .map(|((x0, x1), (y0, y1))| {
            let slope = (y1 - y0) / (x1 - x0);
            (y0 - slope * x0, slope)
          })
          .collect()
      }

      fn eval(&self, t: f64) -> f64 {
        let n = self.x.len();
        if t < self.x[0] || t > self.x[n - 1] {
          return f64::INFINITY;
        }
        let j = self.x.partition_point(|&x| x <= t).clamp(1, n - 1) - 1;
        let (x0, x1) = (self.x[j], self.x[j + 1]);
        self.y[j] + (t - x0) / (x1 - x0) * (self.y[j + 1] - self.y[j])
      }

      fn with_kinks(&self, kinks: &[f64]) -> Self {
        let (lo, hi) = (self.x[0], self.x[self.x.len() - 1]);
        let x: Vec<f64> = self
          .x
          .iter()
          .copied()
          .chain(kinks.iter().copied().filter(|&t| t > lo && t < hi))
          .sorted_by_key(|&t| OrderedFloat(t))
          .dedup()
          .collect();
        let y = x.iter().map(|&t| self.eval(t)).collect();
        Self { x, y }
      }

      fn intersect_kinks(&self, other: &Self) -> Self {
        let lo = self.x[0].max(other.x[0]);
        let hi = self.x[self.x.len() - 1].min(other.x[other.x.len() - 1]);
        let x: Vec<f64> = self
          .x
          .iter()
          .chain(&other.x)
          .copied()
          .chain([lo, hi])
          .filter(|&t| t >= lo && t <= hi)
          .sorted_by_key(|&t| OrderedFloat(t))
          .dedup()
          .collect();
        let y = x.iter().map(|&t| self.eval(t)).collect();
        Self { x, y }
      }

      fn plus(&self, f: impl Fn(f64) -> Result<f64, Report>) -> Result<Self, Report> {
        self.map_values(|t, y| Ok(y + f(t)?))
      }

      fn map_values(&self, f: impl Fn(f64, f64) -> Result<f64, Report>) -> Result<Self, Report> {
        let y = izip!(&self.x, &self.y)
          .map(|(&t, &y)| f(t, y))
          .collect::<Result<_, Report>>()?;
        Ok(Self { x: self.x.clone(), y })
      }
    }

    struct ShiftedProduct<'a> {
      fixed: &'a [f64],
      moving: &'a [f64],
      fixed_lines: Vec<(f64, f64)>,
      moving_lines: Vec<(f64, f64)>,
    }

    impl<'a> ShiftedProduct<'a> {
      fn new(fixed: &'a NegLogPieces, moving: &'a NegLogPieces) -> Self {
        Self {
          fixed: &fixed.x,
          moving: &moving.x,
          fixed_lines: fixed.lines(),
          moving_lines: moving.lines(),
        }
      }

      fn exp_integral(&self, shift: f64) -> f64 {
        let (fixed, moving) = (self.fixed, self.moving);
        let (last_fixed, last_moving) = (fixed.len() - 1, moving.len() - 1);
        let lo = fixed[0].max(moving[0] + shift);
        let hi = fixed[last_fixed].min(moving[last_moving] + shift);
        if lo >= hi {
          return 0.0;
        }
        let mut piece_fixed = fixed.partition_point(|&x| x <= lo).clamp(1, last_fixed) - 1;
        let mut piece_moving = moving.partition_point(|&x| x + shift <= lo).clamp(1, last_moving) - 1;
        let line = |piece_fixed: usize, piece_moving: usize| {
          let (fixed_intercept, fixed_slope) = self.fixed_lines[piece_fixed];
          let (moving_intercept, moving_slope) = self.moving_lines[piece_moving];
          (
            fixed_intercept + moving_intercept - moving_slope * shift,
            fixed_slope + moving_slope,
          )
        };
        let (mut intercept, mut slope) = line(piece_fixed, piece_moving);
        let mut start = lo;
        let mut neglog = intercept + slope * start;
        let mut density = (-neglog).exp();
        let mut total = 0.0;
        for _ in 0..last_fixed + last_moving {
          let next_fixed = fixed[piece_fixed + 1];
          let next_moving = moving[piece_moving + 1] + shift;
          let end = if next_fixed < next_moving {
            next_fixed
          } else {
            next_moving
          };
          let end = if end < hi { end } else { hi };
          let end_neglog = intercept + slope * end;
          let end_density = (-end_neglog).exp();
          total += exp_linear_integral(end - start, (neglog, density), (end_neglog, end_density));
          if end >= hi {
            break;
          }
          let advance_fixed = piece_fixed < last_fixed - 1 && next_fixed <= end;
          let advance_moving = piece_moving < last_moving - 1 && next_moving <= end;
          piece_fixed += usize::from(advance_fixed);
          piece_moving += usize::from(advance_moving);
          if advance_fixed || advance_moving {
            (intercept, slope) = line(piece_fixed, piece_moving);
          }
          (start, neglog, density) = (end, end_neglog, end_density);
        }
        total
      }
    }

    fn exp_linear_integral(width: f64, (y0, density0): (f64, f64), (y1, density1): (f64, f64)) -> f64 {
      let slope = y1 - y0;
      if slope.abs() < 1e-8 {
        width * (-f64::midpoint(y0, y1)).exp() * (1.0 + slope * slope / 24.0)
      } else {
        width * (density0 - density1) / slope
      }
    }
  }
}

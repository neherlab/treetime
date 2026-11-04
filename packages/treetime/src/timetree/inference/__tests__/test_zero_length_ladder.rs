#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::progress::NoopProgress;
  use crate::test_utils::{FullPosterior, full_posterior, parent_edge_key, point_date_constraints};
  use crate::timetree::confidence::CI_FRACTION;
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use crate::timetree::inference::forward_pass::{ForwardPassInputs, propagate_distributions_forward};
  use crate::timetree::inference::result::{BranchLikelihood, NodeBackward, log_marginal_likelihood};
  use crate::timetree::inference::runner::{EPS, GRID_POINTS};
  use crate::timetree::max_posterior_region::{MaxPosteriorRegion, max_posterior_region};
  use argmin::core::{CostFunction, Error, Executor};
  use argmin::solver::brent::BrentRoot;
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::Array1;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use std::ops::RangeInclusive;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, DistributionFunction, MassProfile, NegLog, TailMass, rewindow_to_mass};
  use treetime_graph::graph::Graph;
  use treetime_graph::key_map::{EdgeMap, NodeMap};
  use treetime_graph::node_names::NodeNames;
  use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, GridFn, MaxGridPoints, Side, SoftTailLaw};
  use treetime_io::nwk::nwk_read;
  use treetime_ops::FftPlans;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{make_error, make_report};

  const SITES: f64 = 754.0;
  const CLOCK_RATE: f64 = 0.0014;
  const TIP_DATE: f64 = 2023.0;
  const BRANCH_GRID_END: f64 = 10.0;
  const MOMENT_TAIL_MASS: f64 = 1e-14;
  const MOMENT_POINTS: usize = 40_001;
  const ROOT_TOLERANCE: f64 = 1e-13;
  const ROOT_MAX_ITERS: u64 = 200;
  const REGION_BRACKET_MASS: f64 = 1e-6;
  const SECOND_ORDER_OF_THE_SCHEME: RangeInclusive<f64> = 1.9..=2.1;
  const CI_LOWER_QUANTILE: f64 = 0.5 * (1.0 - CI_FRACTION);
  const CI_UPPER_QUANTILE: f64 = 1.0 - CI_LOWER_QUANTILE;

  #[rustfmt::skip]
  #[rstest]
  #[case::tips_20(  20, [300,  599, 1197])]
  #[case::tips_50(  50, [300,  599, 1197])]
  #[case::tips_200(200, [1197, 2393, 4785])]
  #[trace]
  fn test_zero_length_ladder_log_marginal_likelihood_converges_at_second_order_to_the_exact_value(
    #[case] n_tips: usize,
    #[case] branch_points: [usize; 3],
  ) -> Result<(), Report> {
    let expected = helpers::exact_log_z_from_independent_gaps(n_tips);

    let errors: Vec<f64> = branch_points
      .iter()
      .map(|&points| -> Result<f64, Report> {
        let ladder = helpers::Ladder::new(n_tips, points)?;
        let backward = ladder.backward(MaxGridPoints::default())?;
        let actual = log_marginal_likelihood(&backward, ladder.graph.root_key()?, None)?.expect("the root has a subtree");
        Ok((actual - expected).abs())
      })
      .try_collect()?;

    let reductions = [errors[0] / errors[1], errors[1] / errors[2]];
    assert!(
      reductions.iter().all(|&reduction| reduction >= 3.0),
      "errors {errors:?} fall by {reductions:?} per halving of the branch grid spacing"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tips_20(  20)]
  #[case::tips_50(  50)]
  #[case::tips_200(200)]
  #[trace]
  fn test_zero_length_ladder_root_date_posterior_converges_to_the_hypoexponential_root_age(
    #[case] n_tips: usize,
  ) -> Result<(), Report> {
    let expected = helpers::RootAge::new(n_tips).date_summary()?;

    let actual: Vec<helpers::DateSummary> = [GRID_POINTS, 4 * GRID_POINTS - 3, 16 * GRID_POINTS - 15]
      .iter()
      .map(|&points| -> Result<helpers::DateSummary, Report> {
        let ladder = helpers::Ladder::new(n_tips, points)?;
        let backward = ladder.backward(MaxGridPoints::default())?;
        let posterior = ladder.forward(backward, MaxGridPoints::default())?;
        let root = &posterior[ladder.graph.root_key()?];
        helpers::DateSummary::of(root.distribution.as_deref().expect("the root is dated"))
      })
      .try_collect()?;

    let errors = actual.iter().map(|summary| summary.max_difference(&expected)).collect_vec();
    let finest_pair_order = (errors[1] / errors[2]).log(4.0);
    assert!(
      SECOND_ORDER_OF_THE_SCHEME.contains(&finest_pair_order),
      "root date errors {errors:?} in years against the hypoexponential closed form converge at order {finest_pair_order} on the finest pair of grids, outside the second order {SECOND_ORDER_OF_THE_SCHEME:?} of the scheme: expected {expected:#?}, actual {actual:#?}"
    );
    Ok(())
  }

  #[test]
  fn test_zero_length_ladder_of_400_tips_dates_every_node_within_the_default_grid_point_limit() -> Result<(), Report> {
    let ladder = helpers::Ladder::new(400, GRID_POINTS)?;

    let backward = ladder.backward(MaxGridPoints::default())?;
    let posterior = ladder.forward(backward, MaxGridPoints::default())?;

    let undated = posterior.values().filter(|node| node.time.is_none()).count();
    assert_eq!(0, undated);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) struct Ladder {
      pub(super) graph: Graph,
      names: NodeNames,
      branches: EdgeMap<BranchLikelihood>,
      constraints: DateConstraints,
    }

    impl Ladder {
      pub(super) fn new(n_tips: usize, branch_points: usize) -> Result<Self, Report> {
        let newick = caterpillar_newick(n_tips);
        let parsed = nwk_read(newick.as_bytes())?;
        let names = parsed.names.clone();
        let graph = parsed.graph;
        let branch = zero_mutation_branch(branch_points)?;
        let branches = graph
          .get_nodes()
          .filter(|node| !node.is_root())
          .map(|node| {
            let branch = BranchLikelihood {
              distribution: Some(Arc::clone(&branch)),
              time_length: Some(0.0),
            };
            (parent_edge_key(&graph, node.key()), branch)
          })
          .collect();
        let tip_names = (0..n_tips).map(|i| format!("T{i}")).collect_vec();
        let dates = tip_names.iter().map(|name| (name.as_str(), TIP_DATE)).collect_vec();
        let constraints = point_date_constraints(&graph, &names, &dates);
        Ok(Self {
          graph,
          names,
          branches,
          constraints,
        })
      }

      pub(super) fn backward(&self, max_points: MaxGridPoints) -> Result<NodeMap<NodeBackward>, Report> {
        propagate_distributions_backward(
          &self.graph,
          &self.constraints,
          None,
          &self.bad_branches()?,
          &self.branches,
          max_points,
          &FftPlans::new(),
        )
      }

      pub(super) fn forward(
        &self,
        backward: NodeMap<NodeBackward>,
        max_points: MaxGridPoints,
      ) -> Result<NodeMap<FullPosterior>, Report> {
        let bad_branches = self.bad_branches()?;
        let inputs = ForwardPassInputs {
          constraints: &self.constraints,
          coalescent: None,
          bad_branches: &bad_branches,
          names: &self.names,
          branches: &self.branches,
          max_points,
          fft_plans: &FftPlans::new(),
        };
        propagate_distributions_forward(&self.graph, &inputs, backward, full_posterior, &NoopProgress)
      }

      fn bad_branches(&self) -> Result<NodeMap<bool>, Report> {
        derive_bad_branches(
          &self.graph,
          &self.constraints,
          &bad_leaves(&self.graph, &self.constraints, &BTreeSet::new()),
        )
      }
    }

    #[derive(Clone, Copy, Debug, PartialEq)]
    pub(super) struct DateSummary {
      pub mean: f64,
      pub sd: f64,
      pub quantiles: (f64, f64),
      pub region: (f64, f64),
    }

    impl DateSummary {
      pub(super) fn of(distribution: &Distribution<NegLog>) -> Result<Self, Report> {
        let Distribution::Function(function) = distribution.normalize()? else {
          return make_error!("The root date posterior is not a function");
        };
        let profile = MassProfile::new(&function)?;
        let fractions = profile
          .fractions(TailMass::ZeroWhereUndefined)
          .expect("the root date posterior holds mass");
        let (lower, upper) = (
          fractions.quantile(MOMENT_TAIL_MASS),
          fractions.quantile(1.0 - MOMENT_TAIL_MASS),
        );
        let t = Array1::linspace(lower, upper, MOMENT_POINTS);
        let below = t.mapv(|t| fractions.cdf(t));
        let mean_age = trapezoid(&t, &below);
        let second_age = trapezoid(&t, &(2.0 * (upper - &t) * &below));
        let Some(MaxPosteriorRegion::Grid(region)) = max_posterior_region(distribution, CI_FRACTION)? else {
          return make_error!("The root date posterior has no region on a grid");
        };
        Ok(Self {
          mean: upper - mean_age,
          sd: (second_age - mean_age.powi(2)).sqrt(),
          quantiles: (
            fractions.quantile(CI_LOWER_QUANTILE),
            fractions.quantile(CI_UPPER_QUANTILE),
          ),
          region: (region.lower, region.upper),
        })
      }

      pub(super) fn max_difference(&self, other: &Self) -> f64 {
        [
          self.mean - other.mean,
          self.sd - other.sd,
          self.quantiles.0 - other.quantiles.0,
          self.quantiles.1 - other.quantiles.1,
          self.region.0 - other.region.0,
          self.region.1 - other.region.1,
        ]
        .iter()
        .map(|difference| difference.abs())
        .fold(0.0, f64::max)
      }
    }

    pub(super) struct RootAge {
      n_tips: usize,
      lambda: f64,
    }

    impl RootAge {
      pub(super) fn new(n_tips: usize) -> Self {
        Self {
          n_tips,
          lambda: SITES * CLOCK_RATE,
        }
      }

      pub(super) fn date_summary(&self) -> Result<DateSummary, Report> {
        let rates = (2..=self.n_tips).map(|m| self.lambda * m.to_f64_exact()).collect_vec();
        let mean_age: f64 = rates.iter().map(|rate| rate.recip()).sum();
        let variance: f64 = rates.iter().map(|rate| rate.powi(-2)).sum();
        let (young, old) = self.region()?;
        Ok(DateSummary {
          mean: TIP_DATE - mean_age,
          sd: variance.sqrt(),
          quantiles: (
            TIP_DATE - self.quantile(CI_UPPER_QUANTILE)?,
            TIP_DATE - self.quantile(CI_LOWER_QUANTILE)?,
          ),
          region: (TIP_DATE - old, TIP_DATE - young),
        })
      }

      fn cdf(&self, age: f64) -> f64 {
        let n = self.n_tips.to_f64_exact();
        let u = (-self.lambda * age).exp();
        ((n - 1.0) * (-(-self.lambda * age).exp_m1()).ln() + ((n - 1.0) * u).ln_1p()).exp()
      }

      fn log_density(&self, age: f64) -> f64 {
        let n = self.n_tips.to_f64_exact();
        (self.lambda * n * (n - 1.0)).ln() - 2.0 * self.lambda * age + (n - 2.0) * (-(-self.lambda * age).exp_m1()).ln()
      }

      fn mode(&self) -> f64 {
        (0.5 * self.n_tips.to_f64_exact()).ln() / self.lambda
      }

      fn quantile(&self, q: f64) -> Result<f64, Report> {
        solve_root(|age| self.cdf(age) - q, 0.0, self.far_age())
      }

      fn region(&self) -> Result<(f64, f64), Report> {
        let mode = self.mode();
        let old_end = |young: f64| {
          let level = self.log_density(young);
          solve_root(|age| self.log_density(age) - level, mode, self.far_age())
        };
        let mass =
          |young: f64| -> Result<f64, Report> { Ok(self.cdf(old_end(young)?) - self.cdf(young) - CI_FRACTION) };
        let young = solve_root_fallible(mass, self.quantile(REGION_BRACKET_MASS)?, mode)?;
        Ok((young, old_end(young)?))
      }

      fn far_age(&self) -> f64 {
        (self.mode() + 60.0 / self.lambda).max(1.0)
      }
    }

    fn trapezoid(x: &Array1<f64>, y: &Array1<f64>) -> f64 {
      let h = x[1] - x[0];
      h * (y.sum() - f64::midpoint(y[0], y[y.len() - 1]))
    }

    fn solve_root(f: impl Fn(f64) -> f64, lower: f64, upper: f64) -> Result<f64, Report> {
      solve_root_fallible(|x| Ok(f(x)), lower, upper)
    }

    fn solve_root_fallible(f: impl Fn(f64) -> Result<f64, Report>, lower: f64, upper: f64) -> Result<f64, Report> {
      let result = Executor::new(RootProblem(f), BrentRoot::new(lower, upper, ROOT_TOLERANCE))
        .configure(|cfg| cfg.max_iters(ROOT_MAX_ITERS))
        .run()
        .map_err(|e| make_report!("Root of the closed form on [{lower}, {upper}] not found: {e}"))?;
      result
        .state
        .best_param
        .ok_or_else(|| make_report!("Root solve on [{lower}, {upper}] returned no value"))
    }

    struct RootProblem<F>(F);

    impl<F: Fn(f64) -> Result<f64, Report>> CostFunction for RootProblem<F> {
      type Param = f64;
      type Output = f64;

      fn cost(&self, x: &f64) -> Result<f64, Error> {
        (self.0)(*x).map_err(|report| Error::msg(format!("{report:#}")))
      }
    }

    pub(super) fn exact_log_z_from_independent_gaps(n_tips: usize) -> f64 {
      let lambda = SITES * CLOCK_RATE;
      (0..n_tips - 1).map(|j| -(lambda * (j + 2).to_f64_exact()).ln()).sum()
    }

    fn caterpillar_newick(n_tips: usize) -> String {
      let inner = (1..n_tips).fold(String::from("T0:0"), |subtree, i| {
        format!("({subtree},T{i}:0)N{}:0", n_tips - 1 - i)
      });
      let root = inner
        .strip_suffix(":0")
        .expect("the caterpillar ends with a branch length");
      format!("{root};")
    }

    fn zero_mutation_branch(points: usize) -> Result<Arc<Distribution<NegLog>>, Report> {
      let lambda = SITES * CLOCK_RATE;
      let t = Array1::linspace(0.0, BRANCH_GRID_END, points);
      let f = GridFn::from_range_values((0.0, BRANCH_GRID_END), t.mapv(|t| lambda * t))?;
      let right = SoftTailLaw::fit(&f, Side::Right, DEFAULT_TAIL_FIT_POINTS)?;
      let distribution = DistributionFunction::from_grid_fn(f)
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Linear(right))?;
      let distribution = rewindow_to_mass(
        &Distribution::Function(distribution),
        EPS,
        points,
        MaxGridPoints::default(),
      )?;
      Ok(Arc::new(distribution))
    }
  }
}

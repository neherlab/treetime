#[cfg(test)]
mod tests {
  use crate::coalescent::mean_field::{CoalescentPrior, LineageExpectation};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::GTR;
  use crate::optimize::indel::IndelTerm;
  use crate::partition::optimize;
  use crate::partition::optimize::contribution::OptimizationContribution;
  use crate::progress::NoopProgress;
  use crate::test_utils::{find_node_key_by_name, full_posterior, parent_edge_key, point_date_constraints};
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use crate::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
  use crate::timetree::inference::forward_pass::{ForwardPassInputs, propagate_distributions_forward};
  use crate::timetree::inference::result::{BranchLikelihood, NodePosterior, PosteriorSummaries};
  use crate::timetree::inference::runner::GRID_POINTS;
  use eyre::Report;
  use itertools::{Itertools, izip};
  use ndarray::Array2;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use std::sync::Arc;
  use treetime_distribution::{Distribution, MassProfile, NegLog, TailMass};
  use treetime_grid::MaxGridPoints;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};

  const RATE: f64 = 2.0;
  const SITES: usize = 1000;
  const T_A: f64 = 2012.0;
  const T_B: f64 = 2010.0;
  const T_C: f64 = 2011.0;
  const GRID_START: f64 = 2000.0;
  const GRID_STEP: f64 = 0.005;
  const GRID_LEN: usize = 2401;

  // TODO: kb/issues/M-timetree-time-pass-posteriors-deviate-from-the-exact-posterior.md
  // The CDFs of the time passes differ from the exact posterior by up to 2.9e-3 (node R, both tips at zero).
  #[rustfmt::skip]
  #[rstest]
  #[case::mixed_counts(           [3, 0, 2, 6],  [1e-4, 1e-4])]
  #[case::zero_counts_near_tips(  [1, 0, 0, 2],  [1e-3, 1e-3])]
  #[case::long_branches(          [8, 4, 5, 12], [1e-3, 1e-3])]
  #[case::both_tips_at_zero(      [0, 0, 1, 1],  [1e-4, 1e-2])]
  #[trace]
  fn test_marginal_exact_time_passes_match_the_integrated_posterior_without_a_prior(
    #[case] counts: [usize; 4],
    #[case] max_cdf_difference: [f64; 2],
  ) -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let grid: Vec<f64> = (0..GRID_LEN).map(|i| GRID_START + GRID_STEP * i.to_f64_exact()).collect();
    let (exact_x, exact_r) = helpers::exact_marginals(&gtr, counts, &grid, GRID_STEP);
    let (treetime_x, treetime_r) = helpers::treetime_marginal_cdfs(&gtr, counts, &grid)?;

    let support_end_x = T_A.min(T_B);
    let support_end_r = support_end_x.min(T_C);
    for ((node, exact, support_end, treetime), max_difference) in [
      ("X", exact_x, support_end_x, treetime_x),
      ("R", exact_r, support_end_r, treetime_r),
    ]
    .into_iter()
    .zip(max_cdf_difference)
    {
      let difference = helpers::max_cdf_difference(&grid, &exact, support_end, &treetime);
      assert!(
        difference < max_difference,
        "node {node}: max CDF difference {difference:.4e} to the exact posterior, above {max_difference:e}"
      );
    }
    Ok(())
  }

  const TIPS_4: [(&str, f64); 4] = [("A", 2012.0), ("B", 2010.5), ("C", 2011.5), ("D", 2010.0)];
  const ORACLE_START: f64 = 1995.0;
  const ORACLE_STEP: f64 = 0.00125;
  const BRANCH_GRID_LADDER: [usize; 3] = [GRID_POINTS, 2 * GRID_POINTS - 1, 4 * GRID_POINTS - 3];
  const MIN_ERROR_REDUCTION_PER_HALVING: f64 = 3.0;
  const MAX_ORACLE_ERROR_SHARE: f64 = 0.1;

  #[rustfmt::skip]
  #[rstest]
  #[case::short_internal_branches_strong_prior(   [3, 0, 2, 1, 1, 2],  1.0)]
  #[case::short_internal_branches_weak_prior(     [3, 0, 2, 1, 1, 2], 10.0)]
  #[case::overlapping_clades_strong_prior(        [6, 3, 4, 3, 3, 4],  1.0)]
  #[case::overlapping_clades_weak_prior(          [6, 3, 4, 3, 3, 4], 10.0)]
  #[case::long_internal_branches_strong_prior(    [2, 1, 1, 2, 3, 3],  1.0)]
  #[trace]
  fn test_marginal_exact_mean_field_passes_converge_to_the_exact_mean_field_posterior(
    #[case] counts: [usize; 6],
    #[case] tc: f64,
  ) -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let tree = helpers::MeanFieldTree::new(&gtr, counts, tc)?;

    let nodes = tree.convergence_to_the_exact_mean_field_oracle(&gtr)?;

    let oracle_resolves_the_finest_level = nodes
      .values()
      .all(|node| node.oracle_error <= MAX_ORACLE_ERROR_SHARE * node.errors[BRANCH_GRID_LADDER.len() - 1]);
    assert!(
      oracle_resolves_the_finest_level,
      "the error of the exact mean-field oracle must stay below {MAX_ORACLE_ERROR_SHARE} of the error of the finest branch grid: {nodes:#?}"
    );
    let converged = nodes.values().all(|node| {
      node
        .errors
        .iter()
        .tuple_windows()
        .all(|(coarse, fine)| coarse / fine >= MIN_ERROR_REDUCTION_PER_HALVING)
    });
    assert!(
      converged,
      "the CDF error to the exact mean-field oracle must fall at least {MIN_ERROR_REDUCTION_PER_HALVING} times per halving of the branch grid spacing: {nodes:#?}"
    );
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::clock::date_constraints::DateConstraints;
    use crate::timetree::inference::result::PosteriorView;
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::{EdgeMap, NodeMap};
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::node_names::NodeNames;

    fn transition_probabilities(gtr: &GTR, b: f64) -> (f64, f64) {
      let eigvals = gtr.rate_eigvals();
      let p = |a: usize, c: usize| {
        (0..eigvals.len())
          .map(|j| gtr.v()[[a, j]] * gtr.v_inv()[[j, c]] * (eigvals[j] * b).exp())
          .sum::<f64>()
      };
      (p(0, 0), p(0, 1))
    }

    fn branch_likelihood(gtr: &GTR, tau: f64, k: usize) -> f64 {
      log_branch_likelihood(gtr, tau, k).exp()
    }

    fn log_branch_likelihood(gtr: &GTR, tau: f64, k: usize) -> f64 {
      if tau < 0.0 || (tau == 0.0 && k > 0) {
        return f64::NEG_INFINITY;
      }
      let (same, other) = transition_probabilities(gtr, RATE / SITES.to_f64_exact() * tau);
      let differing = if k == 0 { 0.0 } else { k.to_f64_exact() * other.ln() };
      (SITES - k).to_f64_exact() * same.ln() + differing
    }

    pub(super) fn exact_marginals(
      gtr: &GTR,
      [ka, kb, kx, kc]: [usize; 4],
      grid: &[f64],
      step: f64,
    ) -> (Vec<f64>, Vec<f64>) {
      let n = grid.len();
      let index = |date: f64| f64_to_usize_saturating(((date - grid[0]) / step).round());
      let lag = |steps: usize, k: usize| branch_likelihood(gtr, steps.to_f64_exact() * step, k);
      let (index_a, index_b, end_r) = (index(T_A), index(T_B), index(T_C));
      let end_x = index_a.min(index_b);
      let branch_x: Vec<f64> = (0..=end_x).map(|lag_steps| lag(lag_steps, kx)).collect();
      let tips_x: Vec<f64> = (0..=end_x)
        .map(|i| lag(index_a - i, ka) * lag(index_b - i, kb))
        .collect();
      let tip_r: Vec<f64> = (0..=end_r).map(|j| lag(end_r - j, kc)).collect();
      let mut marginal_x = vec![0.0; n];
      let mut marginal_r = vec![0.0; n];
      for i in 0..=end_x {
        let last_j = i.min(end_r);
        for j in 0..=last_j {
          let value = tips_x[i] * branch_x[i - j] * tip_r[j];
          marginal_x[i] += trapezoid_weight(j, 0, last_j) * value;
          marginal_r[j] += trapezoid_weight(i, j, end_x) * value;
        }
      }
      (marginal_x, marginal_r)
    }

    fn trapezoid_weight(position: usize, first: usize, last: usize) -> f64 {
      if first == last {
        0.0
      } else if position == first || position == last {
        0.5
      } else {
        1.0
      }
    }

    pub(super) fn treetime_marginal_cdfs(
      gtr: &GTR,
      counts: [usize; 4],
      grid: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), Report> {
      let parsed = nwk_read(b"((A:1,B:1)X:1,C:1)R;".as_slice())?;
      let names = parsed.names.clone();
      let graph = parsed.graph;
      let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");
      let branches: EdgeMap<_> = ["A", "B", "X", "C"]
        .iter()
        .zip(counts)
        .map(|(name, k)| {
          let distribution = branch_distribution(gtr, k, GRID_POINTS)?;
          let time_length = distribution.likely_time()?;
          let branch = BranchLikelihood {
            distribution: Some(distribution),
            time_length,
          };
          Ok((parent_edge_key(&graph, key(name)), branch))
        })
        .collect::<Result<_, Report>>()?;
      let constraints = point_date_constraints(&graph, &names, &[("A", T_A), ("B", T_B), ("C", T_C)]);
      let bad_branches = derive_bad_branches(
        &graph,
        &constraints,
        &bad_leaves(&graph, &constraints, &BTreeSet::new()),
      )?;
      let backward = propagate_distributions_backward(
        &graph,
        &constraints,
        None,
        &bad_branches,
        &branches,
        MaxGridPoints::default(),
      )?;
      let forward_inputs = ForwardPassInputs {
        constraints: &constraints,
        coalescent: None,
        bad_branches: &bad_branches,
        names: &names,
        branches: &branches,
        max_points: MaxGridPoints::default(),
      };
      let posterior =
        propagate_distributions_forward(&graph, &forward_inputs, backward, full_posterior, &NoopProgress)?;
      let cdf_of = |name: &str| -> Result<Vec<f64>, Report> {
        let distribution = posterior[&key(name)]
          .distribution
          .as_deref()
          .expect("internal node must have a posterior");
        posterior_cdf(distribution, grid)
      };
      Ok((cdf_of("X")?, cdf_of("R")?))
    }

    fn branch_distribution(gtr: &GTR, k: usize, points: usize) -> Result<Arc<Distribution<NegLog>>, Report> {
      let n_states = gtr.v().nrows();
      let mut coefficients = Array2::<f64>::zeros((SITES, n_states));
      for site in 0..SITES {
        let child = usize::from(site < k);
        for j in 0..n_states {
          coefficients[[site, j]] = gtr.v()[[0, j]] * gtr.v_inv()[[j, child]];
        }
      }
      let contribution =
        OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr.clone()));
      let one_mutation = 1.0 / SITES.to_f64_exact();
      compute_branch_length_distribution(
        &[contribution],
        IndelTerm::new(0, 0.0),
        (k.to_f64_exact() * one_mutation).max(one_mutation),
        one_mutation,
        points,
        RATE * one_mutation,
        MaxGridPoints::default(),
      )
    }

    fn support_ends() -> [f64; 3] {
      let date = |name: &str| {
        TIPS_4
          .iter()
          .find(|(tip, _)| *tip == name)
          .map(|(_, date)| *date)
          .expect("tip date")
      };
      let x = date("A").min(date("B"));
      let y = date("C").min(date("D"));
      [x, y, x.min(y)]
    }

    pub(super) struct MeanFieldTree {
      tree: FourTipTree,
      prior: CoalescentPrior,
    }

    impl MeanFieldTree {
      pub(super) fn new(gtr: &GTR, counts: [usize; 6], tc: f64) -> Result<Self, Report> {
        let tree = FourTipTree::new(counts)?;
        let tc = Distribution::constant(tc);
        let branches = tree.branches(gtr, GRID_POINTS)?;
        let prior_free = tree.passes(&branches, None, survival_summary)?;
        let seed = CoalescentPrior::new(&LineageExpectation::from_posterior(&tree.graph, &prior_free)?, &tc)?;
        let first = tree.passes(&branches, Some(&seed), survival_summary)?;
        let prior = CoalescentPrior::new(&LineageExpectation::from_posterior(&tree.graph, &first)?, &tc)?;
        Ok(Self { tree, prior })
      }

      fn passes_cdfs(&self, gtr: &GTR, branch_points: usize, grid: &[f64]) -> Result<[Vec<f64>; 3], Report> {
        let branches = self.tree.branches(gtr, branch_points)?;
        let posterior = self.tree.passes(&branches, Some(&self.prior), full_posterior)?;
        let cdf_of = |name: &str| -> Result<Vec<f64>, Report> {
          let distribution = posterior[&self.tree.key(name)]
            .distribution
            .as_deref()
            .expect("internal node must have a posterior");
          posterior_cdf(distribution, grid)
        };
        Ok([cdf_of("X")?, cdf_of("Y")?, cdf_of("R")?])
      }

      pub(super) fn convergence_to_the_exact_mean_field_oracle(
        &self,
        gtr: &GTR,
      ) -> Result<BTreeMap<&'static str, NodeConvergence>, Report> {
        let oracle = self.exact_mean_field_oracle(gtr, ORACLE_STEP)?;
        let coarse_oracle = self.exact_mean_field_oracle(gtr, 2.0 * ORACLE_STEP)?;
        let levels: Vec<[Vec<f64>; 3]> = BRANCH_GRID_LADDER
          .iter()
          .map(|&points| self.passes_cdfs(gtr, points, &oracle.grid))
          .try_collect()?;
        Ok(
          izip!(["X", "Y", "R"], &oracle.cdfs, &coarse_oracle.cdfs, 0..)
            .map(|(node, exact, coarse, n)| {
              let convergence = NodeConvergence {
                errors: levels
                  .iter()
                  .map(|level| max_abs_difference(&level[n], exact))
                  .collect(),
                oracle_error: richardson_error(exact, coarse),
              };
              (node, convergence)
            })
            .collect(),
        )
      }

      fn exact_mean_field_oracle(&self, gtr: &GTR, step: f64) -> Result<ExactMarginals, Report> {
        let [ka, kb, kc, kd, kx, ky] = self.tree.counts;
        let ends = support_ends();
        let [end_x, end_y, _] = ends;
        let n = f64_to_usize_saturating(((end_x - ORACLE_START) / step).round()) + 1;
        let grid: Vec<f64> = (0..n).map(|i| ORACLE_START + step * i.to_f64_exact()).collect();
        let index = |date: f64| f64_to_usize_saturating(((date - ORACLE_START) / step).round());
        let (last_x, last_y) = (index(end_x), index(end_y));
        let tip = |name: &str| TIPS_4.iter().find(|(tip, _)| *tip == name).expect("tip must exist").1;
        let prior_of = |name: &str| -> Result<Vec<f64>, Report> {
          let key = self.tree.key(name);
          grid
            .iter()
            .map(|&t| self.prior.node_contribution(key, 2, t).map(|f| -f))
            .collect()
        };
        let (prior_x, prior_y, prior_r) = (prior_of("X")?, prior_of("Y")?, prior_of("R")?);
        let factor_x = rescaled_exp(
          (0..n)
            .map(|i| {
              log_branch_likelihood(gtr, tip("A") - grid[i], ka)
                + log_branch_likelihood(gtr, tip("B") - grid[i], kb)
                + prior_x[i]
            })
            .collect(),
        );
        let factor_y = rescaled_exp(
          (0..n)
            .map(|j| {
              log_branch_likelihood(gtr, tip("C") - grid[j], kc)
                + log_branch_likelihood(gtr, tip("D") - grid[j], kd)
                + prior_y[j]
            })
            .collect(),
        );
        let factor_r = rescaled_exp(prior_r);
        let lag = |k: usize| {
          rescaled_exp(
            (0..n)
              .map(|d| log_branch_likelihood(gtr, step * d.to_f64_exact(), k))
              .collect(),
          )
        };
        let (branch_x, branch_y) = (lag(kx), lag(ky));
        let below = |factor: &[f64], branch: &[f64], last: usize| -> Vec<f64> {
          (0..n)
            .map(|l| {
              if l > last {
                return 0.0;
              }
              (l..=last)
                .map(|i| trapezoid_weight(i, l, last) * factor[i] * branch[i - l])
                .sum()
            })
            .collect()
        };
        let (from_x, from_y) = (below(&factor_x, &branch_x, last_x), below(&factor_y, &branch_y, last_y));
        let above = |factor: &[f64], branch: &[f64], last: usize, sibling: &[f64], last_sibling: usize| -> Vec<f64> {
          (0..n)
            .map(|i| {
              if i > last {
                return 0.0;
              }
              let top = i.min(last_sibling);
              factor[i]
                * (0..=top)
                  .map(|l| trapezoid_weight(l, 0, top) * branch[i - l] * factor_r[l] * sibling[l])
                  .sum::<f64>()
            })
            .collect()
        };
        let density_x = above(&factor_x, &branch_x, last_x, &from_y, last_y);
        let density_y = above(&factor_y, &branch_y, last_y, &from_x, last_x);
        let density_r: Vec<f64> = (0..n).map(|l| factor_r[l] * from_x[l] * from_y[l]).collect();
        let cdfs = [
          cdf(&grid, &density_x, ends[0]),
          cdf(&grid, &density_y, ends[1]),
          cdf(&grid, &density_r, ends[2]),
        ];
        Ok(ExactMarginals { grid, cdfs })
      }
    }

    struct FourTipTree {
      graph: Graph,
      names: NodeNames,
      counts: [usize; 6],
      constraints: DateConstraints,
      bad_branches: NodeMap<bool>,
    }

    impl FourTipTree {
      fn new(counts: [usize; 6]) -> Result<Self, Report> {
        let parsed = nwk_read(b"((A:1,B:1)X:1,(C:1,D:1)Y:1)R;".as_slice())?;
        let names = parsed.names.clone();
        let graph = parsed.graph;
        let constraints = point_date_constraints(&graph, &names, &TIPS_4);
        let bad_branches = derive_bad_branches(
          &graph,
          &constraints,
          &bad_leaves(&graph, &constraints, &BTreeSet::new()),
        )?;
        Ok(Self {
          graph,
          names,
          counts,
          constraints,
          bad_branches,
        })
      }

      fn key(&self, name: &str) -> GraphNodeKey {
        find_node_key_by_name(&self.graph, &self.names, name).expect("fixture node must exist")
      }

      fn branches(&self, gtr: &GTR, points: usize) -> Result<EdgeMap<BranchLikelihood>, Report> {
        ["A", "B", "C", "D", "X", "Y"]
          .iter()
          .zip(self.counts)
          .map(|(name, k)| {
            let distribution = branch_distribution(gtr, k, points)?;
            let time_length = distribution.likely_time()?;
            let branch = BranchLikelihood {
              distribution: Some(distribution),
              time_length,
            };
            Ok((parent_edge_key(&self.graph, self.key(name)), branch))
          })
          .collect()
      }

      fn passes<S: Send + Sync>(
        &self,
        branches: &EdgeMap<BranchLikelihood>,
        prior: Option<&CoalescentPrior>,
        summarize: impl Fn(GraphNodeKey, PosteriorView<'_>) -> Result<S, Report> + Sync,
      ) -> Result<NodeMap<S>, Report> {
        let backward = propagate_distributions_backward(
          &self.graph,
          &self.constraints,
          prior,
          &self.bad_branches,
          branches,
          MaxGridPoints::default(),
        )?;
        let forward_inputs = ForwardPassInputs {
          constraints: &self.constraints,
          coalescent: None,
          bad_branches: &self.bad_branches,
          names: &self.names,
          branches,
          max_points: MaxGridPoints::default(),
        };
        propagate_distributions_forward(&self.graph, &forward_inputs, backward, summarize, &NoopProgress)
      }
    }

    fn survival_summary(key: GraphNodeKey, view: PosteriorView<'_>) -> Result<NodePosterior, Report> {
      let summaries = PosteriorSummaries {
        lineage_survival: true,
        ..PosteriorSummaries::default()
      };
      NodePosterior::summarize(key, view, summaries)
    }

    struct ExactMarginals {
      grid: Vec<f64>,
      cdfs: [Vec<f64>; 3],
    }

    #[derive(Debug)]
    pub(super) struct NodeConvergence {
      pub(super) errors: Vec<f64>,
      pub(super) oracle_error: f64,
    }

    fn rescaled_exp(log_values: Vec<f64>) -> Vec<f64> {
      let max = log_values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
      log_values.into_iter().map(|v| (v - max).exp()).collect()
    }

    fn max_abs_difference(a: &[f64], b: &[f64]) -> f64 {
      a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
    }

    fn richardson_error(fine: &[f64], coarse: &[f64]) -> f64 {
      let on_coarse_points = fine.iter().step_by(2).copied().collect_vec();
      max_abs_difference(&on_coarse_points, coarse) / 3.0
    }

    fn cdf(grid: &[f64], density: &[f64], support_end: f64) -> Vec<f64> {
      let mut cumulative = vec![0.0; grid.len()];
      for i in 1..grid.len() {
        let cell = if grid[i] - support_end < 0.5 * (grid[i] - grid[i - 1]) {
          f64::midpoint(density[i], density[i - 1]) * (grid[i] - grid[i - 1])
        } else {
          0.0
        };
        cumulative[i] = cumulative[i - 1] + cell;
      }
      let total = cumulative[grid.len() - 1];
      cumulative.iter().map(|c| c / total).collect()
    }

    pub(super) fn max_cdf_difference(
      grid: &[f64],
      exact_density: &[f64],
      support_end: f64,
      treetime_cdf: &[f64],
    ) -> f64 {
      cdf(grid, exact_density, support_end)
        .iter()
        .zip(treetime_cdf)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
    }

    fn posterior_cdf(distribution: &Distribution<NegLog>, grid: &[f64]) -> Result<Vec<f64>, Report> {
      let Distribution::Function(f) = distribution else {
        return Err(eyre::eyre!("an internal node posterior must be a function"));
      };
      let profile = MassProfile::new(f)?;
      let fractions = profile
        .fractions(TailMass::Required)
        .ok_or_else(|| eyre::eyre!("an internal node posterior must have a finite mass"))?;
      Ok(grid.iter().map(|&t| fractions.cdf(t)).collect())
    }
  }
}

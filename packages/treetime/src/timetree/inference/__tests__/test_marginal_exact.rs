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
  const GRID_4_START: f64 = 2003.0;
  const GRID_4_STEP: f64 = 0.02;
  const GRID_4_LEN: usize = 376;

  // TODO: kb/issues/M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md
  // The CDFs of the mean-field time passes differ from the exact Kingman posterior by up to 1.6e-2 (node X, overlapping clades, strong prior).
  #[rustfmt::skip]
  #[rstest]
  #[case::short_internal_branches_strong_prior(   [3, 0, 2, 1, 1, 2],  1.0, [1e-2, 1e-1, 1e-2])]
  #[case::short_internal_branches_weak_prior(     [3, 0, 2, 1, 1, 2], 10.0, [1e-3, 1e-2, 1e-2])]
  #[case::overlapping_clades_strong_prior(        [6, 3, 4, 3, 3, 4],  1.0, [1e-1, 1e-1, 1e-2])]
  #[case::overlapping_clades_weak_prior(          [6, 3, 4, 3, 3, 4], 10.0, [1e-2, 1e-2, 1e-2])]
  #[case::long_internal_branches_strong_prior(    [2, 1, 1, 2, 3, 3],  1.0, [1e-2, 1e-2, 1e-2])]
  #[trace]
  fn test_marginal_exact_mean_field_prior_matches_the_kingman_posterior(
    #[case] counts: [usize; 6],
    #[case] tc: f64,
    #[case] max_cdf_difference: [f64; 3],
  ) -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let grid: Vec<f64> = (0..GRID_4_LEN).map(|i| GRID_4_START + GRID_4_STEP * i.to_f64_exact()).collect();
    let exact = helpers::exact_kingman_marginals(&gtr, counts, tc, &grid);
    let treetime = helpers::treetime_mean_field_marginal_cdfs(&gtr, counts, tc, &grid)?;

    let date = |name: &str| TIPS_4.iter().find(|(tip, _)| *tip == name).map(|(_, date)| *date).expect("tip date");
    let support_end_x = date("A").min(date("B"));
    let support_end_y = date("C").min(date("D"));
    let support_ends = [support_end_x, support_end_y, support_end_x.min(support_end_y)];
    for ((((node, exact), support_end), treetime), max_difference) in
      ["X", "Y", "R"].iter().zip(&exact).zip(support_ends).zip(&treetime).zip(max_cdf_difference)
    {
      let difference = helpers::max_cdf_difference(&grid, exact, support_end, treetime);
      assert!(
        difference < max_difference,
        "node {node}: max CDF difference {difference:.4e} to the exact Kingman posterior, above {max_difference:e}"
      );
    }
    Ok(())
  }

  mod helpers {
    use super::*;

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
      if tau < 0.0 || (tau == 0.0 && k > 0) {
        return 0.0;
      }
      let (same, other) = transition_probabilities(gtr, RATE / SITES.to_f64_exact() * tau);
      let differing = if k == 0 { 0.0 } else { k.to_f64_exact() * other.ln() };
      ((SITES - k).to_f64_exact() * same.ln() + differing).exp()
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
      let branches: BTreeMap<_, _> = ["A", "B", "X", "C"]
        .iter()
        .zip(counts)
        .map(|(name, k)| {
          let distribution = branch_distribution(gtr, k)?;
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

    fn branch_distribution(gtr: &GTR, k: usize) -> Result<Arc<Distribution<NegLog>>, Report> {
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
        GRID_POINTS,
        RATE * one_mutation,
        MaxGridPoints::default(),
      )
    }

    pub(super) fn exact_kingman_marginals(gtr: &GTR, counts: [usize; 6], tc: f64, grid: &[f64]) -> [Vec<f64>; 3] {
      let [ka, kb, kc, kd, kx, ky] = counts;
      let n = grid.len();
      let tip = |name: &str| TIPS_4.iter().find(|(tip, _)| *tip == name).expect("tip must exist").1;
      let mut linear = vec![0.0; n];
      for i in 1..n {
        let mid = f64::midpoint(grid[i - 1], grid[i]);
        let sampled_after = TIPS_4.iter().filter(|(_, date)| *date > mid).count().to_f64_exact();
        linear[i] = linear[i - 1] + GRID_4_STEP * (sampled_after - 1.0) / tc;
      }
      let lag = |k: usize| {
        (0..n)
          .map(|i| branch_likelihood(gtr, i.to_f64_exact() * GRID_4_STEP, k))
          .collect::<Vec<_>>()
      };
      let (branch_x, branch_y) = (lag(kx), lag(ky));
      let factor_x: Vec<f64> = (0..n)
        .map(|i| {
          let x = grid[i];
          branch_likelihood(gtr, tip("A") - x, ka) * branch_likelihood(gtr, tip("B") - x, kb) * linear[i].exp()
        })
        .collect();
      let factor_y: Vec<f64> = (0..n)
        .map(|j| {
          let y = grid[j];
          branch_likelihood(gtr, tip("C") - y, kc) * branch_likelihood(gtr, tip("D") - y, kd) * linear[j].exp()
        })
        .collect();
      let factor_r: Vec<f64> = (0..n)
        .map(|l| (linear[l] - 2.0 * (grid[l] - grid[0]) / tc).exp())
        .collect();
      let (mut px, mut py, mut pr) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
      for i in 0..n {
        for j in 0..n {
          let pair = factor_x[i] * factor_y[j] * (-(grid[i.min(j)] - grid[0]) / tc).exp();
          if pair == 0.0 {
            continue;
          }
          for l in 0..i.min(j) {
            let weight = pair * branch_x[i - l] * branch_y[j - l] * factor_r[l];
            px[i] += weight;
            py[j] += weight;
            pr[l] += weight;
          }
        }
      }
      [px, py, pr]
    }

    pub(super) fn treetime_mean_field_marginal_cdfs(
      gtr: &GTR,
      counts: [usize; 6],
      tc: f64,
      grid: &[f64],
    ) -> Result<[Vec<f64>; 3], Report> {
      let parsed = nwk_read(b"((A:1,B:1)X:1,(C:1,D:1)Y:1)R;".as_slice())?;
      let names = parsed.names.clone();
      let graph = parsed.graph;
      let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");
      let branches: BTreeMap<_, _> = ["A", "B", "C", "D", "X", "Y"]
        .iter()
        .zip(counts)
        .map(|(name, k)| {
          let distribution = branch_distribution(gtr, k)?;
          let time_length = distribution.likely_time()?;
          let branch = BranchLikelihood {
            distribution: Some(distribution),
            time_length,
          };
          Ok((parent_edge_key(&graph, key(name)), branch))
        })
        .collect::<Result<_, Report>>()?;
      let constraints = point_date_constraints(&graph, &names, &TIPS_4);
      let bad_branches = derive_bad_branches(
        &graph,
        &constraints,
        &bad_leaves(&graph, &constraints, &BTreeSet::new()),
      )?;
      let backward = |prior: Option<&CoalescentPrior>| {
        propagate_distributions_backward(
          &graph,
          &constraints,
          prior,
          &bad_branches,
          &branches,
          MaxGridPoints::default(),
        )
      };
      let forward_inputs = ForwardPassInputs {
        constraints: &constraints,
        names: &names,
        branches: &branches,
        max_points: MaxGridPoints::default(),
      };
      let summaries = PosteriorSummaries {
        lineage_survival: true,
        ..PosteriorSummaries::default()
      };
      let passes = |prior: Option<&CoalescentPrior>| {
        propagate_distributions_forward(
          &graph,
          &forward_inputs,
          backward(prior)?,
          |key, view| NodePosterior::summarize(key, view, summaries),
          &NoopProgress,
        )
      };
      let tc = Distribution::constant(tc);
      let prior_free = passes(None)?;
      let seed = CoalescentPrior::new(&LineageExpectation::from_posterior(&graph, &prior_free)?, &tc)?;
      let first = passes(Some(&seed))?;
      let rebuilt = CoalescentPrior::new(&LineageExpectation::from_posterior(&graph, &first)?, &tc)?;
      let posterior = propagate_distributions_forward(
        &graph,
        &forward_inputs,
        backward(Some(&rebuilt))?,
        full_posterior,
        &NoopProgress,
      )?;
      let cdf_of = |name: &str| -> Result<Vec<f64>, Report> {
        let distribution = posterior[&key(name)]
          .distribution
          .as_deref()
          .expect("internal node must have a posterior");
        posterior_cdf(distribution, grid)
      };
      Ok([cdf_of("X")?, cdf_of("Y")?, cdf_of("R")?])
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

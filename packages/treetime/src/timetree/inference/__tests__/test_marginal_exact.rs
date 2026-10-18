#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::GTR;
  use crate::optimize::indel::IndelTerm;
  use crate::partition::optimize;
  use crate::partition::optimize::contribution::OptimizationContribution;
  use crate::progress::NoopProgress;
  use crate::test_utils::{find_node_key_by_name, parent_edge_key, point_date_constraints};
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use crate::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
  use crate::timetree::inference::forward_pass::propagate_distributions_forward;
  use crate::timetree::inference::result::BranchLikelihood;
  use crate::timetree::inference::runner::GRID_POINTS;
  use eyre::Report;
  use ndarray::Array2;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_grid::MaxGridPoints;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::cast::ToF64Exact;

  // ((A,B)X,C)R with heterochronous tips: X is bounded by B, the root by no tip directly. Each branch has k of SITES
  // sites differing under JC69 at RATE substitutions per year over all sites. Without a tree prior the joint density
  // of the internal times is the product of the branch likelihoods,
  //   p(x, r) ∝ L_A(t_A - x) L_B(t_B - x) L_X(x - r) L_C(t_C - r),   r < x < min(t_A, t_B),
  // and its marginals, integrated on a uniform grid, are the exact node posteriors that the time passes must match.

  const RATE: f64 = 2.0;
  const SITES: usize = 1000;
  const T_A: f64 = 2012.0;
  const T_B: f64 = 2010.0;
  const T_C: f64 = 2011.0;
  const GRID_START: f64 = 2000.0;
  const GRID_STEP: f64 = 0.005;
  const GRID_LEN: usize = 2401;
  const MAX_CDF_DIFFERENCE: f64 = 0.02;

  #[rustfmt::skip]
  #[rstest]
  #[case::mixed_counts(           [3, 0, 2, 6])]
  #[case::zero_counts_near_tips(  [1, 0, 0, 2])]
  #[case::long_branches(          [8, 4, 5, 12])]
  #[case::both_tips_at_zero(      [0, 0, 1, 1])]
  #[trace]
  fn test_marginal_exact_time_passes_match_the_integrated_posterior_without_a_prior(
    #[case] counts: [usize; 4],
  ) -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let grid: Vec<f64> = (0..GRID_LEN).map(|i| GRID_START + GRID_STEP * i.to_f64_exact()).collect();
    let (exact_x, exact_r) = helpers::exact_marginals(&gtr, counts, &grid);
    let (treetime_x, treetime_r) = helpers::treetime_marginals(&gtr, counts, &grid)?;

    for (node, exact, treetime) in [("X", exact_x, treetime_x), ("R", exact_r, treetime_r)] {
      let difference = helpers::max_cdf_difference(&grid, &exact, &treetime);
      assert!(
        difference < MAX_CDF_DIFFERENCE,
        "node {node}: max CDF difference {difference:.4} to the exact posterior"
      );
    }
    Ok(())
  }

  mod helpers {
    use super::*;

    /// Transition probabilities P(same) and P(a specific other state) of JC69 over branch length b, from the
    /// eigendecomposition the branch likelihood uses: `P(c | a, b) = sum_j v[a, j] v_inv[j, c] exp(lambda_j b)`
    fn transition_probabilities(gtr: &GTR, b: f64) -> (f64, f64) {
      let eigvals = gtr.rate_eigvals();
      let p = |a: usize, c: usize| {
        (0..eigvals.len())
          .map(|j| gtr.v()[[a, j]] * gtr.v_inv()[[j, c]] * (eigvals[j] * b).exp())
          .sum::<f64>()
      };
      (p(0, 0), p(0, 1))
    }

    /// Likelihood of a branch of duration tau with `k` of SITES sites differing; zero for negative durations
    fn branch_likelihood(gtr: &GTR, tau: f64, k: usize) -> f64 {
      if tau < 0.0 || (tau == 0.0 && k > 0) {
        return 0.0;
      }
      let (same, other) = transition_probabilities(gtr, RATE / SITES.to_f64_exact() * tau);
      let differing = if k == 0 { 0.0 } else { k.to_f64_exact() * other.ln() };
      ((SITES - k).to_f64_exact() * same.ln() + differing).exp()
    }

    /// Marginal densities of X and R on `grid`; the X branch likelihood is tabulated on multiples of the grid step
    pub(super) fn exact_marginals(gtr: &GTR, [ka, kb, kx, kc]: [usize; 4], grid: &[f64]) -> (Vec<f64>, Vec<f64>) {
      let n = grid.len();
      let branch_x: Vec<f64> = (0..n)
        .map(|i| branch_likelihood(gtr, i.to_f64_exact() * GRID_STEP, kx))
        .collect();
      let tips_x: Vec<f64> = grid
        .iter()
        .map(|&x| branch_likelihood(gtr, T_A - x, ka) * branch_likelihood(gtr, T_B - x, kb))
        .collect();
      let tip_r: Vec<f64> = grid.iter().map(|&r| branch_likelihood(gtr, T_C - r, kc)).collect();
      let mut marginal_x = vec![0.0; n];
      let mut marginal_r = vec![0.0; n];
      for i in 0..n {
        for j in 0..i {
          let weight = tips_x[i] * branch_x[i - j] * tip_r[j];
          marginal_x[i] += weight;
          marginal_r[j] += weight;
        }
      }
      (marginal_x, marginal_r)
    }

    /// Marginal densities of X and R from the time passes on `grid`, with branch distributions from the production
    /// builder
    pub(super) fn treetime_marginals(
      gtr: &GTR,
      counts: [usize; 4],
      grid: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), Report> {
      let parsed = nwk_read(b"((A:1,B:1)X:1,C:1)R;".as_slice())?;
      let names = parsed.names();
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
      let posterior = propagate_distributions_forward(
        &graph,
        &constraints,
        &names,
        &branches,
        &backward,
        MaxGridPoints::default(),
        &NoopProgress,
      )?;
      let density = |name: &str| -> Result<Vec<f64>, Report> {
        let distribution = posterior[&key(name)]
          .distribution
          .as_deref()
          .expect("internal node must have a posterior")
          .normalize()?;
        let (lower, upper) = distribution.time_bounds().expect("posterior must have bounds");
        grid
          .iter()
          .map(|&t| {
            if t < lower || t > upper {
              Ok(0.0)
            } else {
              distribution.eval(t).map(|neglog| (-neglog).exp())
            }
          })
          .collect()
      };
      Ok((density("X")?, density("R")?))
    }

    /// Branch distribution as the time passes receive it: the JC69 contribution of `k` differing sites, on the grid
    /// of `fn compute_branch_length_distribution()` centred at the maximum-likelihood length
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

    fn cdf(grid: &[f64], density: &[f64]) -> Vec<f64> {
      let mut cumulative = vec![0.0; grid.len()];
      for i in 1..grid.len() {
        cumulative[i] = cumulative[i - 1] + f64::midpoint(density[i], density[i - 1]) * (grid[i] - grid[i - 1]);
      }
      let total = cumulative[grid.len() - 1];
      cumulative.iter().map(|c| c / total).collect()
    }

    pub(super) fn max_cdf_difference(grid: &[f64], a: &[f64], b: &[f64]) -> f64 {
      cdf(grid, a)
        .iter()
        .zip(cdf(grid, b))
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
    }
  }
}

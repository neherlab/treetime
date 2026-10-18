#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::coalescent::mean_field::CoalescentPrior;
  use crate::test_utils::{find_node_key_by_name, parent_edge_key, point_coalescent_prior, unknown_branches};
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::result::BranchLikelihood;
  use crate::timetree::optimization::marginal_clock_rate::branches_at_rate;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::{BoundaryBehavior, Distribution, DistributionFunction, NegLog};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::MaxGridPoints;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::cast::f64_to_usize_saturating;

  // ((A,B)I,C)root with tips at 2015, 2014, 2013 and Poisson branch likelihoods (b L)^k exp(-b L), L = 1000. The log
  // marginal likelihood from the backward pass is compared with a brute-force integral of the same integrand over
  // (t_I, t_root), for the same branch distributions and coalescent factors.

  #[rstest]
  #[case::no_prior_reference_rate(1e-3, None)]
  #[case::no_prior_slower(0.5e-3, None)]
  #[case::coalescent_reference_rate(1e-3, Some(5.0))]
  #[case::coalescent_slower(0.5e-3, Some(5.0))]
  #[case::coalescent_faster_long_tc(2e-3, Some(20.0))]
  #[trace]
  fn test_marginal_likelihood_matches_brute_force(#[case] rate: f64, #[case] tc: Option<f64>) -> Result<(), Report> {
    let setup = helpers::setup()?;
    // branch distributions at `rate`, rescaled from those at the reference rate
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
    let actual = backward
      .log_marginal_likelihood(setup.graph.root_key()?)?
      .expect("root has a distribution");
    let expected = helpers::brute_force_log_z(&setup, branches, model.as_ref())?;
    // without a prior the difference is below 5e-3 (mass windows, grids); the coalescent prior adds the grid resolution
    // of its factors
    assert_abs_diff_eq!(expected, actual, epsilon = 2e-2);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) const REFERENCE_RATE: f64 = 1e-3;
    const SITES: f64 = 1000.0;
    const TIPS: [(&str, f64); 3] = [("A", 2015.0), ("B", 2014.0), ("C", 2013.0)];
    /// substitutions on the branch above each node
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

    /// Peak-1 Poisson likelihood of k substitutions as a function of branch time at the reference rate
    fn poisson_branch(k: f64) -> Result<BranchLikelihood, Report> {
      let mean_time = k / (SITES * REFERENCE_RATE);
      let t = Array1::linspace(1e-4, mean_time * 6.0 + 4.0, 2001);
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

    /// Mean-field prior whose lineage expectation has the root at 2005 and I at 2011
    pub(super) fn coalescent(setup: &Setup, tc: f64) -> Result<CoalescentPrior, Report> {
      let times = [("root", 2005.0), ("I", 2011.0)];
      point_coalescent_prior(
        &setup.graph,
        &setup.names,
        &[TIPS.as_slice(), times.as_slice()].concat(),
        tc,
      )
    }

    /// log of the integral over t_I < 2014 and t_root < min(t_I, 2013) of the product of the branch likelihoods and the
    /// coalescent factors, by the trapezoidal rule
    pub(super) fn brute_force_log_z(
      setup: &Setup,
      branches: &BTreeMap<GraphEdgeKey, BranchLikelihood>,
      model: Option<&CoalescentPrior>,
    ) -> Result<f64, Report> {
      let neglog_branch = |name: &str, dt: f64| -> Result<f64, Report> {
        let branch = &branches[&parent_edge_key(&setup.graph, setup.keys[name])];
        branch.distribution.as_deref().expect("branch distribution").eval(dt)
      };
      let leaf_cost: f64 = match model {
        Some(m) => TIPS
          .iter()
          .map(|&(name, date)| m.node_contribution(setup.keys[name], 0, date))
          .sum::<Result<f64, Report>>()?,
        None => 0.0,
      };
      let n = 1601;
      let t_lo = 1950.0;
      let t_i = Array1::linspace(t_lo, 2014.0, n);
      let mut log_terms = Vec::with_capacity(n * n);
      let dt = t_i[1] - t_i[0];
      for (i, &ti) in t_i.iter().enumerate() {
        let w_i = if i == 0 || i == n - 1 { 0.5 } else { 1.0 };
        let cost_i = neglog_branch("A", 2015.0 - ti)?
          + neglog_branch("B", 2014.0 - ti)?
          + model.map_or(Ok(0.0), |m| m.node_contribution(setup.keys["I"], 2, ti))?;
        let r_hi = ti.min(2013.0);
        let n_r = f64_to_usize_saturating(((r_hi - t_lo) / dt).floor()).max(1) + 1;
        let t_r = Array1::linspace(t_lo, r_hi, n_r);
        let dr = if n_r > 1 { t_r[1] - t_r[0] } else { 0.0 };
        for (j, &tr) in t_r.iter().enumerate() {
          let w_r = if j == 0 || j == n_r - 1 { 0.5 } else { 1.0 };
          let cost = cost_i
            + neglog_branch("I", ti - tr)?
            + neglog_branch("C", 2013.0 - tr)?
            + model.map_or(Ok(0.0), |m| m.node_contribution(setup.keys["root"], 2, tr))?
            + leaf_cost;
          log_terms.push((w_i * w_r * dt * dr).ln() - cost);
        }
      }
      let max = log_terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
      Ok(max + log_terms.iter().map(|v| (v - max).exp()).sum::<f64>().ln())
    }
  }
}

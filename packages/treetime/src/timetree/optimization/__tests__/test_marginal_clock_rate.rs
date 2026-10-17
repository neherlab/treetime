#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::test_utils::{find_node_key_by_name, parent_edge_key, unknown_branches};
  use crate::timetree::inference::result::BranchLikelihood;
  use crate::timetree::optimization::marginal_clock_rate::{
    MarginalRateInputs, marginal_log_lh_at_rate, optimize_marginal_clock_rate,
  };
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::MaxGridPoints;
  use treetime_io::nwk::nwk_read;

  // Two leaves dated t_A = 2015, t_B = 2014 below the root; Gaussian branch likelihoods of precision 4 in branch length
  // with means m_A = 5, m_B = 2. log Z(mu) = -ln(mu) - (mu (t_A - t_B) - (m_A - m_B))^2 + const = -ln(mu) - (mu - 3)^2
  // + const, maximal where 2 mu^2 - 6 mu + 1 = 0.

  #[test]
  fn test_marginal_clock_rate_matches_closed_form_curve() -> Result<(), Report> {
    let setup = helpers::cherry()?;
    let inputs = setup.inputs();
    let at = |mu: f64| marginal_log_lh_at_rate(&inputs, mu).map(|log_lh| log_lh.expect("root has a distribution"));
    let curve = |mu: f64| -mu.ln() - (mu - 3.0).powi(2);
    // near the maximum the two messages overlap within 2 sd; far from it (mu = 0.5: 5 sd apart) their product lies in
    // the messages' outer grid and fitted tails, and the curve is off by about 0.02
    for mu in [2.0, 2.5, 4.0] {
      // differences to mu = 3 drop the constant
      assert_abs_diff_eq!(curve(mu) - curve(3.0), at(mu)? - at(3.0)?, epsilon = 2e-3);
    }
    Ok(())
  }

  #[test]
  fn test_marginal_clock_rate_finds_the_maximum() -> Result<(), Report> {
    let setup = helpers::cherry()?;
    let fit = optimize_marginal_clock_rate(&setup.inputs())?;
    // the search stops at a relative tolerance of 1e-3 on the square root of the rate
    assert_abs_diff_eq!((6.0 + 28.0_f64.sqrt()) / 4.0, fit.clock_rate, epsilon = 6e-3);
    Ok(())
  }

  mod helpers {
    use super::*;

    const PRECISION: f64 = 4.0;

    pub(super) struct Cherry {
      graph: Graph,
      constraints: DateConstraints,
      bad_branches: BTreeMap<GraphNodeKey, bool>,
      branches: BTreeMap<treetime_graph::edge::GraphEdgeKey, BranchLikelihood>,
    }

    impl Cherry {
      pub(super) fn inputs(&self) -> MarginalRateInputs<'_> {
        MarginalRateInputs {
          graph: &self.graph,
          date_constraints: &self.constraints,
          coalescent: None,
          bad_branches: &self.bad_branches,
          branches: &self.branches,
          reference_rate: 1.0,
          max_points: MaxGridPoints::default(),
        }
      }
    }

    /// Branch likelihoods are given at reference rate 1, where time equals branch length
    pub(super) fn cherry() -> Result<Cherry, Report> {
      let parsed = nwk_read(b"(A:5.0,B:2.0)root;".as_slice())?;
      let names = parsed.names();
      let graph = parsed.graph;
      let mut constraints = DateConstraints::default();
      let mut branches = unknown_branches(&graph);
      for (name, date, mean) in [("A", 2015.0, 5.0_f64), ("B", 2014.0, 2.0)] {
        let key = find_node_key_by_name(&graph, &names, name).expect("leaf not found");
        constraints
          .by_node
          .insert(key, Some(Arc::new(Distribution::<NegLog>::point(date, 0.0))));
        let t = Array1::linspace(0.0, 2.0 * mean, 401);
        let y = t.mapv(|t| 0.5 * PRECISION * (t - mean).powi(2));
        let branch = BranchLikelihood {
          distribution: Some(Arc::new(Distribution::function(t, y)?)),
          time_length: Some(mean),
        };
        branches.insert(parent_edge_key(&graph, key), branch);
      }
      let bad_branches = graph.get_nodes().map(|node| (node.key(), false)).collect();
      Ok(Cherry {
        graph,
        constraints,
        bad_branches,
        branches,
      })
    }
  }
}

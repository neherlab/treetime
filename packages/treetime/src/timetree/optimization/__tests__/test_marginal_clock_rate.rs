#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::test_utils::{find_node_key_by_name, parent_edge_key, unknown_branches};
  use crate::timetree::inference::result::BranchLikelihood;
  use crate::timetree::optimization::marginal_clock_rate::{
    MarginalRateInputs, branches_at_rate, marginal_log_lh_at_rate, step_marginal_clock_rate,
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
  fn test_marginal_clock_rate_steps_half_way_to_the_parabola_maximum() -> Result<(), Report> {
    let setup = helpers::cherry()?;
    let reference = 2.5;
    let branches = branches_at_rate(&setup.branches, 1.0, reference)?;
    let step = step_marginal_clock_rate(&setup.inputs_at(&branches, reference))?;

    // parabola in x = ln(mu / 2.5) through the closed-form curve at x = -h, 0, h
    let h = 1.25_f64.ln();
    let curve = |x: f64| {
      let mu = reference * x.exp();
      -mu.ln() - (mu - 3.0).powi(2)
    };
    let (below, at, above) = (curve(-h), curve(0.0), curve(h));
    let target = h * (below - above) / (2.0 * (above + below - 2.0 * at));
    assert_abs_diff_eq!(reference * target.exp(), step.target_rate, epsilon = 1e-2);
    assert_abs_diff_eq!(reference * (0.5 * target).exp(), step.clock_rate, epsilon = 1e-2);
    Ok(())
  }

  #[test]
  fn test_marginal_clock_rate_without_a_maximum_steps_half_the_window_toward_the_better_probe() -> Result<(), Report> {
    // at mu = 1 the curve -ln(mu) - (mu - 3)^2 is convex in ln(mu) over the probes and rises toward faster rates
    let setup = helpers::cherry()?;
    let step = step_marginal_clock_rate(&setup.inputs())?;
    assert_abs_diff_eq!(2.0, step.target_rate, epsilon = 1e-12);
    assert_abs_diff_eq!(2.0_f64.sqrt(), step.clock_rate, epsilon = 1e-12);
    Ok(())
  }

  mod helpers {
    use super::*;

    const PRECISION: f64 = 4.0;

    pub(super) struct Cherry {
      pub graph: Graph,
      pub constraints: DateConstraints,
      pub bad_branches: BTreeMap<GraphNodeKey, bool>,
      pub branches: BTreeMap<treetime_graph::edge::GraphEdgeKey, BranchLikelihood>,
    }

    impl Cherry {
      pub(super) fn inputs(&self) -> MarginalRateInputs<'_> {
        self.inputs_at(&self.branches, 1.0)
      }

      pub(super) fn inputs_at<'a>(
        &'a self,
        branches: &'a BTreeMap<treetime_graph::edge::GraphEdgeKey, BranchLikelihood>,
        reference_rate: f64,
      ) -> MarginalRateInputs<'a> {
        MarginalRateInputs {
          graph: &self.graph,
          date_constraints: &self.constraints,
          coalescent: None,
          bad_branches: &self.bad_branches,
          branches,
          reference_rate,
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

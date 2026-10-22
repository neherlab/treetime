#[cfg(test)]
mod tests {
  use crate::test_utils::arb_simulated_case;
  use helpers::compare_dense_sparse;
  use proptest::prelude::{ProptestConfig, prop_assert, proptest};

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    #[ignore = "kb/issues/H-ancestral-sparse-dense-internal-residues-diverge.md"]
    fn test_prop_marginal_dense_sparse_simulated_agree_on_states(spec in arb_simulated_case()) {
      let comparison = compare_dense_sparse(&spec.simulate().expect("a case of the strategy simulates"))
        .expect("dense and sparse reconstructions run");

      prop_assert!(
        comparison.differing_states.is_empty(),
        "most likely states differ at {} of {} compared node sites (dense posterior gap above {}); {} excluded near-tie sites also differ: {:?}",
        comparison.differing_states.len(),
        comparison.compared_sites,
        helpers::MIN_POSTERIOR_GAP,
        comparison.differing_near_ties,
        comparison.differing_states,
      );
    }

    #[test]
    #[ignore = "kb/issues/M-ancestral-dense-sparse-divergence.md"]
    fn test_prop_marginal_dense_sparse_simulated_agree_on_likelihood(spec in arb_simulated_case()) {
      let comparison = compare_dense_sparse(&spec.simulate().expect("a case of the strategy simulates"))
        .expect("dense and sparse reconstructions run");

      prop_assert!(
        approx::relative_eq!(comparison.log_lh_dense, comparison.log_lh_sparse, max_relative = 1e-5, epsilon = 0.0),
        "log-likelihoods differ: dense {}, sparse {}; {} excluded near-tie sites differ in their most likely state",
        comparison.log_lh_dense,
        comparison.log_lh_sparse,
        comparison.differing_near_ties,
      );
    }
  }

  mod helpers {
    use crate::partition::marginal::posterior_states::most_likely_state;
    use crate::test_utils::{SimulatedCase, gtr_from_model, run_dense_marginal, run_sparse_marginal};
    use eyre::Report;
    use itertools::Itertools;
    use ndarray::{ArrayView1, aview1};
    use ordered_float::OrderedFloat;
    use std::cmp::Reverse;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_simulation::branch_lengths::exact;
    use treetime_utils::collections::container::get_exactly_one;
    use treetime_utils::{make_error, make_report};

    pub(super) const MIN_POSTERIOR_GAP: f64 = 1e-3;

    pub(super) fn compare_dense_sparse(case: &SimulatedCase) -> Result<Comparison, Report> {
      let graph = &case.tree.graph;
      let names = &case.tree.names;
      let lengths: BTreeMap<GraphEdgeKey, Option<f64>> = exact(&case.tree.edge_durations(), &case.model)?
        .into_iter()
        .map(|(edge, length)| (edge, Some(length)))
        .collect();
      let records = case.sequences.leaf_records(graph, names);
      let sequence_lengths = records.iter().map(|record| record.seq.len()).unique().collect_vec();
      let length = *get_exactly_one(&sequence_lengths)?;
      let states = case.model.pi().len();
      let gtr = gtr_from_model(&case.model, 1.0)?;
      let (log_lh_dense, dense) = run_dense_marginal(graph, &lengths, names, &records, gtr.clone())?;
      let (log_lh_sparse, sparse) = run_sparse_marginal(graph, &lengths, names, &records, gtr)?;

      let mut comparison = Comparison {
        log_lh_dense,
        log_lh_sparse,
        compared_sites: 0,
        differing_near_ties: 0,
        differing_states: vec![],
      };
      for node in graph.get_internal_nodes().map(|node| node.key()) {
        let name = &names[node];
        let sparse_state = &sparse.node_states[&node];
        let dense_profile = &dense.node_states[&node].profile.dis;
        if dense_profile.nrows() != length || sparse_state.sequence.len() != length {
          return make_error!(
            "Node {name} has a dense profile of {} sites and a sparse sequence of {} sites, but the alignment has {length}",
            dense_profile.nrows(),
            sparse_state.sequence.len()
          );
        }
        for (pos, dense_posterior) in dense_profile.rows().into_iter().enumerate() {
          let sparse_posterior = if let Some(var) = sparse_state.profile.variable.get(pos) {
            aview1(var.dis)
          } else {
            let fixed_state = sparse_state.sequence[pos];
            aview1(
              sparse_state
                .profile
                .fixed
                .get(fixed_state)
                .ok_or_else(|| {
                  make_report!("Node {name} has no fixed profile of state {fixed_state} at position {pos}")
                })?
                .dis
                .as_slice(),
            )
          };
          check_posterior("dense", &dense_posterior, states, name, pos)?;
          check_posterior("sparse", &sparse_posterior, states, name, pos)?;
          let dense_state = most_likely(&dense_posterior, name, pos)?;
          let sparse_state = most_likely(&sparse_posterior, name, pos)?;
          let gap = top_two_gap(&dense_posterior);
          match (gap > MIN_POSTERIOR_GAP, dense_state == sparse_state) {
            (true, true) => comparison.compared_sites += 1,
            (true, false) => {
              comparison.compared_sites += 1;
              comparison.differing_states.push(format!(
                "node {name} position {pos}: dense posterior {dense_posterior}, sparse posterior {sparse_posterior}"
              ));
            },
            (false, true) => {},
            (false, false) => comparison.differing_near_ties += 1,
          }
        }
      }
      Ok(comparison)
    }

    #[derive(Debug)]
    pub(super) struct Comparison {
      pub log_lh_dense: f64,
      pub log_lh_sparse: f64,
      pub compared_sites: usize,
      pub differing_near_ties: usize,
      pub differing_states: Vec<String>,
    }

    fn check_posterior(
      kind: &str,
      posterior: &ArrayView1<'_, f64>,
      states: usize,
      node: &str,
      pos: usize,
    ) -> Result<(), Report> {
      if posterior.len() != states || !posterior.iter().all(|p| p.is_finite()) {
        return make_error!(
          "The {kind} posterior of node {node} at position {pos} is {posterior}, but it needs {states} finite values"
        );
      }
      Ok(())
    }

    fn most_likely(posterior: &ArrayView1<'_, f64>, node: &str, pos: usize) -> Result<usize, Report> {
      most_likely_state(posterior).ok_or_else(|| make_report!("Posterior of node {node} at position {pos} is empty"))
    }

    fn top_two_gap(posterior: &ArrayView1<'_, f64>) -> f64 {
      let sorted = posterior
        .iter()
        .copied()
        .sorted_by_key(|&p| Reverse(OrderedFloat(p)))
        .collect_vec();
      if let [first, second, ..] = sorted.as_slice() {
        first - second
      } else {
        f64::INFINITY
      }
    }
  }
}

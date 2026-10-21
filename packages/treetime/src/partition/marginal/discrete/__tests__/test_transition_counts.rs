#[cfg(test)]
mod tests {
  use eyre::Report;
  use helpers::{NODES, OracleCase, enumerated_oracle, one_hot_root_state, run_discrete_partition};
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::two_states(  (vec![0.3, 0.7],                  [Some(0), Some(1), Some(0)], [0.1, 0.2, 0.3, 0.4]))]
  #[case::three_states((vec![0.2, 0.3, 0.5],             [Some(0), Some(1), Some(2)], [0.1, 0.2, 0.3, 0.4]))]
  #[case::five_states( (vec![0.1, 0.15, 0.2, 0.25, 0.3], [Some(3), Some(4), None],    [0.1, 0.2, 0.3, 0.4]))]
  #[trace]
  fn test_discrete_transition_counts_match_enumerated_f81_oracle(
    #[case] (pi, leaves, lengths): (Vec<f64>, [Option<usize>; 3], [f64; 4]),
    #[values(false, true)] filter_uninformative_root: bool,
  ) -> Result<(), Report> {
    let case = OracleCase {
      pi: &pi,
      leaves,
      lengths,
    };
    let expected = enumerated_oracle(&case);
    let expected_root_state = one_hot_root_state(&expected.posteriors["root"]);

    let (counts, posteriors) = run_discrete_partition(&case, filter_uninformative_root)?;

    pretty_assert_abs_diff_eq!(expected.nij, counts.nij, epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(expected.dwell_times, counts.Ti, epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(expected_root_state, counts.root_state, epsilon = 1e-14);
    for name in NODES {
      pretty_assert_abs_diff_eq!(expected.posteriors[name], posteriors[name], epsilon = 1e-14, "node {name}");
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::filter_off(false, array![1.0, 0.0])]
  #[case::filter_on( true,  array![0.0, 0.0])]
  #[trace]
  fn test_discrete_transition_counts_near_uniform_root_follows_root_filter(
    #[case] filter_uninformative_root: bool,
    #[case] expected_root_state: Array1<f64>,
  ) -> Result<(), Report> {
    let (t_a, t_b, t_inner) = (0.5, 0.500_000_000_06, 0.2);
    let case = OracleCase {
      pi: &[0.5, 0.5],
      leaves: [Some(0), Some(1), None],
      lengths: [t_a, t_b, t_inner, 0.3],
    };
    let expected = enumerated_oracle(&case);
    let (decay_a, decay_b): (f64, f64) = ((-t_a).exp(), (-t_b).exp());
    let root_excess_over_uniform_by_f81_pruning =
      (-t_inner).exp() * (decay_a - decay_b) / (2.0 * (1.0 - decay_a * decay_b));
    let expected_root_posterior = array![
      0.5 + root_excess_over_uniform_by_f81_pruning,
      0.5 - root_excess_over_uniform_by_f81_pruning
    ];

    let (counts, posteriors) = run_discrete_partition(&case, filter_uninformative_root)?;

    pretty_assert_abs_diff_eq!(expected_root_posterior, posteriors["root"], epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(expected.nij, counts.nij, epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(expected.dwell_times, counts.Ti, epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(expected_root_state, counts.root_state, epsilon = 1e-14);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::state_alphabet::StateAlphabet;
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
    use crate::gtr::gtr::{GTR, Mu};
    use crate::gtr::infer_gtr::common::MutationCounts;
    use crate::partition::marginal::discrete::partition::{PartitionMarginalDiscrete, StateObservation};
    use crate::partition::marginal::posterior_states::most_likely_state;
    use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate, TransitionCounts};
    use crate::test_utils::find_node_key_by_name;
    use eyre::Report;
    use itertools::{iproduct, izip};
    use ndarray::{Array1, Array2, Array3, Axis, aview1};
    use std::collections::BTreeMap;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read;
    use treetime_utils::cast::ToF64Exact;

    pub(super) const NODES: [&str; 5] = ["root", "inner", "A", "B", "C"];
    const LEAVES: [&str; 3] = ["A", "B", "C"];

    pub(super) struct OracleCase<'a> {
      pub(super) pi: &'a [f64],
      pub(super) leaves: [Option<usize>; 3],
      pub(super) lengths: [f64; 4],
    }

    pub(super) struct Oracle {
      pub(super) nij: Array2<f64>,
      pub(super) dwell_times: Array1<f64>,
      pub(super) posteriors: BTreeMap<&'static str, Array1<f64>>,
    }

    pub(super) fn enumerated_oracle(case: &OracleCase<'_>) -> Oracle {
      let n = case.pi.len();
      let [t_a, t_b, t_inner, t_c] = case.lengths;
      let [leaf_a, leaf_b, leaf_c] = case.leaves;
      let edges = [(0, 1, t_inner), (0, 4, t_c), (1, 2, t_a), (1, 3, t_b)];
      let transition = |child: usize, parent: usize, t: f64| f81_transition_probability(case.pi, child, parent, t);

      let mut node_marginals = Array2::<f64>::zeros((NODES.len(), n));
      let mut edge_child_parent_marginals = Array3::<f64>::zeros((edges.len(), n, n));
      let mut total = 0.0;
      for (root, inner, a, b, c) in iproduct!(0..n, 0..n, 0..n, 0..n, 0..n) {
        let states = [root, inner, a, b, c];
        let joint = case.pi[root]
          * transition(inner, root, t_inner)
          * transition(c, root, t_c)
          * transition(a, inner, t_a)
          * transition(b, inner, t_b)
          * observation(leaf_a, a)
          * observation(leaf_b, b)
          * observation(leaf_c, c);
        total += joint;
        for (node, state) in states.iter().enumerate() {
          node_marginals[[node, *state]] += joint;
        }
        for (edge, (parent, child, _)) in edges.iter().enumerate() {
          edge_child_parent_marginals[[edge, states[*child], states[*parent]]] += joint;
        }
      }
      node_marginals /= total;
      edge_child_parent_marginals /= total;

      let mut nij = edge_child_parent_marginals.sum_axis(Axis(0));
      nij.diag_mut().fill(0.0);
      let dwell_times = edges.iter().fold(Array1::zeros(n), |dwell, &(parent, child, t)| {
        dwell + 0.5 * t * (&node_marginals.row(parent) + &node_marginals.row(child))
      });
      let posteriors = izip!(NODES, node_marginals.outer_iter())
        .map(|(name, row)| (name, row.to_owned()))
        .collect();
      Oracle {
        nij,
        dwell_times,
        posteriors,
      }
    }

    pub(super) fn one_hot_root_state(root_posterior: &Array1<f64>) -> Array1<f64> {
      let n = root_posterior.len();
      let most_likely = most_likely_state(root_posterior).expect("root posterior is not empty");
      assert!(
        root_posterior[most_likely] > 1.0 / n.to_f64_exact() + 1e-10,
        "case must have an informative root posterior: {root_posterior}"
      );
      Array1::from_shape_fn(n, |state| if state == most_likely { 1.0 } else { 0.0 })
    }

    pub(super) fn run_discrete_partition(
      case: &OracleCase<'_>,
      filter_uninformative_root: bool,
    ) -> Result<(MutationCounts, BTreeMap<&'static str, Array1<f64>>), Report> {
      let n = case.pi.len();
      let [t_a, t_b, t_inner, t_c] = case.lengths;
      let parsed = nwk_read(format!("((A:{t_a},B:{t_b})inner:{t_inner},C:{t_c})root;").as_bytes())?;
      let names = parsed.names.clone();
      let graph = parsed.graph;
      let key = |name: &str| -> GraphNodeKey {
        find_node_key_by_name(&graph, &names, name).unwrap_or_else(|| panic!("missing node '{name}'"))
      };

      let states = StateAlphabet::new((0..n).map(|state| format!("state {state}")))?;
      let observations = izip!(LEAVES, case.leaves)
        .map(|(name, leaf)| {
          (
            key(name),
            leaf.map_or(StateObservation::Unknown, StateObservation::Known),
          )
        })
        .collect();
      let partition = PartitionMarginalDiscrete::new(
        states,
        &graph,
        &observations,
        MIN_BRANCH_LENGTH_FRACTION,
        filter_uninformative_root,
      )?;
      let gtr = GTR::builder()
        .n_states(n)
        .mu(Mu::FromMatrix(1.0))
        .pi(aview1(case.pi).to_owned())
        .build()?;
      let branch_lengths = branch_lengths_or_zero(&parsed.branch_lengths);

      let MarginalUpdate { node_states, edges, .. } = partition.marginal_update(&gtr, &graph, &branch_lengths, &())?;
      let counts = partition.count_transitions(
        &gtr,
        &graph,
        &branch_lengths,
        &node_states,
        &edges.backward,
        &edges.forward,
      )?;
      let state_posteriors = partition.state_posteriors(&graph, &node_states)?;
      let posteriors = NODES
        .into_iter()
        .map(|name| (name, state_posteriors[&key(name)].profile().clone()))
        .collect();
      Ok((counts, posteriors))
    }

    fn f81_transition_probability(pi: &[f64], child: usize, parent: usize, t: f64) -> f64 {
      let stays = if child == parent { 1.0 } else { 0.0 };
      pi[child] + (stays - pi[child]) * (-t).exp()
    }

    fn observation(leaf: Option<usize>, state: usize) -> f64 {
      match leaf {
        Some(known) if known != state => 0.0,
        _ => 1.0,
      }
    }
  }
}

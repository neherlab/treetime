#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::{marginal_update, profile_branch_lengths};
  use crate::o;
  use crate::partition::marginal::shared::pass::{marginal_process_backward_indexed, marginal_process_forward_indexed};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::{NwkParse, nwk_read_str};
  use treetime_utils::assert_error;

  #[test]
  fn test_discrete_marginal_attach_traits_maps_observed_and_missing_profiles() -> Result<(), Report> {
    let NwkParse { graph, names, .. } = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let graph: Graph = graph;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("?"),
    };

    partition.attach_traits(&graph, &traits, &names)?;

    let node_a_profile = helpers::get_node_profile(&graph, &names, &partition, "A");
    assert_abs_diff_eq!(node_a_profile[0], 0.0, epsilon = 1e-10);
    assert_abs_diff_eq!(node_a_profile[1], 1.0, epsilon = 1e-10);

    let node_b_profile = helpers::get_node_profile(&graph, &names, &partition, "B");
    assert_abs_diff_eq!(node_b_profile[0], 0.5, epsilon = 1e-10);
    assert_abs_diff_eq!(node_b_profile[1], 0.5, epsilon = 1e-10);

    assert_eq!(graph.get_edges().len(), partition.data.edges.len());

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_attach_traits_rejects_tree_leaf_missing_from_metadata() -> Result<(), Report> {
    let NwkParse { graph, names, .. } = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let graph: Graph = graph;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
    };

    let result = partition.attach_traits(&graph, &traits, &names);
    assert_error!(result, "Mugration: tree leaves missing from metadata: B");

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_attach_traits_accepts_metadata_name_missing_from_tree() -> Result<(), Report> {
    let NwkParse { graph, names, .. } = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let graph: Graph = graph;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
      o!("C") => o!("usa"),
    };

    // "C" has no matching tree leaf. Attachment warns and proceeds; matched leaves are unaffected.
    partition.attach_traits(&graph, &traits, &names)?;

    let node_a_profile = helpers::get_node_profile(&graph, &names, &partition, "A");
    assert_abs_diff_eq!(node_a_profile[0], 0.0, epsilon = 1e-10);
    assert_abs_diff_eq!(node_a_profile[1], 1.0, epsilon = 1e-10);

    let node_b_profile = helpers::get_node_profile(&graph, &names, &partition, "B");
    assert_abs_diff_eq!(node_b_profile[0], 1.0, epsilon = 1e-10);
    assert_abs_diff_eq!(node_b_profile[1], 0.0, epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_passes_normalize_backward_and_forward_profiles() -> Result<(), Report> {
    let (graph, names, raw_branch_lengths) = helpers::make_fixture_graph()?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = helpers::make_fixture_traits();

    partition.attach_traits(&graph, &traits, &names)?;

    let branch_lengths = profile_branch_lengths(&raw_branch_lengths);
    marginal_process_backward_indexed(&mut partition, &graph, &branch_lengths)?;

    let root_profile = helpers::get_node_profile(&graph, &names, &partition, "root");
    helpers::assert_profile_normalized(&root_profile);

    let inner_to_root_msg = helpers::get_edge_msg_from_child(&graph, &names, &partition, "root", "inner");
    helpers::assert_profile_normalized(&inner_to_root_msg);

    let leaf_to_inner_msg = helpers::get_edge_msg_from_child(&graph, &names, &partition, "inner", "A");
    helpers::assert_profile_normalized(&leaf_to_inner_msg);

    marginal_process_forward_indexed(&mut partition, &graph, &branch_lengths)?;

    let root_profile = helpers::get_node_profile(&graph, &names, &partition, "root");
    helpers::assert_profile_normalized(&root_profile);

    let inner_profile = helpers::get_node_profile(&graph, &names, &partition, "inner");
    helpers::assert_profile_normalized(&inner_profile);

    let root_to_c_msg = helpers::get_edge_msg_to_child(&graph, &names, &partition, "root", "C");
    helpers::assert_profile_normalized(&root_to_c_msg);

    let inner_to_b_msg = helpers::get_edge_msg_to_child(&graph, &names, &partition, "inner", "B");
    helpers::assert_profile_normalized(&inner_to_b_msg);

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_run_returns_finite_log_lh_and_reconstructs_internal_trait() -> Result<(), Report> {
    let (graph, names, raw_branch_lengths) = helpers::make_fixture_graph()?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = helpers::make_fixture_traits();

    partition.attach_traits(&graph, &traits, &names)?;

    let mut partition = partition;
    let actual_log_lh = marginal_update(
      &graph,
      &profile_branch_lengths(&raw_branch_lengths),
      std::slice::from_mut(&mut partition),
    )?
    .value();

    assert!(actual_log_lh.is_finite());
    assert!(
      actual_log_lh <= 0.0,
      "Log-likelihood must be non-positive: {actual_log_lh}"
    );

    let partition = &partition;
    let inner_profile = helpers::get_node_profile(&graph, &names, partition, "inner");
    helpers::assert_profile_normalized(&inner_profile);

    let inner_key = helpers::get_node_key(&graph, &names, "inner");
    let expected_trait = Some(o!("usa"));
    let actual_trait = partition.get_reconstructed_trait(inner_key);
    assert_eq!(expected_trait, actual_trait);

    Ok(())
  }

  mod helpers {
    use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
    use crate::gtr::gtr::{GTR, GTRParams};
    use crate::o;
    use crate::partition::marginal::discrete::partition::PartitionMarginalDiscrete;
    use crate::partition::storage::discrete::DiscreteStates;
    use crate::test_utils::{find_edge_key, find_node_key_by_name};
    use eyre::Report;
    use maplit::btreemap;
    use ndarray::Array1;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::{NwkParse, nwk_read_str};
    use treetime_utils::pretty_assert_abs_diff_eq;

    pub(super) fn make_partition(states: [&str; 2]) -> Result<PartitionMarginalDiscrete, Report> {
      let discrete_states = DiscreteStates::from_values(states.into_iter(), "?");
      let n_states = discrete_states.len();
      let gtr = GTR::new(GTRParams {
        n_states,
        mu: 1.0,
        W: None,
        pi: Array1::from_elem(n_states, 1.0 / n_states as f64),
      })?;

      Ok(PartitionMarginalDiscrete::new(
        gtr,
        discrete_states,
        MIN_BRANCH_LENGTH_FRACTION,
        false,
      ))
    }

    pub(super) fn make_fixture_graph() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let NwkParse {
        graph,
        names,
        branch_lengths,
        ..
      } = nwk_read_str("((A:0.01,B:0.01)inner:0.01,C:0.25)root;")?;
      Ok((graph, names, branch_lengths))
    }

    pub(super) fn make_fixture_traits() -> BTreeMap<String, String> {
      btreemap! {
        o!("A") => o!("usa"),
        o!("B") => o!("usa"),
        o!("C") => o!("germany"),
      }
    }

    pub(super) fn assert_profile_normalized(profile: &Array1<f64>) {
      assert!(profile.iter().all(|value| value.is_finite()));
      let expected_sum = 1.0;
      let actual_sum = profile.sum();
      pretty_assert_abs_diff_eq!(expected_sum, actual_sum, epsilon = 1e-12);
    }

    pub(super) fn get_node_key(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      name: &str,
    ) -> GraphNodeKey {
      find_node_key_by_name(graph, names, name).unwrap_or_else(|| panic!("Missing test node '{name}'"))
    }

    pub(super) fn get_node_profile(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      partition: &PartitionMarginalDiscrete,
      name: &str,
    ) -> Array1<f64> {
      let node_key = get_node_key(graph, names, name);
      let node = partition
        .data
        .nodes
        .get(&node_key)
        .unwrap_or_else(|| panic!("Missing discrete node data for '{name}'"));
      node.profile.dis.row(0).to_owned()
    }

    pub(super) fn get_edge_msg_from_child(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      partition: &PartitionMarginalDiscrete,
      source_name: &str,
      target_name: &str,
    ) -> Array1<f64> {
      let edge_key = find_edge_key(graph, names, source_name, target_name)
        .unwrap_or_else(|| panic!("Missing test edge '{source_name}' -> '{target_name}'"));
      let edge = partition
        .data
        .edges
        .get(&edge_key)
        .unwrap_or_else(|| panic!("Missing discrete edge data for '{source_name}' -> '{target_name}'"));
      edge.msg_from_child.dis.row(0).to_owned()
    }

    pub(super) fn get_edge_msg_to_child(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      partition: &PartitionMarginalDiscrete,
      source_name: &str,
      target_name: &str,
    ) -> Array1<f64> {
      let edge_key = find_edge_key(graph, names, source_name, target_name)
        .unwrap_or_else(|| panic!("Missing test edge '{source_name}' -> '{target_name}'"));
      let edge = partition
        .data
        .edges
        .get(&edge_key)
        .unwrap_or_else(|| panic!("Missing discrete edge data for '{source_name}' -> '{target_name}'"));
      edge.msg_to_child.dis.row(0).to_owned()
    }
  }
}

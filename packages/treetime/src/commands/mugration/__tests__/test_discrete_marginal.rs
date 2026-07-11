#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::update_marginal;
  use crate::o;
  use crate::partition::traits::PartitionMarginalPasses;
  use crate::payload::ancestral::GraphAncestral;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;

  #[test]
  fn test_discrete_marginal_attach_traits_maps_observed_and_missing_profiles() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("?"),
    };

    partition.attach_traits(&graph, &traits)?;

    let node_a_profile = helpers::get_node_profile(&graph, &partition, "A");
    assert_abs_diff_eq!(node_a_profile[0], 0.0, epsilon = 1e-10);
    assert_abs_diff_eq!(node_a_profile[1], 1.0, epsilon = 1e-10);

    let node_b_profile = helpers::get_node_profile(&graph, &partition, "B");
    assert_abs_diff_eq!(node_b_profile[0], 0.5, epsilon = 1e-10);
    assert_abs_diff_eq!(node_b_profile[1], 0.5, epsilon = 1e-10);

    assert_eq!(graph.get_edges().len(), partition.data.edges.len());

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_attach_traits_rejects_tree_leaf_missing_from_metadata() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
    };

    let result = partition.attach_traits(&graph, &traits);
    assert_error!(result, "Mugration: tree leaves missing from metadata: B");

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_attach_traits_rejects_metadata_name_missing_from_tree() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
      o!("C") => o!("usa"),
    };

    let result = partition.attach_traits(&graph, &traits);
    assert_error!(result, "Mugration: metadata names missing from tree leaves: C");

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_passes_normalize_backward_and_forward_profiles() -> Result<(), Report> {
    let graph: GraphAncestral = helpers::make_fixture_graph()?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = helpers::make_fixture_traits();

    partition.attach_traits(&graph, &traits)?;

    partition.process_backward_pass(&graph)?;

    let root_profile = helpers::get_node_profile(&graph, &partition, "root");
    helpers::assert_profile_normalized(&root_profile);

    let inner_to_root_msg = helpers::get_edge_msg_from_child(&graph, &partition, "root", "inner");
    helpers::assert_profile_normalized(&inner_to_root_msg);

    let leaf_to_inner_msg = helpers::get_edge_msg_from_child(&graph, &partition, "inner", "A");
    helpers::assert_profile_normalized(&leaf_to_inner_msg);

    partition.process_forward_pass(&graph)?;

    let root_profile = helpers::get_node_profile(&graph, &partition, "root");
    helpers::assert_profile_normalized(&root_profile);

    let inner_profile = helpers::get_node_profile(&graph, &partition, "inner");
    helpers::assert_profile_normalized(&inner_profile);

    let root_to_c_msg = helpers::get_edge_msg_to_child(&graph, &partition, "root", "C");
    helpers::assert_profile_normalized(&root_to_c_msg);

    let inner_to_b_msg = helpers::get_edge_msg_to_child(&graph, &partition, "inner", "B");
    helpers::assert_profile_normalized(&inner_to_b_msg);

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_run_returns_finite_log_lh_and_reconstructs_internal_trait() -> Result<(), Report> {
    let graph: GraphAncestral = helpers::make_fixture_graph()?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = helpers::make_fixture_traits();

    partition.attach_traits(&graph, &traits)?;

    let partition = Arc::new(RwLock::new(partition));
    let actual_log_lh = update_marginal(&graph, std::slice::from_ref(&partition))?;

    assert!(actual_log_lh.is_finite());
    assert!(
      actual_log_lh <= 0.0,
      "Log-likelihood must be non-positive: {actual_log_lh}"
    );

    let partition = partition.read_arc();
    let inner_profile = helpers::get_node_profile(&graph, &partition, "inner");
    helpers::assert_profile_normalized(&inner_profile);

    let inner_key = helpers::get_node_key(&graph, "inner");
    let expected_trait = Some(o!("usa"));
    let actual_trait = partition.get_reconstructed_trait(inner_key);
    assert_eq!(expected_trait, actual_trait);

    Ok(())
  }

  mod helpers {
    use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
    use crate::gtr::gtr::{GTR, GTRParams};
    use crate::o;
    use crate::partition::discrete_states::DiscreteStates;
    use crate::partition::marginal_discrete::PartitionMarginalDiscrete;
    use crate::payload::ancestral::GraphAncestral;
    use crate::test_utils::{find_edge_key, find_node_key_by_name};
    use eyre::Report;
    use maplit::btreemap;
    use ndarray::Array1;
    use std::collections::BTreeMap;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read_str;
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

    pub(super) fn make_fixture_graph() -> Result<GraphAncestral, Report> {
      nwk_read_str("((A:0.01,B:0.01)inner:0.01,C:0.25)root;")
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

    pub(super) fn get_node_key(graph: &GraphAncestral, name: &str) -> GraphNodeKey {
      find_node_key_by_name(graph, name).unwrap_or_else(|| panic!("Missing test node '{name}'"))
    }

    pub(super) fn get_node_profile(
      graph: &GraphAncestral,
      partition: &PartitionMarginalDiscrete,
      name: &str,
    ) -> Array1<f64> {
      let node_key = get_node_key(graph, name);
      let node = partition
        .data
        .nodes
        .get(&node_key)
        .unwrap_or_else(|| panic!("Missing discrete node data for '{name}'"));
      node.profile.dis.row(0).to_owned()
    }

    pub(super) fn get_edge_msg_from_child(
      graph: &GraphAncestral,
      partition: &PartitionMarginalDiscrete,
      source_name: &str,
      target_name: &str,
    ) -> Array1<f64> {
      let edge_key = find_edge_key(graph, source_name, target_name)
        .unwrap_or_else(|| panic!("Missing test edge '{source_name}' -> '{target_name}'"));
      let edge = partition
        .data
        .edges
        .get(&edge_key)
        .unwrap_or_else(|| panic!("Missing discrete edge data for '{source_name}' -> '{target_name}'"));
      edge.msg_from_child.dis.row(0).to_owned()
    }

    pub(super) fn get_edge_msg_to_child(
      graph: &GraphAncestral,
      partition: &PartitionMarginalDiscrete,
      source_name: &str,
      target_name: &str,
    ) -> Array1<f64> {
      let edge_key = find_edge_key(graph, source_name, target_name)
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

#[cfg(test)]
mod tests {
  use crate::commands::mugration::discrete_marginal::{attach_traits, run_discrete_marginal};
  use crate::o;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  #[test]
  fn test_discrete_marginal_attach_traits_maps_observed_and_missing_profiles() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("?"),
    };

    attach_traits(&mut partition, &graph, &traits)?;

    let node_a = helpers::get_node_data(&graph, &partition, "A");
    let actual_node_a = json_write_str(node_a, JsonPretty(true))?;
    let expected_node_a = o!(indoc! {r#"{
      "observed": 1,
      "profile": [
        0.0,
        1.0
      ],
      "log_lh": 0.0
    }"#});
    assert_eq!(expected_node_a, actual_node_a);

    let node_b = helpers::get_node_data(&graph, &partition, "B");
    let actual_node_b = json_write_str(node_b, JsonPretty(true))?;
    let expected_node_b = o!(indoc! {r#"{
      "observed": null,
      "profile": [
        0.5,
        0.5
      ],
      "log_lh": 0.0
    }"#});
    assert_eq!(expected_node_b, actual_node_b);

    assert_eq!(graph.get_edges().len(), partition.edges.len());

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_attach_traits_rejects_tree_leaf_missing_from_metadata() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
    };

    let result = attach_traits(&mut partition, &graph, &traits);
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

    let result = attach_traits(&mut partition, &graph, &traits);
    assert_error!(result, "Mugration: metadata names missing from tree leaves: C");

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_passes_normalize_backward_and_forward_profiles() -> Result<(), Report> {
    let graph: GraphAncestral = helpers::make_fixture_graph()?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = helpers::make_fixture_traits();

    attach_traits(&mut partition, &graph, &traits)?;

    graph.iter_breadth_first_reverse(|node| {
      partition.process_node_backward(&node).unwrap();
    });

    let root_after_backward = helpers::get_node_data(&graph, &partition, "root");
    helpers::assert_profile_normalized(&root_after_backward.profile);

    let inner_to_root_after_backward = helpers::get_edge_data(&graph, &partition, "root", "inner");
    helpers::assert_profile_normalized(&inner_to_root_after_backward.msg_from_child);

    let leaf_to_inner_after_backward = helpers::get_edge_data(&graph, &partition, "inner", "A");
    helpers::assert_profile_normalized(&leaf_to_inner_after_backward.msg_from_child);

    graph.iter_breadth_first_forward(|node| {
      partition.process_node_forward(&graph, &node).unwrap();
    });

    let root_after_forward = helpers::get_node_data(&graph, &partition, "root");
    helpers::assert_profile_normalized(&root_after_forward.profile);

    let inner_after_forward = helpers::get_node_data(&graph, &partition, "inner");
    helpers::assert_profile_normalized(&inner_after_forward.profile);

    let root_to_c_after_forward = helpers::get_edge_data(&graph, &partition, "root", "C");
    helpers::assert_profile_normalized(&root_to_c_after_forward.msg_to_child);

    let inner_to_b_after_forward = helpers::get_edge_data(&graph, &partition, "inner", "B");
    helpers::assert_profile_normalized(&inner_to_b_after_forward.msg_to_child);

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_run_returns_finite_log_lh_and_reconstructs_internal_trait() -> Result<(), Report> {
    let graph: GraphAncestral = helpers::make_fixture_graph()?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = helpers::make_fixture_traits();

    attach_traits(&mut partition, &graph, &traits)?;

    let actual_log_lh = run_discrete_marginal(&graph, &mut partition)?;

    assert!(actual_log_lh.is_finite());

    let inner_node = helpers::get_node_data(&graph, &partition, "inner");
    helpers::assert_profile_normalized(&inner_node.profile);

    let inner_key = helpers::get_node_key(&graph, "inner");
    let expected_trait = Some(o!("usa"));
    let actual_trait = partition.get_reconstructed_trait(inner_key);
    assert_eq!(expected_trait, actual_trait);

    Ok(())
  }

  mod helpers {
    use crate::gtr::gtr::{GTR, GTRParams};
    use crate::o;
    use crate::representation::discrete_states::DiscreteStates;
    use crate::representation::partition::discrete::PartitionDiscrete;
    use crate::representation::payload::ancestral::GraphAncestral;
    use crate::representation::payload::discrete::{DiscreteEdgeData, DiscreteNodeData};
    use crate::test_utils::{find_edge_key, find_node_key_by_name};
    use eyre::Report;
    use maplit::btreemap;
    use ndarray::Array1;
    use std::collections::BTreeMap;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read_str;
    use treetime_utils::pretty_assert_abs_diff_eq;

    pub(super) fn make_partition(states: [&str; 2]) -> Result<PartitionDiscrete, Report> {
      let discrete_states = DiscreteStates::from_values(states.into_iter(), "?");
      let n_states = discrete_states.len();
      let gtr = GTR::new(GTRParams {
        n_states,
        mu: 1.0,
        W: None,
        pi: Array1::from_elem(n_states, 1.0 / n_states as f64),
      })?;

      Ok(PartitionDiscrete::new(0, gtr, discrete_states))
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

    pub(super) fn get_node_data<'a>(
      graph: &GraphAncestral,
      partition: &'a PartitionDiscrete,
      name: &str,
    ) -> &'a DiscreteNodeData {
      let node_key = get_node_key(graph, name);
      partition
        .nodes
        .get(&node_key)
        .unwrap_or_else(|| panic!("Missing discrete node data for '{name}'"))
    }

    pub(super) fn get_edge_data<'a>(
      graph: &GraphAncestral,
      partition: &'a PartitionDiscrete,
      source_name: &str,
      target_name: &str,
    ) -> &'a DiscreteEdgeData {
      let edge_key = find_edge_key(graph, source_name, target_name)
        .unwrap_or_else(|| panic!("Missing test edge '{source_name}' -> '{target_name}'"));
      partition
        .edges
        .get(&edge_key)
        .unwrap_or_else(|| panic!("Missing discrete edge data for '{source_name}' -> '{target_name}'"))
    }
  }
}

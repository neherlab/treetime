#[cfg(test)]
mod tests {
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::gtr::{GTR, Mu};
  use crate::o;
  use crate::partition::marginal::shared::update::{MarginalBackward, MarginalForward, MarginalPasses, MarginalUpdate};
  use crate::progress::NoopProgress;
  use crate::test_utils::RecordingLog;
  use eyre::Report;
  use itertools::{Itertools, izip};
  use maplit::btreemap;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::array::ndarray::argmax_first;
  use treetime_utils::assert_error;

  #[test]
  fn test_discrete_marginal_new_maps_observed_and_missing_data_profiles_without_warning() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2)root;".as_slice())?;
    let names = nwk_parsed.names();
    let graph: Graph = nwk_parsed.graph;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("?"),
    };
    let log = RecordingLog::default();

    let (partition, _) = helpers::make_partition(&["usa", "germany"], &graph, &traits, &names, &log)?;

    let expected = vec![(o!("A"), array![0.0, 1.0]), (o!("B"), array![1.0, 1.0])];
    let actual = helpers::leaf_observations_by_name(&graph, &names, &partition, &["A", "B"]);
    assert_eq!(expected, actual);
    assert_eq!(Vec::<String>::new(), log.warnings());

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_new_fills_and_warns_tree_leaf_missing_from_metadata() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.3)root;".as_slice())?;
    let names = nwk_parsed.names();
    let graph: Graph = nwk_parsed.graph;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("C") => o!("germany"),
    };
    let log = RecordingLog::default();

    let (partition, _) = helpers::make_partition(&["usa", "germany"], &graph, &traits, &names, &log)?;

    let expected = vec![
      (o!("A"), array![0.0, 1.0]),
      (o!("B"), array![1.0, 1.0]),
      (o!("C"), array![1.0, 0.0]),
    ];
    let actual = helpers::leaf_observations_by_name(&graph, &names, &partition, &["A", "B", "C"]);
    assert_eq!(expected, actual);
    let expected_warnings = vec![
      o!("No trait value found for leaf 'B' in metadata; treating it as fully ambiguous (missing data)."),
      o!("1 of 3 tips have no matching trait value in metadata and are treated as fully ambiguous."),
    ];
    assert_eq!(expected_warnings, log.warnings());

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_new_accepts_metadata_name_missing_from_tree() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2)root;".as_slice())?;
    let names = nwk_parsed.names();
    let graph: Graph = nwk_parsed.graph;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
      o!("C") => o!("usa"),
    };
    let log = RecordingLog::default();

    let (partition, _) = helpers::make_partition(&["usa", "germany"], &graph, &traits, &names, &log)?;

    let expected = vec![(o!("A"), array![0.0, 1.0]), (o!("B"), array![1.0, 0.0])];
    let actual = helpers::leaf_observations_by_name(&graph, &names, &partition, &["A", "B"]);
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_new_maps_many_named_states_to_their_profile_columns() -> Result<(), Report> {
    let state_names = (0..140)
      .map(|i| format!("division {i:03}"))
      .chain([o!("Côte d'Ivoire"), o!("北京"), o!("São Paulo")])
      .collect_vec();
    let leaf_names = (0..state_names.len()).map(|i| format!("leaf{i}")).collect_vec();
    let nwk = format!(
      "({})root;",
      leaf_names.iter().map(|name| format!("{name}:0.1")).join(",")
    );
    let nwk_parsed = nwk_read(nwk.as_bytes())?;
    let names = nwk_parsed.names();
    let graph: Graph = nwk_parsed.graph;
    let traits: BTreeMap<String, String> =
      izip!(leaf_names.iter().cloned(), state_names.iter().rev().cloned()).collect();
    let state_refs = state_names.iter().map(String::as_str).collect_vec();
    let log = RecordingLog::default();

    let (partition, _) = helpers::make_partition(&state_refs, &graph, &traits, &names, &log)?;

    assert_eq!(state_names.len(), partition.n_states());
    let sorted_states = state_names.iter().sorted().collect_vec();
    let expected = traits
      .iter()
      .map(|(leaf, state)| {
        let column = sorted_states
          .iter()
          .position(|name| *name == state)
          .expect("state is listed");
        let profile = Array1::from_shape_fn(state_names.len(), |index| if index == column { 1.0 } else { 0.0 });
        (leaf.clone(), Some(state.clone()), profile)
      })
      .collect_vec();
    let leaf_refs = leaf_names.iter().map(String::as_str).collect_vec();
    let actual = helpers::leaf_observations_by_name(&graph, &names, &partition, &leaf_refs)
      .into_iter()
      .map(|(leaf, profile)| {
        let state = argmax_first(&profile.view())
          .and_then(|index| partition.states.name(index))
          .map(ToOwned::to_owned);
        (leaf, state, profile)
      })
      .sorted_by(|left, right| left.0.cmp(&right.0))
      .collect_vec();
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_passes_normalize_backward_and_forward_profiles() -> Result<(), Report> {
    let (graph, names, raw_branch_lengths) = helpers::make_fixture_graph()?;
    let traits = helpers::make_fixture_traits();
    let (partition, gtr) = helpers::make_partition(&["usa", "germany"], &graph, &traits, &names, &NoopProgress)?;

    let branch_lengths = branch_lengths_or_zero(&raw_branch_lengths);
    let MarginalBackward { node_states, backward } = partition.marginal_backward(&gtr, &graph, &branch_lengths, &())?;

    let root_profile = helpers::get_node_profile(&graph, &names, &node_states, "root");
    helpers::assert_profile_normalized(&root_profile);

    let inner_to_root_msg = helpers::get_edge_msg_from_child(&graph, &names, &backward, "root", "inner");
    helpers::assert_profile_normalized(&inner_to_root_msg);

    let leaf_to_inner_msg = helpers::get_edge_msg_from_child(&graph, &names, &backward, "inner", "A");
    helpers::assert_profile_normalized(&leaf_to_inner_msg);

    let MarginalForward {
      node_states, forward, ..
    } = partition.marginal_forward(&gtr, &graph, &branch_lengths, node_states, &backward)?;

    let root_profile = helpers::get_node_profile(&graph, &names, &node_states, "root");
    helpers::assert_profile_normalized(&root_profile);

    let inner_profile = helpers::get_node_profile(&graph, &names, &node_states, "inner");
    helpers::assert_profile_normalized(&inner_profile);

    let root_to_c_msg = helpers::get_edge_msg_to_child(&graph, &names, &forward, "root", "C");
    helpers::assert_profile_normalized(&root_to_c_msg);

    let inner_to_b_msg = helpers::get_edge_msg_to_child(&graph, &names, &forward, "inner", "B");
    helpers::assert_profile_normalized(&inner_to_b_msg);

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_passes_reject_site_rates_for_more_than_one_site() -> Result<(), Report> {
    let (graph, names, raw_branch_lengths) = helpers::make_fixture_graph()?;
    let traits = helpers::make_fixture_traits();
    let (partition, _) = helpers::make_partition(&["usa", "germany"], &graph, &traits, &names, &NoopProgress)?;
    let gtr = GTR::builder()
      .n_states(2)
      .mu(Mu::Normalized(1.0))
      .pi(array![0.5, 0.5])
      .site_rates(array![1.0, 1.0])
      .build()?;

    let branch_lengths = branch_lengths_or_zero(&raw_branch_lengths);
    assert_error!(
      partition.marginal_backward(&gtr, &graph, &branch_lengths, &()),
      "Number of site rates (2) differs from the number of sites in the partition (1)"
    );

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_run_returns_finite_log_lh_and_reconstructs_internal_trait() -> Result<(), Report> {
    let (graph, names, raw_branch_lengths) = helpers::make_fixture_graph()?;
    let traits = helpers::make_fixture_traits();
    let (partition, gtr) = helpers::make_partition(&["usa", "germany"], &graph, &traits, &names, &NoopProgress)?;

    let MarginalUpdate {
      node_states, log_lh, ..
    } = partition.marginal_update(&gtr, &graph, &branch_lengths_or_zero(&raw_branch_lengths), &())?;
    let actual_log_lh = log_lh.value();

    assert!(actual_log_lh.is_finite());
    assert!(
      actual_log_lh <= 0.0,
      "Log-likelihood must be non-positive: {actual_log_lh}"
    );

    let inner_profile = helpers::get_node_profile(&graph, &names, &node_states, "inner");
    helpers::assert_profile_normalized(&inner_profile);

    let inner_key = helpers::get_node_key(&graph, &names, "inner");
    let expected_trait = Some(o!("usa"));
    let actual_trait = partition.get_reconstructed_trait(&node_states, inner_key);
    assert_eq!(expected_trait, actual_trait);

    Ok(())
  }

  mod helpers {
    use crate::alphabet::state_alphabet::StateAlphabet;
    use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
    use crate::gtr::gtr::{GTR, Mu};
    use crate::mugration::leaf_observations::leaf_observations;
    use crate::o;
    use crate::partition::marginal::discrete::partition::PartitionMarginalDiscrete;
    use crate::partition::storage::dense::{DenseEdgeBackward, DenseEdgeForward, DenseNodeState};
    use crate::progress::LogSink;
    use crate::test_utils::traits_by_node;
    use crate::test_utils::{find_edge_key, find_node_key_by_name};
    use eyre::Report;
    use maplit::btreemap;
    use ndarray::Array1;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read;
    use treetime_utils::cast::ToF64Exact;
    use treetime_utils::pretty_assert_abs_diff_eq;

    pub(super) fn make_partition(
      states: &[&str],
      graph: &Graph,
      traits: &BTreeMap<String, String>,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      log: &dyn LogSink,
    ) -> Result<(PartitionMarginalDiscrete, GTR), Report> {
      let states = StateAlphabet::new(states.iter().copied())?;
      let n_states = states.len();
      let gtr = GTR::builder()
        .n_states(n_states)
        .mu(Mu::FromMatrix(1.0))
        .pi(Array1::from_elem(n_states, 1.0 / n_states.to_f64_exact()))
        .build()?;
      let traits = traits_by_node(traits.clone(), graph, names).traits;
      let observations = leaf_observations(graph, &traits, names, &states, "?", "metadata", log)?;
      let partition = PartitionMarginalDiscrete::new(states, graph, &observations, MIN_BRANCH_LENGTH_FRACTION, false)?;
      Ok((partition, gtr))
    }

    pub(super) fn leaf_observations_by_name(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      partition: &PartitionMarginalDiscrete,
      leaf_names: &[&str],
    ) -> Vec<(String, Array1<f64>)> {
      leaf_names
        .iter()
        .map(|name| {
          let profile = partition.obs_leaves[&get_node_key(graph, names, name)]
            .row(0)
            .to_owned();
          ((*name).to_owned(), profile)
        })
        .collect()
    }

    pub(super) fn make_fixture_graph() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read(b"((A:0.01,B:0.01)inner:0.01,C:0.25)root;".as_slice())?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
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
      node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
      name: &str,
    ) -> Array1<f64> {
      let node_key = get_node_key(graph, names, name);
      let node = node_states
        .get(&node_key)
        .unwrap_or_else(|| panic!("Missing discrete node data for '{name}'"));
      node.profile.dis.row(0).to_owned()
    }

    pub(super) fn get_edge_msg_from_child(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
      source_name: &str,
      target_name: &str,
    ) -> Array1<f64> {
      let edge_key = find_edge_key(graph, names, source_name, target_name)
        .unwrap_or_else(|| panic!("Missing test edge '{source_name}' -> '{target_name}'"));
      let edge = backward
        .get(&edge_key)
        .unwrap_or_else(|| panic!("Missing discrete edge data for '{source_name}' -> '{target_name}'"));
      edge.msg_from_child.dis.row(0).to_owned()
    }

    pub(super) fn get_edge_msg_to_child(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      forward: &BTreeMap<GraphEdgeKey, DenseEdgeForward>,
      source_name: &str,
      target_name: &str,
    ) -> Array1<f64> {
      let edge_key = find_edge_key(graph, names, source_name, target_name)
        .unwrap_or_else(|| panic!("Missing test edge '{source_name}' -> '{target_name}'"));
      let edge = forward
        .get(&edge_key)
        .unwrap_or_else(|| panic!("Missing discrete edge data for '{source_name}' -> '{target_name}'"));
      edge.msg_to_child.dis.row(0).to_owned()
    }
  }
}

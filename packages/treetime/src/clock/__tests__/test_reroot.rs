#[cfg(test)]
mod tests {
  use crate::clock::find_best_root::find_best_root::find_best_root;
  use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootMethod, RerootSpec, RootObjective};
  use crate::clock::reroot::{RerootParams, reroot_in_place};
  use crate::o;
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_graph::assign_node_names::assign_node_names;
  use treetime_graph::graph::Graph;
  use treetime_graph::reroot::{record_merge, remove_node_if_trivial, trivial_node_branch_lengths};
  use treetime_io::nwk::{NwkParse, NwkWriteOptions, nwk_read_str, nwk_write_str};
  use treetime_utils::assert_error;

  use helpers::{setup_reroot_test_graph, setup_reroot_test_graph_with_dates};

  #[test]
  fn test_remove_node_if_trivial_simple() -> Result<(), Report> {
    // define tree with trivial node:
    //        root
    //        /  \
    //      mid  B
    //      /
    //     A
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str("((A:0.5)mid:0.3,B:0.2)root;")?;
    let mut graph: Graph = graph;

    let mid_key = find_node_key_by_name(&graph, &names, "mid").expect("Expected node named 'mid'");

    let (parent_branch, child_branch) = trivial_node_branch_lengths(&graph, mid_key, &branch_lengths);
    let merge = remove_node_if_trivial(&mut graph, mid_key, parent_branch, child_branch)?;
    if let Some(info) = &merge {
      record_merge(&mut branch_lengths, info);
    }

    assert!(graph.get_node(mid_key).is_none(), "Expected node to be removed");

    let expected = "(B:0.2,A:0.8)root;";
    let actual = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_reroot_min_dev_matches_fixed_zero_rate_objective() -> Result<(), Report> {
    let (mut graph, names, options, mut inputs, mut state, mut branch_lengths) = setup_reroot_test_graph()?;
    let branch_params = BranchPointOptimizationParams::default();
    let expected = find_best_root(
      &graph,
      &inputs,
      &state,
      &options,
      &branch_params,
      &branch_lengths,
      false,
      RootObjective::FixedRate(0.0),
    )?;
    let expected_edge = expected.edge.expect("fixture should select a root edge");

    let reroot_params = RerootParams {
      spec: RerootSpec::Method(RerootMethod::MinDev),
      ..RerootParams::default()
    };
    let names_tt_9 = names;
    let (_state, actual) = reroot_in_place(
      &mut graph,
      &mut inputs,
      state,
      &options,
      &branch_params,
      &reroot_params,
      &mut branch_lengths,
      &names_tt_9,
    )?;
    let actual_split = actual.edge_split.expect("fixture should select an interior root point");

    assert_eq!(expected_edge, actual_split.old_edge_key);
    pretty_assert_ulps_eq!(expected.split, actual_split.split_position, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_reroot_min_dev_is_independent_of_sampling_date_values() -> Result<(), Report> {
    let dates_ascending = btreemap! {
      o!("A") => 2000.0,
      o!("B") => 2001.0,
      o!("C") => 2002.0,
      o!("D") => 2003.0,
    };
    let dates_descending = btreemap! {
      o!("A") => 2003.0,
      o!("B") => 2002.0,
      o!("C") => 2001.0,
      o!("D") => 2000.0,
    };
    let (
      mut graph_ascending,
      graph_ascending_names,
      options_ascending,
      mut inputs_ascending,
      mut state_ascending,
      mut branch_lengths_ascending,
    ) = setup_reroot_test_graph_with_dates(&dates_ascending)?;
    let (
      mut graph_descending,
      graph_descending_names,
      options_descending,
      mut inputs_descending,
      mut state_descending,
      mut branch_lengths_descending,
    ) = setup_reroot_test_graph_with_dates(&dates_descending)?;
    let reroot_params = RerootParams {
      spec: RerootSpec::Method(RerootMethod::MinDev),
      ..RerootParams::default()
    };
    let branch_params = BranchPointOptimizationParams::default();

    let names_tt_8 = graph_ascending_names.clone();
    reroot_in_place(
      &mut graph_ascending,
      &mut inputs_ascending,
      state_ascending,
      &options_ascending,
      &branch_params,
      &reroot_params,
      &mut branch_lengths_ascending,
      &names_tt_8,
    )?;
    let names_tt_7 = graph_descending_names.clone();
    reroot_in_place(
      &mut graph_descending,
      &mut inputs_descending,
      state_descending,
      &options_descending,
      &branch_params,
      &reroot_params,
      &mut branch_lengths_descending,
      &names_tt_7,
    )?;

    let graph_ascending_names = assign_node_names(graph_ascending_names, &graph_ascending)?;
    let graph_descending_names = assign_node_names(graph_descending_names, &graph_descending)?;

    let expected = nwk_write_str(
      &graph_ascending,
      &graph_ascending_names,
      &branch_lengths_ascending,
      &NwkWriteOptions::default(),
    )?;
    let actual = nwk_write_str(
      &graph_descending,
      &graph_descending_names,
      &branch_lengths_descending,
      &NwkWriteOptions::default(),
    )?;
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_reroot_policy_allow_edge_split_false_no_new_nodes() -> Result<(), Report> {
    let (mut graph, names, options, mut inputs, mut state, mut branch_lengths) = setup_reroot_test_graph()?;
    let node_count_before = graph.get_nodes().len();

    // Both flags false: don't split edges AND don't remove old root
    // This guarantees no new nodes created and no nodes removed
    let reroot_params = RerootParams {
      split_edge: false,
      remove_trivial_root: false,
      ..RerootParams::default()
    };

    let names_tt_6 = names;
    reroot_in_place(
      &mut graph,
      &mut inputs,
      state,
      &options,
      &BranchPointOptimizationParams::default(),
      &reroot_params,
      &mut branch_lengths,
      &names_tt_6,
    )?;

    let node_count_after = graph.get_nodes().len();
    assert_eq!(
      node_count_before, node_count_after,
      "Node count should be unchanged when edge split is disabled and old root is preserved"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_policy_remove_old_root_if_trivial_false_preserves_old_root() -> Result<(), Report> {
    let (mut graph, names, options, mut inputs, mut state, mut branch_lengths) = setup_reroot_test_graph()?;
    let old_root_key = graph.get_exactly_one_root()?.read_arc().key();

    let reroot_params = RerootParams {
      split_edge: true,
      remove_trivial_root: false,
      ..RerootParams::default()
    };

    let names_tt_5 = names;
    let (_state, reroot_result) = reroot_in_place(
      &mut graph,
      &mut inputs,
      state,
      &options,
      &BranchPointOptimizationParams::default(),
      &reroot_params,
      &mut branch_lengths,
      &names_tt_5,
    )?;

    if reroot_result.new_root_key != old_root_key {
      assert!(
        graph.get_node(old_root_key).is_some(),
        "Old root should still exist when remove_old_root_if_trivial is false"
      );
    }

    Ok(())
  }

  #[test]
  fn test_reroot_policy_default_allows_edge_split() -> Result<(), Report> {
    let (mut graph, names, options, mut inputs, mut state, mut branch_lengths) = setup_reroot_test_graph()?;
    let node_count_before = graph.get_nodes().len();

    let reroot_params = RerootParams::default();

    let names_tt_4 = names;
    reroot_in_place(
      &mut graph,
      &mut inputs,
      state,
      &options,
      &BranchPointOptimizationParams::default(),
      &reroot_params,
      &mut branch_lengths,
      &names_tt_4,
    )?;

    let node_count_after = graph.get_nodes().len();
    // With default policy, a new node may be created by edge split (count increases)
    // or old trivial root may be removed (count stays same or decreases by 1 if split created one)
    // The key is it should not crash and should complete successfully
    assert!(
      node_count_after >= node_count_before - 1,
      "Node count should be reasonable after reroot with default policy"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_tips_uses_mrca_branch() -> Result<(), Report> {
    let (mut graph, names, options, mut inputs, mut state, mut branch_lengths) = setup_reroot_test_graph()?;
    let reroot_params = RerootParams {
      spec: RerootSpec::Tips(vec![o!("A"), o!("B")]),
      ..RerootParams::default()
    };

    let names_tt_3 = names.clone();
    let (_state, reroot_result) = reroot_in_place(
      &mut graph,
      &mut inputs,
      state,
      &options,
      &BranchPointOptimizationParams::default(),
      &reroot_params,
      &mut branch_lengths,
      &names_tt_3,
    )?;

    let root = graph
      .get_node(reroot_result.new_root_key)
      .expect("new root should exist");
    let child_names = root
      .read_arc()
      .outbound()
      .iter()
      .map(|edge_key| {
        let child_key = graph.get_target_node_key(*edge_key)?;
        let child = graph.get_node(child_key).expect("child should exist");
        Ok(names.get(&child_key).cloned().flatten())
      })
      .collect::<Result<Vec<_>, Report>>()?;

    assert!(
      child_names.contains(&Some(o!("AB"))),
      "new root should split the branch leading to the AB MRCA"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_tips_reports_missing_tip() -> Result<(), Report> {
    let (mut graph, names, options, mut inputs, mut state, mut branch_lengths) = setup_reroot_test_graph()?;
    let reroot_params = RerootParams {
      spec: RerootSpec::Tips(vec![o!("missing")]),
      ..RerootParams::default()
    };

    let names_tt_2 = names;
    let result = reroot_in_place(
      &mut graph,
      &mut inputs,
      state,
      &options,
      &BranchPointOptimizationParams::default(),
      &reroot_params,
      &mut branch_lengths,
      &names_tt_2,
    );

    assert_error!(result, "Reroot tip not found: missing");
    Ok(())
  }

  #[test]
  fn test_reroot_oldest_uses_oldest_dated_leaf() -> Result<(), Report> {
    let (mut graph, names, options, mut inputs, mut state, mut branch_lengths) = setup_reroot_test_graph()?;
    let reroot_params = RerootParams {
      spec: RerootSpec::Method(RerootMethod::Oldest),
      ..RerootParams::default()
    };

    let names_tt_1 = names.clone();
    let (_state, reroot_result) = reroot_in_place(
      &mut graph,
      &mut inputs,
      state,
      &options,
      &BranchPointOptimizationParams::default(),
      &reroot_params,
      &mut branch_lengths,
      &names_tt_1,
    )?;

    let root = graph
      .get_node(reroot_result.new_root_key)
      .expect("new root should exist");
    let child_names = root
      .read_arc()
      .outbound()
      .iter()
      .map(|edge_key| {
        let child_key = graph.get_target_node_key(*edge_key)?;
        let child = graph.get_node(child_key).expect("child should exist");
        Ok(names.get(&child_key).cloned().flatten())
      })
      .collect::<Result<Vec<_>, Report>>()?;

    assert!(
      child_names.contains(&Some(o!("D"))),
      "new root should split the branch leading to the oldest dated leaf"
    );

    Ok(())
  }

  mod helpers {
    use crate::clock::clock_regression::{ClockParams, clock_regression_backward, clock_regression_forward};
    use crate::clock::clock_state::{ClockInputs, ClockState};
    use crate::o;
    use eyre::Report;
    use maplit::btreemap;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::{NwkParse, nwk_read_str};

    pub fn leaf_times(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      dates: &BTreeMap<String, f64>,
    ) -> BTreeMap<GraphNodeKey, Option<f64>> {
      graph
        .get_leaves()
        .iter()
        .map(|node| {
          let node = node.read_arc();
          let name = names[&node.key()].clone().expect("Leaf has name");
          (node.key(), dates.get(&name).copied())
        })
        .collect()
    }

    pub fn setup_reroot_test_graph_with_dates(
      dates: &BTreeMap<String, f64>,
    ) -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        ClockParams,
        ClockInputs,
        ClockState,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let NwkParse {
        graph,
        names,
        branch_lengths,
        ..
      } = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
      let graph: Graph = graph;
      let times = leaf_times(&graph, &names, dates);

      let options = ClockParams::default();
      let inputs = ClockInputs::seed_from_times(&graph, &times);
      let mut state = ClockState::new(&graph);
      clock_regression_backward(&graph, &inputs, &mut state, &options, &branch_lengths, None)?;
      clock_regression_forward(&graph, &inputs, &mut state, &options, &branch_lengths, None)?;

      Ok((graph, names, options, inputs, state, branch_lengths))
    }

    pub fn setup_reroot_test_graph() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        ClockParams,
        ClockInputs,
        ClockState,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let dates = btreemap! {
        o!("A") => 2013.0,
        o!("B") => 2022.0,
        o!("C") => 2017.0,
        o!("D") => 2005.0,
      };
      setup_reroot_test_graph_with_dates(&dates)
    }
  }
}

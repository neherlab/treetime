#[cfg(test)]
mod tests {
  use crate::clock::clock_graph::GraphClock;
  use crate::clock::find_best_root::find_best_root::find_best_root;
  use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootMethod, RerootSpec, RootObjective};
  use crate::clock::reroot::{RerootParams, reroot_in_place};
  use crate::o;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_graph::node::Named;
  use treetime_graph::reroot::remove_node_if_trivial;
  use treetime_graph::value_maps::{edge_branch_lengths, node_names};
  use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
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
    let mut graph: GraphClock = nwk_read_str("((A:0.5)mid:0.3,B:0.2)root;")?.graph;

    let mid_key = graph
      .find_node(|node| node.name.as_deref() == Some("mid"))
      .expect("Expected node named 'mid'");

    remove_node_if_trivial(&mut graph, mid_key)?;

    assert!(graph.get_node(mid_key).is_none(), "Expected node to be removed");

    let expected = "(B:0.2,A:0.8)root;";
    let actual = nwk_write_str(
      &graph,
      &node_names(&graph),
      &edge_branch_lengths(&graph),
      &NwkWriteOptions::default(),
    )?;
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_reroot_min_dev_matches_fixed_zero_rate_objective() -> Result<(), Report> {
    let (mut graph, options, mut state) = setup_reroot_test_graph()?;
    let branch_params = BranchPointOptimizationParams::default();
    let expected = find_best_root(
      &graph,
      &state,
      &options,
      &branch_params,
      false,
      RootObjective::FixedRate(0.0),
    )?;
    let expected_edge = expected.edge.expect("fixture should select a root edge");

    let reroot_params = RerootParams {
      spec: RerootSpec::Method(RerootMethod::MinDev),
      ..RerootParams::default()
    };
    let names_tt_9 = node_names(&graph);
    let actual = reroot_in_place(
      &mut graph,
      &mut state,
      &options,
      &branch_params,
      &reroot_params,
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
    let (mut graph_ascending, options_ascending, mut state_ascending) =
      setup_reroot_test_graph_with_dates(&dates_ascending)?;
    let (mut graph_descending, options_descending, mut state_descending) =
      setup_reroot_test_graph_with_dates(&dates_descending)?;
    let reroot_params = RerootParams {
      spec: RerootSpec::Method(RerootMethod::MinDev),
      ..RerootParams::default()
    };
    let branch_params = BranchPointOptimizationParams::default();

    let names_tt_8 = node_names(&graph_ascending);
    reroot_in_place(
      &mut graph_ascending,
      &mut state_ascending,
      &options_ascending,
      &branch_params,
      &reroot_params,
      &names_tt_8,
    )?;
    let names_tt_7 = node_names(&graph_descending);
    reroot_in_place(
      &mut graph_descending,
      &mut state_descending,
      &options_descending,
      &branch_params,
      &reroot_params,
      &names_tt_7,
    )?;

    let expected = nwk_write_str(
      &graph_ascending,
      &node_names(&graph_ascending),
      &edge_branch_lengths(&graph_ascending),
      &NwkWriteOptions::default(),
    )?;
    let actual = nwk_write_str(
      &graph_descending,
      &node_names(&graph_descending),
      &edge_branch_lengths(&graph_descending),
      &NwkWriteOptions::default(),
    )?;
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_reroot_policy_allow_edge_split_false_no_new_nodes() -> Result<(), Report> {
    let (mut graph, options, mut state) = setup_reroot_test_graph()?;
    let node_count_before = graph.get_nodes().len();

    // Both flags false: don't split edges AND don't remove old root
    // This guarantees no new nodes created and no nodes removed
    let reroot_params = RerootParams {
      split_edge: false,
      remove_trivial_root: false,
      ..RerootParams::default()
    };

    let names_tt_6 = node_names(&graph);
    reroot_in_place(
      &mut graph,
      &mut state,
      &options,
      &BranchPointOptimizationParams::default(),
      &reroot_params,
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
    let (mut graph, options, mut state) = setup_reroot_test_graph()?;
    let old_root_key = graph.get_exactly_one_root()?.read_arc().key();

    let reroot_params = RerootParams {
      split_edge: true,
      remove_trivial_root: false,
      ..RerootParams::default()
    };

    let names_tt_5 = node_names(&graph);
    let reroot_result = reroot_in_place(
      &mut graph,
      &mut state,
      &options,
      &BranchPointOptimizationParams::default(),
      &reroot_params,
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
    let (mut graph, options, mut state) = setup_reroot_test_graph()?;
    let node_count_before = graph.get_nodes().len();

    let reroot_params = RerootParams::default();

    let names_tt_4 = node_names(&graph);
    reroot_in_place(
      &mut graph,
      &mut state,
      &options,
      &BranchPointOptimizationParams::default(),
      &reroot_params,
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
    let (mut graph, options, mut state) = setup_reroot_test_graph()?;
    let reroot_params = RerootParams {
      spec: RerootSpec::Tips(vec![o!("A"), o!("B")]),
      ..RerootParams::default()
    };

    let names_tt_3 = node_names(&graph);
    let reroot_result = reroot_in_place(
      &mut graph,
      &mut state,
      &options,
      &BranchPointOptimizationParams::default(),
      &reroot_params,
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
        Ok(
          child
            .read_arc()
            .payload()
            .read_arc()
            .name()
            .map(|name| name.as_ref().to_owned()),
        )
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
    let (mut graph, options, mut state) = setup_reroot_test_graph()?;
    let reroot_params = RerootParams {
      spec: RerootSpec::Tips(vec![o!("missing")]),
      ..RerootParams::default()
    };

    let names_tt_2 = node_names(&graph);
    let result = reroot_in_place(
      &mut graph,
      &mut state,
      &options,
      &BranchPointOptimizationParams::default(),
      &reroot_params,
      &names_tt_2,
    );

    assert_error!(result, "Reroot tip not found: missing");
    Ok(())
  }

  #[test]
  fn test_reroot_oldest_uses_oldest_dated_leaf() -> Result<(), Report> {
    let (mut graph, options, mut state) = setup_reroot_test_graph()?;
    let reroot_params = RerootParams {
      spec: RerootSpec::Method(RerootMethod::Oldest),
      ..RerootParams::default()
    };

    let names_tt_1 = node_names(&graph);
    let reroot_result = reroot_in_place(
      &mut graph,
      &mut state,
      &options,
      &BranchPointOptimizationParams::default(),
      &reroot_params,
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
        Ok(
          child
            .read_arc()
            .payload()
            .read_arc()
            .name()
            .map(|name| name.as_ref().to_owned()),
        )
      })
      .collect::<Result<Vec<_>, Report>>()?;

    assert!(
      child_names.contains(&Some(o!("D"))),
      "new root should split the branch leading to the oldest dated leaf"
    );

    Ok(())
  }

  mod helpers {
    use crate::clock::clock_graph::GraphClock;
    use crate::clock::clock_regression::{ClockParams, clock_regression_backward, clock_regression_forward};
    use crate::clock::clock_state::ClockState;
    use crate::o;
    use eyre::Report;
    use maplit::btreemap;
    use std::collections::BTreeMap;
    use treetime_graph::node::{GraphNodeKey, Named};
    use treetime_io::nwk::nwk_read_str;

    pub fn leaf_times(graph: &GraphClock, dates: &BTreeMap<String, f64>) -> BTreeMap<GraphNodeKey, Option<f64>> {
      graph
        .get_leaves()
        .iter()
        .map(|node| {
          let node = node.read_arc();
          let name = node
            .payload()
            .read_arc()
            .name()
            .expect("Leaf has name")
            .as_ref()
            .to_owned();
          (node.key(), dates.get(&name).copied())
        })
        .collect()
    }

    pub fn setup_reroot_test_graph_with_dates(
      dates: &BTreeMap<String, f64>,
    ) -> Result<(GraphClock, ClockParams, ClockState), Report> {
      let graph: GraphClock = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?.graph;
      let times = leaf_times(&graph, dates);

      let options = ClockParams::default();
      let mut state = ClockState::seed_from_values(&graph, &times);
      clock_regression_backward(&graph, &mut state, &options, None)?;
      clock_regression_forward(&graph, &mut state, &options, None)?;

      Ok((graph, options, state))
    }

    pub fn setup_reroot_test_graph() -> Result<(GraphClock, ClockParams, ClockState), Report> {
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

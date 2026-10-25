#[cfg(test)]
mod tests {
  use crate::o;
  use crate::reroot::orchestrate::RerootTopologyParams;
  use crate::reroot::params::{RerootMethod, RerootParams, RerootSpec};
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use helpers::setup_reroot_test_graph;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_graph::reroot::{record_merge, remove_node_if_trivial, trivial_node_branch_lengths};
  use treetime_graph::schedule::GraphSchedule;
  use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};
  use treetime_utils::assert_error;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  #[test]
  fn test_remove_node_if_trivial_simple() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.5)mid:0.3,B:0.2)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let mid_key = find_node_key_by_name(&graph, &names, "mid").expect("Expected node named 'mid'");

    let (parent_branch, child_branch) = trivial_node_branch_lengths(&graph, mid_key, &branch_lengths);
    let merge = remove_node_if_trivial(&mut graph, mid_key, parent_branch, child_branch)?;
    if let Some(info) = &merge {
      record_merge(&mut branch_lengths, info);
    }

    assert!(graph.get_node(mid_key).is_none(), "Expected node to be removed");

    let expected = "(B:0.2,A:0.8)root;";
    let actual = nwk_write_str(
      &GraphSchedule::new(&graph)?.into_tree()?,
      &names,
      &branch_lengths,
      &NwkWriteOptions::default(),
      &NwkNodeComments::new(),
    )?;
    assert_eq!(expected, actual);

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
    let reroot_params = RerootParams {
      spec: RerootSpec::Method(RerootMethod::MinDev),
      ..RerootParams::default()
    };

    let expected = helpers::rerooted_newick(&dates_ascending, &reroot_params)?;
    let actual = helpers::rerooted_newick(&dates_descending, &reroot_params)?;
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_reroot_policy_allow_edge_split_false_no_new_nodes() -> Result<(), Report> {
    let fixture = setup_reroot_test_graph()?;
    let node_count_before = fixture.tree.graph.get_nodes().count();

    let reroot_params = RerootParams {
      topology: RerootTopologyParams {
        split_edge: false,
        remove_trivial_root: false,
      },
      ..RerootParams::default()
    };
    let (tree, names, _) = helpers::reroot(fixture, &reroot_params)?;

    let node_count_after = tree.graph.get_nodes().count();
    assert_eq!(
      node_count_before, node_count_after,
      "Node count should be unchanged when edge split is disabled and old root is preserved"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_policy_remove_old_root_if_trivial_false_preserves_old_root() -> Result<(), Report> {
    let fixture = setup_reroot_test_graph()?;
    let old_root_key = fixture.tree.graph.get_exactly_one_root()?.key();

    let reroot_params = RerootParams {
      topology: RerootTopologyParams {
        split_edge: true,
        remove_trivial_root: false,
      },
      ..RerootParams::default()
    };
    let (tree, names, new_root_key) = helpers::reroot(fixture, &reroot_params)?;

    if new_root_key != old_root_key {
      assert!(
        tree.graph.get_node(old_root_key).is_some(),
        "Old root should still exist when remove_old_root_if_trivial is false"
      );
    }

    Ok(())
  }

  #[test]
  fn test_reroot_policy_default_allows_edge_split() -> Result<(), Report> {
    let fixture = setup_reroot_test_graph()?;
    let node_count_before = fixture.tree.graph.get_nodes().count();

    let (tree, names, _) = helpers::reroot(fixture, &RerootParams::default())?;

    let node_count_after = tree.graph.get_nodes().count();
    assert!(
      node_count_after >= node_count_before - 1,
      "Node count should be reasonable after reroot with default policy"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_tips_uses_mrca_branch() -> Result<(), Report> {
    let fixture = setup_reroot_test_graph()?;
    let reroot_params = RerootParams {
      spec: RerootSpec::Tips(helpers::tip_keys(&fixture.tree.graph, &fixture.names, &["A", "B"])),
      ..RerootParams::default()
    };

    let (tree, names, new_root_key) = helpers::reroot(fixture, &reroot_params)?;

    let child_names = helpers::child_names(&tree.graph, &names, new_root_key)?;
    assert!(
      child_names.contains(&o!("AB")),
      "new root should split the branch leading to the AB MRCA"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_tips_rejects_an_empty_tip_list() -> Result<(), Report> {
    let reroot_params = RerootParams {
      spec: RerootSpec::Tips(vec![]),
      ..RerootParams::default()
    };

    let result = helpers::reroot(setup_reroot_test_graph()?, &reroot_params);

    assert_error!(result, "--outgroup requires at least one tip name");
    Ok(())
  }

  #[test]
  fn test_reroot_oldest_uses_oldest_dated_leaf() -> Result<(), Report> {
    let fixture = setup_reroot_test_graph()?;
    let reroot_params = RerootParams {
      spec: RerootSpec::Method(RerootMethod::Oldest),
      ..RerootParams::default()
    };

    let (tree, names, new_root_key) = helpers::reroot(fixture, &reroot_params)?;

    let child_names = helpers::child_names(&tree.graph, &names, new_root_key)?;
    assert!(
      child_names.contains(&o!("D")),
      "new root should split the branch leading to the oldest dated leaf"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_min_dev_places_the_root_at_the_minimum_variance_point() -> Result<(), Report> {
    let fixture = setup_reroot_test_graph()?;
    let reroot_params = RerootParams {
      spec: RerootSpec::Method(RerootMethod::MinDev),
      ..RerootParams::default()
    };

    let (tree, names, result) = helpers::reroot_with_result(fixture, &reroot_params)?;

    let expected = btreemap! { o!("AB") => 0.08, o!("CD") => 0.07 };
    let actual = helpers::root_child_lengths(&tree, &names, result.new_root_key)?;
    pretty_assert_map_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_reroot_merged_old_root_leaves_no_clock_inputs_behind() -> Result<(), Report> {
    let fixture = setup_reroot_test_graph()?;
    let old_root_key = fixture.tree.graph.get_exactly_one_root()?.key();

    let (tree, names, result) = helpers::reroot_with_result(fixture, &RerootParams::default())?;

    let merge = result.edge_merge.expect("the old root becomes trivial and is merged");
    assert_eq!(old_root_key, merge.removed_node_key);
    assert!(tree.graph.get_node(old_root_key).is_none());
    helpers::assert_inputs_match_graph(&tree);
    Ok(())
  }

  #[test]
  fn test_reroot_removes_an_undated_stem_root_when_the_root_moves() -> Result<(), Report> {
    let fixture = helpers::setup_stem_graph(&helpers::leaf_dates())?;
    let names = fixture.names.clone();
    let stem_key = find_node_key_by_name(&fixture.tree.graph, &names, "STEM").expect("STEM exists");
    let reroot_params = RerootParams {
      spec: RerootSpec::Tips(helpers::tip_keys(&fixture.tree.graph, &fixture.names, &["A", "B"])),
      ..RerootParams::default()
    };

    let (tree, names, result) = helpers::reroot_with_result(fixture, &reroot_params)?;

    let stem = result.stem_removal.expect("the undated stem is removed");
    assert_eq!(stem_key, stem.removed_node_key);
    assert!(tree.graph.get_node(stem_key).is_none());
    assert!(!tree.branch_lengths.contains_key(stem.removed_edge_key));
    helpers::assert_inputs_match_graph(&tree);
    assert_eq!(
      vec![o!("A"), o!("B"), o!("C"), o!("D")],
      helpers::leaf_names(&tree, &names)
    );
    let middle_of_the_unrooted_ab_cd_branch = btreemap! { o!("AB") => 0.075, o!("CD") => 0.075 };
    let actual = helpers::root_child_lengths(&tree, &names, result.new_root_key)?;
    pretty_assert_map_abs_diff_eq!(middle_of_the_unrooted_ab_cd_branch, actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_reroot_least_squares_root_ignores_an_undated_stem_root() -> Result<(), Report> {
    let without_stem = helpers::setup_unstemmed_graph(&helpers::leaf_dates())?;
    let with_stem = helpers::setup_stem_graph(&helpers::leaf_dates())?;

    let (tree_without_stem, names_without_stem, result_without_stem) =
      helpers::reroot_with_result(without_stem, &RerootParams::default())?;
    let (tree_with_stem, names_with_stem, result_with_stem) =
      helpers::reroot_with_result(with_stem, &RerootParams::default())?;

    assert!(result_with_stem.stem_removal.is_some());
    let expected = helpers::root_child_lengths(
      &tree_without_stem,
      &names_without_stem,
      result_without_stem.new_root_key,
    )?;
    let actual = helpers::root_child_lengths(&tree_with_stem, &names_with_stem, result_with_stem.new_root_key)?;
    pretty_assert_map_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_reroot_keeps_a_dated_stem_root_as_a_dated_leaf() -> Result<(), Report> {
    let mut dates = helpers::leaf_dates();
    dates.insert(o!("STEM"), 1990.0);
    let fixture = helpers::setup_stem_graph(&dates)?;
    let reroot_params = RerootParams {
      spec: RerootSpec::Tips(helpers::tip_keys(&fixture.tree.graph, &fixture.names, &["A", "B"])),
      ..RerootParams::default()
    };

    let (tree, names, result) = helpers::reroot_with_result(fixture, &reroot_params)?;

    assert!(result.stem_removal.is_none());
    assert_eq!(
      vec![o!("A"), o!("B"), o!("C"), o!("D"), o!("STEM")],
      helpers::leaf_names(&tree, &names)
    );
    let stem_key = find_node_key_by_name(&tree.graph, &names, "STEM").expect("STEM stays");
    assert_eq!(Some(1990.0), tree.inputs.likely_time(stem_key));
    helpers::assert_inputs_match_graph(&tree);
    Ok(())
  }

  #[test]
  fn test_reroot_past_a_tip_roots_at_a_new_node_with_a_zero_length_branch_to_the_tip() -> Result<(), Report> {
    let fixture = helpers::setup_tip_past_optimum_graph()?;

    let (tree, names, result) = helpers::reroot_at_fixed_rate(fixture, &RerootParams::default())?;

    let a_key = find_node_key_by_name(&tree.graph, &names, "A").expect("A exists");
    assert!(
      tree.graph.is_leaf(a_key),
      "v0 roots a split at the tip end with Biopython root_with_outgroup(tip): a new root with a 0-length branch to \
       the tip (packages/legacy/treetime/treetime/treeregression.py#L443-L463)"
    );
    assert_eq!("NODE_0000000", &names[result.new_root_key]);
    helpers::assert_inputs_match_graph(&tree);
    let expected = btreemap! { o!("A") => 0.0, o!("N") => 0.002 };
    let actual = helpers::root_child_lengths(&tree, &names, result.new_root_key)?;
    pretty_assert_map_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_reroot_past_a_tip_without_edge_split_keeps_the_parent_end() -> Result<(), Report> {
    let fixture = helpers::setup_tip_past_optimum_graph()?;
    let reroot_params = RerootParams {
      topology: RerootTopologyParams {
        split_edge: false,
        ..RerootTopologyParams::default()
      },
      ..RerootParams::default()
    };

    let (tree, names, result) = helpers::reroot_at_fixed_rate(fixture, &reroot_params)?;

    assert_eq!("R", &names[result.new_root_key]);
    let expected = btreemap! { o!("A") => 0.001, o!("N") => 0.001 };
    let actual = helpers::root_child_lengths(&tree, &names, result.new_root_key)?;
    pretty_assert_map_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  mod helpers {
    use crate::clock::clock_regression::ClockTree;
    use crate::clock::clock_state::ClockInputs;
    use crate::o;
    use crate::reroot::params::RerootParams;
    use crate::test_utils::{find_node_key_by_name, reroot_and_fit_clock};
    use eyre::Report;
    use itertools::Itertools;
    use maplit::btreemap;
    use pretty_assertions::assert_eq;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::node_names::NodeNames;
    use treetime_graph::reroot::RerootResult;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};

    const TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    const STEM_TREE: &str = "(((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)R:0.001)STEM;";

    const UNSTEMMED_TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)R;";

    const TIP_PAST_OPTIMUM_TREE: &str = "(A:0.001,(B:0.002,C:0.003)N:0.001)R;";

    const TIP_PAST_OPTIMUM_RATE: f64 = 0.01;

    pub(super) struct RerootFixture {
      pub tree: ClockTree,
      pub names: NodeNames,
    }

    pub(super) fn tip_keys(graph: &Graph, names: &NodeNames, tips: &[&str]) -> Vec<GraphNodeKey> {
      tips
        .iter()
        .map(|tip| find_node_key_by_name(graph, names, tip).expect("the tip exists"))
        .collect()
    }

    pub(super) fn leaf_dates() -> BTreeMap<String, f64> {
      btreemap! {
        o!("A") => 2013.0,
        o!("B") => 2022.0,
        o!("C") => 2017.0,
        o!("D") => 2005.0,
      }
    }

    pub(super) fn setup_reroot_test_graph_with_dates(dates: &BTreeMap<String, f64>) -> Result<RerootFixture, Report> {
      setup_graph(TREE, dates)
    }

    pub(super) fn setup_reroot_test_graph() -> Result<RerootFixture, Report> {
      setup_reroot_test_graph_with_dates(&leaf_dates())
    }

    pub(super) fn setup_stem_graph(dates: &BTreeMap<String, f64>) -> Result<RerootFixture, Report> {
      setup_graph(STEM_TREE, dates)
    }

    pub(super) fn setup_unstemmed_graph(dates: &BTreeMap<String, f64>) -> Result<RerootFixture, Report> {
      setup_graph(UNSTEMMED_TREE, dates)
    }

    pub(super) fn setup_tip_past_optimum_graph() -> Result<RerootFixture, Report> {
      let dates = btreemap! {
        o!("A") => 2000.0,
        o!("B") => 2001.0,
        o!("C") => 2001.0,
      };
      setup_graph(TIP_PAST_OPTIMUM_TREE, &dates)
    }

    pub(super) fn reroot_at_fixed_rate(
      fixture: RerootFixture,
      reroot_params: &RerootParams,
    ) -> Result<(ClockTree, NodeNames, RerootResult), Report> {
      let RerootFixture { mut tree, mut names } = fixture;
      let (result, _) = reroot_and_fit_clock(
        &mut tree,
        &BTreeSet::new(),
        Some(TIP_PAST_OPTIMUM_RATE),
        reroot_params,
        &mut names,
      )?;
      Ok((tree, names, result))
    }

    pub(super) fn reroot_with_result(
      fixture: RerootFixture,
      reroot_params: &RerootParams,
    ) -> Result<(ClockTree, NodeNames, RerootResult), Report> {
      let RerootFixture { mut tree, mut names } = fixture;
      let (result, _) = reroot_and_fit_clock(&mut tree, &BTreeSet::new(), None, reroot_params, &mut names)?;
      Ok((tree, names, result))
    }

    pub(super) fn reroot(
      fixture: RerootFixture,
      reroot_params: &RerootParams,
    ) -> Result<(ClockTree, NodeNames, GraphNodeKey), Report> {
      let (tree, names, result) = reroot_with_result(fixture, reroot_params)?;
      Ok((tree, names, result.new_root_key))
    }

    pub(super) fn rerooted_newick(
      dates: &BTreeMap<String, f64>,
      reroot_params: &RerootParams,
    ) -> Result<String, Report> {
      let fixture = setup_reroot_test_graph_with_dates(dates)?;
      let (tree, names, _) = reroot(fixture, reroot_params)?;
      nwk_write_str(
        &GraphSchedule::new(&tree.graph)?.into_tree()?,
        &names,
        &tree.branch_lengths,
        &NwkWriteOptions::default(),
        &NwkNodeComments::new(),
      )
    }

    pub(super) fn child_names(graph: &Graph, names: &NodeNames, node_key: GraphNodeKey) -> Result<Vec<String>, Report> {
      let node = graph.get_node(node_key).expect("node should exist");
      node
        .outbound()
        .iter()
        .map(|edge_key| {
          let child_key = graph.get_target_node_key(*edge_key)?;
          Ok(names[child_key].to_owned())
        })
        .collect()
    }

    pub(super) fn root_child_lengths(
      tree: &ClockTree,
      names: &NodeNames,
      root_key: GraphNodeKey,
    ) -> Result<BTreeMap<String, f64>, Report> {
      let root = tree.graph.get_node(root_key).expect("root should exist");
      root
        .outbound()
        .iter()
        .map(|edge_key| {
          let child_key = tree.graph.get_target_node_key(*edge_key)?;
          let name = names[child_key].to_owned();
          let length = tree.branch_lengths[edge_key].expect("root edges have a length");
          Ok((name, length))
        })
        .collect()
    }

    pub(super) fn leaf_names(tree: &ClockTree, names: &NodeNames) -> Vec<String> {
      tree
        .graph
        .get_leaves()
        .map(|leaf| names[leaf.key()].to_owned())
        .sorted()
        .collect()
    }

    pub(super) fn assert_inputs_match_graph(tree: &ClockTree) {
      let node_keys: BTreeSet<_> = tree.graph.get_nodes().map(|node| node.key()).collect();
      let edge_keys: BTreeSet<_> = tree.graph.get_edges().map(|edge| edge.key()).collect();
      assert_eq!(node_keys, tree.inputs.nodes.keys().collect());
      assert_eq!(edge_keys, tree.inputs.edges.keys().collect());
      assert_eq!(edge_keys, tree.branch_lengths.keys().collect());
    }

    fn setup_graph(newick: &str, dates: &BTreeMap<String, f64>) -> Result<RerootFixture, Report> {
      let nwk_parsed = nwk_read(newick.as_bytes())?;
      let names = nwk_parsed.names;
      let graph = nwk_parsed.graph;
      let times = graph
        .get_nodes()
        .map(|node| (node.key(), dates.get(&names[node.key()]).copied()))
        .collect();
      let inputs = ClockInputs::from_times(&graph, &times, &EdgeMap::new());
      Ok(RerootFixture {
        tree: ClockTree {
          graph,
          branch_lengths: nwk_parsed.branch_lengths,
          inputs,
        },
        names,
      })
    }
  }
}

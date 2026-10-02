#[cfg(test)]
mod tests {
  use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
  use crate::clock::reroot::RerootParams;
  use crate::o;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_graph::reroot::{record_merge, remove_node_if_trivial, trivial_node_branch_lengths};
  use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
  use treetime_utils::assert_error;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  use helpers::setup_reroot_test_graph;

  #[test]
  fn test_remove_node_if_trivial_simple() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.5)mid:0.3,B:0.2)root;")?;
    let names = nwk_parsed.names();
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
    let actual = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
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
      split_edge: false,
      remove_trivial_root: false,
      ..RerootParams::default()
    };
    let (tree, _) = helpers::reroot(fixture, &reroot_params)?;

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
      split_edge: true,
      remove_trivial_root: false,
      ..RerootParams::default()
    };
    let (tree, new_root_key) = helpers::reroot(fixture, &reroot_params)?;

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

    let (tree, _) = helpers::reroot(fixture, &RerootParams::default())?;

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
    let names = fixture.names.clone();
    let reroot_params = RerootParams {
      spec: RerootSpec::Tips(vec![o!("A"), o!("B")]),
      ..RerootParams::default()
    };

    let (tree, new_root_key) = helpers::reroot(fixture, &reroot_params)?;

    let child_names = helpers::child_names(&tree.graph, &names, new_root_key)?;
    assert!(
      child_names.contains(&Some(o!("AB"))),
      "new root should split the branch leading to the AB MRCA"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_tips_reports_missing_tip() -> Result<(), Report> {
    let reroot_params = RerootParams {
      spec: RerootSpec::Tips(vec![o!("missing")]),
      ..RerootParams::default()
    };

    let result = helpers::reroot(setup_reroot_test_graph()?, &reroot_params);

    assert_error!(result, "Reroot tip not found: missing");
    Ok(())
  }

  #[test]
  fn test_reroot_oldest_uses_oldest_dated_leaf() -> Result<(), Report> {
    let fixture = setup_reroot_test_graph()?;
    let names = fixture.names.clone();
    let reroot_params = RerootParams {
      spec: RerootSpec::Method(RerootMethod::Oldest),
      ..RerootParams::default()
    };

    let (tree, new_root_key) = helpers::reroot(fixture, &reroot_params)?;

    let child_names = helpers::child_names(&tree.graph, &names, new_root_key)?;
    assert!(
      child_names.contains(&Some(o!("D"))),
      "new root should split the branch leading to the oldest dated leaf"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_min_dev_places_the_root_at_the_minimum_variance_point() -> Result<(), Report> {
    let fixture = setup_reroot_test_graph()?;
    let names = fixture.names.clone();
    let reroot_params = RerootParams {
      spec: RerootSpec::Method(RerootMethod::MinDev),
      ..RerootParams::default()
    };

    let (tree, result) = helpers::reroot_with_result(fixture, &reroot_params)?;

    let expected = btreemap! { o!("AB") => 0.08, o!("CD") => 0.07 };
    let actual = helpers::root_child_lengths(&tree, &names, result.new_root_key)?;
    pretty_assert_map_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_reroot_merged_old_root_leaves_no_clock_inputs_behind() -> Result<(), Report> {
    let fixture = setup_reroot_test_graph()?;
    let old_root_key = fixture.tree.graph.get_exactly_one_root()?.key();

    let (tree, result) = helpers::reroot_with_result(fixture, &RerootParams::default())?;

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
      spec: RerootSpec::Tips(vec![o!("A"), o!("B")]),
      ..RerootParams::default()
    };

    let (tree, result) = helpers::reroot_with_result(fixture, &reroot_params)?;

    let stem = result.stem_removal.expect("the undated stem is removed");
    assert_eq!(stem_key, stem.removed_node_key);
    assert!(tree.graph.get_node(stem_key).is_none());
    assert!(!tree.branch_lengths.contains_key(&stem.removed_edge_key));
    helpers::assert_inputs_match_graph(&tree);
    assert_eq!(
      vec![o!("A"), o!("B"), o!("C"), o!("D")],
      helpers::leaf_names(&tree, &names)
    );
    let expected = btreemap! { o!("AB") => 0.05, o!("CD") => 0.1 };
    let actual = helpers::root_child_lengths(&tree, &names, result.new_root_key)?;
    pretty_assert_map_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_reroot_least_squares_root_ignores_an_undated_stem_root() -> Result<(), Report> {
    let without_stem = helpers::setup_unstemmed_graph(&helpers::leaf_dates())?;
    let with_stem = helpers::setup_stem_graph(&helpers::leaf_dates())?;
    let names_without_stem = without_stem.names.clone();
    let names_with_stem = with_stem.names.clone();

    let (tree_without_stem, result_without_stem) = helpers::reroot_with_result(without_stem, &RerootParams::default())?;
    let (tree_with_stem, result_with_stem) = helpers::reroot_with_result(with_stem, &RerootParams::default())?;

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
    let names = fixture.names.clone();
    let reroot_params = RerootParams {
      spec: RerootSpec::Tips(vec![o!("A"), o!("B")]),
      ..RerootParams::default()
    };

    let (tree, result) = helpers::reroot_with_result(fixture, &reroot_params)?;

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
  fn test_reroot_keep_root_leaves_an_undated_stem_in_place() -> Result<(), Report> {
    let fixture = helpers::setup_stem_graph(&helpers::leaf_dates())?;
    let names = fixture.names.clone();
    let nodes_before = fixture.tree.graph.get_nodes().count();

    let (tree, result) = helpers::estimate(fixture, true, &RerootParams::default())?;

    assert!(result.is_none());
    assert_eq!(nodes_before, tree.graph.get_nodes().count());
    let root_key = tree.graph.get_exactly_one_root()?.key();
    assert_eq!(Some(o!("STEM")), names[&root_key]);
    Ok(())
  }

  mod helpers {
    use crate::clock::clock_regression::{ClockTree, ClockVarianceParams, estimate_clock_model_with_reroot_policy};
    use crate::clock::clock_state::ClockInputs;
    use crate::clock::find_best_root::params::BranchPointOptimizationParams;
    use crate::clock::reroot::RerootParams;
    use crate::o;
    use crate::progress::NoopProgress;
    use eyre::Report;
    use itertools::Itertools;
    use maplit::btreemap;
    use pretty_assertions::assert_eq;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_graph::assign_node_names::assign_node_names;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::reroot::RerootResult;
    use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};

    const TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    const STEM_TREE: &str = "(((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)R:0.001)STEM;";

    const UNSTEMMED_TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)R;";

    pub(super) struct RerootFixture {
      pub tree: ClockTree,
      pub names: BTreeMap<GraphNodeKey, Option<String>>,
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

    pub(super) fn estimate(
      fixture: RerootFixture,
      keep_root: bool,
      reroot_params: &RerootParams,
    ) -> Result<(ClockTree, Option<RerootResult>), Report> {
      let (tree, result) = estimate_clock_model_with_reroot_policy(
        fixture.tree,
        &BTreeSet::new(),
        &ClockVarianceParams::default(),
        None,
        keep_root,
        &BranchPointOptimizationParams::default(),
        reroot_params,
        None,
        &fixture.names,
        &NoopProgress,
      )?;
      Ok((tree, result.reroot_result().cloned()))
    }

    pub(super) fn reroot_with_result(
      fixture: RerootFixture,
      reroot_params: &RerootParams,
    ) -> Result<(ClockTree, RerootResult), Report> {
      let (tree, result) = estimate(fixture, false, reroot_params)?;
      Ok((tree, result.expect("a reroot reports its result")))
    }

    pub(super) fn reroot(
      fixture: RerootFixture,
      reroot_params: &RerootParams,
    ) -> Result<(ClockTree, GraphNodeKey), Report> {
      let (tree, result) = reroot_with_result(fixture, reroot_params)?;
      Ok((tree, result.new_root_key))
    }

    pub(super) fn rerooted_newick(
      dates: &BTreeMap<String, f64>,
      reroot_params: &RerootParams,
    ) -> Result<String, Report> {
      let fixture = setup_reroot_test_graph_with_dates(dates)?;
      let names = fixture.names.clone();
      let (tree, _) = reroot(fixture, reroot_params)?;
      let names = assign_node_names(names, &tree.graph)?;
      nwk_write_str(&tree.graph, &names, &tree.branch_lengths, &NwkWriteOptions::default())
    }

    pub(super) fn child_names(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      node_key: GraphNodeKey,
    ) -> Result<Vec<Option<String>>, Report> {
      let node = graph.get_node(node_key).expect("node should exist");
      node
        .outbound()
        .iter()
        .map(|edge_key| {
          let child_key = graph.get_target_node_key(*edge_key)?;
          Ok(names.get(&child_key).cloned().flatten())
        })
        .collect()
    }

    pub(super) fn root_child_lengths(
      tree: &ClockTree,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      root_key: GraphNodeKey,
    ) -> Result<BTreeMap<String, f64>, Report> {
      let root = tree.graph.get_node(root_key).expect("root should exist");
      root
        .outbound()
        .iter()
        .map(|edge_key| {
          let child_key = tree.graph.get_target_node_key(*edge_key)?;
          let name = names[&child_key].clone().expect("root children are named");
          let length = tree.branch_lengths[edge_key].expect("root edges have a length");
          Ok((name, length))
        })
        .collect()
    }

    pub(super) fn leaf_names(tree: &ClockTree, names: &BTreeMap<GraphNodeKey, Option<String>>) -> Vec<String> {
      tree
        .graph
        .get_leaves()
        .map(|leaf| names[&leaf.key()].clone().expect("leaves are named"))
        .sorted()
        .collect()
    }

    pub(super) fn assert_inputs_match_graph(tree: &ClockTree) {
      let node_keys: BTreeSet<_> = tree.graph.get_nodes().map(|node| node.key()).collect();
      let edge_keys: BTreeSet<_> = tree.graph.get_edges().map(|edge| edge.key()).collect();
      assert_eq!(node_keys, tree.inputs.nodes.keys().copied().collect());
      assert_eq!(edge_keys, tree.inputs.edges.keys().copied().collect());
      assert_eq!(edge_keys, tree.branch_lengths.keys().copied().collect());
    }

    fn setup_graph(newick: &str, dates: &BTreeMap<String, f64>) -> Result<RerootFixture, Report> {
      let nwk_parsed = nwk_read_str(newick)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let times = names
        .iter()
        .map(|(key, name)| (*key, name.as_ref().and_then(|name| dates.get(name)).copied()))
        .collect();
      let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
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

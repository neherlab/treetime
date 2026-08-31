#[cfg(test)]
mod tests {
  use crate::partition::timetree::partition::{GraphTimetree, PartitionTimetreeRef};
  use crate::payload::clock_set::ClockSet;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::optimization::polytomy::{prepare_tree_after_topology_change, resolve_polytomies};
  use eyre::Report;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rand::RngCore;
  use std::collections::BTreeSet;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::assign_node_names::assign_node_names;
  use treetime_graph::edge::HasBranchLength;
  use treetime_graph::node::{GraphNodeKey, Named};
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::sync::random::get_random_number_generator;
  use treetime_utils::{make_report, pretty_assert_abs_diff_eq};

  /// Expected substitutions per unit time over the whole alignment.
  const TEST_MUTATION_RATE: f64 = 0.1;

  /// Per-branch coalescent merger rate. Chosen so that the fixtures' time windows comfortably
  /// admit the mergers a full resolution needs, without making them a foregone conclusion.
  const TEST_MERGER_RATE: f64 = 0.15;

  fn set_time(graph: &GraphTimetree, name: &str, time: f64) -> Result<GraphNodeKey, Report> {
    let key = find_node_key_by_name(graph, name).ok_or_else(|| make_report!("{name} not found"))?;
    graph
      .get_node(key)
      .expect("Node must exist")
      .write_arc()
      .payload()
      .write_arc()
      .time = Some(time);
    Ok(key)
  }

  /// `((A,B,C)ABC)root` with a wide time window above the polytomy.
  fn polytomy_tree() -> Result<GraphTimetree, Report> {
    let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2,C:0.15)ABC:0.05)root;")?;
    for (name, time) in [
      ("A", 2020.0),
      ("B", 2015.0),
      ("C", 2018.0),
      ("ABC", 1990.0),
      ("root", 1980.0),
    ] {
      set_time(&graph, name, time)?;
    }
    for edge in graph.get_edges() {
      let mut payload = edge.read_arc().payload().write_arc();
      payload.time_length = payload.branch_length();
      payload.set_branch_length(Some(0.0));
    }
    Ok(graph)
  }

  /// A 6-way polytomy, closer to what the sweep is meant for.
  fn wide_polytomy_tree() -> Result<GraphTimetree, Report> {
    let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.1,C:0.1,D:0.1,E:0.1,F:0.1)P:0.05)root;")?;
    for (name, time) in [
      ("A", 2020.0),
      ("B", 2019.0),
      ("C", 2018.0),
      ("D", 2017.0),
      ("E", 2016.0),
      ("F", 2015.0),
      ("P", 1980.0),
      ("root", 1970.0),
    ] {
      set_time(&graph, name, time)?;
    }
    for edge in graph.get_edges() {
      let mut payload = edge.read_arc().payload().write_arc();
      payload.time_length = payload.branch_length();
      payload.set_branch_length(Some(0.0));
    }
    Ok(graph)
  }

  fn binary_tree() -> Result<GraphTimetree, Report> {
    let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2)AB:0.05,(C:0.15,D:0.1)CD:0.08)root;")?;
    for (name, time) in [
      ("A", 2020.0),
      ("B", 2015.0),
      ("C", 2018.0),
      ("D", 2012.0),
      ("AB", 2000.0),
      ("CD", 2000.0),
      ("root", 1990.0),
    ] {
      set_time(&graph, name, time)?;
    }
    Ok(graph)
  }

  fn no_partitions() -> Vec<PartitionTimetreeRef> {
    vec![]
  }

  /// Resolve under a constant merger rate.
  ///
  /// These fixtures carry no partitions, so every branch has zero reconstructed substitutions
  /// and the total alignment length never enters: the sampled history is shaped by the merger
  /// rate alone.
  fn resolve(graph: &mut GraphTimetree, rng: &mut dyn RngCore) -> Result<usize, Report> {
    let merger_rate = PiecewiseConstantFn::new(array![], array![TEST_MERGER_RATE]);
    resolve_polytomies(graph, &no_partitions(), TEST_MUTATION_RATE, 0, &merger_rate, rng)
  }

  /// Names of the leaves reachable from `node_key`.
  fn leaf_names_under(graph: &GraphTimetree, node_key: GraphNodeKey) -> BTreeSet<String> {
    let mut names = BTreeSet::new();
    let mut stack = vec![node_key];
    while let Some(key) = stack.pop() {
      let node = graph.get_node(key).expect("Node must exist");
      let node = node.read_arc();
      if node.is_leaf() {
        if let Some(name) = node.payload().read_arc().name() {
          names.insert(name.as_ref().to_owned());
        }
        continue;
      }
      for &edge_key in node.outbound() {
        stack.push(graph.get_edge(edge_key).expect("Edge must exist").read_arc().target());
      }
    }
    names
  }

  #[test]
  fn test_resolve_polytomies_leaves_a_binary_tree_alone() -> Result<(), Report> {
    let mut graph = binary_tree()?;
    let mut rng = get_random_number_generator(Some(1));

    let created = resolve(&mut graph, &mut rng)?;

    assert_eq!(created, 0, "a binary tree has no polytomy to resolve");
    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_resolves_a_three_way_polytomy() -> Result<(), Report> {
    let mut graph = polytomy_tree()?;
    let abc_key = find_node_key_by_name(&graph, "ABC").ok_or_else(|| make_report!("ABC not found"))?;
    let mut rng = get_random_number_generator(Some(11));

    let created = resolve(&mut graph, &mut rng)?;

    assert_eq!(created, 1, "a 3-way polytomy needs one merger to become a bifurcation");
    let degree = graph
      .get_node(abc_key)
      .expect("Node must exist")
      .read_arc()
      .degree_out();
    assert_eq!(degree, 2);
    Ok(())
  }

  proptest! {
    #[test]
    fn test_prop_resolve_polytomies_preserves_every_leaf(seed in any::<u64>()) {
      let mut graph = wide_polytomy_tree().unwrap();
      let parent_key = find_node_key_by_name(&graph, "P").expect("P must exist");
      let before = leaf_names_under(&graph, parent_key);
      let mut rng = get_random_number_generator(Some(seed));

      resolve(&mut graph, &mut rng).unwrap();

      let after = leaf_names_under(&graph, parent_key);
      prop_assert_eq!(before, after);
    }

    #[test]
    fn test_prop_resolve_polytomies_leaves_no_single_child_nodes(seed in any::<u64>()) {
      let mut graph = wide_polytomy_tree().unwrap();
      let mut rng = get_random_number_generator(Some(seed));
      resolve(&mut graph, &mut rng).unwrap();

      let has_single_child_node = graph.get_nodes().into_iter().any(|node| {
        let node = node.read_arc();
        node.inbound().len() == 1 && node.outbound().len() == 1
      });
      prop_assert!(!has_single_child_node);
    }
  }

  #[test]
  fn test_resolve_polytomies_is_reproducible_under_the_same_seed() -> Result<(), Report> {
    // Compare the resolved topology by the set of leaf-name clusters it induces, which is
    // independent of node keys and traversal order.
    let clusters = |seed: u64| -> Result<BTreeSet<Vec<String>>, Report> {
      let mut graph = wide_polytomy_tree()?;
      let mut rng = get_random_number_generator(Some(seed));
      resolve(&mut graph, &mut rng)?;
      Ok(
        graph
          .get_nodes()
          .into_iter()
          .map(|node| {
            let key = node.read_arc().key();
            leaf_names_under(&graph, key).into_iter().collect::<Vec<_>>()
          })
          .collect(),
      )
    };

    assert_eq!(
      clusters(4)?,
      clusters(4)?,
      "the same seed must produce the same topology"
    );
    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_different_seeds_can_differ() -> Result<(), Report> {
    let clusters = |seed: u64| -> Result<BTreeSet<Vec<String>>, Report> {
      let mut graph = wide_polytomy_tree()?;
      let mut rng = get_random_number_generator(Some(seed));
      resolve(&mut graph, &mut rng)?;
      Ok(
        graph
          .get_nodes()
          .into_iter()
          .map(|node| {
            let key = node.read_arc().key();
            leaf_names_under(&graph, key).into_iter().collect::<Vec<_>>()
          })
          .collect(),
      )
    };

    let distinct: BTreeSet<BTreeSet<Vec<String>>> = (0..20).map(clusters).collect::<Result<_, _>>()?;
    assert!(
      distinct.len() > 1,
      "resolution is stochastic, so 20 seeds should not all agree"
    );
    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_without_a_time_window_is_a_noop() -> Result<(), Report> {
    // The polytomy sits at the same time as its children, so no merger fits above them.
    let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2,C:0.15)ABC:0.05)root;")?;
    for (name, time) in [
      ("A", 2010.0),
      ("B", 2010.0),
      ("C", 2010.0),
      ("ABC", 2010.0),
      ("root", 2000.0),
    ] {
      set_time(&graph, name, time)?;
    }
    let mut graph = graph;
    let mut rng = get_random_number_generator(Some(1));

    let created = resolve(&mut graph, &mut rng)?;

    assert_eq!(created, 0, "no window above the polytomy means no resolution");
    let abc_key = find_node_key_by_name(&graph, "ABC").ok_or_else(|| make_report!("ABC not found"))?;
    let degree = graph
      .get_node(abc_key)
      .expect("Node must exist")
      .read_arc()
      .degree_out();
    assert_eq!(degree, 3, "the multifurcation must survive intact");
    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_dates_new_nodes_between_parent_and_children() -> Result<(), Report> {
    let mut graph = wide_polytomy_tree()?;
    let parent_key = find_node_key_by_name(&graph, "P").ok_or_else(|| make_report!("P not found"))?;
    let parent_time = 1980.0;
    let mut rng = get_random_number_generator(Some(9));

    resolve(&mut graph, &mut rng)?;

    for node in graph.get_nodes() {
      let node = node.read_arc();
      if node.is_leaf() || node.payload().read_arc().name().is_some() {
        continue;
      }
      let time = node.payload().read_arc().time.expect("new nodes must be dated");
      assert!(
        time > parent_time,
        "new node at {time} must be more recent than the polytomy at {parent_time}"
      );
    }

    // Every edge in the resolved subtree runs forward in time.
    let mut stack = vec![parent_key];
    while let Some(key) = stack.pop() {
      let node = graph.get_node(key).expect("Node must exist");
      let node = node.read_arc();
      let time = node.payload().read_arc().time.expect("node must be dated");
      for &edge_key in node.outbound() {
        let edge = graph.get_edge(edge_key).expect("Edge must exist");
        let target = edge.read_arc().target();
        let time_length = edge.read_arc().payload().read_arc().time_length;
        let child_time = graph
          .get_node(target)
          .expect("Node must exist")
          .read_arc()
          .payload()
          .read_arc()
          .time
          .expect("node must be dated");
        assert!(child_time > time, "edge must run forward in time");
        pretty_assert_abs_diff_eq!(
          time_length.expect("time_length must be set"),
          child_time - time,
          epsilon = 1e-9
        );
        stack.push(target);
      }
    }

    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_names_new_nodes() -> Result<(), Report> {
    let mut graph = polytomy_tree()?;
    let mut rng = get_random_number_generator(Some(11));

    let created = resolve(&mut graph, &mut rng)?;
    assert_eq!(created, 1);

    assign_node_names(&graph)?;

    let mut names: Vec<String> = graph
      .get_nodes()
      .iter()
      .filter_map(|node| {
        node
          .read_arc()
          .payload()
          .read_arc()
          .name()
          .map(|n| n.as_ref().to_owned())
      })
      .collect();
    names.sort();

    assert_eq!(names, vec!["A", "ABC", "B", "C", "NODE_0000000", "root"]);
    Ok(())
  }

  #[test]
  fn test_prepare_tree_after_topology_change_resets_derived_state_preserves_inputs() -> Result<(), Report> {
    let graph = polytomy_tree()?;

    for edge in graph.get_edges() {
      let mut payload = edge.read_arc().payload().write_arc();
      payload.msg_to_parent = Some(Arc::new(Distribution::point(1.0, 1.0)));
      payload.gamma = 2.0;
      payload.clock_to_parent = ClockSet::leaf_contribution(Some(2020.0));
      payload.clock_to_child = ClockSet::leaf_contribution(Some(2021.0));
      payload.clock_from_child = ClockSet::leaf_contribution(Some(2022.0));
      payload.set_branch_length(Some(0.25));
      payload.time_length = Some(3.0);
    }

    let date_dist = Arc::new(Distribution::point(2020.0, 1.0));
    for name in ["A", "B", "C"] {
      let key = find_node_key_by_name(&graph, name).ok_or_else(|| make_report!("{name} not found"))?;
      graph
        .get_node(key)
        .expect("Node must exist")
        .read_arc()
        .payload()
        .write_arc()
        .time_distribution = Some(Arc::clone(&date_dist));
    }

    let leaf_b_key = find_node_key_by_name(&graph, "B").ok_or_else(|| make_report!("B not found"))?;
    graph
      .get_node(leaf_b_key)
      .expect("Node must exist")
      .read_arc()
      .payload()
      .write_arc()
      .bad_branch = true;

    let abc_key = find_node_key_by_name(&graph, "ABC").ok_or_else(|| make_report!("ABC not found"))?;
    {
      let node = graph.get_node(abc_key).expect("Node must exist");
      let mut payload = node.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(2010.0, 1.0)));
      payload.bad_branch = true;
    }

    prepare_tree_after_topology_change(&graph)?;

    for name in ["A", "B", "C"] {
      let key = find_node_key_by_name(&graph, name).ok_or_else(|| make_report!("{name} not found"))?;
      let node = graph.get_node(key).expect("Node must exist");
      let payload = node.read_arc().payload().read_arc();
      assert!(
        payload.time_distribution.is_some(),
        "leaf time_distribution must survive topology change"
      );
    }

    {
      let node = graph.get_node(leaf_b_key).expect("Node must exist");
      assert!(
        node.read_arc().payload().read_arc().bad_branch,
        "leaf bad_branch flag must survive topology change"
      );
    }

    {
      let node = graph.get_node(abc_key).expect("Node must exist");
      let payload = node.read_arc().payload().read_arc();
      pretty_assert_abs_diff_eq!(
        payload
          .time_distribution
          .as_ref()
          .and_then(|distribution| distribution.likely_time())
          .expect("internal node time distribution must be rebuilt"),
        1990.0,
        epsilon = 1e-10
      );
      assert!(payload.bad_branch, "internal node bad_branch must be preserved");
    }

    for edge in graph.get_edges() {
      let payload = edge.read_arc().payload().read_arc();
      assert!(
        payload.branch_length_distribution.is_none(),
        "edge branch_length_distribution must be cleared"
      );
      assert!(payload.msg_to_parent.is_none(), "edge msg_to_parent must be cleared");
      pretty_assert_abs_diff_eq!(payload.gamma, 1.0, epsilon = 1e-10);
      pretty_assert_abs_diff_eq!(
        payload.branch_length().expect("branch length must be preserved"),
        0.25,
        epsilon = 1e-10
      );
      pretty_assert_abs_diff_eq!(
        payload.time_length.expect("time length must be preserved"),
        3.0,
        epsilon = 1e-10
      );
    }

    Ok(())
  }
}

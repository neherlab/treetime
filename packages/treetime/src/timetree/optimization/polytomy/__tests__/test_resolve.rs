#[cfg(test)]
mod tests {
  use crate::test_utils::{find_node_key_by_name, point_date_constraints};
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use crate::timetree::inference::result::NodeTimes;
  use crate::timetree::optimization::polytomy::resolve::{
    PolytomyResolution, require_internal_node_times, resolve_polytomies,
  };
  use eyre::Report;
  use maplit::btreeset;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rand::RngCore;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_graph::assign_node_names::assign_node_names;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::sync::random::get_random_number_generator;
  use treetime_utils::{assert_error, make_report};

  const TEST_MUTATION_RATE: f64 = 0.1;

  const TEST_MERGER_RATE: f64 = 0.15;

  const MERGER_FLAG_SEEDS: u64 = 20;

  #[test]
  fn test_resolve_polytomies_leaves_a_binary_tree_alone() -> Result<(), Report> {
    let (graph, names, mut node_times, branch_lengths) = binary_tree()?;
    let mut rng = get_random_number_generator(1);

    let (graph, created) = resolve(graph, branch_lengths, &mut node_times, &mut rng)?;

    assert_eq!(created, 0, "a binary tree has no polytomy to resolve");
    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_resolves_a_three_way_polytomy() -> Result<(), Report> {
    let (graph, names, mut node_times, branch_lengths) = polytomy_tree()?;
    let abc_key = find_node_key_by_name(&graph, &names, "ABC").ok_or_else(|| make_report!("ABC not found"))?;
    let mut rng = get_random_number_generator(11);

    let (graph, created) = resolve(graph, branch_lengths, &mut node_times, &mut rng)?;

    assert_eq!(created, 1, "a 3-way polytomy needs one merger to become a bifurcation");
    let degree = graph.get_node(abc_key).expect("Node must exist").degree_out();
    assert_eq!(degree, 2);
    Ok(())
  }

  proptest! {
    #[test]
    fn test_prop_resolve_polytomies_preserves_every_leaf(seed in any::<u64>()) {
      let (graph, names, mut node_times, branch_lengths) = wide_polytomy_tree().unwrap();
      let parent_key = find_node_key_by_name(&graph, &names, "P").expect("P must exist");
      let before = leaf_names_under(&graph, &names, parent_key);
      let mut rng = get_random_number_generator(seed);

      let (graph, _) = resolve(graph, branch_lengths, &mut node_times, &mut rng).unwrap();

      let after = leaf_names_under(&graph, &names, parent_key);
      prop_assert_eq!(before, after);
    }

    #[test]
    fn test_prop_resolve_polytomies_leaves_no_single_child_nodes(seed in any::<u64>()) {
      let (graph, names, mut node_times, branch_lengths) = wide_polytomy_tree().unwrap();
      let mut rng = get_random_number_generator(seed);
      let (graph, _) = resolve(graph, branch_lengths, &mut node_times, &mut rng).unwrap();

      let has_single_child_node = graph.get_nodes().into_iter().any(|node| {
        node.inbound().len() == 1 && node.outbound().len() == 1
      });
      prop_assert!(!has_single_child_node);
    }
  }

  #[test]
  fn test_resolve_polytomies_is_reproducible_under_the_same_seed() -> Result<(), Report> {
    let clusters = |seed: u64| -> Result<BTreeSet<Vec<String>>, Report> {
      let (graph, names, mut node_times, branch_lengths) = wide_polytomy_tree()?;
      let mut rng = get_random_number_generator(seed);
      let (graph, _) = resolve(graph, branch_lengths, &mut node_times, &mut rng)?;
      Ok(
        graph
          .get_nodes()
          .map(|node| {
            let key = node.key();
            leaf_names_under(&graph, &names, key).into_iter().collect::<Vec<_>>()
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
      let (graph, names, mut node_times, branch_lengths) = wide_polytomy_tree()?;
      let mut rng = get_random_number_generator(seed);
      let (graph, _) = resolve(graph, branch_lengths, &mut node_times, &mut rng)?;
      Ok(
        graph
          .get_nodes()
          .map(|node| {
            let key = node.key();
            leaf_names_under(&graph, &names, key).into_iter().collect::<Vec<_>>()
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
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2,C:0.15)ABC:0.05)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut node_times = undated(&graph);
    for (name, time) in [
      ("A", 2010.0),
      ("B", 2010.0),
      ("C", 2010.0),
      ("ABC", 2010.0),
      ("root", 2000.0),
    ] {
      set_time(&graph, &names, &mut node_times, name, time)?;
    }
    let mut rng = get_random_number_generator(1);

    let (graph, created) = resolve(graph, branch_lengths, &mut node_times, &mut rng)?;

    assert_eq!(created, 0, "no window above the polytomy means no resolution");
    let abc_key = find_node_key_by_name(&graph, &names, "ABC").ok_or_else(|| make_report!("ABC not found"))?;
    let degree = graph.get_node(abc_key).expect("Node must exist").degree_out();
    assert_eq!(degree, 3, "the multifurcation must survive intact");
    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_dates_new_nodes_between_parent_and_children() -> Result<(), Report> {
    let (graph, names, mut node_times, branch_lengths) = wide_polytomy_tree()?;
    let parent_key = find_node_key_by_name(&graph, &names, "P").ok_or_else(|| make_report!("P not found"))?;
    let parent_time = 1980.0;
    let mut rng = get_random_number_generator(9);

    let (graph, _) = resolve(graph, branch_lengths, &mut node_times, &mut rng)?;

    for node in graph.get_nodes() {
      if node.is_leaf() || names.get(&node.key()).and_then(|x| x.as_ref()).is_some() {
        continue;
      }
      let time = node_times[&node.key()].expect("new nodes must be dated");
      assert!(
        time > parent_time,
        "new node at {time} must be more recent than the polytomy at {parent_time}"
      );
    }

    let mut stack = vec![parent_key];
    while let Some(key) = stack.pop() {
      let node = graph.get_node(key).expect("Node must exist");
      let time = node_times[&key].expect("node must be dated");
      for &edge_key in node.outbound() {
        let edge = graph.get_edge(edge_key).expect("Edge must exist");
        let target = edge.target();
        let child_time = node_times[&target].expect("node must be dated");
        assert!(child_time > time, "edge must run forward in time");
        stack.push(target);
      }
    }

    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_names_new_nodes() -> Result<(), Report> {
    let (graph, names, mut node_times, branch_lengths) = polytomy_tree()?;
    let mut rng = get_random_number_generator(11);

    let (graph, created) = resolve(graph, branch_lengths, &mut node_times, &mut rng)?;
    assert_eq!(created, 1);

    let names = assign_node_names(names, &graph)?;

    let mut name_list: Vec<String> = names.values().filter_map(Clone::clone).collect();
    name_list.sort();

    assert_eq!(name_list, vec!["A", "ABC", "B", "C", "NODE_0000000", "root"]);
    Ok(())
  }

  #[test]
  fn test_require_internal_node_times_accepts_a_fully_dated_tree() -> Result<(), Report> {
    let (graph, _names, node_times, _branch_lengths) = polytomy_tree()?;
    require_internal_node_times(&graph, &node_times)
  }

  #[test]
  fn test_require_internal_node_times_rejects_an_undated_internal_node() -> Result<(), Report> {
    let (graph, names, mut node_times, _branch_lengths) = polytomy_tree()?;
    let abc_key = find_node_key_by_name(&graph, &names, "ABC").ok_or_else(|| make_report!("ABC not found"))?;
    node_times.insert(abc_key, None);

    assert_error!(
      require_internal_node_times(&graph, &node_times),
      format!("Polytomy resolution requires an inferred time for every internal node, but node {abc_key:?} has none")
    );
    Ok(())
  }

  #[test]
  fn test_require_internal_node_times_ignores_undated_leaves() -> Result<(), Report> {
    let (graph, names, mut node_times, _branch_lengths) = polytomy_tree()?;
    let leaf_key = find_node_key_by_name(&graph, &names, "B").ok_or_else(|| make_report!("B not found"))?;
    node_times.insert(leaf_key, None);

    require_internal_node_times(&graph, &node_times)
  }

  #[test]
  fn test_resolve_polytomies_merger_is_a_bad_branch_iff_every_leaf_below_it_is_undated() -> Result<(), Report> {
    let undated_leaves: BTreeSet<String> = ["A", "B", "C"].map(str::to_owned).into();
    let mut observed_flags = BTreeSet::new();
    for seed in 0..MERGER_FLAG_SEEDS {
      let (graph, names, mut node_times, branch_lengths) = wide_polytomy_tree()?;
      let constraints = point_date_constraints(&graph, &names, &[("D", 2017.0), ("E", 2016.0), ("F", 2015.0)]);
      let leaf_bad_branches = bad_leaves(&graph, &constraints, &BTreeSet::new());
      let mut rng = get_random_number_generator(seed);
      let (graph, _) = resolve(graph, branch_lengths, &mut node_times, &mut rng)?;

      let bad_branches = derive_bad_branches(&graph, &constraints, &leaf_bad_branches)?;

      for merger in graph.get_nodes().filter(|node| !names.contains_key(&node.key())) {
        let expected = leaf_names_under(&graph, &names, merger.key()).is_subset(&undated_leaves);
        assert_eq!(
          expected,
          bad_branches[&merger.key()],
          "merger {:?} under seed {seed}",
          merger.key()
        );
        observed_flags.insert(expected);
      }
    }
    assert_eq!(
      btreeset! { false, true },
      observed_flags,
      "the seeds must produce both a bad and a good merger"
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn set_time(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      times: &mut NodeTimes,
      name: &str,
      time: f64,
    ) -> Result<GraphNodeKey, Report> {
      let key = find_node_key_by_name(graph, names, name).ok_or_else(|| make_report!("{name} not found"))?;
      times.insert(key, Some(time));
      Ok(key)
    }

    pub(super) fn undated(graph: &Graph) -> NodeTimes {
      graph.get_nodes().map(|node| (node.key(), None)).collect()
    }

    pub(super) fn polytomy_tree() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        NodeTimes,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str("((A:0.1,B:0.2,C:0.15)ABC:0.05)root;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let mut branch_lengths = nwk_parsed.branch_lengths;
      let mut node_times = undated(&graph);
      for (name, time) in [
        ("A", 2020.0),
        ("B", 2015.0),
        ("C", 2018.0),
        ("ABC", 1990.0),
        ("root", 1980.0),
      ] {
        set_time(&graph, &names, &mut node_times, name, time)?;
      }
      for edge in graph.get_edges() {
        branch_lengths.insert(edge.key(), Some(0.0));
      }
      Ok((graph, names, node_times, branch_lengths))
    }

    pub(super) fn wide_polytomy_tree() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        NodeTimes,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str("((A:0.1,B:0.1,C:0.1,D:0.1,E:0.1,F:0.1)P:0.05)root;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let mut branch_lengths = nwk_parsed.branch_lengths;
      let mut node_times = undated(&graph);
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
        set_time(&graph, &names, &mut node_times, name, time)?;
      }
      for edge in graph.get_edges() {
        branch_lengths.insert(edge.key(), Some(0.0));
      }
      Ok((graph, names, node_times, branch_lengths))
    }

    pub(super) fn binary_tree() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        NodeTimes,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.05,(C:0.15,D:0.1)CD:0.08)root;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let mut node_times = undated(&graph);
      for (name, time) in [
        ("A", 2020.0),
        ("B", 2015.0),
        ("C", 2018.0),
        ("D", 2012.0),
        ("AB", 2000.0),
        ("CD", 2000.0),
        ("root", 1990.0),
      ] {
        set_time(&graph, &names, &mut node_times, name, time)?;
      }
      Ok((graph, names, node_times, branch_lengths))
    }

    pub(super) fn resolve(
      graph: Graph,
      branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
      times: &mut NodeTimes,
      rng: &mut dyn RngCore,
    ) -> Result<(Graph, usize), Report> {
      let merger_rate = PiecewiseConstantFn::new(array![], array![TEST_MERGER_RATE]);
      let PolytomyResolution {
        graph, merger_times, ..
      } = resolve_polytomies(
        graph,
        branch_lengths,
        &BranchModel::Input,
        TEST_MUTATION_RATE,
        0,
        &merger_rate,
        rng,
        times,
      )?;
      let created = merger_times.len();
      times.extend(merger_times.into_iter().map(|(key, time)| (key, Some(time))));
      Ok((graph, created))
    }

    pub(super) fn leaf_names_under(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      node_key: GraphNodeKey,
    ) -> BTreeSet<String> {
      let mut result = BTreeSet::new();
      let mut stack = vec![node_key];
      while let Some(key) = stack.pop() {
        let node = graph.get_node(key).expect("Node must exist");
        if node.is_leaf() {
          if let Some(name) = names.get(&node.key()).cloned().flatten() {
            result.insert(name);
          }
          continue;
        }
        for &edge_key in node.outbound() {
          stack.push(graph.get_edge(edge_key).expect("Edge must exist").target());
        }
      }
      result
    }
  }

  use helpers::*;
}

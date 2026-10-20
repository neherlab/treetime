#[cfg(test)]
mod tests {
  use self::helpers::{
    BRAID_LEVELS, BRAID_WIDTH, braid, braid_levels, caterpillar_parents, inputs_of, run_in_pool,
    run_in_small_stack_pool, subtree_sizes,
  };
  use super::super::test_graph_support::tests::helpers::{NamedGraph, graph_from_dag, graph_from_parents};
  use crate::edge::GraphEdgeKey;
  use crate::graph::Graph;
  use crate::node::GraphNodeKey;
  use crate::schedule::GraphSchedule;
  use crate::schedule_scan::{NodeEdgeMaps, PerEdge, ScanOutputs};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use smallvec::smallvec;
  use std::collections::{BTreeMap, BTreeSet};
  use std::sync::Arc;
  use treetime_utils::{assert_error, make_report};

  #[test]
  fn test_schedule_scan_up_sums_subtrees() -> Result<(), Report> {
    let tree = NamedGraph::new(
      &["root", "AB", "A", "B", "C"],
      &[("root", "AB"), ("root", "C"), ("AB", "A"), ("AB", "B")],
    )?;
    let schedule = GraphSchedule::new(&tree.graph)?;
    let weights = btreemap! {
      tree.key("root") => 100, tree.key("AB") => 10, tree.key("A") => 1, tree.key("B") => 2, tree.key("C") => 3,
    };

    let sums = schedule
      .scan_up(weights, |_, weight, children: PerEdge<usize>| {
        let sum = weight + children.iter().sum::<usize>();
        Ok((sum, sum))
      })?
      .into_map();

    let expected = btreemap! {
      tree.key("root") => 116, tree.key("AB") => 13, tree.key("A") => 1, tree.key("B") => 2, tree.key("C") => 3,
    };
    assert_eq!(expected, sums);
    Ok(())
  }

  #[test]
  fn test_schedule_scan_up_delivers_children_in_tree_order() -> Result<(), Report> {
    let tree = NamedGraph::new(
      &["root", "AB", "A", "B", "C"],
      &[("root", "C"), ("root", "AB"), ("AB", "B"), ("AB", "A")],
    )?;
    let schedule = GraphSchedule::new(&tree.graph)?;

    let received = schedule
      .scan_up(
        inputs_of(&tree.graph, ()),
        |key, (), children: PerEdge<&'static str>| Ok((children.into_vec(), tree.name(key))),
      )?
      .into_map();

    let expected = btreemap! {
      tree.key("root") => vec!["C", "AB"], tree.key("AB") => vec!["B", "A"],
      tree.key("A") => vec![], tree.key("B") => vec![], tree.key("C") => vec![],
    };
    assert_eq!(expected, received);
    Ok(())
  }

  #[test]
  fn test_schedule_scan_down_accumulates_from_the_root() -> Result<(), Report> {
    let tree = NamedGraph::new(
      &["root", "AB", "A", "B", "C"],
      &[("root", "AB"), ("root", "C"), ("AB", "A"), ("AB", "B")],
    )?;
    let schedule = GraphSchedule::new(&tree.graph)?;
    let weights = btreemap! {
      tree.key("root") => 100, tree.key("AB") => 10, tree.key("A") => 1, tree.key("B") => 2, tree.key("C") => 3,
    };

    let sums = schedule
      .scan_down(weights, |_, weight, parents: PerEdge<usize>| {
        let sum = weight + parents.iter().sum::<usize>();
        Ok((sum, sum))
      })?
      .into_map();

    let expected = btreemap! {
      tree.key("root") => 100, tree.key("AB") => 110, tree.key("A") => 111, tree.key("B") => 112, tree.key("C") => 103,
    };
    assert_eq!(expected, sums);
    Ok(())
  }

  #[test]
  fn test_schedule_scan_down_single_node_receives_no_message() -> Result<(), Report> {
    let tree = NamedGraph::new(&["root"], &[])?;
    let schedule = GraphSchedule::new(&tree.graph)?;

    let received = schedule
      .scan_down(inputs_of(&tree.graph, ()), |_, (), parents: PerEdge<u8>| {
        Ok((parents.into_vec(), 0))
      })?
      .into_map();

    assert_eq!(btreemap! { tree.key("root") => vec![] }, received);
    Ok(())
  }

  #[test]
  fn test_schedule_scan_up_single_node_receives_no_messages() -> Result<(), Report> {
    let tree = NamedGraph::new(&["root"], &[])?;
    let schedule = GraphSchedule::new(&tree.graph)?;

    let received = schedule
      .scan_up(inputs_of(&tree.graph, ()), |_, (), children: PerEdge<u8>| {
        Ok((children.into_vec(), 0))
      })?
      .into_map();

    assert_eq!(btreemap! { tree.key("root") => vec![] }, received);
    Ok(())
  }

  #[rstest]
  #[case::one_thread(1)]
  #[case::four_threads(4)]
  #[trace]
  fn test_schedule_scan_up_star_delivers_every_leaf_in_order(#[case] threads: usize) -> Result<(), Report> {
    let graph = graph_from_parents(&vec![0; 1000])?;
    let schedule = GraphSchedule::new(&graph)?;

    let received = run_in_pool(threads, || {
      schedule
        .scan_up(inputs_of(&graph, ()), |key, (), children: PerEdge<GraphNodeKey>| {
          Ok((children.into_vec(), key))
        })
        .map(ScanOutputs::into_map)
    })?;

    let leaves = schedule
      .children(schedule.roots()[0])
      .iter()
      .map(|&(child, _)| child)
      .collect_vec();
    assert_eq!(999, leaves.len());
    assert_eq!(leaves, received[&schedule.roots()[0]]);
    Ok(())
  }

  #[rstest]
  #[case::one_thread(1)]
  #[case::four_threads(4)]
  #[trace]
  fn test_schedule_scan_down_star_sends_the_root_message_to_every_leaf(#[case] threads: usize) -> Result<(), Report> {
    let graph = graph_from_parents(&vec![0; 1000])?;
    let schedule = GraphSchedule::new(&graph)?;

    let received = run_in_pool(threads, || {
      schedule
        .scan_down(inputs_of(&graph, ()), |key, (), parents: PerEdge<GraphNodeKey>| {
          Ok((parents.into_vec(), key))
        })
        .map(ScanOutputs::into_map)
    })?;

    let root = schedule.roots()[0];
    let expected = graph
      .get_nodes()
      .map(|node| (node.key(), if node.key() == root { vec![] } else { vec![root] }))
      .collect::<BTreeMap<_, _>>();
    assert_eq!(expected, received);
    Ok(())
  }

  #[rstest]
  #[case::leaf_first_one_thread(false, 1)]
  #[case::leaf_first_sixteen_threads(false, 16)]
  #[case::spine_first_one_thread(true, 1)]
  #[case::spine_first_sixteen_threads(true, 16)]
  #[trace]
  fn test_schedule_scan_up_caterpillar_of_100000_nodes_fits_small_worker_stacks(
    #[case] spine_first: bool,
    #[case] threads: usize,
  ) -> Result<(), Report> {
    let graph = graph_from_parents(&caterpillar_parents(100_000, spine_first))?;
    let schedule = GraphSchedule::new(&graph)?;

    let sizes = run_in_small_stack_pool(threads, || {
      schedule
        .scan_up(inputs_of(&graph, ()), |_, (), children: PerEdge<usize>| {
          let size = 1 + children.iter().sum::<usize>();
          Ok((size, size))
        })
        .map(ScanOutputs::into_map)
    })?;

    assert_eq!(subtree_sizes(&schedule), sizes);
    assert_eq!(100_000, sizes[&schedule.roots()[0]]);
    Ok(())
  }

  #[rstest]
  #[case::leaf_first_one_thread(false, 1)]
  #[case::leaf_first_sixteen_threads(false, 16)]
  #[case::spine_first_one_thread(true, 1)]
  #[case::spine_first_sixteen_threads(true, 16)]
  #[trace]
  fn test_schedule_scan_down_caterpillar_of_100000_nodes_fits_small_worker_stacks(
    #[case] spine_first: bool,
    #[case] threads: usize,
  ) -> Result<(), Report> {
    let graph = graph_from_parents(&caterpillar_parents(100_000, spine_first))?;
    let schedule = GraphSchedule::new(&graph)?;

    let depths = run_in_small_stack_pool(threads, || {
      schedule
        .scan_down(inputs_of(&graph, ()), |_, (), parents: PerEdge<usize>| {
          let depth = parents.first().map_or(0, |depth| depth + 1);
          Ok((depth, depth))
        })
        .map(ScanOutputs::into_map)
    })?;

    let max_depth = depths.values().copied().max();
    assert_eq!((100_000, Some(50_000)), (depths.len(), max_depth));
    Ok(())
  }

  #[test]
  fn test_schedule_scan_up_braid_of_100000_nodes_runs_every_node_after_its_children_on_small_worker_stacks()
  -> Result<(), Report> {
    let graph = graph_from_dag(&braid(BRAID_WIDTH, BRAID_LEVELS))?;
    let schedule = GraphSchedule::new(&graph)?;

    let heights = |threads| {
      run_in_small_stack_pool(threads, || {
        schedule
          .scan_up(
            inputs_of(&graph, ()),
            |key, (), children: PerEdge<(GraphNodeKey, usize)>| {
              let senders = children.iter().map(|&(sender, _)| sender).collect_vec();
              let expected = schedule.children(key).iter().map(|&(child, _)| child).collect_vec();
              if senders != expected {
                return Err(make_report!("node {key} received {senders:?} instead of {expected:?}"));
              }
              let height = children.iter().map(|&(_, height)| height + 1).max().unwrap_or(0);
              Ok((height, (key, height)))
            },
          )
          .map(ScanOutputs::into_map)
      })
    };

    let expected = braid_levels()
      .into_iter()
      .map(|(key, level)| (key, BRAID_LEVELS - 1 - level))
      .collect::<BTreeMap<_, _>>();
    assert_eq!((expected.clone(), expected), (heights(1)?, heights(16)?));
    Ok(())
  }

  #[test]
  fn test_schedule_scan_down_braid_of_100000_nodes_runs_every_node_after_its_parents_on_small_worker_stacks()
  -> Result<(), Report> {
    let graph = graph_from_dag(&braid(BRAID_WIDTH, BRAID_LEVELS))?;
    let schedule = GraphSchedule::new(&graph)?;

    let depths = |threads| {
      run_in_small_stack_pool(threads, || {
        schedule
          .scan_down(
            inputs_of(&graph, ()),
            |key, (), parents: PerEdge<(GraphNodeKey, usize)>| {
              let senders = parents.iter().map(|&(sender, _)| sender).collect_vec();
              let expected = schedule.parents(key).iter().map(|&(parent, _)| parent).collect_vec();
              if senders != expected {
                return Err(make_report!("node {key} received {senders:?} instead of {expected:?}"));
              }
              let depth = parents.iter().map(|&(_, depth)| depth + 1).max().unwrap_or(0);
              Ok((depth, (key, depth)))
            },
          )
          .map(ScanOutputs::into_map)
      })
    };

    let expected = braid_levels();
    assert_eq!((expected.clone(), expected), (depths(1)?, depths(16)?));
    Ok(())
  }

  #[test]
  fn test_schedule_scan_up_skips_deleted_key_gaps() -> Result<(), Report> {
    let mut graph = Graph::new();
    let root = graph.add_node();
    let removed = graph.add_node();
    let child = graph.add_node();
    graph.add_edge(root, child)?;
    graph.remove_node(removed)?;
    graph.build()?;
    let schedule = GraphSchedule::new(&graph)?;

    let sums = schedule
      .scan_up(
        btreemap! { root => 100, child => 1 },
        |_, weight, children: PerEdge<usize>| {
          let sum = weight + children.iter().sum::<usize>();
          Ok((sum, sum))
        },
      )?
      .into_map();

    assert_eq!(btreemap! { root => 101, child => 1 }, sums);
    Ok(())
  }

  #[test]
  fn test_schedule_scan_rejects_missing_input() -> Result<(), Report> {
    let tree = NamedGraph::new(&["root", "A", "B"], &[("root", "A"), ("root", "B")])?;
    let schedule = GraphSchedule::new(&tree.graph)?;
    let visited = Mutex::new(0_usize);
    let inputs = btreemap! { tree.key("root") => (), tree.key("A") => () };

    let result = schedule.scan_up(inputs, |_, (), _: PerEdge<()>| {
      *visited.lock() += 1;
      Ok(((), ()))
    });

    assert_error!(
      result,
      format!(
        "Scan input is missing node {}. This is an internal error. Please report it to developers.",
        tree.key("B")
      )
    );
    assert_eq!(0, *visited.lock());
    Ok(())
  }

  #[test]
  fn test_schedule_scan_rejects_input_of_unknown_node() -> Result<(), Report> {
    let tree = NamedGraph::new(&["root", "A"], &[("root", "A")])?;
    let schedule = GraphSchedule::new(&tree.graph)?;
    let unknown = GraphNodeKey(7);
    let inputs = btreemap! { tree.key("root") => (), tree.key("A") => (), unknown => () };

    let result = schedule.scan_down(inputs, |_, (), _: PerEdge<()>| Ok(((), ())));

    assert_error!(
      result,
      "Scan input names node 7, which is not in the graph. This is an internal error. Please report it to developers."
    );
    Ok(())
  }

  #[test]
  fn test_schedule_scan_rejects_input_of_deleted_node() -> Result<(), Report> {
    let mut graph = Graph::new();
    let root = graph.add_node();
    let removed = graph.add_node();
    let child = graph.add_node();
    graph.add_edge(root, child)?;
    graph.remove_node(removed)?;
    graph.build()?;
    let schedule = GraphSchedule::new(&graph)?;
    let inputs = btreemap! { root => (), removed => (), child => () };

    let result = schedule.scan_up(inputs, |_, (), _: PerEdge<()>| Ok(((), ())));

    assert_error!(
      result,
      format!(
        "Scan input names node {removed}, which is not in the graph. This is an internal error. Please report it to \
         developers."
      )
    );
    Ok(())
  }

  #[test]
  fn test_schedule_scan_rejects_input_named_twice() -> Result<(), Report> {
    let tree = NamedGraph::new(&["root", "A"], &[("root", "A")])?;
    let schedule = GraphSchedule::new(&tree.graph)?;
    let inputs = [(tree.key("root"), ()), (tree.key("A"), ()), (tree.key("A"), ())];

    let result = schedule.scan_up(inputs, |_, (), _: PerEdge<()>| Ok(((), ())));

    assert_error!(
      result,
      format!(
        "Scan input names node {} twice. This is an internal error. Please report it to developers.",
        tree.key("A")
      )
    );
    Ok(())
  }

  #[test]
  fn test_schedule_scan_up_failed_node_stops_its_ancestors() -> Result<(), Report> {
    let tree = NamedGraph::new(
      &["root", "AB", "A", "B", "C"],
      &[("root", "AB"), ("root", "C"), ("AB", "A"), ("AB", "B")],
    )?;
    let schedule = GraphSchedule::new(&tree.graph)?;
    let visited = Mutex::new(BTreeSet::new());

    let result = schedule.scan_up(inputs_of(&tree.graph, ()), |key, (), _: PerEdge<()>| {
      visited.lock().insert(tree.name(key));
      if tree.name(key) == "B" {
        return Err(make_report!("injected failure of node B"));
      }
      Ok(((), ()))
    });

    assert_error!(result, "injected failure of node B");
    let visited = visited.into_inner();
    assert_eq!(
      (false, false),
      (visited.contains("AB"), visited.contains("root")),
      "the ancestors of the failed node must not run"
    );
    Ok(())
  }

  #[test]
  fn test_schedule_scan_down_failed_node_stops_its_descendants() -> Result<(), Report> {
    let tree = NamedGraph::new(
      &["root", "AB", "A", "B", "C"],
      &[("root", "AB"), ("root", "C"), ("AB", "A"), ("AB", "B")],
    )?;
    let schedule = GraphSchedule::new(&tree.graph)?;
    let visited = Mutex::new(BTreeSet::new());

    let result = schedule.scan_down(inputs_of(&tree.graph, ()), |key, (), _: PerEdge<()>| {
      visited.lock().insert(tree.name(key));
      if tree.name(key) == "AB" {
        return Err(make_report!("injected failure of node AB"));
      }
      Ok(((), ()))
    });

    assert_error!(result, "injected failure of node AB");
    let visited = visited.into_inner();
    assert_eq!(
      (false, false),
      (visited.contains("A"), visited.contains("B")),
      "the children of the failed node must not run"
    );
    Ok(())
  }

  #[test]
  fn test_schedule_scan_up_returns_the_first_error_and_stops_later_visits() -> Result<(), Report> {
    let graph = graph_from_parents(&[0, 0, 0, 0])?;
    let schedule = GraphSchedule::new(&graph)?;
    let visited = Mutex::new(vec![]);

    let result = run_in_pool(1, || {
      schedule.scan_up(inputs_of(&graph, ()), |key, (), _: PerEdge<()>| {
        visited.lock().push(key);
        Err::<((), ()), _>(make_report!("injected failure of node {key}"))
      })
    });

    assert_error!(result, "injected failure of node 1");
    assert_eq!(vec![GraphNodeKey(1)], visited.into_inner());
    Ok(())
  }

  #[test]
  fn test_schedule_scan_up_after_a_failed_scan_succeeds_on_the_same_view() -> Result<(), Report> {
    let tree = NamedGraph::new(
      &["root", "AB", "A", "B", "C"],
      &[("root", "AB"), ("root", "C"), ("AB", "A"), ("AB", "B")],
    )?;
    let schedule = GraphSchedule::new(&tree.graph)?;
    let failed = schedule.scan_up(inputs_of(&tree.graph, 1), |_, _, _: PerEdge<usize>| {
      Err::<(usize, usize), _>(make_report!("injected failure"))
    });
    assert_error!(failed, "injected failure");

    let sizes = schedule
      .scan_up(inputs_of(&tree.graph, 1), |_, weight, children: PerEdge<usize>| {
        let sum = weight + children.iter().sum::<usize>();
        Ok((sum, sum))
      })?
      .into_map();

    let expected = btreemap! {
      tree.key("root") => 5, tree.key("AB") => 3, tree.key("A") => 1, tree.key("B") => 1, tree.key("C") => 1,
    };
    assert_eq!(expected, sizes);
    Ok(())
  }

  #[test]
  fn test_schedule_scan_join_and_split_parent_edges_roundtrip() -> Result<(), Report> {
    let tree = NamedGraph::new(&["root", "A", "B"], &[("root", "A"), ("root", "B")])?;
    let schedule = GraphSchedule::new(&tree.graph)?;
    let (root, a, b) = (tree.key("root"), tree.key("A"), tree.key("B"));
    let (root_a, root_b) = (tree.edge("root", "A")?, tree.edge("root", "B")?);
    let nodes = btreemap! { root => "root", a => "A", b => "B" };
    let edges = btreemap! { root_a => "root-A", root_b => "root-B" };

    let joined = schedule
      .join_parent_edges(nodes.clone(), edges.clone())?
      .collect::<BTreeMap<_, _>>();
    let split = schedule
      .scan_down(joined.clone(), |_, joined, _: PerEdge<()>| Ok((joined, ())))?
      .split_parent_edges()?;

    let expected_joined = btreemap! {
      root => ("root", smallvec![]), a => ("A", smallvec!["root-A"]), b => ("B", smallvec!["root-B"]),
    };
    assert_eq!(expected_joined, joined);
    assert_eq!(NodeEdgeMaps { nodes, edges }, split);
    Ok(())
  }

  #[test]
  fn test_schedule_scan_join_parent_edges_rejects_missing_edge() -> Result<(), Report> {
    let tree = NamedGraph::new(&["root", "A", "B"], &[("root", "A"), ("root", "B")])?;
    let schedule = GraphSchedule::new(&tree.graph)?;
    let nodes = inputs_of(&tree.graph, ());
    let edges = btreemap! { tree.edge("root", "A")? => () };

    let result = schedule.join_parent_edges(nodes, edges);

    assert_error!(
      result,
      format!(
        "Scan input is missing edge {} above node {}. This is an internal error. Please report it to developers.",
        tree.edge("root", "B")?,
        tree.key("B")
      )
    );
    Ok(())
  }

  #[test]
  fn test_schedule_scan_join_parent_edges_rejects_edge_outside_the_tree() -> Result<(), Report> {
    let tree = NamedGraph::new(&["root", "A"], &[("root", "A")])?;
    let schedule = GraphSchedule::new(&tree.graph)?;
    let nodes = inputs_of(&tree.graph, ());
    let edges = btreemap! { tree.edge("root", "A")? => (), GraphEdgeKey(7) => () };

    let result = schedule.join_parent_edges(nodes, edges);

    assert_error!(
      result,
      "Scan input names edge 7, which is not in the graph. This is an internal error. Please report it to developers."
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::value_at_the_root(        "r", 1, "holds 1 parent edge value for 0 parent edges")]
  #[case::no_value_at_a_child(      "a", 0, "holds 0 parent edge values for 1 parent edge")]
  #[case::one_value_for_two_parents("s", 1, "holds 1 parent edge value for 2 parent edges")]
  #[trace]
  fn test_schedule_scan_split_parent_edges_rejects_wrong_value_count(
    #[case] wrong: &str,
    #[case] count: usize,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let dag = NamedGraph::new(&["r", "a", "s"], &[("r", "a"), ("r", "s"), ("a", "s")])?;
    let schedule = GraphSchedule::new(&dag.graph)?;
    let outputs = schedule.scan_up(inputs_of(&dag.graph, ()), |key, (), _: PerEdge<()>| {
      let values = if dag.name(key) == wrong {
        count
      } else {
        schedule.parents(key).len()
      };
      Ok((((), smallvec![(); values]), ()))
    })?;

    let result = outputs.split_parent_edges();

    assert_error!(
      result,
      format!(
        "Scan output of node {} {expected}. This is an internal error. Please report it to developers.",
        dag.key(wrong)
      )
    );
    Ok(())
  }

  #[test]
  fn test_schedule_scan_unwrap_shared_rejects_output_with_another_holder() -> Result<(), Report> {
    let tree = NamedGraph::new(&["root", "A"], &[("root", "A")])?;
    let schedule = GraphSchedule::new(&tree.graph)?;
    let shared = Arc::new(2);
    let outputs = schedule.scan_up(inputs_of(&tree.graph, ()), |key, (), _: PerEdge<()>| {
      let output = if key == tree.key("A") {
        Arc::clone(&shared)
      } else {
        Arc::new(1)
      };
      Ok((output, ()))
    })?;

    let result = outputs.unwrap_shared().into_map();

    assert_error!(
      result,
      format!(
        "Scan output of node {} is still shared after the scan. This is an internal error. Please report it to developers.",
        tree.key("A")
      )
    );
    Ok(())
  }

  #[rstest]
  #[case::one_thread(1)]
  #[case::four_threads(4)]
  #[trace]
  fn test_schedule_scan_down_delivers_every_parent_message_in_edge_order(#[case] threads: usize) -> Result<(), Report> {
    let dag = NamedGraph::new(&["p", "q", "x", "s"], &[("p", "x"), ("q", "s"), ("p", "s"), ("x", "s")])?;
    let schedule = GraphSchedule::new(&dag.graph)?;

    let received = run_in_pool(threads, || {
      schedule
        .scan_down(inputs_of(&dag.graph, ()), |key, (), parents: PerEdge<&'static str>| {
          Ok((parents.into_vec(), dag.name(key)))
        })
        .map(ScanOutputs::into_map)
    })?;

    let expected = btreemap! {
      dag.key("p") => vec![], dag.key("q") => vec![], dag.key("x") => vec!["p"], dag.key("s") => vec!["q", "p", "x"],
    };
    assert_eq!(expected, received);
    Ok(())
  }

  #[rstest]
  #[case::one_thread(1)]
  #[case::four_threads(4)]
  #[trace]
  fn test_schedule_scan_up_sends_the_message_of_a_shared_child_to_every_parent(
    #[case] threads: usize,
  ) -> Result<(), Report> {
    let dag = NamedGraph::new(&["p", "q", "x", "s"], &[("p", "s"), ("q", "s"), ("p", "x"), ("x", "s")])?;
    let schedule = GraphSchedule::new(&dag.graph)?;

    let received = run_in_pool(threads, || {
      schedule
        .scan_up(inputs_of(&dag.graph, ()), |key, (), children: PerEdge<&'static str>| {
          Ok((children.into_vec(), dag.name(key)))
        })
        .map(ScanOutputs::into_map)
    })?;

    let expected = btreemap! {
      dag.key("p") => vec!["s", "x"], dag.key("q") => vec!["s"], dag.key("x") => vec!["s"], dag.key("s") => vec![],
    };
    assert_eq!(expected, received);
    Ok(())
  }

  #[test]
  fn test_schedule_scan_join_and_split_parent_edges_roundtrip_with_two_parents() -> Result<(), Report> {
    let dag = NamedGraph::new(&["p", "q", "s"], &[("q", "s"), ("p", "s")])?;
    let schedule = GraphSchedule::new(&dag.graph)?;
    let (p, q, s) = (dag.key("p"), dag.key("q"), dag.key("s"));
    let (q_s, p_s) = (dag.edge("q", "s")?, dag.edge("p", "s")?);
    let nodes = btreemap! { p => "p", q => "q", s => "s" };
    let edges = btreemap! { q_s => "q-s", p_s => "p-s" };

    let joined = schedule
      .join_parent_edges(nodes.clone(), edges.clone())?
      .collect::<BTreeMap<_, _>>();
    let split = schedule
      .scan_down(joined.clone(), |_, joined, _: PerEdge<()>| Ok((joined, ())))?
      .split_parent_edges()?;

    let expected_joined =
      btreemap! { p => ("p", smallvec![]), q => ("q", smallvec![]), s => ("s", smallvec!["q-s", "p-s"]) };
    assert_eq!(expected_joined, joined);
    assert_eq!(NodeEdgeMaps { nodes, edges }, split);
    Ok(())
  }

  #[test]
  fn test_schedule_scan_join_parent_edges_rejects_missing_edge_of_a_second_parent() -> Result<(), Report> {
    let dag = NamedGraph::new(&["p", "q", "s"], &[("q", "s"), ("p", "s")])?;
    let schedule = GraphSchedule::new(&dag.graph)?;
    let edges = btreemap! { dag.edge("q", "s")? => () };

    let result = schedule.join_parent_edges(inputs_of(&dag.graph, ()), edges);

    assert_error!(
      result,
      format!(
        "Scan input is missing edge {} above node {}. This is an internal error. Please report it to developers.",
        dag.edge("p", "s")?,
        dag.key("s")
      )
    );
    Ok(())
  }

  mod helpers {
    use super::super::super::test_graph_support::tests::helpers::DagSpec;
    use crate::graph::Graph;
    use crate::node::GraphNodeKey;
    use crate::schedule::GraphSchedule;
    use eyre::Report;
    use rayon::ThreadPoolBuilder;
    use std::collections::BTreeMap;

    const SMALL_WORKER_STACK_BYTES: usize = 32 * 1024;
    pub(super) const BRAID_WIDTH: usize = 8;
    pub(super) const BRAID_LEVELS: usize = 12_500;

    pub(super) fn inputs_of<N: Clone>(graph: &Graph, value: N) -> BTreeMap<GraphNodeKey, N> {
      graph.get_nodes().map(|node| (node.key(), value.clone())).collect()
    }

    pub(super) fn run_in_pool<T: Send>(
      threads: usize,
      run: impl FnOnce() -> Result<T, Report> + Send,
    ) -> Result<T, Report> {
      ThreadPoolBuilder::new().num_threads(threads).build()?.install(run)
    }

    pub(super) fn run_in_small_stack_pool<T: Send>(
      threads: usize,
      run: impl FnOnce() -> Result<T, Report> + Send,
    ) -> Result<T, Report> {
      ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(SMALL_WORKER_STACK_BYTES)
        .build()?
        .install(run)
    }

    pub(super) fn caterpillar_parents(count: usize, spine_first: bool) -> Vec<usize> {
      (0..count)
        .map(|index| {
          let previous = index.saturating_sub(1);
          let pair_start = previous - previous % 2;
          if spine_first {
            pair_start.saturating_sub(1)
          } else {
            pair_start
          }
        })
        .collect()
    }

    pub(super) fn braid(width: usize, levels: usize) -> DagSpec {
      let edges = (width..width * levels).flat_map(|child| {
        let (level_start, column) = (child - child % width - width, child % width);
        [
          (level_start + column, child),
          (level_start + (column + 1) % width, child),
        ]
      });
      DagSpec {
        keys: (0..width * levels).collect(),
        edges: edges.collect(),
      }
    }

    pub(super) fn braid_levels() -> BTreeMap<GraphNodeKey, usize> {
      (0..BRAID_LEVELS)
        .flat_map(|level| (0..BRAID_WIDTH).map(move |column| (GraphNodeKey(level * BRAID_WIDTH + column), level)))
        .collect()
    }

    pub(super) fn subtree_sizes(schedule: &GraphSchedule<'_>) -> BTreeMap<GraphNodeKey, usize> {
      let mut sizes = BTreeMap::new();
      for &key in schedule.preorder().iter().rev() {
        let size = 1
          + schedule
            .children(key)
            .iter()
            .map(|(child, _)| sizes[child])
            .sum::<usize>();
        sizes.insert(key, size);
      }
      sizes
    }
  }
}

#![allow(
  clippy::collection_is_never_read,
  reason = "test and benchmark code: index and expected-value casts, and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::gtr::jc_distance::jukes_cantor_distance;
  use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
  use crate::partition::marginal::reconstruction::SparseReconstruction;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
  use crate::seq::mutation::Sub;
  use crate::test_utils::sparse_topology_partitions;
  use crate::test_utils::{deletion, find_edge_key, find_node_key_by_name, fitch_edge_obs};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::keymap;

  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_graph::node_names::NodeNames;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::{AsciiChar, Seq};
  use treetime_utils::pretty_assert_relative_eq;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
    Sub::new(c(reff), pos, c(qry)).unwrap()
  }

  fn make_partition(
    graph: &Graph,
    names: &NodeNames,
    length: usize,
    edge_mutations: &[(&str, &str, Vec<Sub>)],
  ) -> Result<PartitionMarginalSparse, Report> {
    let alphabet = Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?;

    let mut ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
    for (_, _, subs) in edge_mutations {
      for s in subs {
        if s.pos() < length {
          ref_seq[s.pos()] = s.reff();
        }
      }
    }

    let mut obs_nodes = keymap! {};
    let mut node_states = btreemap! {};
    for node in graph.get_nodes() {
      let key = node.key();
      obs_nodes.insert(key, SparseNodeObs::new(&ref_seq, &alphabet));
      node_states.insert(key, SparseNodeState::leaf(&ref_seq));
    }

    let mut obs_edges = keymap! {};
    for (source, target, subs) in edge_mutations {
      let edge_key = find_edge_key(graph, names, source, target)
        .unwrap_or_else(|| panic!("edge {source}->{target} not found in graph"));
      obs_edges.insert(edge_key, fitch_edge_obs(subs.clone()));
    }

    let partition = PartitionMarginalSparse {
      alphabet,
      length,
      total_length: length,
      root_sequence: ref_seq,
      obs_nodes,
      obs_edges,
      outside_at_no_data: false,
    };

    Ok(partition)
  }

  fn input_name(names: &NodeNames, key: GraphNodeKey) -> Option<&str> {
    let name = &names[key];
    (!name.starts_with("NODE_")).then_some(name)
  }

  fn find_generated_internal_nodes(graph: &Graph, names: &NodeNames) -> Vec<GraphNodeKey> {
    graph
      .get_nodes()
      .filter(|node| !node.is_leaf() && !node.is_root() && names[node.key()].starts_with("NODE_"))
      .map(|node| node.key())
      .collect()
  }

  #[test]
  fn test_merge_no_polytomy() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)internal:0.3)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "internal", vec![sub(b'A', 0, b'T')]),
        ("internal", "A", vec![sub(b'A', 0, b'T')]),
        ("internal", "B", vec![sub(b'A', 0, b'T')]),
      ],
    )?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 0);
    assert_eq!(graph.get_nodes().count(), 4);
    Ok(())
  }

  #[test]
  fn test_merge_polytomy_no_shared_mutations() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.3)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T')]),
        ("root", "B", vec![sub(b'G', 5, b'C')]),
        ("root", "C", vec![sub(b'T', 10, b'A')]),
      ],
    )?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 0);
    assert_eq!(graph.get_nodes().count(), 4);
    Ok(())
  }

  #[test]
  fn test_merge_polytomy_two_siblings_share_all_mutations() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.3)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let shared = vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')];
    let partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "A", shared.clone()),
        ("root", "B", shared),
        ("root", "C", vec![sub(b'T', 10, b'A')]),
      ],
    )?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 1);

    graph.build()?;
    assert_eq!(graph.get_nodes().count(), 5);
    assert_eq!(graph.get_edges().count(), 4);

    let created = find_generated_internal_nodes(&graph, &names);
    assert_eq!(created.len(), 1);

    let p = &partitions[0].partition;
    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = input_name(&names, target.key());
      if let Some(edge_data) = p.obs_edges.get(edge.key()) {
        match target_name {
          Some("A" | "B") => assert_eq!(
            edge_data.fitch_subs().len(),
            0,
            "child should have no remaining mutations"
          ),
          Some("C") => assert_eq!(edge_data.fitch_subs().len(), 1),
          None => assert_eq!(
            edge_data.fitch_subs().len(),
            2,
            "new internal edge should carry shared mutations"
          ),
          _ => {},
        }
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_polytomy_partial_overlap() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.3)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')]),
        (
          "root",
          "B",
          vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C'), sub(b'T', 10, b'A')],
        ),
        ("root", "C", vec![sub(b'C', 20, b'G')]),
      ],
    )?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 1);

    graph.build()?;
    let p = &partitions[0].partition;

    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = input_name(&names, target.key());
      if let Some(edge_data) = p.obs_edges.get(edge.key()) {
        match target_name {
          Some("A") => assert_eq!(edge_data.fitch_subs().len(), 0),
          Some("B") => {
            assert_eq!(edge_data.fitch_subs().len(), 1);
            assert_eq!(edge_data.fitch_subs()[0], sub(b'T', 10, b'A'));
          },
          Some("C") => assert_eq!(edge_data.fitch_subs().len(), 1),
          None => assert_eq!(
            edge_data.fitch_subs().len(),
            2,
            "internal edge carries shared mutations"
          ),
          _ => {},
        }
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_greedy_picks_best_pair() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.3,D:0.4)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let partition = make_partition(
      &graph,
      &names,
      100,
      &[
        (
          "root",
          "A",
          vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C'), sub(b'T', 10, b'A')],
        ),
        (
          "root",
          "B",
          vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C'), sub(b'T', 10, b'A')],
        ),
        ("root", "C", vec![sub(b'C', 20, b'G')]),
        ("root", "D", vec![sub(b'T', 30, b'A')]),
      ],
    )?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 1);

    graph.build()?;
    assert_eq!(graph.get_nodes().count(), 6);

    let p = &partitions[0].partition;
    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = input_name(&names, target.key());
      if target_name.is_none() {
        if let Some(edge_data) = p.obs_edges.get(edge.key()) {
          assert_eq!(edge_data.fitch_subs().len(), 3);
        }
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_branch_length_adjustment() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.3)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let shared = vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')];
    let partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "A", shared.clone()),
        ("root", "B", shared),
        ("root", "C", vec![sub(b'T', 10, b'A')]),
      ],
    )?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    graph.build()?;

    let d = -0.75 * f64::ln(1.0 - 4.0 * 0.02 / 3.0);
    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = input_name(&names, target.key());
      let bl = branch_lengths[&edge.key()];

      match target_name {
        None => pretty_assert_relative_eq!(bl.unwrap(), d, epsilon = 1e-15),
        Some("A") => pretty_assert_relative_eq!(bl.unwrap(), 0.0, epsilon = 1e-15),
        Some("B") => pretty_assert_relative_eq!(bl.unwrap(), 0.0, epsilon = 1e-15),
        Some("C") => pretty_assert_relative_eq!(bl.unwrap(), 0.3, epsilon = 1e-6),
        _ => {},
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_child_bl_zero_when_all_shared() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.05,B:0.2,C:0.3)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let shared: Vec<Sub> = (0..10).map(|i| sub(b'A', i, b'T')).collect();
    let partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "A", shared.clone()),
        ("root", "B", shared),
        ("root", "C", vec![sub(b'T', 50, b'A')]),
      ],
    )?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = input_name(&names, target.key());
      let bl = branch_lengths[&edge.key()];

      if target_name == Some("A") {
        pretty_assert_relative_eq!(bl.unwrap(), 0.0, epsilon = 1e-15);
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_multiple_partitions() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.3)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let p1 = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T')]),
        ("root", "B", vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')]),
        ("root", "C", vec![sub(b'T', 10, b'A')]),
      ],
    )?;

    let p2_alphabet = Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?;
    let mut p2_ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(200).collect();
    p2_ref_seq[50] = c(b'C');

    let mut p2_obs_nodes = keymap! {};
    let mut p2_node_states = btreemap! {};
    for node in graph.get_nodes() {
      let key = node.key();
      p2_obs_nodes.insert(key, SparseNodeObs::empty(&p2_alphabet));
      p2_node_states.insert(key, SparseNodeState::leaf(&p2_ref_seq));
    }

    let edge_a = find_edge_key(&graph, &names, "root", "A").unwrap();
    let edge_b = find_edge_key(&graph, &names, "root", "B").unwrap();
    let edge_c = find_edge_key(&graph, &names, "root", "C").unwrap();
    let mut p2_obs_edges = keymap! {};
    p2_obs_edges.insert(edge_a, fitch_edge_obs(vec![sub(b'C', 50, b'G')]));
    p2_obs_edges.insert(edge_b, fitch_edge_obs(vec![sub(b'C', 50, b'G')]));
    p2_obs_edges.insert(edge_c, SparseEdgeObs::default());

    let p2 = PartitionMarginalSparse {
      alphabet: p2_alphabet,
      length: 200,
      total_length: 200,
      root_sequence: p2_ref_seq,
      obs_nodes: p2_obs_nodes,
      obs_edges: p2_obs_edges,
      outside_at_no_data: false,
    };

    let mut partitions = sparse_topology_partitions(vec![p1, p2]);

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 1);
    graph.build()?;

    let p_pooled = 2.0 / 300.0;
    let d_expected = -0.75 * f64::ln(1.0 - 4.0 * p_pooled / 3.0);
    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = input_name(&names, target.key());
      let bl = branch_lengths[&edge.key()];

      if target_name.is_none() {
        pretty_assert_relative_eq!(bl.unwrap(), d_expected, epsilon = 1e-15);
      }
    }

    let p1 = &partitions[0].partition;
    let p2 = &partitions[1].partition;
    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = input_name(&names, target.key());
      if target_name == Some("B") {
        assert_eq!(p1.obs_edges[&edge.key()].fitch_subs().len(), 1);
        assert_eq!(p1.obs_edges[&edge.key()].fitch_subs()[0], sub(b'G', 5, b'C'));
        assert_eq!(p2.obs_edges[&edge.key()].fitch_subs().len(), 0);
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_repeated_until_exhausted() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.1,C:0.1,D:0.1,E:0.1)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T')]),
        ("root", "B", vec![sub(b'A', 0, b'T')]),
        ("root", "C", vec![sub(b'G', 5, b'C')]),
        ("root", "D", vec![sub(b'G', 5, b'C')]),
        ("root", "E", vec![sub(b'T', 10, b'A')]),
      ],
    )?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 2);

    graph.build()?;
    assert_eq!(graph.get_nodes().count(), 8);
    assert_eq!(graph.get_edges().count(), 7);

    Ok(())
  }

  #[test]
  fn test_merge_empty_partitions() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.3)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut partitions: Vec<SparseReconstruction> = vec![];

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 0);
    Ok(())
  }

  #[test]
  fn test_merge_preserves_tree_structure_for_non_polytomies() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.1,C:0.1)internal1:0.1,D:0.2)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "internal1", vec![sub(b'C', 20, b'G')]),
        ("internal1", "A", vec![sub(b'A', 0, b'T')]),
        ("internal1", "B", vec![sub(b'A', 0, b'T')]),
        ("internal1", "C", vec![sub(b'G', 5, b'C')]),
        ("root", "D", vec![sub(b'T', 10, b'A')]),
      ],
    )?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 1);

    graph.build()?;
    assert_eq!(graph.get_nodes().count(), 7);

    assert!(find_node_key_by_name(&graph, &names, "D").is_some());
    assert!(find_edge_key(&graph, &names, "root", "D").is_some());

    Ok(())
  }

  #[test]
  fn test_merge_single_mutation_shared() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.05,B:0.05,C:0.05)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let partition = make_partition(
      &graph,
      &names,
      1000,
      &[
        ("root", "A", vec![sub(b'A', 42, b'T')]),
        ("root", "B", vec![sub(b'A', 42, b'T')]),
        ("root", "C", vec![]),
      ],
    )?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 1);
    graph.build()?;

    let d_expected = -0.75 * f64::ln(1.0 - 4.0 * 0.001 / 3.0);
    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = input_name(&names, target.key());
      let bl = branch_lengths[&edge.key()];

      if target_name.is_none() {
        pretty_assert_relative_eq!(bl.unwrap(), d_expected, epsilon = 1e-15);
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_branch_length_jc_correction_differs_from_raw() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.5,B:0.5,C:0.5)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let shared: Vec<Sub> = (0..10).map(|i| sub(b'A', i, b'T')).collect();
    let partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "A", shared.clone()),
        ("root", "B", shared),
        ("root", "C", vec![sub(b'T', 50, b'A')]),
      ],
    )?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    graph.build()?;

    let p = 0.10;
    let d = -0.75 * f64::ln(1.0 - 4.0 * p / 3.0);
    assert!(d > p * 1.05, "JC correction must exceed raw p by >5%: d={d} p={p}");
    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = input_name(&names, target.key());
      let bl = branch_lengths[&edge.key()];

      match target_name {
        None => pretty_assert_relative_eq!(bl.unwrap(), d, epsilon = 1e-15),
        Some("A" | "B") => pretty_assert_relative_eq!(bl.unwrap(), 0.0, epsilon = 1e-15),
        _ => {},
      }
    }

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::basic_remaining(       (2, 0, 1),  100, 0.2)]
  #[case::zero_remaining(        (2, 0, 0),  100, 0.2)]
  #[case::asymmetric_remaining(  (3, 2, 5), 1000, 0.5)]
  #[case::newick_bl_independent( (2, 0, 1),  100, 99.0)]
  #[trace]
  fn test_merge_child_bl(
    #[case] (n_shared, n_unique_a, n_unique_b): (usize, usize, usize),
    #[case] length: usize,
    #[case] newick_bl: f64,
  ) -> Result<(), Report> {
    let newick = format!("(A:{newick_bl},B:{newick_bl},C:{newick_bl})root;");
    let nwk_parsed = nwk_read(newick.as_bytes())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let edge_subs = helpers::build_shared_unique_subs(n_shared, n_unique_a, n_unique_b);
    let partition = make_partition(&graph, &names, length, &edge_subs)?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    graph.build()?;

    let bls = helpers::extract_branch_lengths(&names, &graph, &branch_lengths);
    let expected_bl_a = helpers::expected_jc_bl(n_unique_a, length);
    let expected_bl_b = helpers::expected_jc_bl(n_unique_b, length);
    pretty_assert_relative_eq!(bls["A"], expected_bl_a, epsilon = 1e-15);
    pretty_assert_relative_eq!(bls["B"], expected_bl_b, epsilon = 1e-15);

    Ok(())
  }

  #[test]
  fn test_merge_child_bl_includes_indels_in_remaining() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.5,B:0.5,C:0.5)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let shared = vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')];
    let mut partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "A", shared.clone()),
        ("root", "B", shared),
        ("root", "C", vec![sub(b'T', 20, b'A')]),
      ],
    )?;

    {
      let p = &mut partition;
      let edge_a = find_edge_key(&graph, &names, "root", "A").expect("edge root->A");
      p.obs_edges.get_mut(edge_a).expect("partition edge A").indels =
        vec![deletion((10, 13), Seq::try_from_str("GTA").unwrap())];
    }

    let mut partitions = sparse_topology_partitions(vec![partition]);
    let mut branch_lengths = branch_lengths;
    merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    graph.build()?;

    let bls = helpers::extract_branch_lengths(&names, &graph, &branch_lengths);
    pretty_assert_relative_eq!(bls["A"], helpers::expected_jc_bl(1, 100), epsilon = 1e-15);
    pretty_assert_relative_eq!(bls["B"], helpers::expected_jc_bl(0, 100), epsilon = 1e-15);

    Ok(())
  }

  #[test]
  fn test_merge_child_bl_across_partitions() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.5,B:0.5,C:0.5)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let p1 = make_partition(
      &graph,
      &names,
      100,
      &[
        (
          "root",
          "A",
          vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C'), sub(b'T', 10, b'A')],
        ),
        ("root", "B", vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')]),
        ("root", "C", vec![sub(b'T', 20, b'A')]),
      ],
    )?;

    let p2 = helpers::make_second_partition(
      &graph,
      &names,
      200,
      &[
        (
          ("root", "A"),
          vec![sub(b'C', 50, b'G'), sub(b'G', 60, b'T'), sub(b'T', 70, b'A')],
        ),
        (("root", "B"), vec![sub(b'C', 50, b'G')]),
        (("root", "C"), vec![]),
      ],
    )?;

    let mut partitions = sparse_topology_partitions(vec![p1, p2]);
    let mut branch_lengths = branch_lengths;
    merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    graph.build()?;

    let bls = helpers::extract_branch_lengths(&names, &graph, &branch_lengths);
    pretty_assert_relative_eq!(bls["A"], helpers::expected_jc_bl(3, 300), epsilon = 1e-15);
    pretty_assert_relative_eq!(bls["B"], helpers::expected_jc_bl(0, 300), epsilon = 1e-15);

    Ok(())
  }

  #[test]
  fn test_merge_group_three_siblings_same_mutation() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.1,C:0.1,D:0.1)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T')]),
        ("root", "B", vec![sub(b'A', 0, b'T')]),
        ("root", "C", vec![sub(b'A', 0, b'T')]),
        ("root", "D", vec![sub(b'G', 5, b'C')]),
      ],
    )?;
    let mut partitions = sparse_topology_partitions(vec![partition]);

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 1);
    graph.build()?;

    let created = find_generated_internal_nodes(&graph, &names);
    assert_eq!(created.len(), 1);

    let new_node = graph.get_node(created[0]).expect("new internal node");
    assert_eq!(new_node.degree_out(), 3);

    Ok(())
  }

  #[test]
  fn test_merge_shared_indels_only() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.1,C:0.1)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "A", vec![]),
        ("root", "B", vec![]),
        ("root", "C", vec![sub(b'T', 20, b'A')]),
      ],
    )?;

    let shared_indel = deletion((5, 8), Seq::try_from_str("GTA").unwrap());
    {
      let p = &mut partition;
      let edge_a = find_edge_key(&graph, &names, "root", "A").expect("edge root->A");
      let edge_b = find_edge_key(&graph, &names, "root", "B").expect("edge root->B");
      p.obs_edges.get_mut(edge_a).expect("partition edge A").indels = vec![shared_indel.clone()];
      p.obs_edges.get_mut(edge_b).expect("partition edge B").indels = vec![shared_indel];
    }

    let mut partitions = sparse_topology_partitions(vec![partition]);
    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 1);

    Ok(())
  }

  #[test]
  fn test_merge_shared_subs_and_indels_split_correctly() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.1,C:0.1)root;".as_slice())?;
    let mut names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut partition = make_partition(
      &graph,
      &names,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T'), sub(b'G', 15, b'C')]),
        ("root", "B", vec![sub(b'A', 0, b'T')]),
        ("root", "C", vec![sub(b'T', 20, b'A')]),
      ],
    )?;

    let shared_indel = deletion((5, 8), Seq::try_from_str("GTA").unwrap());
    {
      let p = &mut partition;
      let edge_a = find_edge_key(&graph, &names, "root", "A").expect("edge root->A");
      let edge_b = find_edge_key(&graph, &names, "root", "B").expect("edge root->B");
      p.obs_edges.get_mut(edge_a).expect("partition edge A").indels = vec![shared_indel.clone()];
      p.obs_edges.get_mut(edge_b).expect("partition edge B").indels = vec![shared_indel];
    }

    let mut partitions = sparse_topology_partitions(vec![partition]);
    let mut branch_lengths = branch_lengths;
    merge_shared_mutation_branches(&mut graph, &mut names, &mut partitions, &mut branch_lengths)?;
    graph.build()?;

    let edge_data = helpers::extract_edge_mutation_counts(&names, &graph, &partitions[0].partition);
    assert_eq!(edge_data[&None], (1, 1), "parent edge: 1 shared sub, 1 shared indel");
    assert_eq!(edge_data[&Some("A")], (1, 0), "A: 1 unique sub, no indels");
    assert_eq!(edge_data[&Some("B")], (0, 0), "B: no remaining mutations");

    Ok(())
  }

  mod helpers {
    use super::*;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::keymap;
    use treetime_utils::cast::ToF64Exact;

    pub(super) fn expected_jc_bl(count: usize, length: usize) -> f64 {
      if length == 0 || count == 0 {
        return 0.0;
      }
      jukes_cantor_distance(count.to_f64_exact() / length.to_f64_exact(), 4)
    }

    pub(super) fn extract_branch_lengths(
      names: &NodeNames,
      graph: &Graph,
      branch_lengths: &EdgeMap<Option<f64>>,
    ) -> BTreeMap<String, f64> {
      graph
        .get_edges()
        .filter_map(|edge| {
          let target = graph.get_node(edge.target())?;
          let name = input_name(names, target.key())?.to_owned();
          let bl = branch_lengths[&edge.key()]?;
          Some((name, bl))
        })
        .collect()
    }

    pub(super) fn extract_edge_mutation_counts<'a>(
      names: &NodeNames,
      graph: &Graph,
      partition: &PartitionMarginalSparse,
    ) -> BTreeMap<Option<&'a str>, (usize, usize)> {
      let p = &partition;
      let mut result = BTreeMap::new();
      for edge in graph.get_edges() {
        let target = graph.get_node(edge.target()).expect("target node");
        let target_name = input_name(names, target.key());
        if let Some(edge_data) = p.obs_edges.get(edge.key()) {
          let key: Option<&'a str> = match target_name {
            Some("A") => Some("A"),
            Some("B") => Some("B"),
            Some("C") => Some("C"),
            Some("D") => Some("D"),
            None => None,
            _ => continue,
          };
          result.insert(key, (edge_data.fitch_subs().len(), edge_data.indels.len()));
        }
      }
      result
    }

    pub(super) fn build_shared_unique_subs<'a>(
      n_shared: usize,
      n_unique_a: usize,
      n_unique_b: usize,
    ) -> Vec<(&'a str, &'a str, Vec<Sub>)> {
      let nucs = *b"ACGT";
      let shared: Vec<Sub> = (0..n_shared).map(|i| sub(nucs[i % 4], i, nucs[(i + 1) % 4])).collect();

      let mut subs_a = shared.clone();
      for i in 0..n_unique_a {
        let pos = n_shared + i;
        subs_a.push(sub(nucs[pos % 4], pos, nucs[(pos + 1) % 4]));
      }
      subs_a.sort();

      let mut subs_b = shared;
      for i in 0..n_unique_b {
        let pos = n_shared + n_unique_a + i;
        subs_b.push(sub(nucs[pos % 4], pos, nucs[(pos + 1) % 4]));
      }
      subs_b.sort();

      let c_pos = n_shared + n_unique_a + n_unique_b;
      vec![
        ("root", "A", subs_a),
        ("root", "B", subs_b),
        ("root", "C", vec![sub(nucs[c_pos % 4], c_pos, nucs[(c_pos + 1) % 4])]),
      ]
    }

    pub(super) fn make_second_partition(
      graph: &Graph,
      names: &NodeNames,
      length: usize,
      edge_subs: &[((&str, &str), Vec<Sub>)],
    ) -> Result<PartitionMarginalSparse, Report> {
      let alphabet = Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?;

      let mut ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
      for (_, subs) in edge_subs {
        for s in subs {
          if s.pos() < length {
            ref_seq[s.pos()] = s.reff();
          }
        }
      }

      let mut obs_nodes = keymap! {};
      let mut node_states = btreemap! {};
      for node in graph.get_nodes() {
        let key = node.key();
        obs_nodes.insert(key, SparseNodeObs::empty(&alphabet));
        node_states.insert(key, SparseNodeState::leaf(&ref_seq));
      }

      let mut obs_edges = keymap! {};
      for ((source, target), subs) in edge_subs {
        let edge_key = find_edge_key(graph, names, source, target)
          .unwrap_or_else(|| panic!("edge {source}->{target} not found in graph"));
        obs_edges.insert(edge_key, fitch_edge_obs(subs.clone()));
      }

      let partition = PartitionMarginalSparse {
        alphabet,
        length,
        total_length: length,
        root_sequence: ref_seq,
        obs_nodes,
        obs_edges,
        outside_at_no_data: false,
      };

      Ok(partition)
    }
  }
}

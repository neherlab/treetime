#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::prune::run::{collapse_sparse_edges_from_leaf_recursive, get_edge_num_muts, prune_nodes};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
  use crate::representation::payload::sparse::{SparseEdgePartition, SparseNodePartition};
  use crate::seq::indel::InDel;
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use approx::{assert_relative_eq, assert_ulps_eq};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::{btreemap, btreeset};
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeSet;
  use std::sync::Arc;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
  use treetime_primitives::AsciiChar;
  use treetime_utils::make_report;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  /// Populate partition node entries for all graph nodes with dummy reference sequences.
  /// Required for `edge_subs()` to work (it accesses node data to reconstruct states).
  fn populate_test_nodes(partition: &mut PartitionMarginalSparse, graph: &GraphAncestral) {
    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      partition.nodes.entry(key).or_insert_with(|| {
        let mut node_part = SparseNodePartition::empty(&partition.alphabet);
        node_part.seq.sequence = treetime_primitives::Seq::from_iter((0..partition.length).map(|_| c(b'A')));
        node_part
      });
    }
  }

  fn create_test_graph_with_partitions(
    nwk: &str,
    edge_mutations: &[(usize, Option<usize>)], // (edge_index, num_mutations)
  ) -> Result<(GraphAncestral, Vec<Arc<RwLock<PartitionMarginalSparse>>>), Report> {
    let graph: GraphAncestral = nwk_read_str(nwk)?;

    let partitions = if edge_mutations.is_empty() {
      vec![]
    } else {
      let mut partition = PartitionMarginalSparse {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
        length: 100, // dummy length
        nodes: btreemap! {},
        edges: btreemap! {},
      };

      populate_test_nodes(&mut partition, &graph);

      for (edge_index, num_muts) in edge_mutations {
        if let Some(edge) = graph.get_edges().get(*edge_index) {
          let edge_key = edge.read_arc().key();
          if let Some(num_muts) = num_muts {
            partition.edges.insert(
              edge_key,
              SparseEdgePartition {
                subs: (0..*num_muts)
                  .map(|i| Sub::new(c(b'A'), i, c(b'T')).unwrap())
                  .collect_vec(),
                ..SparseEdgePartition::default()
              },
            );
          }
        }
      }

      vec![Arc::new(RwLock::new(partition))]
    };

    Ok((graph, partitions))
  }

  /// Create test graph with edge mutations specified by node names.
  /// edge_mutations: (source_name, target_name, num_mutations) - None means default (no mutations entry)
  fn create_test_graph_with_named_edge_mutations(
    nwk: &str,
    edge_mutations: &[(&str, &str, Option<usize>)],
  ) -> Result<(GraphAncestral, Vec<Arc<RwLock<PartitionMarginalSparse>>>), Report> {
    let graph: GraphAncestral = nwk_read_str(nwk)?;

    let partitions = if edge_mutations.is_empty() {
      vec![]
    } else {
      let mut partition = PartitionMarginalSparse {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
        length: 100,
        nodes: btreemap! {},
        edges: btreemap! {},
      };

      populate_test_nodes(&mut partition, &graph);

      for (source_name, target_name, num_muts) in edge_mutations {
        if let Some(edge_key) = find_edge_key(&graph, source_name, target_name) {
          match num_muts {
            Some(n) => {
              partition.edges.insert(
                edge_key,
                SparseEdgePartition {
                  subs: (0..*n).map(|i| Sub::new(c(b'A'), i, c(b'T')).unwrap()).collect_vec(),
                  ..SparseEdgePartition::default()
                },
              );
            },
            None => {
              partition.edges.insert(edge_key, SparseEdgePartition::default());
            },
          }
        }
      }

      vec![Arc::new(RwLock::new(partition))]
    };

    Ok((graph, partitions))
  }

  #[test]
  fn test_prune_nodes_basic() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.0,B:0.1)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.0), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0,B:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_with_threshold() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.01,B:0.02,C:0.1)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.05), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.01,B:0.02,C:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_preserves_large_edges() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.1,B:0.2)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.01), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.1,B:0.2)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_empty_graph() -> Result<(), Report> {
    let mut graph: GraphAncestral = Graph::new();
    let partitions = vec![];
    prune_nodes(&mut graph, &partitions, Some(0.0), false, &btreeset! {})?;
    assert!(graph.get_nodes().is_empty());
    Ok(())
  }

  #[test]
  fn test_prune_nodes_handles_none_weights() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.0,B:0.1)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.0), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0,B:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_preserves_terminal_nodes() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.00001,B:0.1)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.001), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:1.0e-5,B:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_complex_tree() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions(
      "(((A:0,B:0.1)internal1:0.00002,(C:0.00003,D:0.1)internal2:0.1)internal3:0.00004,E:0.00005)root;",
      &[],
    )?;
    prune_nodes(&mut graph, &partitions, Some(0.01), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(
      output_nwk,
      "(E:5.0e-5,(C:3.0e-5,D:0.1)internal2:0.1,A:6.0e-5,B:0.1)root;"
    );
    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_preserves_leaves() -> Result<(), Report> {
    // Edge to A has no mutations, edge to B has 2 mutations
    let (mut graph, partitions) = create_test_graph_with_named_edge_mutations(
      "(A:0.1,B:0.1)root;",
      &[("root", "A", None), ("root", "B", Some(2))],
    )?;

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // Both leaves should be preserved even if edge to A has no mutations
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2); // root->A, root->B

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_internal_nodes() -> Result<(), Report> {
    // Structure: ((A:0.1,B:0.1)internal:0.1)root;
    // root->internal edge: no mutations (empty)
    // internal->A edge: 1 mutation
    // internal->B edge: 2 mutations
    let (mut graph, partitions) = create_test_graph_with_named_edge_mutations(
      "((A:0.1,B:0.1)internal:0.1)root;",
      &[
        ("root", "internal", None),
        ("internal", "A", Some(1)),
        ("internal", "B", Some(2)),
      ],
    )?;

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // Internal node should be collapsed, but leaves preserved
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2); // root->A, root->B

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_none_mutations() -> Result<(), Report> {
    // Only add edge data for internal->A edge, leaving root->internal edge unknown (None mutations)
    let (mut graph, partitions) = create_test_graph_with_named_edge_mutations(
      "((A:0.1)internal:0.1)root;",
      &[("internal", "A", Some(1))], // root->internal has no entry (unknown mutations)
    )?;

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // Internal node should be preserved when mutations is None (unknown)
    assert_eq!(graph.get_nodes().len(), 3); // root, internal, A
    assert_eq!(graph.get_edges().len(), 2); // root->internal, internal->A

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_simple_leaf_case() -> Result<(), Report> {
    // Single leaf edge with no mutations should be preserved
    let (mut graph, partitions) = create_test_graph_with_named_edge_mutations("(A:0.1)root;", &[("root", "A", None)])?;

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // Leaf should be preserved even with no mutations
    assert_eq!(graph.get_nodes().len(), 2); // root and A
    assert_eq!(graph.get_edges().len(), 1); // root->A

    Ok(())
  }

  #[test]
  fn test_prune_nodes_combined_prune_short_and_empty() -> Result<(), Report> {
    // Tree: root -> internal1 (short edge 0.001, has muts) -> A (leaf)
    //            -> internal2 (normal edge 0.1, no muts)   -> B (leaf)
    // Topology-only test: mutation counts trigger pruning but content is not verified
    let (mut graph, partitions) = create_test_graph_with_named_edge_mutations(
      "((A:0.1)internal1:0.001,(B:0.1)internal2:0.1)root;",
      &[
        ("root", "internal1", Some(1)), // short edge with mutations
        ("root", "internal2", None),    // empty edge (no mutations)
        ("internal1", "A", None),       // leaf, no overlapping subs with parent
        ("internal2", "B", Some(2)),    // leaf with mutations
      ],
    )?;

    prune_nodes(&mut graph, &partitions, Some(0.01), true, &btreeset! {})?;

    // Both internal nodes should be collapsed, leaves preserved
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2); // root->A, root->B

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_short_threshold_exact() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.05,B:0.05,C:0.051)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.05), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.05,B:0.05,C:0.051)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_short_threshold_below() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.049,B:0.05,C:0.051)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.05), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    // All edges remain because A, B, C are leaves and leaves are never collapsed
    assert_eq!(output_nwk, "(A:0.049,B:0.05,C:0.051)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_complex_tree() -> Result<(), Report> {
    // Complex tree structure: root -> internal1 (with muts) -> A (leaf)
    //                              -> internal1 -> internal3 (no muts) -> C,D (leaves)
    //                         root -> internal2 (no muts) -> B (leaf)
    let (mut graph, partitions) = create_test_graph_with_named_edge_mutations(
      "(((C:0.1,D:0.1)internal3:0.1,A:0.1)internal1:0.1,(B:0.1)internal2:0.1)root;",
      &[
        ("root", "internal1", Some(2)),      // has muts
        ("root", "internal2", Some(0)),      // no muts, to internal
        ("internal1", "A", Some(1)),         // leaf, has muts
        ("internal1", "internal3", Some(0)), // no muts, to internal
        ("internal2", "B", Some(0)),         // leaf, no muts (preserved)
        ("internal3", "C", Some(1)),         // leaf, has muts
        ("internal3", "D", Some(2)),         // leaf, has muts
      ],
    )?;

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // internal2 and internal3 should be collapsed (empty internal edges), leaves preserved
    // Result: root -> internal1 -> A, C, D and root -> B
    assert_eq!(graph.get_nodes().len(), 6); // root, internal1, A, B, C, D
    assert_eq!(graph.get_edges().len(), 5); // root->internal1, internal1->A, internal1->C, internal1->D, root->B

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_both_disabled() -> Result<(), Report> {
    // Very short edge (root->internal) with no mutations - should be preserved when both pruning options disabled
    let (mut graph, partitions) = create_test_graph_with_named_edge_mutations(
      "((A:0.1)internal:0.0001)root;",
      &[("root", "internal", Some(0)), ("internal", "A", Some(1))],
    )?;

    prune_nodes(&mut graph, &partitions, None, false, &btreeset! {})?;

    // Nothing should be collapsed
    assert_eq!(graph.get_nodes().len(), 3); // root, internal, A
    assert_eq!(graph.get_edges().len(), 2); // root->internal, internal->A

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_basic() -> Result<(), Report> {
    // Tree structure: root -> internal1 -> internal2 -> A (the path to prune)
    //                      -> B (to keep)
    //                      -> C (to keep)
    let mut graph: GraphAncestral = nwk_read_str("(((A:0.1)internal2:0.1)internal1:0.1,B:0.2,C:0.3)root;")?;

    let partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];

    // Find the edge leading to leaf A
    let a_inbound_edge =
      find_edge_key(&graph, "internal2", "A").ok_or_else(|| make_report!("Edge internal2->A not found"))?;

    // Recursively prune leaf A and its childless ancestors
    collapse_sparse_edges_from_leaf_recursive(&mut graph, &partitions, a_inbound_edge)?;

    // The result should be: root -> B, root -> C (internal1 and internal2 should be removed)
    assert_eq!(graph.get_nodes().len(), 3); // root, B, C
    assert_eq!(graph.get_edges().len(), 2); // root->B, root->C

    // Verify the remaining nodes by name
    assert!(find_node_key_by_name(&graph, "root").is_some());
    assert!(find_node_key_by_name(&graph, "B").is_some());
    assert!(find_node_key_by_name(&graph, "C").is_some());

    // Verify removed nodes
    assert!(find_node_key_by_name(&graph, "A").is_none());
    assert!(find_node_key_by_name(&graph, "internal1").is_none());
    assert!(find_node_key_by_name(&graph, "internal2").is_none());

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_stops_at_node_with_children() -> Result<(), Report> {
    // Tree structure: root -> internal1 -> A (to prune)
    //                               -> B (to keep)
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)internal1:0.1)root;")?;

    let partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];

    // Find the edge leading to leaf A
    let a_inbound_edge =
      find_edge_key(&graph, "internal1", "A").ok_or_else(|| make_report!("Edge internal1->A not found"))?;

    // Recursively prune leaf A; internal1 becomes unary and should be collapsed upward
    collapse_sparse_edges_from_leaf_recursive(&mut graph, &partitions, a_inbound_edge)?;

    // The result should be: root -> B (internal1 collapsed)
    assert_eq!(graph.get_nodes().len(), 2); // root, B
    assert_eq!(graph.get_edges().len(), 1); // root->B

    // Verify the remaining nodes by name
    assert!(find_node_key_by_name(&graph, "root").is_some());
    assert!(find_node_key_by_name(&graph, "B").is_some());

    // Verify removed nodes
    assert!(find_node_key_by_name(&graph, "A").is_none());
    assert!(find_node_key_by_name(&graph, "internal1").is_none());

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_stops_at_root() -> Result<(), Report> {
    // Tree: root -> A (only child)
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let partitions = vec![];

    // Collapse the path starting at leaf A; should stop at root
    let a_inbound_edge = find_edge_key(&graph, "root", "A").unwrap();
    collapse_sparse_edges_from_leaf_recursive(&mut graph, &partitions, a_inbound_edge)?;

    // Only root should remain
    assert_eq!(graph.get_nodes().len(), 1);
    assert!(find_node_key_by_name(&graph, "root").is_some());
    assert!(find_node_key_by_name(&graph, "A").is_none());
    assert_eq!(graph.get_edges().len(), 0);

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_invalid_edge_key_errors() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });

    graph.add_edge(
      root,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;
    graph.build()?;

    let partitions = vec![];

    // Use a non-existent edge key to ensure we surface an error
    let bogus = GraphEdgeKey(usize::MAX);
    let res = collapse_sparse_edges_from_leaf_recursive(&mut graph, &partitions, bogus);
    assert!(res.is_err());

    Ok(())
  }

  #[test]
  fn test_create_test_edge_num_muts_none_vs_some_zero() -> Result<(), Report> {
    // Test that we can distinguish between unknown mutations (None) and zero mutations (Some(0))
    let (graph, partitions) = create_test_graph_with_partitions("(A:0.1)root;", &[(0, None)])?;
    let edge_unknown_key = graph.get_edges()[0].read_arc().key();

    let (graph2, partitions2) = create_test_graph_with_partitions("(A:0.1)root;", &[(0, Some(0))])?;
    let edge_zero_key = graph2.get_edges()[0].read_arc().key();

    assert_eq!(get_edge_num_muts(&graph, &partitions, edge_unknown_key)?, None);
    assert_eq!(get_edge_num_muts(&graph2, &partitions2, edge_zero_key)?, Some(0));

    Ok(())
  }

  #[test]
  fn test_prune_nodes_single_named_leaf() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.1,B:0.2,C:0.3)root;", &[])?;
    prune_nodes(&mut graph, &partitions, None, false, &btreeset! { "B".to_owned() })?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.1,C:0.3)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_multiple_named_leaves() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.1,B:0.2,C:0.3,D:0.4)root;", &[])?;
    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "A".to_owned(), "C".to_owned() },
    )?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(B:0.2,D:0.4)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_nonexistent_name_is_noop() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.1,B:0.2)root;", &[])?;
    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "nonexistent".to_owned() },
    )?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.1,B:0.2)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_named_internal_node() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("((A:0.1,B:0.2)internal:0.3,C:0.4)root;", &[])?;
    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
    )?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    // Internal node collapsed, its children moved to root
    assert_eq!(output_nwk, "(C:0.4,A:0.4,B:0.5)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_mixed_internal_and_leaf_names() -> Result<(), Report> {
    let (mut graph, partitions) =
      create_test_graph_with_partitions("((A:0.1,B:0.2)internal:0.3,C:0.4,D:0.5)root;", &[])?;
    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "internal".to_owned(), "D".to_owned() },
    )?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    // Internal collapsed, D removed
    assert_eq!(output_nwk, "(C:0.4,A:0.4,B:0.5)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_empty_names_set_is_noop() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.1,B:0.2,C:0.3)root;", &[])?;
    prune_nodes(&mut graph, &partitions, None, false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.1,B:0.2,C:0.3)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_all_leaves_preserves_root() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.1,B:0.2)root;", &[])?;
    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "A".to_owned(), "B".to_owned() },
    )?;
    // Root should remain
    assert_eq!(graph.get_nodes().len(), 1);
    Ok(())
  }

  #[test]
  fn test_prune_nodes_deep_nested_leaf_removal() -> Result<(), Report> {
    let (mut graph, partitions) =
      create_test_graph_with_partitions("(((A:0.1,B:0.2)i1:0.3,C:0.4)i2:0.5,D:0.6)root;", &[])?;
    prune_nodes(&mut graph, &partitions, None, false, &btreeset! { "A".to_owned() })?;

    // A removed, i1 becomes unary and collapses
    // Verify structure: root has 2 children (i2, D), i2 has 2 children (B, C)
    assert_eq!(graph.get_nodes().len(), 5); // root, i2, B, C, D

    // Verify the leaf names are preserved
    let leaf_names: BTreeSet<_> = graph
      .get_nodes()
      .iter()
      .filter(|n| graph.is_leaf(n.read_arc().key()))
      .filter_map(|n| n.read_arc().payload().read_arc().name.clone())
      .collect();
    assert_eq!(leaf_names, btreeset! { "B".to_owned(), "C".to_owned(), "D".to_owned() });

    Ok(())
  }

  #[test]
  fn test_collapse_edge_compose_non_overlapping() -> Result<(), Report> {
    // Non-overlapping positions: all subs kept from both edges
    // Tree: root -> internal (subs at pos 0,1) -> A (subs at pos 2,3), B
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)internal:0.1)root;")?;

    let root_internal_edge_key = find_edge_key(&graph, "root", "internal").unwrap();
    let internal_a_edge_key = find_edge_key(&graph, "internal", "A").unwrap();
    let internal_b_edge_key = find_edge_key(&graph, "internal", "B").unwrap();

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    // Root -> internal has mutations at positions 0, 1
    partition.edges.insert(
      root_internal_edge_key,
      SparseEdgePartition {
        subs: vec![
          Sub::new(c(b'A'), 0_usize, c(b'T'))?,
          Sub::new(c(b'A'), 1_usize, c(b'C'))?,
        ],
        ..SparseEdgePartition::default()
      },
    );
    // Internal -> A has mutations at positions 2, 3
    partition.edges.insert(
      internal_a_edge_key,
      SparseEdgePartition {
        subs: vec![
          Sub::new(c(b'G'), 2_usize, c(b'T'))?,
          Sub::new(c(b'C'), 3_usize, c(b'A'))?,
        ],
        ..SparseEdgePartition::default()
      },
    );
    // Internal -> B has no mutations
    partition
      .edges
      .insert(internal_b_edge_key, SparseEdgePartition::default());

    let partitions = vec![Arc::new(RwLock::new(partition))];

    // Collapse internal node by name
    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
    )?;

    // After collapse: root -> A should have all 4 mutations, root -> B should have 2 mutations
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2);

    // Find new edge keys after collapse
    let partition = partitions[0].read_arc();
    let mut a_muts_count = 0;
    let mut b_muts_count = 0;
    for edge in graph.get_edges() {
      let edge_key = edge.read_arc().key();
      let target_key = edge.read_arc().target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| n.read_arc().payload().read_arc().name.clone());

      if let Some(edge_partition) = partition.edges.get(&edge_key) {
        if target_name.as_deref() == Some("A") {
          a_muts_count = edge_partition.subs.len();
        } else if target_name.as_deref() == Some("B") {
          b_muts_count = edge_partition.subs.len();
        }
      }
    }

    assert_eq!(a_muts_count, 4); // 2 from root->internal + 2 from internal->A
    assert_eq!(b_muts_count, 2); // 2 from root->internal + 0 from internal->B

    Ok(())
  }

  #[test]
  fn test_collapse_edge_compose_chain() -> Result<(), Report> {
    // Chain composition: parent A->G + child G->T = net A->T at same position
    // Tree: root -> internal -> A, B
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)internal:0.1)root;")?;

    let root_internal_edge_key = find_edge_key(&graph, "root", "internal").unwrap();
    let internal_a_edge_key = find_edge_key(&graph, "internal", "A").unwrap();
    let internal_b_edge_key = find_edge_key(&graph, "internal", "B").unwrap();

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    // Parent edge: A->G at pos 0
    partition.edges.insert(
      root_internal_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new(c(b'A'), 0_usize, c(b'G'))?],
        ..SparseEdgePartition::default()
      },
    );
    // Child edge: G->T at pos 0 (ref=G matches parent qry)
    partition.edges.insert(
      internal_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new(c(b'G'), 0_usize, c(b'T'))?],
        ..SparseEdgePartition::default()
      },
    );
    partition
      .edges
      .insert(internal_b_edge_key, SparseEdgePartition::default());

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
    )?;

    // After collapse: net A->T at pos 0 (1 composed sub)
    let partition = partitions[0].read_arc();
    for edge in graph.get_edges() {
      let edge_key = edge.read_arc().key();
      let target_key = edge.read_arc().target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| n.read_arc().payload().read_arc().name.clone());

      if target_name.as_deref() == Some("A") {
        let edge_partition = &partition.edges[&edge_key];
        assert_eq!(edge_partition.subs.len(), 1);
        assert_eq!(edge_partition.subs[0].reff(), c(b'A'));
        assert_eq!(edge_partition.subs[0].qry(), c(b'T'));
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_compose_cancellation() -> Result<(), Report> {
    // Cancellation: parent A->G + child G->A = no net change at same position
    // Tree: root -> internal -> A, B
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)internal:0.1)root;")?;

    let root_internal_edge_key = find_edge_key(&graph, "root", "internal").unwrap();
    let internal_a_edge_key = find_edge_key(&graph, "internal", "A").unwrap();
    let internal_b_edge_key = find_edge_key(&graph, "internal", "B").unwrap();

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    // Parent edge: A->G at pos 0
    partition.edges.insert(
      root_internal_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new(c(b'A'), 0_usize, c(b'G'))?],
        ..SparseEdgePartition::default()
      },
    );
    // Child edge: G->A at pos 0 (reverts back to original state)
    partition.edges.insert(
      internal_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new(c(b'G'), 0_usize, c(b'A'))?],
        ..SparseEdgePartition::default()
      },
    );
    partition
      .edges
      .insert(internal_b_edge_key, SparseEdgePartition::default());

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
    )?;

    // After collapse: mutations cancel, edge to A has 0 subs
    let partition = partitions[0].read_arc();
    for edge in graph.get_edges() {
      let edge_key = edge.read_arc().key();
      let target_key = edge.read_arc().target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| n.read_arc().payload().read_arc().name.clone());

      if target_name.as_deref() == Some("A") {
        let edge_partition = &partition.edges[&edge_key];
        assert_eq!(edge_partition.subs.len(), 0);
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_compose_multiple_partitions() -> Result<(), Report> {
    // Substitutions composed independently per partition
    // Tree: root -> internal -> A, B
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)internal:0.1)root;")?;

    let root_internal_edge_key = find_edge_key(&graph, "root", "internal").unwrap();
    let internal_a_edge_key = find_edge_key(&graph, "internal", "A").unwrap();
    let internal_b_edge_key = find_edge_key(&graph, "internal", "B").unwrap();

    // Create two partitions with different mutations
    let mut partition1 = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    partition1.edges.insert(
      root_internal_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new(c(b'A'), 0_usize, c(b'T'))?],
        ..SparseEdgePartition::default()
      },
    );
    partition1.edges.insert(
      internal_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new(c(b'G'), 1_usize, c(b'C'))?],
        ..SparseEdgePartition::default()
      },
    );
    partition1
      .edges
      .insert(internal_b_edge_key, SparseEdgePartition::default());

    let mut partition2 = PartitionMarginalSparse {
      index: 1,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    partition2.edges.insert(
      root_internal_edge_key,
      SparseEdgePartition {
        subs: vec![
          Sub::new(c(b'C'), 10_usize, c(b'A'))?,
          Sub::new(c(b'T'), 11_usize, c(b'G'))?,
        ],
        ..SparseEdgePartition::default()
      },
    );
    partition2.edges.insert(
      internal_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new(c(b'A'), 12_usize, c(b'T'))?],
        ..SparseEdgePartition::default()
      },
    );
    partition2
      .edges
      .insert(internal_b_edge_key, SparseEdgePartition::default());

    let partitions = vec![Arc::new(RwLock::new(partition1)), Arc::new(RwLock::new(partition2))];

    // Collapse internal node by name
    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
    )?;

    // Check partition 1: edge to A should have 2 mutations
    let p1 = partitions[0].read_arc();
    let p2 = partitions[1].read_arc();

    for edge in graph.get_edges() {
      let edge_key = edge.read_arc().key();
      let target_key = edge.read_arc().target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| n.read_arc().payload().read_arc().name.clone());

      if target_name.as_deref() == Some("A") {
        assert_eq!(p1.edges[&edge_key].subs.len(), 2); // 1 + 1
        assert_eq!(p2.edges[&edge_key].subs.len(), 3); // 2 + 1
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_branch_length_sum_both_some() -> Result<(), Report> {
    // When both edges have branch lengths, they should be summed
    // Tree: root -> internal:0.3 -> A:0.2, B:0.1
    let mut graph: GraphAncestral = nwk_read_str("((A:0.2,B:0.1)internal:0.3)root;")?;

    let partitions = vec![];

    // Mark internal node for pruning by name
    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
    )?;

    // After collapse: root -> A should have 0.3 + 0.2 = 0.5
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| n.read_arc().payload().read_arc().name.clone());
      let branch_length = edge.payload().read_arc().branch_length;

      if target_name.as_deref() == Some("A") {
        // 0.3 + 0.2 = 0.5
        assert_relative_eq!(branch_length.unwrap(), 0.5, epsilon = 1e-6);
      } else if target_name.as_deref() == Some("B") {
        // 0.3 + 0.1 = 0.4
        assert_relative_eq!(branch_length.unwrap(), 0.4, epsilon = 1e-6);
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_branch_length_sum_precision() -> Result<(), Report> {
    // Verify precision is preserved when summing small branch lengths
    // Tree: ((A:2e-10,B:3e-10)internal:1e-10)root;
    let mut graph: GraphAncestral = nwk_read_str("((A:2e-10,B:3e-10)internal:1e-10)root;")?;

    let partitions = vec![];

    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
    )?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| n.read_arc().payload().read_arc().name.clone());
      let branch_length = edge.payload().read_arc().branch_length;

      if target_name.as_deref() == Some("A") {
        // 1e-10 + 2e-10 = 3e-10
        assert_ulps_eq!(branch_length.unwrap(), 3e-10, max_ulps = 4);
      } else if target_name.as_deref() == Some("B") {
        // 1e-10 + 3e-10 = 4e-10
        assert_ulps_eq!(branch_length.unwrap(), 4e-10, max_ulps = 4);
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_branch_length_none_plus_some() -> Result<(), Report> {
    // When removed edge has Some and new edge has None, result should stay None
    // (based on collapse_sparse_edge logic: only sums when both are Some)
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });
    let internal = graph.add_node(NodeAncestral {
      name: Some("internal".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });

    // Root -> internal has branch length, internal -> A has None
    graph.add_edge(
      root,
      internal,
      EdgeAncestral {
        branch_length: Some(0.5),
      },
    )?;
    graph.add_edge(internal, a, EdgeAncestral { branch_length: None })?;
    graph.add_edge(
      internal,
      b,
      EdgeAncestral {
        branch_length: Some(0.2),
      },
    )?;

    graph.build()?;

    let partitions = vec![];

    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
    )?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| n.read_arc().payload().read_arc().name.clone());
      let branch_length = edge.payload().read_arc().branch_length;

      if target_name.as_deref() == Some("A") {
        // One was None, so result stays None (no sum performed)
        assert!(branch_length.is_none());
      } else if target_name.as_deref() == Some("B") {
        // Both Some: 0.5 + 0.2 = 0.7
        assert_ulps_eq!(branch_length.unwrap(), 0.7, max_ulps = 4);
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_branch_length_some_plus_none() -> Result<(), Report> {
    // When removed edge has None and new edge has Some, result should stay as-is (Some)
    // because condition requires both Some
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });
    let internal = graph.add_node(NodeAncestral {
      name: Some("internal".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });

    // Root -> internal has None, internal -> A has Some
    graph.add_edge(root, internal, EdgeAncestral { branch_length: None })?;
    graph.add_edge(
      internal,
      a,
      EdgeAncestral {
        branch_length: Some(0.3),
      },
    )?;
    graph.add_edge(
      internal,
      b,
      EdgeAncestral {
        branch_length: Some(0.2),
      },
    )?;

    graph.build()?;

    let partitions = vec![];

    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
    )?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| n.read_arc().payload().read_arc().name.clone());
      let branch_length = edge.payload().read_arc().branch_length;

      if target_name.as_deref() == Some("A") {
        // Removed edge had None, so no sum - original Some(0.3) preserved
        assert_ulps_eq!(branch_length.unwrap(), 0.3, max_ulps = 4);
      } else if target_name.as_deref() == Some("B") {
        // Removed edge had None, so no sum - original Some(0.2) preserved
        assert_ulps_eq!(branch_length.unwrap(), 0.2, max_ulps = 4);
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_branch_length_both_none() -> Result<(), Report> {
    // When both edges have None, result should stay None
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });
    let internal = graph.add_node(NodeAncestral {
      name: Some("internal".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      desc: None,
      ..NodeAncestral::default()
    });

    graph.add_edge(root, internal, EdgeAncestral { branch_length: None })?;
    graph.add_edge(internal, a, EdgeAncestral { branch_length: None })?;
    graph.add_edge(internal, b, EdgeAncestral { branch_length: None })?;

    graph.build()?;

    let partitions = vec![];

    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
    )?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let branch_length = edge.payload().read_arc().branch_length;
      // Both None -> stays None
      assert!(branch_length.is_none());
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_indel_preservation() -> Result<(), Report> {
    // Indels from both removed (parent) and retained (child) edges must be preserved
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)internal:0.1)root;")?;

    let root_internal_edge_key = find_edge_key(&graph, "root", "internal").unwrap();
    let internal_a_edge_key = find_edge_key(&graph, "internal", "A").unwrap();
    let internal_b_edge_key = find_edge_key(&graph, "internal", "B").unwrap();

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let parent_indel = InDel::del((10, 15), [c(b'A'), c(b'C'), c(b'G'), c(b'T'), c(b'A')].as_slice());
    let child_indel = InDel::ins((20, 23), [c(b'G'), c(b'G'), c(b'C')].as_slice());

    partition.edges.insert(
      root_internal_edge_key,
      SparseEdgePartition {
        indels: vec![parent_indel],
        ..SparseEdgePartition::default()
      },
    );
    partition.edges.insert(
      internal_a_edge_key,
      SparseEdgePartition {
        indels: vec![child_indel],
        ..SparseEdgePartition::default()
      },
    );
    partition
      .edges
      .insert(internal_b_edge_key, SparseEdgePartition::default());

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(
      &mut graph,
      &partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
    )?;

    let partition = partitions[0].read_arc();
    for edge in graph.get_edges() {
      let edge_key = edge.read_arc().key();
      let target_key = edge.read_arc().target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| n.read_arc().payload().read_arc().name.clone());

      if target_name.as_deref() == Some("A") {
        let edge_partition = &partition.edges[&edge_key];
        // Parent indel first, then child indel
        assert_eq!(edge_partition.indels.len(), 2);
        assert_eq!(edge_partition.indels[0].range, (10, 15));
        assert!(edge_partition.indels[0].deletion);
        assert_eq!(edge_partition.indels[1].range, (20, 23));
        assert!(!edge_partition.indels[1].deletion);
      } else if target_name.as_deref() == Some("B") {
        let edge_partition = &partition.edges[&edge_key];
        // Only parent indel (child had none)
        assert_eq!(edge_partition.indels.len(), 1);
        assert_eq!(edge_partition.indels[0].range, (10, 15));
      }
    }

    Ok(())
  }
}

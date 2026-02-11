#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::prune::run::{collapse_sparse_edges_from_leaf_recursive, get_edge_num_muts, prune_nodes};
  use crate::graph::edge::GraphEdgeKey;
  use crate::graph::graph::Graph;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
  use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
  use crate::representation::graph_sparse::SparseEdgePartition;
  use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
  use crate::seq::mutation::Sub;
  use eyre::Report;
  use itertools::Itertools;
  use maplit::{btreemap, btreeset};
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;

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
        alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
        length: 100, // dummy length
        nodes: btreemap! {},
        edges: btreemap! {},
      };

      for (edge_index, num_muts) in edge_mutations {
        if let Some(edge) = graph.get_edges().get(*edge_index) {
          let edge_key = edge.read_arc().key();
          if let Some(num_muts) = num_muts {
            partition.edges.insert(
              edge_key,
              SparseEdgePartition {
                subs: (0..*num_muts).map(|i| Sub::new('A', i, 'T').unwrap()).collect_vec(),
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
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      desc: None,
    });

    // Edge with no mutations to leaf A (should be preserved)
    graph.add_edge(
      root,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;
    // Edge with mutations to leaf B
    graph.add_edge(
      root,
      b,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let root_a_edge_key = graph.get_edges()[0].read_arc().key();
    let root_b_edge_key = graph.get_edges()[1].read_arc().key();

    // Edge to A has no mutations
    partition.edges.insert(root_a_edge_key, SparseEdgePartition::default());
    // Edge to B has 2 mutations
    partition.edges.insert(
      root_b_edge_key,
      SparseEdgePartition {
        subs: (0_usize..2).map(|i| Sub::new('A', i, 'T').unwrap()).collect_vec(),
        ..SparseEdgePartition::default()
      },
    );

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // Both leaves should be preserved even if edge to A has no mutations
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2); // root->A, root->B

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_internal_nodes() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
    });
    let internal = graph.add_node(NodeAncestral {
      name: Some("internal".to_owned()),
      desc: None,
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      desc: None,
    });

    // Internal edge with no mutations
    graph.add_edge(
      root,
      internal,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;
    // Leaf edges with mutations
    graph.add_edge(
      internal,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;
    graph.add_edge(
      internal,
      b,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let root_internal_edge_key = graph.get_edges()[0].read_arc().key();
    let internal_a_edge_key = graph.get_edges()[1].read_arc().key();
    let internal_b_edge_key = graph.get_edges()[2].read_arc().key();

    // Internal edge has no mutations
    partition
      .edges
      .insert(root_internal_edge_key, SparseEdgePartition::default());
    // Leaf edges have mutations
    partition.edges.insert(
      internal_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..SparseEdgePartition::default()
      },
    );
    partition.edges.insert(
      internal_b_edge_key,
      SparseEdgePartition {
        subs: (0_usize..2).map(|i| Sub::new('A', i, 'T').unwrap()).collect_vec(),
        ..SparseEdgePartition::default()
      },
    );

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // Internal node should be collapsed, but leaves preserved
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2); // root->A, root->B

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_none_mutations() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
    });
    let internal = graph.add_node(NodeAncestral {
      name: Some("internal".to_owned()),
      desc: None,
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
    });

    // Edge with None mutations (unknown number of mutations - should not be pruned)
    graph.add_edge(
      root,
      internal,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;
    graph.add_edge(
      internal,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let internal_a_edge_key = graph.get_edges()[1].read_arc().key();

    // Only add edge data for internal->A edge, leaving root->internal edge unknown (None mutations)
    partition.edges.insert(
      internal_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..SparseEdgePartition::default()
      },
    );

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // Internal node should be preserved when mutations is None (unknown)
    assert_eq!(graph.get_nodes().len(), 3); // root, internal, A
    assert_eq!(graph.get_edges().len(), 2); // root->internal, internal->A

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_simple_leaf_case() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
    });

    // Single leaf edge with no mutations should be preserved
    graph.add_edge(
      root,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let root_a_edge_key = graph.get_edges()[0].read_arc().key();
    partition.edges.insert(root_a_edge_key, SparseEdgePartition::default());

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // Leaf should be preserved even with no mutations
    assert_eq!(graph.get_nodes().len(), 2); // root and A
    assert_eq!(graph.get_edges().len(), 1); // root->A

    Ok(())
  }

  #[test]
  fn test_prune_nodes_combined_prune_short_and_empty() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
    });
    let internal1 = graph.add_node(NodeAncestral {
      name: Some("internal1".to_owned()),
      desc: None,
    });
    let internal2 = graph.add_node(NodeAncestral {
      name: Some("internal2".to_owned()),
      desc: None,
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      desc: None,
    });

    // Short edge (should be pruned by threshold)
    graph.add_edge(
      root,
      internal1,
      EdgeAncestral {
        branch_length: Some(0.001),
      },
    )?;
    // Empty edge (should be pruned by empty check)
    graph.add_edge(
      root,
      internal2,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;
    // Leaf edges (should be preserved)
    graph.add_edge(
      internal1,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;
    graph.add_edge(
      internal2,
      b,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let root_internal1_edge_key = graph.get_edges()[0].read_arc().key();
    let root_internal2_edge_key = graph.get_edges()[1].read_arc().key();
    let internal1_a_edge_key = graph.get_edges()[2].read_arc().key();
    let internal2_b_edge_key = graph.get_edges()[3].read_arc().key();

    // Short edge with mutations
    partition.edges.insert(
      root_internal1_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..SparseEdgePartition::default()
      },
    );
    // Empty edge (no mutations)
    partition
      .edges
      .insert(root_internal2_edge_key, SparseEdgePartition::default());
    // Leaf edges
    partition.edges.insert(
      internal1_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..SparseEdgePartition::default()
      },
    );
    partition.edges.insert(
      internal2_b_edge_key,
      SparseEdgePartition {
        subs: (0_usize..2).map(|i| Sub::new('A', i, 'T').unwrap()).collect_vec(),
        ..SparseEdgePartition::default()
      },
    );

    let partitions = vec![Arc::new(RwLock::new(partition))];

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
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
    });
    let internal1 = graph.add_node(NodeAncestral {
      name: Some("internal1".to_owned()),
      desc: None,
    });
    let internal2 = graph.add_node(NodeAncestral {
      name: Some("internal2".to_owned()),
      desc: None,
    });
    let internal3 = graph.add_node(NodeAncestral {
      name: Some("internal3".to_owned()),
      desc: None,
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      desc: None,
    });
    let c = graph.add_node(NodeAncestral {
      name: Some("C".to_owned()),
      desc: None,
    });
    let d = graph.add_node(NodeAncestral {
      name: Some("D".to_owned()),
      desc: None,
    });

    // Complex tree structure: root -> internal1 (with muts) -> A (leaf)
    //                              -> internal1 -> internal3 (no muts) -> C,D (leaves)
    //                         root -> internal2 (no muts) -> B (leaf)
    graph.add_edge(
      root,
      internal1,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?; // has muts
    graph.add_edge(
      root,
      internal2,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?; // no muts, to internal
    graph.add_edge(
      internal1,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?; // leaf, has muts
    graph.add_edge(
      internal1,
      internal3,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?; // no muts, to internal
    graph.add_edge(
      internal2,
      b,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?; // leaf, no muts (preserved)
    graph.add_edge(
      internal3,
      c,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?; // leaf, has muts
    graph.add_edge(
      internal3,
      d,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?; // leaf, has muts

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let root_internal1_edge_key = graph.get_edges()[0].read_arc().key();
    let root_internal2_edge_key = graph.get_edges()[1].read_arc().key();
    let internal1_a_edge_key = graph.get_edges()[2].read_arc().key();
    let internal1_internal3_edge_key = graph.get_edges()[3].read_arc().key();
    let internal2_b_edge_key = graph.get_edges()[4].read_arc().key();
    let internal3_c_edge_key = graph.get_edges()[5].read_arc().key();
    let internal3_d_edge_key = graph.get_edges()[6].read_arc().key();

    // Set up mutations for each edge
    partition.edges.insert(
      root_internal1_edge_key,
      SparseEdgePartition {
        subs: (0_usize..2).map(|i| Sub::new('A', i, 'T').unwrap()).collect_vec(),
        ..SparseEdgePartition::default()
      },
    ); // has muts
    partition
      .edges
      .insert(root_internal2_edge_key, SparseEdgePartition::default()); // no muts, to internal
    partition.edges.insert(
      internal1_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..SparseEdgePartition::default()
      },
    ); // leaf, has muts
    partition
      .edges
      .insert(internal1_internal3_edge_key, SparseEdgePartition::default()); // no muts, to internal
    partition
      .edges
      .insert(internal2_b_edge_key, SparseEdgePartition::default()); // leaf, no muts (preserved)
    partition.edges.insert(
      internal3_c_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..SparseEdgePartition::default()
      },
    ); // leaf, has muts
    partition.edges.insert(
      internal3_d_edge_key,
      SparseEdgePartition {
        subs: (0_usize..2).map(|i| Sub::new('A', i, 'T').unwrap()).collect_vec(),
        ..SparseEdgePartition::default()
      },
    ); // leaf, has muts

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // internal2 and internal3 should be collapsed (empty internal edges), leaves preserved
    // Result: root -> internal1 -> A, C, D and root -> B
    assert_eq!(graph.get_nodes().len(), 6); // root, internal1, A, B, C, D
    assert_eq!(graph.get_edges().len(), 5); // root->internal1, internal1->A, internal1->C, internal1->D, root->B

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_both_disabled() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
    });
    let internal = graph.add_node(NodeAncestral {
      name: Some("internal".to_owned()),
      desc: None,
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
    });

    // Very short edge with no mutations - should be preserved when both pruning options disabled
    graph.add_edge(
      root,
      internal,
      EdgeAncestral {
        branch_length: Some(0.0001),
      },
    )?;
    graph.add_edge(
      internal,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let root_internal_edge_key = graph.get_edges()[0].read_arc().key();
    let internal_a_edge_key = graph.get_edges()[1].read_arc().key();

    partition
      .edges
      .insert(root_internal_edge_key, SparseEdgePartition::default());
    partition.edges.insert(
      internal_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..SparseEdgePartition::default()
      },
    );

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, None, false, &btreeset! {})?;

    // Nothing should be collapsed
    assert_eq!(graph.get_nodes().len(), 3); // root, internal, A
    assert_eq!(graph.get_edges().len(), 2); // root->internal, internal->A

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_basic() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
    });
    let internal1 = graph.add_node(NodeAncestral {
      name: Some("internal1".to_owned()),
      desc: None,
    });
    let internal2 = graph.add_node(NodeAncestral {
      name: Some("internal2".to_owned()),
      desc: None,
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      desc: None,
    });
    let c = graph.add_node(NodeAncestral {
      name: Some("C".to_owned()),
      desc: None,
    });

    // Tree structure: root -> internal1 -> internal2 -> A (the path to prune)
    //                     -> B (to keep)
    //                     -> C (to keep)
    graph.add_edge(
      root,
      internal1,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;
    graph.add_edge(
      root,
      b,
      EdgeAncestral {
        branch_length: Some(0.2),
      },
    )?;
    graph.add_edge(
      root,
      c,
      EdgeAncestral {
        branch_length: Some(0.3),
      },
    )?;
    graph.add_edge(
      internal1,
      internal2,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;
    graph.add_edge(
      internal2,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;

    graph.build()?;

    let partitions = vec![];

    // Find the edge leading to leaf A
    let a_inbound_edge = {
      let a_node = graph.get_node(a).unwrap();
      a_node.read_arc().inbound()[0]
    };

    // Recursively prune leaf A and its childless ancestors
    collapse_sparse_edges_from_leaf_recursive(&mut graph, &partitions, a_inbound_edge)?;

    // The result should be: root -> B, root -> C (internal1 and internal2 should be removed)
    assert_eq!(graph.get_nodes().len(), 3); // root, B, C
    assert_eq!(graph.get_edges().len(), 2); // root->B, root->C

    // Verify the remaining nodes
    assert!(graph.get_node(root).is_some());
    assert!(graph.get_node(b).is_some());
    assert!(graph.get_node(c).is_some());

    // Verify removed nodes
    assert!(graph.get_node(a).is_none());
    assert!(graph.get_node(internal1).is_none());
    assert!(graph.get_node(internal2).is_none());

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_stops_at_node_with_children() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
    });
    let internal1 = graph.add_node(NodeAncestral {
      name: Some("internal1".to_owned()),
      desc: None,
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      desc: None,
    });

    // Tree structure: root -> internal1 -> A (to prune)
    //                               -> B (to keep)
    graph.add_edge(
      root,
      internal1,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;
    graph.add_edge(
      internal1,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;
    graph.add_edge(
      internal1,
      b,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;

    graph.build()?;

    let partitions = vec![];

    // Find the edge leading to leaf A
    let a_inbound_edge = {
      let a_node = graph.get_node(a).unwrap();
      a_node.read_arc().inbound()[0]
    };

    // Recursively prune leaf A; internal1 becomes unary and should be collapsed upward
    collapse_sparse_edges_from_leaf_recursive(&mut graph, &partitions, a_inbound_edge)?;

    // The result should be: root -> B (internal1 collapsed)
    assert_eq!(graph.get_nodes().len(), 2); // root, B
    assert_eq!(graph.get_edges().len(), 1); // root->B

    // Verify the remaining nodes
    assert!(graph.get_node(root).is_some());
    assert!(graph.get_node(b).is_some());

    // Verify removed node
    assert!(graph.get_node(a).is_none());
    assert!(graph.get_node(internal1).is_none());

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_stops_at_root() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
    });

    // Tree: root -> A (only child)
    graph.add_edge(
      root,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
      },
    )?;
    graph.build()?;

    let partitions = vec![];

    // Collapse the path starting at leaf A; should stop at root
    let a_inbound_edge = {
      let a_node = graph.get_node(a).unwrap();
      a_node.read_arc().inbound()[0]
    };
    collapse_sparse_edges_from_leaf_recursive(&mut graph, &partitions, a_inbound_edge)?;

    // Only root should remain
    assert_eq!(graph.get_nodes().len(), 1);
    assert!(graph.get_node(root).is_some());
    assert!(graph.get_node(a).is_none());
    assert_eq!(graph.get_edges().len(), 0);

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_invalid_edge_key_errors() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      desc: None,
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      desc: None,
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

    assert_eq!(get_edge_num_muts(&partitions, edge_unknown_key), None);
    assert_eq!(get_edge_num_muts(&partitions2, edge_zero_key), Some(0));

    Ok(())
  }
}

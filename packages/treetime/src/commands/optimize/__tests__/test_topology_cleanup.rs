#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::update_marginal;
  use crate::commands::optimize::optimize_unified::{initial_guess_mixed, run_optimize_mixed};
  use crate::commands::optimize::run::{
    apply_damping, collapse_edge_for_optimize, collect_optimize_partitions, find_zero_optimal_internal_edges,
    prune_and_merge_in_loop, save_branch_lengths,
  };
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::representation::payload::sparse::{SparseEdgePartition, SparseNodePartition};
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AsciiChar;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
    Sub::new(c(reff), pos, c(qry)).unwrap()
  }

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

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_empty_graph() -> Result<(), Report> {
    let graph = GraphAncestral::new();
    let edges = find_zero_optimal_internal_edges(&graph);
    assert_eq!(edges.len(), 0);
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_no_zero_edges() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let edges = find_zero_optimal_internal_edges(&graph);
    assert_eq!(edges.len(), 0);
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_skips_leaves() -> Result<(), Report> {
    // A has bl=0.0 but is a leaf: should NOT be collected
    let graph: GraphAncestral = nwk_read_str("(A:0.0,B:0.2)root;")?;
    let edges = find_zero_optimal_internal_edges(&graph);
    assert_eq!(edges.len(), 0);
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_collects_internal() -> Result<(), Report> {
    // I has bl=0.0 and is internal: should be collected
    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)I:0.0,C:0.3)root;")?;
    let edges = find_zero_optimal_internal_edges(&graph);
    assert_eq!(edges.len(), 1);
    let edge_key = edges[0];
    let edge = graph.get_edge(edge_key).unwrap();
    let target = edge.read_arc().target();
    let target_name = graph
      .get_node(target)
      .and_then(|n| n.read_arc().payload().read_arc().name.clone());
    assert_eq!(target_name.as_deref(), Some("I"));
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_multiple() -> Result<(), Report> {
    // Both internal nodes have bl=0.0
    let graph: GraphAncestral = nwk_read_str("(((A:0.1,B:0.1)I1:0.0,C:0.1)I2:0.0,D:0.1)root;")?;
    let edges = find_zero_optimal_internal_edges(&graph);
    assert_eq!(edges.len(), 2);
    Ok(())
  }

  #[test]
  fn test_optimize_collapse_edge_sparse_composes_subs() -> Result<(), Report> {
    // Tree: root -> I (bl=0.0) -> A, B
    // I has sub A0T; A has sub G5C; B has no subs
    // After collapse: root -> A has {A0T, G5C}, root -> B has {A0T}
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;

    let ri_key = find_edge_key(&graph, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, "I", "B").unwrap();

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(AlphabetName::Nuc)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    populate_test_nodes(&mut partition, &graph);

    partition.edges.insert(
      ri_key,
      SparseEdgePartition {
        subs: vec![sub(b'A', 0, b'T')],
        ..SparseEdgePartition::default()
      },
    );
    partition.edges.insert(
      ia_key,
      SparseEdgePartition {
        subs: vec![sub(b'G', 5, b'C')],
        ..SparseEdgePartition::default()
      },
    );
    partition.edges.insert(ib_key, SparseEdgePartition::default());

    let sparse = vec![Arc::new(RwLock::new(partition))];
    let dense: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];

    // Save node key before collapse (node will be removed from graph)
    let i_node_key = find_node_key_by_name(&graph, "I").unwrap();

    collapse_edge_for_optimize(&mut graph, &sparse, &dense, ri_key)?;
    graph.build()?;

    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2);

    // Check substitution composition
    let p = sparse[0].read_arc();
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      if let Some(edge_data) = p.edges.get(&edge.key()) {
        match target_name.as_deref() {
          Some("A") => assert_eq!(edge_data.subs.len(), 2, "A: parent sub + own sub"),
          Some("B") => assert_eq!(edge_data.subs.len(), 1, "B: parent sub only"),
          _ => {},
        }
      }
    }

    // Stale entries removed
    assert!(!p.nodes.contains_key(&i_node_key));
    assert!(!p.edges.contains_key(&ri_key));

    Ok(())
  }

  #[test]
  fn test_optimize_collapse_edge_dense_cleanup() -> Result<(), Report> {
    // Dense partition: stale node/edge entries should be removed after collapse
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;

    let ri_key = find_edge_key(&graph, "root", "I").unwrap();
    let i_key = find_node_key_by_name(&graph, "I").unwrap();

    let mut dense_partition = PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(AlphabetName::Nuc)?,
      length: 10,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    // Add dummy entries for all nodes/edges to verify cleanup
    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      dense_partition.nodes.insert(
        key,
        crate::representation::payload::dense::DenseNodePartition {
          seq: crate::representation::payload::dense::DenseSeqInfo::default(),
          profile: crate::representation::payload::dense::DenseSeqDis::default(),
        },
      );
    }
    for edge in graph.get_edges() {
      let key = edge.read_arc().key();
      dense_partition.edges.insert(
        key,
        crate::representation::payload::dense::DenseEdgePartition::default(),
      );
    }

    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let dense = vec![Arc::new(RwLock::new(dense_partition))];

    collapse_edge_for_optimize(&mut graph, &sparse, &dense, ri_key)?;
    graph.build()?;

    let d = dense[0].read_arc();
    assert!(!d.nodes.contains_key(&i_key), "removed node should be cleaned up");
    assert!(!d.edges.contains_key(&ri_key), "removed edge should be cleaned up");

    Ok(())
  }

  #[test]
  fn test_optimize_collapse_edge_branch_length_sum() -> Result<(), Report> {
    // Collapsed edge has bl=0.0, child edge has bl=0.1
    // After collapse: child edge bl = 0.0 + 0.1 = 0.1 (unchanged)
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)I:0.0)root;")?;

    let ri_key = find_edge_key(&graph, "root", "I").unwrap();

    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let dense: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];

    collapse_edge_for_optimize(&mut graph, &sparse, &dense, ri_key)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let bl = edge.payload().read_arc().branch_length();

      match target_name.as_deref() {
        Some("A") => assert_abs_diff_eq!(bl.unwrap(), 0.1, epsilon = 1e-8),
        Some("B") => assert_abs_diff_eq!(bl.unwrap(), 0.2, epsilon = 1e-8),
        _ => {},
      }
    }

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_empty_list() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let dense: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];

    let changed = prune_and_merge_in_loop(&mut graph, &sparse, &dense, &[])?;
    assert!(!changed);
    assert_eq!(graph.get_nodes().len(), 4);
    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_collapses_and_merges() -> Result<(), Report> {
    // Tree: root -> I (bl=0.0) -> A (subs: A0T), B (subs: A0T)
    //       root -> C (subs: A0T)
    //       root -> D (subs: G5C)
    //
    // After collapse of I: root has 4 children (A, B, C, D) - polytomy
    // A, B, C share sub A0T -> merge creates new internal node
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)I:0.0,C:0.1,D:0.1)root;")?;

    let ri_key = find_edge_key(&graph, "root", "I").unwrap();

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(AlphabetName::Nuc)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    populate_test_nodes(&mut partition, &graph);

    partition.edges.insert(ri_key, SparseEdgePartition::default());

    let ia_key = find_edge_key(&graph, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, "I", "B").unwrap();
    let rc_key = find_edge_key(&graph, "root", "C").unwrap();
    let rd_key = find_edge_key(&graph, "root", "D").unwrap();

    partition.edges.insert(
      ia_key,
      SparseEdgePartition {
        subs: vec![sub(b'A', 0, b'T')],
        ..SparseEdgePartition::default()
      },
    );
    partition.edges.insert(
      ib_key,
      SparseEdgePartition {
        subs: vec![sub(b'A', 0, b'T')],
        ..SparseEdgePartition::default()
      },
    );
    partition.edges.insert(
      rc_key,
      SparseEdgePartition {
        subs: vec![sub(b'A', 0, b'T')],
        ..SparseEdgePartition::default()
      },
    );
    partition.edges.insert(
      rd_key,
      SparseEdgePartition {
        subs: vec![sub(b'G', 5, b'C')],
        ..SparseEdgePartition::default()
      },
    );

    let sparse = vec![Arc::new(RwLock::new(partition))];
    let dense: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];

    // Set bl to 0.0 and pass damped value (simulating damping override in prune_and_merge_in_loop)
    if let Some(edge) = graph.get_edge(ri_key) {
      edge.write_arc().payload().write_arc().set_branch_length(Some(0.0));
    }

    let changed = prune_and_merge_in_loop(&mut graph, &sparse, &dense, &[ri_key])?;
    assert!(changed);

    // I should be gone
    assert!(find_node_key_by_name(&graph, "I").is_none());

    // D remains directly under root
    assert!(find_node_key_by_name(&graph, "D").is_some());

    // Root should have 2 children after merging A, B, C into a new subtree
    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let root_node = graph.get_node(root_key).unwrap();
    assert_eq!(root_node.read_arc().degree_out(), 2);

    Ok(())
  }

  #[test]
  fn test_optimize_loop_with_topology_cleanup_sparse() -> Result<(), Report> {
    // Tree with a zero-length internal branch that should be collapsed during optimization.
    // After Fitch compression + marginal, the optimizer should detect the zero-optimal
    // branch and the loop should collapse it.
    let nuc = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        ACGTACGTACGT
        >B
        ACGTACGTACGT
        >C
        ACGTACGTACGG
        >D
        TCGTACGTACGT
      "#},
      &nuc,
    )?;

    // A and B are identical: the internal edge AB should be optimized to zero
    let mut graph: GraphAncestral = nwk_read_str("((A:0.01,B:0.01)AB:0.01,(C:0.01,D:0.01)CD:0.01)root:0.0;")?;

    let sparse_partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: nuc,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &sparse_partitions, &aln)?;
    update_marginal(&graph, &sparse_partitions)?;

    let dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];
    let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);

    initial_guess_mixed(&graph, &mixed_partitions)?;

    let initial_node_count = graph.get_nodes().len();

    // Run optimize loop with topology cleanup
    let mut lh_prev = f64::MIN;
    for i in 0..10 {
      let sparse_lh = update_marginal(&graph, &sparse_partitions)?;
      let total_lh = sparse_lh;

      if (total_lh - lh_prev).abs() < 1e-2 {
        break;
      }

      let old_branch_lengths = save_branch_lengths(&graph);
      run_optimize_mixed(&graph, &mixed_partitions)?;

      let zero_optimal_edges = find_zero_optimal_internal_edges(&graph);

      apply_damping(&graph, &old_branch_lengths, 0.75, i);

      prune_and_merge_in_loop(&mut graph, &sparse_partitions, &dense_partitions, &zero_optimal_edges)?;

      lh_prev = total_lh;
    }

    // A and B are identical sequences: the AB internal edge should have been
    // collapsed, reducing the node count
    let final_node_count = graph.get_nodes().len();
    assert!(
      final_node_count < initial_node_count,
      "Expected topology simplification: {initial_node_count} nodes -> {final_node_count} nodes"
    );

    // All remaining branch lengths should be non-negative
    for edge in graph.get_edges() {
      let bl = edge.read_arc().payload().read_arc().branch_length().unwrap_or(0.0);
      assert!(bl >= 0.0, "Negative branch length after optimization: {bl}");
    }

    Ok(())
  }

  #[test]
  fn test_optimize_loop_no_collapse_when_branches_nonzero() -> Result<(), Report> {
    // All branches have genuine signal: no edges should be collapsed
    let nuc = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        ACGTACGTACGT
        >B
        TCGTACGTACGT
        >C
        ACGTACGTACGG
        >D
        ACGTACGTACGA
      "#},
      &nuc,
    )?;

    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;")?;

    let sparse_partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: nuc,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &sparse_partitions, &aln)?;
    update_marginal(&graph, &sparse_partitions)?;

    let dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];
    let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);

    initial_guess_mixed(&graph, &mixed_partitions)?;

    let initial_node_count = graph.get_nodes().len();

    let mut lh_prev = f64::MIN;
    for i in 0..10 {
      let total_lh = update_marginal(&graph, &sparse_partitions)?;
      if (total_lh - lh_prev).abs() < 1e-2 {
        break;
      }

      let old_branch_lengths = save_branch_lengths(&graph);
      run_optimize_mixed(&graph, &mixed_partitions)?;

      let zero_optimal_edges = find_zero_optimal_internal_edges(&graph);
      apply_damping(&graph, &old_branch_lengths, 0.75, i);
      prune_and_merge_in_loop(&mut graph, &sparse_partitions, &dense_partitions, &zero_optimal_edges)?;

      lh_prev = total_lh;
    }

    // No edges should have been collapsed - all branches carry genuine signal
    assert_eq!(graph.get_nodes().len(), initial_node_count);

    Ok(())
  }
}

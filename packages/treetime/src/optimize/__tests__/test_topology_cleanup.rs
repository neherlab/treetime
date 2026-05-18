#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::iteration::{apply_damping, save_branch_lengths};
  use crate::optimize::run_loop::{find_zero_optimal_internal_edges, prune_and_merge_in_loop};
  use crate::optimize::dispatch::{initial_guess_mixed, run_optimize_mixed};
  use crate::optimize::run_loop::collect_optimize_partitions;
  use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::partition::payload::ancestral::GraphAncestral;
  use crate::partition::payload::sparse::{SparseEdgePartition, SparseNodePartition};
  use crate::seq::alignment::get_common_length;
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AsciiChar;
  use treetime_primitives::seq;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
    Sub::new(c(reff), pos, c(qry)).unwrap()
  }

  fn populate_test_nodes(partition: &mut PartitionMarginalSparse, graph: &GraphAncestral) {
    let ref_seq: treetime_primitives::Seq = std::iter::repeat_with(|| c(b'A')).take(partition.length).collect();
    if partition.root_sequence.is_empty() {
      partition.root_sequence = ref_seq.clone();
    }
    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      partition.nodes.entry(key).or_insert_with(|| {
        let mut node_part = SparseNodePartition::empty(&partition.alphabet);
        node_part.seq.sequence = ref_seq.clone();
        node_part
      });
    }
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_empty_graph() -> Result<(), Report> {
    let graph = GraphAncestral::new();
    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse);
    assert_eq!(edges.len(), 0);
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_no_zero_edges() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse);
    assert_eq!(edges.len(), 0);
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_skips_leaves() -> Result<(), Report> {
    // A has bl=0.0 but is a leaf: should NOT be collected
    let graph: GraphAncestral = nwk_read_str("(A:0.0,B:0.2)root;")?;
    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse);
    assert_eq!(edges.len(), 0);
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_collects_internal() -> Result<(), Report> {
    // I has bl=0.0 and is internal: should be collected
    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)I:0.0,C:0.3)root;")?;
    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse);
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
    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse);
    assert_eq!(edges.len(), 2);
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
      root_sequence: seq![],
    };

    populate_test_nodes(&mut partition, &graph);

    partition.edges.insert(ri_key, SparseEdgePartition::default());

    let ia_key = find_edge_key(&graph, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, "I", "B").unwrap();
    let rc_key = find_edge_key(&graph, "root", "C").unwrap();
    let rd_key = find_edge_key(&graph, "root", "D").unwrap();

    partition
      .edges
      .insert(ia_key, SparseEdgePartition::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .edges
      .insert(ib_key, SparseEdgePartition::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .edges
      .insert(rc_key, SparseEdgePartition::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .edges
      .insert(rd_key, SparseEdgePartition::with_fitch_subs(vec![sub(b'G', 5, b'C')]));

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

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_optimize_loop_with_topology_cleanup_sparse(#[case] method: BranchOptMethod) -> Result<(), Report> {
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

    let fitch = create_fitch_partition(&graph, 0, nuc, &aln)?;
    let sparse_partitions = vec![Arc::new(RwLock::new(fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?))];
    update_marginal(&graph, &sparse_partitions)?;

    let dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];
    let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);

    initial_guess_mixed(&graph, &mixed_partitions, true)?;

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
      run_optimize_mixed(&graph, &mixed_partitions, method)?;

      let zero_optimal_edges = find_zero_optimal_internal_edges(&graph, &sparse_partitions);

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

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_optimize_loop_no_collapse_when_branches_nonzero(#[case] method: BranchOptMethod) -> Result<(), Report> {
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

    let fitch = create_fitch_partition(&graph, 0, nuc, &aln)?;
    let sparse_partitions = vec![Arc::new(RwLock::new(fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?))];
    update_marginal(&graph, &sparse_partitions)?;

    let dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];
    let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);

    initial_guess_mixed(&graph, &mixed_partitions, true)?;

    let initial_node_count = graph.get_nodes().len();

    let mut lh_prev = f64::MIN;
    for i in 0..10 {
      let total_lh = update_marginal(&graph, &sparse_partitions)?;
      if (total_lh - lh_prev).abs() < 1e-2 {
        break;
      }

      let old_branch_lengths = save_branch_lengths(&graph);
      run_optimize_mixed(&graph, &mixed_partitions, method)?;

      let zero_optimal_edges = find_zero_optimal_internal_edges(&graph, &sparse_partitions);
      apply_damping(&graph, &old_branch_lengths, 0.75, i);
      prune_and_merge_in_loop(&mut graph, &sparse_partitions, &dense_partitions, &zero_optimal_edges)?;

      lh_prev = total_lh;
    }

    // No edges should have been collapsed - all branches carry genuine signal
    assert_eq!(graph.get_nodes().len(), initial_node_count);

    Ok(())
  }

  #[test]
  fn test_optimize_merge_then_marginal_finite_likelihood() -> Result<(), Report> {
    // After merge creates new internal nodes, update_marginal must produce
    // finite log-likelihood. This exercises the composition propagation fix:
    // merge-created nodes inherit the parent's composition so the backward
    // pass computes correct fixed-site contributions.
    //
    // Tree: star polytomy with 5 leaves. A and B share a derived state (T at
    // pos 0) while C, D, E retain the ancestral state (A at pos 0). With 3-vs-2
    // majority, the root's marginal MAP resolves to A. Edges to A and B carry
    // the shared mutation A->T, triggering merge.
    let nuc = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        TCGTACGTACGTACGT
        >B
        TCGTACGTACGTACGT
        >C
        ACGTACGTACGTACGT
        >D
        ACGTACGTACGTACGT
        >E
        ACGTACGTACGTACGT
      "#},
      &nuc,
    )?;

    let mut graph: GraphAncestral = nwk_read_str("(A:0.001,B:0.001,C:0.001,D:0.001,E:0.001)root:0.0;")?;

    let fitch = create_fitch_partition(&graph, 0, nuc, &aln)?;
    let sparse_partitions = vec![Arc::new(RwLock::new(
      fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?,
    ))];
    update_marginal(&graph, &sparse_partitions)?;

    let initial_node_count = graph.get_nodes().len();

    // A and B share mutation A->T at pos 0 (root MAP = A due to 3-vs-2 majority)
    let merged = merge_shared_mutation_branches(&mut graph, &sparse_partitions)?;
    assert!(merged > 0, "A and B should share mutation A->T, triggering merge");
    graph.build()?;

    assert!(
      graph.get_nodes().len() > initial_node_count,
      "merge should have created new internal nodes"
    );

    // The critical test: update_marginal after merge must produce finite log-likelihood.
    // Before the composition fix, the merge-created node had zero composition,
    // causing the backward pass to produce incorrect values.
    let lh = update_marginal(&graph, &sparse_partitions)?;
    assert!(lh.is_finite(), "log-likelihood must be finite after merge: {lh}");
    assert!(lh < 0.0, "log-likelihood must be negative: {lh}");

    Ok(())
  }

  #[test]
  fn test_optimize_cascading_collapse_parent_child_both_zero() -> Result<(), Report> {
    // Parent and child internal edges are both zero-optimal.
    // Tree: root -> I1 (bl=0.0) -> I2 (bl=0.0) -> A, B
    //       root -> C
    //
    // Both I1 and I2 should be collapsed. The guard for already-removed edges
    // must handle the case where collapsing I1 removes I2's inbound edge.
    let mut graph: GraphAncestral = nwk_read_str("(((A:0.1,B:0.1)I2:0.0)I1:0.0,C:0.1)root;")?;

    let ri1_key = find_edge_key(&graph, "root", "I1").unwrap();
    let i1i2_key = find_edge_key(&graph, "I1", "I2").unwrap();

    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let dense: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];

    let changed = prune_and_merge_in_loop(&mut graph, &sparse, &dense, &[ri1_key, i1i2_key])?;
    assert!(changed);

    // Both I1 and I2 should be gone. A, B become children of root.
    assert!(find_node_key_by_name(&graph, "I1").is_none());
    assert!(find_node_key_by_name(&graph, "I2").is_none());

    // root should have 3 children: A, B, C
    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let root_node = graph.get_node(root_key).unwrap();
    assert_eq!(root_node.read_arc().degree_out(), 3);

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_optimize_loop_with_topology_cleanup_dense(#[case] method: BranchOptMethod) -> Result<(), Report> {
    // Dense-mode integration test: identical sequences A and B should cause
    // the AB internal edge to be collapsed during optimization.
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

    let mut graph: GraphAncestral = nwk_read_str("((A:0.01,B:0.01)AB:0.01,(C:0.01,D:0.01)CD:0.01)root:0.0;")?;

    let dense_partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: nuc,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    initialize_marginal(&graph, &dense_partitions, &aln)?;
    update_marginal(&graph, &dense_partitions)?;

    let sparse_partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);

    initial_guess_mixed(&graph, &mixed_partitions, true)?;

    let initial_node_count = graph.get_nodes().len();

    let mut lh_prev = f64::MIN;
    for i in 0..10 {
      let dense_lh = update_marginal(&graph, &dense_partitions)?;
      if (dense_lh - lh_prev).abs() < 1e-2 {
        break;
      }

      let old_branch_lengths = save_branch_lengths(&graph);
      run_optimize_mixed(&graph, &mixed_partitions, method)?;

      let zero_optimal_edges = find_zero_optimal_internal_edges(&graph, &sparse_partitions);
      apply_damping(&graph, &old_branch_lengths, 0.75, i);
      prune_and_merge_in_loop(&mut graph, &sparse_partitions, &dense_partitions, &zero_optimal_edges)?;

      lh_prev = dense_lh;
    }

    // A and B are identical: AB edge should have been collapsed
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
}

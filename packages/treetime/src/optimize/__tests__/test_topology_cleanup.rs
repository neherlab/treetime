#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::{initial_guess_mixed, run_optimize_mixed};
  use crate::optimize::iteration::apply_damping;
  use crate::optimize::params::{BranchOptMethod, TopologyOps};
  use crate::optimize::run_loop::{
    OptimizeReadouts, find_zero_optimal_internal_edges, marginal_update_dense, marginal_update_sparse,
    prune_and_merge_in_loop, run_optimize_loop,
  };
  use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
  use crate::seq::alignment::get_common_length;
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::{NwkParse, nwk_read_str};
  use treetime_primitives::AsciiChar;
  use treetime_primitives::seq;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
    Sub::new(c(reff), pos, c(qry)).unwrap()
  }

  fn empty_sparse_recon() -> Result<SparseReconstruction, Report> {
    Ok(SparseReconstruction {
      partition: PartitionMarginalSparse {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet: Alphabet::new(AlphabetName::Nuc)?,
        length: 100,
        root_sequence: seq![],
        obs_nodes: btreemap! {},
        obs_edges: btreemap! {},
      },
      node_states: btreemap! {},
      backward: btreemap! {},
      forward: btreemap! {},
      estimates: btreemap! {},
    })
  }

  fn populate_test_nodes(recon: &mut SparseReconstruction, graph: &Graph) {
    let ref_seq: treetime_primitives::Seq = std::iter::repeat_with(|| c(b'A'))
      .take(recon.partition.length)
      .collect();
    if recon.partition.root_sequence.is_empty() {
      recon.partition.root_sequence = ref_seq.clone();
    }
    let alphabet = recon.partition.alphabet.clone();
    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      recon
        .partition
        .obs_nodes
        .entry(key)
        .or_insert_with(|| SparseNodeObs::empty(&alphabet));
      recon
        .node_states
        .entry(key)
        .or_insert_with(|| SparseNodeState::leaf(&ref_seq));
    }
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_empty_graph() -> Result<(), Report> {
    let graph = Graph::new();
    let sparse: Vec<SparseReconstruction> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse, &BTreeMap::new());
    assert_eq!(edges.len(), 0);
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_no_zero_edges() -> Result<(), Report> {
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let graph: Graph = graph;
    let sparse: Vec<SparseReconstruction> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse, &branch_lengths);
    assert_eq!(edges.len(), 0);
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_skips_leaves() -> Result<(), Report> {
    // A has bl=0.0 but is a leaf: should NOT be collected
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("(A:0.0,B:0.2)root;")?;
    let graph: Graph = graph;
    let sparse: Vec<SparseReconstruction> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse, &branch_lengths);
    assert_eq!(edges.len(), 0);
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_collects_internal() -> Result<(), Report> {
    // I has bl=0.0 and is internal: should be collected
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.2)I:0.0,C:0.3)root;")?;
    let graph: Graph = graph;
    let sparse: Vec<SparseReconstruction> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse, &branch_lengths);
    assert_eq!(edges.len(), 1);
    let edge_key = edges[0];
    let edge = graph.get_edge(edge_key).unwrap();
    let target = edge.read_arc().target();
    let target_name = graph
      .get_node(target)
      .and_then(|n| names.get(&n.read_arc().key()).cloned().flatten());
    assert_eq!(target_name.as_deref(), Some("I"));
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_multiple() -> Result<(), Report> {
    // Both internal nodes have bl=0.0
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("(((A:0.1,B:0.1)I1:0.0,C:0.1)I2:0.0,D:0.1)root;")?;
    let graph: Graph = graph;
    let sparse: Vec<SparseReconstruction> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse, &branch_lengths);
    assert_eq!(edges.len(), 2);
    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_empty_list() -> Result<(), Report> {
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let mut graph: Graph = graph;
    let mut sparse: Vec<SparseReconstruction> = vec![];
    let mut dense: Vec<DenseReconstruction> = vec![];

    let mut names_tt_13 = names;
    let changed = prune_and_merge_in_loop(
      &mut graph,
      &mut sparse,
      &mut dense,
      &[],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_13,
    )?;
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
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.1)I:0.0,C:0.1,D:0.1)root;")?;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();

    let mut partition = empty_sparse_recon()?;

    populate_test_nodes(&mut partition, &graph);

    partition.partition.obs_edges.insert(ri_key, SparseEdgeObs::default());

    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();
    let rc_key = find_edge_key(&graph, &names, "root", "C").unwrap();
    let rd_key = find_edge_key(&graph, &names, "root", "D").unwrap();

    partition
      .partition
      .obs_edges
      .insert(ia_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(ib_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(rc_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(rd_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'G', 5, b'C')]));

    let mut sparse = vec![partition];
    let mut dense: Vec<DenseReconstruction> = vec![];

    // Set bl to 0.0 and pass damped value (simulating damping override in prune_and_merge_in_loop)
    branch_lengths.insert(ri_key, Some(0.0));

    let mut names_tt_12 = names.clone();
    let changed = prune_and_merge_in_loop(
      &mut graph,
      &mut sparse,
      &mut dense,
      &[ri_key],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_12,
    )?;
    assert!(changed);

    // I should be gone
    assert!(find_node_key_by_name(&graph, &names, "I").is_none());

    // D remains directly under root
    assert!(find_node_key_by_name(&graph, &names, "D").is_some());

    // Root should have 2 children after merging A, B, C into a new subtree
    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
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
    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str("((A:0.01,B:0.01)AB:0.01,(C:0.01,D:0.01)CD:0.01)root:0.0;")?;
    let mut graph: Graph = graph;

    let fitch = create_fitch_partition(&graph, 0, nuc, &aln, &names)?;
    let (sp_partition, sp_node_states) = fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?;
    let mut sparse_partitions = vec![SparseReconstruction { partition: sp_partition, node_states: sp_node_states, backward: btreemap!{}, forward: btreemap!{}, estimates: btreemap!{} }];
    marginal_update_sparse(&graph, &profile_branch_lengths(&branch_lengths), &mut sparse_partitions)?.value();

    let mut dense_partitions: Vec<DenseReconstruction> = vec![];

    initial_guess_mixed(&graph, &OptimizeReadouts::new(&dense_partitions, &sparse_partitions).view(), true, false, &mut branch_lengths)?;

    let initial_node_count = graph.get_nodes().len();

    // Run optimize loop with topology cleanup
    let mut lh_prev = f64::MIN;
    for i in 0..10 {
      let sparse_lh = marginal_update_sparse(&graph, &profile_branch_lengths(&branch_lengths), &mut sparse_partitions)?.value();
      let total_lh = sparse_lh;

      if (total_lh - lh_prev).abs() < 1e-2 {
        break;
      }

      let old_branch_lengths = branch_lengths.clone();
      run_optimize_mixed(&graph, &OptimizeReadouts::new(&dense_partitions, &sparse_partitions).view(), method, &mut branch_lengths)?;
        let zero_optimal_edges = find_zero_optimal_internal_edges(&graph, &sparse_partitions, &branch_lengths);
      apply_damping(&mut branch_lengths, &old_branch_lengths, 0.75, i);
      let mut names_tt_11 = names.clone();
      prune_and_merge_in_loop(&mut graph, &mut sparse_partitions, &mut dense_partitions, &zero_optimal_edges, TopologyOps::default(), &mut branch_lengths, &mut names_tt_11)?;

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
      let bl = branch_lengths.get(&edge.read_arc().key()).copied().flatten().unwrap_or(0.0);
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

    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str("((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;")?;

    let mut graph: Graph = graph;

    let fitch = create_fitch_partition(&graph, 0, nuc, &aln, &names)?;
    let (sp_partition, sp_node_states) = fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?;
    let mut sparse_partitions = vec![SparseReconstruction { partition: sp_partition, node_states: sp_node_states, backward: btreemap!{}, forward: btreemap!{}, estimates: btreemap!{} }];
    marginal_update_sparse(&graph, &profile_branch_lengths(&branch_lengths), &mut sparse_partitions)?.value();

    let mut dense_partitions: Vec<DenseReconstruction> = vec![];

    initial_guess_mixed(&graph, &OptimizeReadouts::new(&dense_partitions, &sparse_partitions).view(), true, false, &mut branch_lengths)?;

    let initial_node_count = graph.get_nodes().len();

    let mut lh_prev = f64::MIN;
    for i in 0..10 {
      let total_lh = marginal_update_sparse(&graph, &profile_branch_lengths(&branch_lengths), &mut sparse_partitions)?.value();
      if (total_lh - lh_prev).abs() < 1e-2 {
        break;
      }

      let old_branch_lengths = branch_lengths.clone();
      run_optimize_mixed(&graph, &OptimizeReadouts::new(&dense_partitions, &sparse_partitions).view(), method, &mut branch_lengths)?;
        let zero_optimal_edges = find_zero_optimal_internal_edges(&graph, &sparse_partitions, &branch_lengths);
      apply_damping(&mut branch_lengths, &old_branch_lengths, 0.75, i);
      let mut names_tt_10 = names.clone();
      prune_and_merge_in_loop(&mut graph, &mut sparse_partitions, &mut dense_partitions, &zero_optimal_edges, TopologyOps::default(), &mut branch_lengths, &mut names_tt_10)?;

      lh_prev = total_lh;
    }

    // No edges should have been collapsed - all branches carry genuine signal
    assert_eq!(graph.get_nodes().len(), initial_node_count);

    Ok(())
  }

  #[test]
  fn test_optimize_merge_then_marginal_finite_likelihood() -> Result<(), Report> {
    // After merge creates new internal nodes, marginal_update must produce
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

    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str("(A:0.001,B:0.001,C:0.001,D:0.001,E:0.001)root:0.0;")?;

    let mut graph: Graph = graph;

    let fitch = create_fitch_partition(&graph, 0, nuc, &aln, &names)?;
    let (sp_partition, sp_node_states) = fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?;
    let mut sparse_partitions = vec![SparseReconstruction {
      partition: sp_partition,
      node_states: sp_node_states,
      backward: btreemap! {},
      forward: btreemap! {},
      estimates: btreemap! {},
    }];
    marginal_update_sparse(&graph, &profile_branch_lengths(&branch_lengths), &mut sparse_partitions)?.value();

    let initial_node_count = graph.get_nodes().len();

    // A and B share mutation A->T at pos 0 (root MAP = A due to 3-vs-2 majority)
    let merged = merge_shared_mutation_branches(&mut graph, &mut sparse_partitions, &mut branch_lengths)?;
    assert!(merged > 0, "A and B should share mutation A->T, triggering merge");
    graph.build()?;

    assert!(
      graph.get_nodes().len() > initial_node_count,
      "merge should have created new internal nodes"
    );

    // The critical test: marginal_update after merge must produce finite log-likelihood.
    // Before the composition fix, the merge-created node had zero composition,
    // causing the backward pass to produce incorrect values.
    let lh = marginal_update_sparse(&graph, &profile_branch_lengths(&branch_lengths), &mut sparse_partitions)?.value();
    assert!(lh.is_finite(), "log-likelihood must be finite after merge: {lh}");
    assert!(lh < 0.0, "log-likelihood must be negative: {lh}");

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_hoists_reversion_without_collapse() -> Result<(), Report> {
    // Reversion polytomy with no zero-optimal edge to collapse. The loop must still resolve
    // it: merge C1+C2 (shared reversion), hoist the reverting group, retire the helper.
    // Tree: root -> U -> V -> {C1, C2, C3}. U->V carries {A0T, C5G}; C1 and C2 revert A0T,
    // C3 keeps it. Parsimony optimum is 2 mutations.
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str("(((C1:0.1,C2:0.1,C3:0.1)V:0.2)U:0.1)root:0.0;")?;
    let mut graph: Graph = graph;

    let mut partition = empty_sparse_recon()?;
    populate_test_nodes(&mut partition, &graph);

    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let vc1 = find_edge_key(&graph, &names, "V", "C1").unwrap();
    let vc2 = find_edge_key(&graph, &names, "V", "C2").unwrap();
    let vc3 = find_edge_key(&graph, &names, "V", "C3").unwrap();
    partition.partition.obs_edges.insert(
      uv,
      SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T'), sub(b'C', 5, b'G')]),
    );
    partition
      .partition
      .obs_edges
      .insert(vc1, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));
    partition
      .partition
      .obs_edges
      .insert(vc2, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));
    partition.partition.obs_edges.insert(vc3, SparseEdgeObs::default());

    let mut sparse = vec![partition];
    let mut dense: Vec<DenseReconstruction> = vec![];

    // Empty zero-optimal list: the old loop was a no-op here. The hoist must still fire.
    let mut names_tt_9 = names;
    let changed = prune_and_merge_in_loop(
      &mut graph,
      &mut sparse,
      &mut dense,
      &[],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_9,
    )?;
    assert!(changed, "reversion polytomy must be resolved even without a collapse");

    let p = &sparse[0];
    let total_subs: usize = graph
      .get_edges()
      .iter()
      .filter_map(|e| p.partition.obs_edges.get(&e.read_arc().key()))
      .map(|e| e.fitch_subs().len())
      .sum();
    assert_eq!(total_subs, 2, "reaches the parsimony optimum");

    let reversion_remains = graph
      .get_edges()
      .iter()
      .filter_map(|e| p.partition.obs_edges.get(&e.read_arc().key()))
      .any(|e| e.fitch_subs().contains(&sub(b'T', 0, b'A')));
    assert!(!reversion_remains, "reversion must be removed");

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
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str("(((A:0.1,B:0.1)I2:0.0)I1:0.0,C:0.1)root;")?;
    let mut graph: Graph = graph;

    let ri1_key = find_edge_key(&graph, &names, "root", "I1").unwrap();
    let i1i2_key = find_edge_key(&graph, &names, "I1", "I2").unwrap();

    let mut sparse: Vec<SparseReconstruction> = vec![];
    let mut dense: Vec<DenseReconstruction> = vec![];

    let mut names_tt_8 = names.clone();
    let changed = prune_and_merge_in_loop(
      &mut graph,
      &mut sparse,
      &mut dense,
      &[ri1_key, i1i2_key],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_8,
    )?;
    assert!(changed);

    // Both I1 and I2 should be gone. A, B become children of root.
    assert!(find_node_key_by_name(&graph, &names, "I1").is_none());
    assert!(find_node_key_by_name(&graph, &names, "I2").is_none());

    // root should have 3 children: A, B, C
    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
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

    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str("((A:0.01,B:0.01)AB:0.01,(C:0.01,D:0.01)CD:0.01)root:0.0;")?;

    let mut graph: Graph = graph;

    let dense_partition = PartitionMarginalDense::new(0, jc69(JC69Params::default())?, nuc, get_common_length(&aln)?);
    let dense_node_states = dense_partition.attach_sequences(&graph, &aln, &names)?;
    let mut dense_partitions = vec![DenseReconstruction {
      partition: dense_partition,
      node_states: dense_node_states,
      backward: btreemap! {},
      forward: btreemap! {},
      estimates: btreemap! {},
    }];

    marginal_update_dense(&graph, &profile_branch_lengths(&branch_lengths), &mut dense_partitions)?.value();

    let mut sparse_partitions: Vec<SparseReconstruction> = vec![];

    initial_guess_mixed(&graph, &OptimizeReadouts::new(&dense_partitions, &sparse_partitions).view(), true, false, &mut branch_lengths)?;

    let initial_node_count = graph.get_nodes().len();

    let mut lh_prev = f64::MIN;
    for i in 0..10 {
      let dense_lh = marginal_update_dense(&graph, &profile_branch_lengths(&branch_lengths), &mut dense_partitions)?.value();
      if (dense_lh - lh_prev).abs() < 1e-2 {
        break;
      }

      let old_branch_lengths = branch_lengths.clone();
      run_optimize_mixed(&graph, &OptimizeReadouts::new(&dense_partitions, &sparse_partitions).view(), method, &mut branch_lengths)?;
        let zero_optimal_edges = find_zero_optimal_internal_edges(&graph, &sparse_partitions, &branch_lengths);
      apply_damping(&mut branch_lengths, &old_branch_lengths, 0.75, i);
      let mut names_tt_7 = names.clone();
      prune_and_merge_in_loop(&mut graph, &mut sparse_partitions, &mut dense_partitions, &zero_optimal_edges, TopologyOps::default(), &mut branch_lengths, &mut names_tt_7)?;

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
      let bl = branch_lengths.get(&edge.read_arc().key()).copied().flatten().unwrap_or(0.0);
      assert!(bl >= 0.0, "Negative branch length after optimization: {bl}");
    }

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_names_new_nodes() -> Result<(), Report> {
    // Collapse zero-length I, then merge A+B+C (shared sub A0T) under a new
    // internal node. The new node must receive a NODE_NNNNNNN name.
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.1)I:0.0,C:0.1,D:0.1)root;")?;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();

    let mut partition = empty_sparse_recon()?;

    populate_test_nodes(&mut partition, &graph);

    partition.partition.obs_edges.insert(ri_key, SparseEdgeObs::default());

    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();
    let rc_key = find_edge_key(&graph, &names, "root", "C").unwrap();
    let rd_key = find_edge_key(&graph, &names, "root", "D").unwrap();

    partition
      .partition
      .obs_edges
      .insert(ia_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(ib_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(rc_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(rd_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'G', 5, b'C')]));

    let mut sparse = vec![partition];
    let mut dense: Vec<DenseReconstruction> = vec![];

    branch_lengths.insert(ri_key, Some(0.0));

    let mut names_tt_6 = names;
    let changed = prune_and_merge_in_loop(
      &mut graph,
      &mut sparse,
      &mut dense,
      &[ri_key],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_6,
    )?;
    assert!(changed);

    let mut names: Vec<String> = graph
      .get_nodes()
      .iter()
      .filter_map(|n| names_tt_6.get(&n.read_arc().key()).cloned().flatten())
      .collect();
    names.sort();

    // I collapsed, A+B+C merged under a new NODE_0000000
    assert_eq!(names, vec!["A", "B", "C", "D", "NODE_0000000", "root"]);

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_merge_disabled_keeps_polytomy() -> Result<(), Report> {
    // Same setup as the collapse+merge test, but with merge-siblings disabled. Collapsing the
    // zero-length internal edge still forms the polytomy; without merge, the shared-mutation
    // siblings A, B, C stay as direct children of root rather than being grouped under a node.
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.1)I:0.0,C:0.1,D:0.1)root;")?;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();

    let mut partition = empty_sparse_recon()?;
    populate_test_nodes(&mut partition, &graph);
    partition.partition.obs_edges.insert(ri_key, SparseEdgeObs::default());

    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();
    let rc_key = find_edge_key(&graph, &names, "root", "C").unwrap();
    let rd_key = find_edge_key(&graph, &names, "root", "D").unwrap();
    partition
      .partition
      .obs_edges
      .insert(ia_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(ib_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(rc_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(rd_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'G', 5, b'C')]));

    let mut sparse = vec![partition];
    let mut dense: Vec<DenseReconstruction> = vec![];

    branch_lengths.insert(ri_key, Some(0.0));

    let ops = TopologyOps {
      merge_siblings: false,
      ..TopologyOps::default()
    };
    let mut names_tt_5 = names.clone();
    let changed = prune_and_merge_in_loop(
      &mut graph,
      &mut sparse,
      &mut dense,
      &[ri_key],
      ops,
      &mut branch_lengths,
      &mut names_tt_5,
    )?;
    assert!(changed, "collapse still fires even with merge disabled");

    // I collapsed away.
    assert!(find_node_key_by_name(&graph, &names, "I").is_none());

    // Without merge, root keeps all four children A, B, C, D (no grouping under a new node).
    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let root_node = graph.get_node(root_key).unwrap();
    assert_eq!(root_node.read_arc().degree_out(), 4);

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_flip_disabled_keeps_reversion() -> Result<(), Report> {
    // Reversion polytomy. With flip-parent-child disabled, merge still groups the two reverting
    // children, but the reverting mutation is not hoisted away, so it remains in the tree.
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str("(((C1:0.1,C2:0.1,C3:0.1)V:0.2)U:0.1)root:0.0;")?;
    let mut graph: Graph = graph;

    let mut partition = empty_sparse_recon()?;
    populate_test_nodes(&mut partition, &graph);

    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let vc1 = find_edge_key(&graph, &names, "V", "C1").unwrap();
    let vc2 = find_edge_key(&graph, &names, "V", "C2").unwrap();
    let vc3 = find_edge_key(&graph, &names, "V", "C3").unwrap();
    partition.partition.obs_edges.insert(
      uv,
      SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T'), sub(b'C', 5, b'G')]),
    );
    partition
      .partition
      .obs_edges
      .insert(vc1, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));
    partition
      .partition
      .obs_edges
      .insert(vc2, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));
    partition.partition.obs_edges.insert(vc3, SparseEdgeObs::default());

    let mut sparse = vec![partition];
    let mut dense: Vec<DenseReconstruction> = vec![];

    let ops = TopologyOps {
      flip_parent_child: false,
      ..TopologyOps::default()
    };
    let mut names_tt_4 = names;
    let changed = prune_and_merge_in_loop(
      &mut graph,
      &mut sparse,
      &mut dense,
      &[],
      ops,
      &mut branch_lengths,
      &mut names_tt_4,
    )?;
    assert!(changed, "merge still groups the reverting siblings");

    let p = &sparse[0];
    let reversion_remains = graph
      .get_edges()
      .iter()
      .filter_map(|e| p.partition.obs_edges.get(&e.read_arc().key()))
      .any(|e| e.fitch_subs().contains(&sub(b'T', 0, b'A')));
    assert!(
      reversion_remains,
      "reversion is kept when flip-parent-child is disabled"
    );

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_all_ops_disabled_is_noop() -> Result<(), Report> {
    // With every topology step disabled, the reversion polytomy is left untouched: no collapse,
    // no merge, no hoist. The tree shape and its mutation content are unchanged.
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str("(((C1:0.1,C2:0.1,C3:0.1)V:0.2)U:0.1)root:0.0;")?;
    let mut graph: Graph = graph;

    let mut partition = empty_sparse_recon()?;
    populate_test_nodes(&mut partition, &graph);

    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let vc1 = find_edge_key(&graph, &names, "V", "C1").unwrap();
    let vc2 = find_edge_key(&graph, &names, "V", "C2").unwrap();
    let vc3 = find_edge_key(&graph, &names, "V", "C3").unwrap();
    partition.partition.obs_edges.insert(
      uv,
      SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T'), sub(b'C', 5, b'G')]),
    );
    partition
      .partition
      .obs_edges
      .insert(vc1, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));
    partition
      .partition
      .obs_edges
      .insert(vc2, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));
    partition.partition.obs_edges.insert(vc3, SparseEdgeObs::default());

    let mut sparse = vec![partition];
    let mut dense: Vec<DenseReconstruction> = vec![];

    let node_count_before = graph.get_nodes().len();
    let ops = TopologyOps {
      collapse_short_branches: false,
      merge_siblings: false,
      flip_parent_child: false,
    };
    let mut names_tt_3 = names;
    let changed = prune_and_merge_in_loop(
      &mut graph,
      &mut sparse,
      &mut dense,
      &[],
      ops,
      &mut branch_lengths,
      &mut names_tt_3,
    )?;
    assert!(!changed, "no topology step runs when all are disabled");
    assert_eq!(graph.get_nodes().len(), node_count_before);

    let p = &sparse[0];
    let total_subs: usize = graph
      .get_edges()
      .iter()
      .filter_map(|e| p.partition.obs_edges.get(&e.read_arc().key()))
      .map(|e| e.fitch_subs().len())
      .sum();
    assert_eq!(total_subs, 4, "mutation content is unchanged");

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::brent_sqrt(BranchOptMethod::BrentSqrt)]
  #[case::newton(    BranchOptMethod::Newton)]
  #[trace]
  fn test_run_optimize_loop_collapse_disabled_keeps_zero_edge(#[case] method: BranchOptMethod) -> Result<(), Report> {
    // A and B are identical, so the AB internal edge optimizes toward zero. With collapse
    // disabled, run_optimize_loop must leave that edge in place: the node count is unchanged.
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

    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str("((A:0.01,B:0.01)AB:0.01,(C:0.01,D:0.01)CD:0.01)root:0.0;")?;

    let mut graph: Graph = graph;

    let fitch = create_fitch_partition(&graph, 0, nuc, &aln, &names)?;
    let (sp_partition, sp_node_states) = fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?;
    let mut sparse_partitions = vec![SparseReconstruction { partition: sp_partition, node_states: sp_node_states, backward: btreemap!{}, forward: btreemap!{}, estimates: btreemap!{} }];
    marginal_update_sparse(&graph, &profile_branch_lengths(&branch_lengths), &mut sparse_partitions)?.value();

    let mut dense_partitions: Vec<DenseReconstruction> = vec![];
    initial_guess_mixed(&graph, &OptimizeReadouts::new(&dense_partitions, &sparse_partitions).view(), true, false, &mut branch_lengths)?;

    let initial_node_count = graph.get_nodes().len();

    let ops = TopologyOps {
      collapse_short_branches: false,
      ..TopologyOps::default()
    };
    let names_tt_2 = names;
    run_optimize_loop(
      &mut graph,
      &mut sparse_partitions,
      &mut dense_partitions,
      10,
      1e-2,
      0.75,
      method,
      false,
      ops,
      branch_lengths,
      &names_tt_2,
    )?;

    assert_eq!(
      graph.get_nodes().len(),
      initial_node_count,
      "no collapse when collapse_short_branches is disabled"
    );

    Ok(())
  }

  // Rollback validation (T1.4 gate). When a topology change fires, the internal best-branch-length
  // map is discarded (best is reset to `None`, best_lh to IMPOSSIBLE), so no rollback can restore
  // lengths keyed to the superseded tree. The observable guarantee is that the returned map is
  // keyed exactly by the post-change edge set, with no stale keys and no missing edges.
  // Oracle: the topology-change branch in `run_optimize_loop` resets the best and the topology
  // producers keep the map in step with the current edge set.
  #[test]
  fn test_run_optimize_loop_topology_change_map_matches_edge_set() -> Result<(), Report> {
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

    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str("((A:0.01,B:0.01)AB:0.01,(C:0.01,D:0.01)CD:0.01)root:0.0;")?;

    let mut graph: Graph = graph;

    let fitch = create_fitch_partition(&graph, 0, nuc, &aln, &names)?;
    let (sp_partition, sp_node_states) = fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?;
    let mut sparse_partitions = vec![SparseReconstruction {
      partition: sp_partition,
      node_states: sp_node_states,
      backward: btreemap! {},
      forward: btreemap! {},
      estimates: btreemap! {},
    }];
    marginal_update_sparse(&graph, &profile_branch_lengths(&branch_lengths), &mut sparse_partitions)?.value();

    let mut dense_partitions: Vec<DenseReconstruction> = vec![];
    initial_guess_mixed(
      &graph,
      &OptimizeReadouts::new(&dense_partitions, &sparse_partitions).view(),
      true,
      false,
      &mut branch_lengths,
    )?;

    let initial_node_count = graph.get_nodes().len();

    let names_tt_1 = names;
    let result = run_optimize_loop(
      &mut graph,
      &mut sparse_partitions,
      &mut dense_partitions,
      10,
      1e-2,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_1,
    )?;

    // A and B are identical, so the AB internal edge collapses: topology changed.
    assert!(
      graph.get_nodes().len() < initial_node_count,
      "expected a collapse (topology change) with collapse enabled"
    );

    // The returned map is keyed exactly by the post-change edge set: no stale keys, no gaps.
    let map_keys: Vec<_> = result.branch_lengths.keys().copied().collect();
    let mut edge_keys: Vec<_> = graph.get_edges().iter().map(|edge| edge.read_arc().key()).collect();
    edge_keys.sort_unstable();
    assert_eq!(map_keys, edge_keys);

    for bl in result.branch_lengths.values().flatten() {
      assert!(
        bl.is_finite() && *bl >= 0.0,
        "branch length must be finite and non-negative: {bl}"
      );
    }
    Ok(())
  }
}

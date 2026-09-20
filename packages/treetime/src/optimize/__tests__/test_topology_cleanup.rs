#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::{initial_guess_mixed, run_optimize_mixed};
  use crate::optimize::gather::{
    gather_edge_contributions, gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts,
    total_sequence_length,
  };
  use crate::optimize::iteration::apply_damping;
  use crate::optimize::params::{BranchOptMethod, TopologyOps};
  use crate::optimize::run_loop::{
    find_zero_optimal_internal_edges, marginal_update_dense, marginal_update_sparse, prune_and_merge_in_loop,
    run_optimize_loop,
  };
  use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::shared::reconcile::{live_node_keys, reconcile_node_states};
  use crate::partition::marginal::shared::update::MarginalEdges;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use eyre::Report;
  use indoc::indoc;
  use itertools::{Itertools, izip};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;
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
        alphabet: Alphabet::new(AlphabetName::Nuc)?,
        length: 100,
        root_sequence: seq![],
        obs_nodes: btreemap! {},
        obs_edges: btreemap! {},
      },
      gtr: jc69(JC69Params::default())?,
      node_states: btreemap! {},
      edges: MarginalEdges::default(),
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
      let key = node.key();
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
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let sparse: Vec<SparseReconstruction> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse, &branch_lengths);
    assert_eq!(edges.len(), 0);
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_skips_leaves() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.0,B:0.2)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let sparse: Vec<SparseReconstruction> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse, &branch_lengths);
    assert_eq!(edges.len(), 0);
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_collects_internal() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)I:0.0,C:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let sparse: Vec<SparseReconstruction> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse, &branch_lengths);
    assert_eq!(edges.len(), 1);
    let edge_key = edges[0];
    let edge = graph.get_edge(edge_key).unwrap();
    let target = edge.target();
    let target_name = graph
      .get_node(target)
      .and_then(|n| names.get(&n.key()).cloned().flatten());
    assert_eq!(target_name.as_deref(), Some("I"));
    Ok(())
  }

  #[test]
  fn test_optimize_find_zero_optimal_internal_edges_multiple() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(((A:0.1,B:0.1)I1:0.0,C:0.1)I2:0.0,D:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let sparse: Vec<SparseReconstruction> = vec![];
    let edges = find_zero_optimal_internal_edges(&graph, &sparse, &branch_lengths);
    assert_eq!(edges.len(), 2);
    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_empty_list() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let sparse: Vec<SparseReconstruction> = vec![];
    let dense: Vec<DenseReconstruction> = vec![];

    let mut names_tt_13 = names;
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_13,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(!changed);
    assert_eq!(graph.get_nodes().count(), 4);
    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_collapses_and_merges() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)I:0.0,C:0.1,D:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
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

    let sparse = vec![partition];
    let dense: Vec<DenseReconstruction> = vec![];

    branch_lengths.insert(ri_key, Some(0.0));

    let mut names_tt_12 = names.clone();
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[ri_key],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_12,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(changed);

    assert!(find_node_key_by_name(&graph, &names, "I").is_none());

    assert!(find_node_key_by_name(&graph, &names, "D").is_some());

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let root_node = graph.get_node(root_key).unwrap();
    assert_eq!(root_node.degree_out(), 2);

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
    let nuc = Alphabet::new(AlphabetName::Nuc)?;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.01,B:0.01)AB:0.01,(C:0.01,D:0.01)CD:0.01)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let fitch = create_fitch_partition(&graph, 0, nuc, &node_seq_inputs(&graph, &names, aln))?;
    let (sp_partition, sp_node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partitions = vec![SparseReconstruction::seeded(sp_partition, jc69(JC69Params::default())?, sp_node_states)];
    let (mut sparse_partitions, _) = marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;

    let mut dense_partitions: Vec<DenseReconstruction> = vec![];

    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);
    let sub_counts = gather_edge_sub_counts(&graph, &dense_partitions, &sparse_partitions)?;
    let effective_lengths = gather_edge_effective_lengths(&graph, &dense_partitions, &sparse_partitions)?;
    initial_guess_mixed(&graph, total_length, &indel_counts, &sub_counts, &effective_lengths, true, false, &mut branch_lengths)?;

    let initial_node_count = graph.get_nodes().count();

    let mut lh_prev = f64::MIN;
    for i in 0..10 {
      let sparse_lh;
      (sparse_partitions, sparse_lh) = marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;
      let sparse_lh = sparse_lh.value();
      let total_lh = sparse_lh;

      if (total_lh - lh_prev).abs() < 1e-2 {
        break;
      }

      let old_branch_lengths = branch_lengths.clone();
      let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
      let contributions = gather_edge_contributions(&graph, &dense_partitions, &sparse_partitions)?;
      let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
        let zero_optimal_edges = find_zero_optimal_internal_edges(&graph, &sparse_partitions, &branch_lengths);
      apply_damping(&mut branch_lengths, &old_branch_lengths, 0.75, i);
      let mut names_tt_11 = names.clone();
      let cleanup = prune_and_merge_in_loop(&mut graph, sparse_partitions, dense_partitions, &zero_optimal_edges, TopologyOps::default(), &mut branch_lengths, &mut names_tt_11)?;
      sparse_partitions = cleanup.sparse_partitions;
      dense_partitions = cleanup.dense_partitions;

      lh_prev = total_lh;
    }

    let final_node_count = graph.get_nodes().count();
    assert!(
      final_node_count < initial_node_count,
      "Expected topology simplification: {initial_node_count} nodes -> {final_node_count} nodes"
    );

    for edge in graph.get_edges() {
      let bl = branch_lengths.get(&edge.key()).copied().flatten().unwrap_or(0.0);
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
    let nuc = Alphabet::new(AlphabetName::Nuc)?;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let mut graph: Graph = graph;

    let fitch = create_fitch_partition(&graph, 0, nuc, &node_seq_inputs(&graph, &names, aln))?;
    let (sp_partition, sp_node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partitions = vec![SparseReconstruction::seeded(sp_partition, jc69(JC69Params::default())?, sp_node_states)];
    let (mut sparse_partitions, _) = marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;

    let mut dense_partitions: Vec<DenseReconstruction> = vec![];

    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);
    let sub_counts = gather_edge_sub_counts(&graph, &dense_partitions, &sparse_partitions)?;
    let effective_lengths = gather_edge_effective_lengths(&graph, &dense_partitions, &sparse_partitions)?;
    initial_guess_mixed(&graph, total_length, &indel_counts, &sub_counts, &effective_lengths, true, false, &mut branch_lengths)?;

    let initial_node_count = graph.get_nodes().count();

    let mut lh_prev = f64::MIN;
    for i in 0..10 {
      let total_lh;
      (sparse_partitions, total_lh) = marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;
      let total_lh = total_lh.value();
      if (total_lh - lh_prev).abs() < 1e-2 {
        break;
      }

      let old_branch_lengths = branch_lengths.clone();
      let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
      let contributions = gather_edge_contributions(&graph, &dense_partitions, &sparse_partitions)?;
      let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
        let zero_optimal_edges = find_zero_optimal_internal_edges(&graph, &sparse_partitions, &branch_lengths);
      apply_damping(&mut branch_lengths, &old_branch_lengths, 0.75, i);
      let mut names_tt_10 = names.clone();
      let cleanup = prune_and_merge_in_loop(&mut graph, sparse_partitions, dense_partitions, &zero_optimal_edges, TopologyOps::default(), &mut branch_lengths, &mut names_tt_10)?;
      sparse_partitions = cleanup.sparse_partitions;
      dense_partitions = cleanup.dense_partitions;

      lh_prev = total_lh;
    }

    assert_eq!(graph.get_nodes().count(), initial_node_count);

    Ok(())
  }

  #[test]
  fn test_optimize_merge_then_marginal_finite_likelihood() -> Result<(), Report> {
    let nuc = Alphabet::new(AlphabetName::Nuc)?;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("(A:0.001,B:0.001,C:0.001,D:0.001,E:0.001)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let mut graph: Graph = graph;

    let fitch = create_fitch_partition(&graph, 0, nuc, &node_seq_inputs(&graph, &names, aln))?;
    let (sp_partition, sp_node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partitions = vec![SparseReconstruction::seeded(
      sp_partition,
      jc69(JC69Params::default())?,
      sp_node_states,
    )];
    let (sparse_partitions, _) =
      marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;

    let initial_node_count = graph.get_nodes().count();

    let (mut sparse_obs, sparse_gtrs, sparse_node_states): (Vec<_>, Vec<_>, Vec<_>) = sparse_partitions
      .into_iter()
      .map(|family| (family.partition, family.gtr, family.node_states))
      .multiunzip();
    let merged = merge_shared_mutation_branches(&mut graph, &mut sparse_obs, &mut branch_lengths)?;
    assert!(merged > 0, "A and B should share mutation A->T, triggering merge");
    graph.build()?;

    assert!(
      graph.get_nodes().count() > initial_node_count,
      "merge should have created new internal nodes"
    );

    let live_nodes = live_node_keys(&graph);
    let sparse_partitions = izip!(sparse_obs, sparse_gtrs, sparse_node_states)
      .map(|(partition, gtr, node_states)| {
        SparseReconstruction::seeded(
          partition,
          gtr,
          reconcile_node_states(node_states, &live_nodes, SparseNodeState::empty),
        )
      })
      .collect_vec();

    let (sparse_partitions, lh) =
      marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;
    let lh = lh.value();
    assert!(lh.is_finite(), "log-likelihood must be finite after merge: {lh}");
    assert!(lh < 0.0, "log-likelihood must be negative: {lh}");

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_hoists_reversion_without_collapse() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(((C1:0.1,C2:0.1,C3:0.1)V:0.2)U:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
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

    let sparse = vec![partition];
    let dense: Vec<DenseReconstruction> = vec![];

    let mut names_tt_9 = names;
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_9,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(changed, "reversion polytomy must be resolved even without a collapse");

    let p = &sparse[0];
    let total_subs: usize = graph
      .get_edges()
      .filter_map(|e| p.partition.obs_edges.get(&e.key()))
      .map(|e| e.fitch_subs().len())
      .sum();
    assert_eq!(total_subs, 2, "reaches the parsimony optimum");

    let reversion_remains = graph
      .get_edges()
      .filter_map(|e| p.partition.obs_edges.get(&e.key()))
      .any(|e| e.fitch_subs().contains(&sub(b'T', 0, b'A')));
    assert!(!reversion_remains, "reversion must be removed");

    Ok(())
  }

  #[test]
  fn test_optimize_cascading_collapse_parent_child_both_zero() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(((A:0.1,B:0.1)I2:0.0)I1:0.0,C:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let ri1_key = find_edge_key(&graph, &names, "root", "I1").unwrap();
    let i1i2_key = find_edge_key(&graph, &names, "I1", "I2").unwrap();

    let sparse: Vec<SparseReconstruction> = vec![];
    let dense: Vec<DenseReconstruction> = vec![];

    let mut names_tt_8 = names.clone();
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[ri1_key, i1i2_key],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_8,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(changed);

    assert!(find_node_key_by_name(&graph, &names, "I1").is_none());
    assert!(find_node_key_by_name(&graph, &names, "I2").is_none());

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let root_node = graph.get_node(root_key).unwrap();
    assert_eq!(root_node.degree_out(), 3);

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
    let nuc = Alphabet::new(AlphabetName::Nuc)?;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.01,B:0.01)AB:0.01,(C:0.01,D:0.01)CD:0.01)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let mut graph: Graph = graph;

    let dense_partition = PartitionMarginalDense::new(0, nuc, get_common_length(&aln)?);
    let dense_node_states = dense_partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, aln))?;
    let dense_partitions = vec![DenseReconstruction::seeded(dense_partition, jc69(JC69Params::default())?, dense_node_states)];

    let (mut dense_partitions, _) = marginal_update_dense(&graph, &branch_lengths_or_zero(&branch_lengths), dense_partitions)?;

    let mut sparse_partitions: Vec<SparseReconstruction> = vec![];

    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);
    let sub_counts = gather_edge_sub_counts(&graph, &dense_partitions, &sparse_partitions)?;
    let effective_lengths = gather_edge_effective_lengths(&graph, &dense_partitions, &sparse_partitions)?;
    initial_guess_mixed(&graph, total_length, &indel_counts, &sub_counts, &effective_lengths, true, false, &mut branch_lengths)?;

    let initial_node_count = graph.get_nodes().count();

    let mut lh_prev = f64::MIN;
    for i in 0..10 {
      let dense_lh;
      (dense_partitions, dense_lh) = marginal_update_dense(&graph, &branch_lengths_or_zero(&branch_lengths), dense_partitions)?;
      let dense_lh = dense_lh.value();
      if (dense_lh - lh_prev).abs() < 1e-2 {
        break;
      }

      let old_branch_lengths = branch_lengths.clone();
      let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
      let contributions = gather_edge_contributions(&graph, &dense_partitions, &sparse_partitions)?;
      let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
        let zero_optimal_edges = find_zero_optimal_internal_edges(&graph, &sparse_partitions, &branch_lengths);
      apply_damping(&mut branch_lengths, &old_branch_lengths, 0.75, i);
      let mut names_tt_7 = names.clone();
      let cleanup = prune_and_merge_in_loop(&mut graph, sparse_partitions, dense_partitions, &zero_optimal_edges, TopologyOps::default(), &mut branch_lengths, &mut names_tt_7)?;
      sparse_partitions = cleanup.sparse_partitions;
      dense_partitions = cleanup.dense_partitions;

      lh_prev = dense_lh;
    }

    let final_node_count = graph.get_nodes().count();
    assert!(
      final_node_count < initial_node_count,
      "Expected topology simplification: {initial_node_count} nodes -> {final_node_count} nodes"
    );

    for edge in graph.get_edges() {
      let bl = branch_lengths.get(&edge.key()).copied().flatten().unwrap_or(0.0);
      assert!(bl >= 0.0, "Negative branch length after optimization: {bl}");
    }

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_names_new_nodes() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)I:0.0,C:0.1,D:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
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

    let sparse = vec![partition];
    let dense: Vec<DenseReconstruction> = vec![];

    branch_lengths.insert(ri_key, Some(0.0));

    let mut names_tt_6 = names;
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[ri_key],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_6,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(changed);

    let mut names: Vec<String> = graph
      .get_nodes()
      .filter_map(|n| names_tt_6.get(&n.key()).cloned().flatten())
      .collect();
    names.sort();

    assert_eq!(names, vec!["A", "B", "C", "D", "NODE_0000000", "root"]);

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_merge_disabled_keeps_polytomy() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)I:0.0,C:0.1,D:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
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

    let sparse = vec![partition];
    let dense: Vec<DenseReconstruction> = vec![];

    branch_lengths.insert(ri_key, Some(0.0));

    let ops = TopologyOps {
      merge_siblings: false,
      ..TopologyOps::default()
    };
    let mut names_tt_5 = names.clone();
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[ri_key],
      ops,
      &mut branch_lengths,
      &mut names_tt_5,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(changed, "collapse still fires even with merge disabled");

    assert!(find_node_key_by_name(&graph, &names, "I").is_none());

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let root_node = graph.get_node(root_key).unwrap();
    assert_eq!(root_node.degree_out(), 4);

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_flip_disabled_keeps_reversion() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(((C1:0.1,C2:0.1,C3:0.1)V:0.2)U:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
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

    let sparse = vec![partition];
    let dense: Vec<DenseReconstruction> = vec![];

    let ops = TopologyOps {
      flip_parent_child: false,
      ..TopologyOps::default()
    };
    let mut names_tt_4 = names;
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[],
      ops,
      &mut branch_lengths,
      &mut names_tt_4,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(changed, "merge still groups the reverting siblings");

    let p = &sparse[0];
    let reversion_remains = graph
      .get_edges()
      .filter_map(|e| p.partition.obs_edges.get(&e.key()))
      .any(|e| e.fitch_subs().contains(&sub(b'T', 0, b'A')));
    assert!(
      reversion_remains,
      "reversion is kept when flip-parent-child is disabled"
    );

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_all_ops_disabled_is_noop() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(((C1:0.1,C2:0.1,C3:0.1)V:0.2)U:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
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

    let sparse = vec![partition];
    let dense: Vec<DenseReconstruction> = vec![];

    let node_count_before = graph.get_nodes().count();
    let ops = TopologyOps {
      collapse_short_branches: false,
      merge_siblings: false,
      flip_parent_child: false,
    };
    let mut names_tt_3 = names;
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[],
      ops,
      &mut branch_lengths,
      &mut names_tt_3,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(!changed, "no topology step runs when all are disabled");
    assert_eq!(graph.get_nodes().count(), node_count_before);

    let p = &sparse[0];
    let total_subs: usize = graph
      .get_edges()
      .filter_map(|e| p.partition.obs_edges.get(&e.key()))
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
    let nuc = Alphabet::new(AlphabetName::Nuc)?;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.01,B:0.01)AB:0.01,(C:0.01,D:0.01)CD:0.01)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let mut graph: Graph = graph;

    let fitch = create_fitch_partition(&graph, 0, nuc, &node_seq_inputs(&graph, &names, aln))?;
    let (sp_partition, sp_node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partitions = vec![SparseReconstruction::seeded(sp_partition, jc69(JC69Params::default())?, sp_node_states)];
    let (sparse_partitions, _) = marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;

    let dense_partitions: Vec<DenseReconstruction> = vec![];
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);
    let sub_counts = gather_edge_sub_counts(&graph, &dense_partitions, &sparse_partitions)?;
    let effective_lengths = gather_edge_effective_lengths(&graph, &dense_partitions, &sparse_partitions)?;
    initial_guess_mixed(&graph, total_length, &indel_counts, &sub_counts, &effective_lengths, true, false, &mut branch_lengths)?;

    let initial_node_count = graph.get_nodes().count();

    let ops = TopologyOps {
      collapse_short_branches: false,
      ..TopologyOps::default()
    };
    let names_tt_2 = names;
    run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
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
      graph.get_nodes().count(),
      initial_node_count,
      "no collapse when collapse_short_branches is disabled"
    );

    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_topology_change_map_matches_edge_set() -> Result<(), Report> {
    let nuc = Alphabet::new(AlphabetName::Nuc)?;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.01,B:0.01)AB:0.01,(C:0.01,D:0.01)CD:0.01)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let mut graph: Graph = graph;

    let fitch = create_fitch_partition(&graph, 0, nuc, &node_seq_inputs(&graph, &names, aln))?;
    let (sp_partition, sp_node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partitions = vec![SparseReconstruction::seeded(
      sp_partition,
      jc69(JC69Params::default())?,
      sp_node_states,
    )];
    let (sparse_partitions, _) =
      marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;

    let dense_partitions: Vec<DenseReconstruction> = vec![];
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);
    let sub_counts = gather_edge_sub_counts(&graph, &dense_partitions, &sparse_partitions)?;
    let effective_lengths = gather_edge_effective_lengths(&graph, &dense_partitions, &sparse_partitions)?;
    initial_guess_mixed(
      &graph,
      total_length,
      &indel_counts,
      &sub_counts,
      &effective_lengths,
      true,
      false,
      &mut branch_lengths,
    )?;

    let initial_node_count = graph.get_nodes().count();

    let names_tt_1 = names;
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      10,
      1e-2,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_1,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    assert!(
      graph.get_nodes().count() < initial_node_count,
      "expected a collapse (topology change) with collapse enabled"
    );

    let map_keys: Vec<_> = result.branch_lengths.keys().copied().collect();
    let mut edge_keys: Vec<_> = graph.get_edges().map(|edge| edge.key()).collect();
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

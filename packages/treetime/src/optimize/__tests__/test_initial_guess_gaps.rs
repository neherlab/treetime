#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::gather::{
    gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts, total_sequence_length,
  };
  use crate::optimize::run_loop::{marginal_update_dense, marginal_update_sparse};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_primitives::AlignmentRecord;

  use pretty_assertions::assert_eq;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn gap_free_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        ACGTACGTACGTACGT
        >B
        ACGTACGTACGTACGA
        >C
        ACGTACGTACGTACGG
        >D
        ACGTACGTACGTACGC
      "#},
      &alphabet,
    )
  }

  /// All four sequences share gaps at positions 4-7. Remaining 12 positions
  /// are identical to `gap_free_alignment`.
  fn gappy_alignment_shared() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        ACGT----ACGTACGT
        >B
        ACGT----ACGTACGA
        >C
        ACGT----ACGTACGG
        >D
        ACGT----ACGTACGC
      "#},
      &alphabet,
    )
  }

  /// B alone has gaps at positions 4-7. Other sequences have ACGT there.
  fn gappy_alignment_one_leaf() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        ACGTACGTACGTACGT
        >B
        ACGT----ACGTACGA
        >C
        ACGTACGTACGTACGG
        >D
        ACGTACGTACGTACGC
      "#},
      &alphabet,
    )
  }

  fn setup_sparse(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<Vec<SparseReconstruction>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let fitch = create_fitch_partition(graph, 0, alphabet, &node_seq_inputs(graph, names, aln.to_vec()))?;
    let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
    let partitions = vec![SparseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
      node_states,
    )];
    let (partitions, _) = marginal_update_sparse(graph, &branch_lengths_or_zero(branch_lengths), partitions)?;

    Ok(partitions)
  }

  fn setup_dense(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<Vec<DenseReconstruction>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partition = PartitionMarginalDense::new(0, alphabet, get_common_length(aln)?);
    let node_states = partition.attach_sequences(graph, &node_seq_inputs(graph, names, aln.to_vec()))?;
    let partitions = vec![DenseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
      node_states,
    )];

    let (partitions, _) = marginal_update_dense(graph, &branch_lengths_or_zero(branch_lengths), partitions)?;

    Ok(partitions)
  }

  fn get_branch_lengths(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Vec<f64> {
    graph
      .get_edges()
      .map(|edge| branch_lengths[&edge.key()].unwrap_or(0.0))
      .collect()
  }

  #[test]
  fn test_sparse_effective_length_no_gaps() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = gap_free_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let partitions = setup_sparse(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let effective = partitions[0].edge_effective_length(&graph, edge_key)?;
      assert_eq!(16, effective);
    }

    Ok(())
  }

  #[test]
  fn test_dense_effective_length_no_gaps() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = gap_free_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let partitions = setup_dense(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let effective = partitions[0].edge_effective_length(&graph, edge_key)?;
      assert_eq!(16, effective);
    }

    Ok(())
  }

  #[test]
  fn test_sparse_effective_length_shared_gaps() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = gappy_alignment_shared()?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let partitions = setup_sparse(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let effective = partitions[0].edge_effective_length(&graph, edge_key)?;
      // All nodes share gaps at positions 4-7, so effective = 16 - 4 = 12
      assert_eq!(12, effective);
    }

    Ok(())
  }

  #[test]
  fn test_dense_effective_length_shared_gaps() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = gappy_alignment_shared()?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let partitions = setup_dense(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let effective = partitions[0].edge_effective_length(&graph, edge_key)?;
      // All nodes share gaps at positions 4-7, so effective = 16 - 4 = 12
      assert_eq!(12, effective);
    }

    Ok(())
  }

  #[test]
  fn test_sparse_effective_length_one_leaf_gapped() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = gappy_alignment_one_leaf()?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let partitions = setup_sparse(&graph, &names, &aln, &branch_lengths)?;

    let mut found_reduced = false;
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let effective = partitions[0].edge_effective_length(&graph, edge_key)?;
      // At least one edge (B→AB) should have reduced effective length
      if effective < 16 {
        found_reduced = true;
      }
    }
    assert!(found_reduced);

    Ok(())
  }

  #[test]
  fn test_dense_edge_subs_excludes_gap_positions() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = gappy_alignment_one_leaf()?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let partitions = setup_dense(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let subs = partitions[0].edge_subs(&graph, edge_key)?;
      // No substitution should involve a gap position (4-7)
      for sub in &subs {
        assert!(sub.pos() < 4 || sub.pos() >= 8);
      }
    }

    Ok(())
  }

  /// With shared gaps at positions 4-7 (which have identical nucleotides in
  /// the gap-free version), the initial guess should produce the same
  /// substitution rate per informative site. Branch lengths should be
  /// proportionally adjusted: subs/12 for gappy vs subs/16 for gap-free.
  #[test]
  fn test_initial_guess_sparse_gap_adjusted_rate() -> Result<(), Report> {
    let aln_clean: Vec<AlignmentRecord> = gap_free_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_clean_names = nwk_parsed.names();
    let graph_clean = nwk_parsed.graph;
    let mut branch_lengths_clean = nwk_parsed.branch_lengths;
    let partitions_clean = setup_sparse(&graph_clean, &graph_clean_names, &aln_clean, &branch_lengths_clean)?;
    {
      let total_length = total_sequence_length(&[], &partitions_clean);
      let indel_counts = gather_edge_indel_counts(&graph_clean, &[], &partitions_clean);
      let sub_counts = gather_edge_sub_counts(&graph_clean, &[], &partitions_clean)?;
      let effective_lengths = gather_edge_effective_lengths(&graph_clean, &[], &partitions_clean)?;
      initial_guess_mixed(
        &graph_clean,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        true,
        false,
        &mut branch_lengths_clean,
      )?;
    }
    let bl_clean = get_branch_lengths(&graph_clean, &branch_lengths_clean);

    let aln_gappy: Vec<AlignmentRecord> = gappy_alignment_shared()?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_gappy_names = nwk_parsed.names();
    let graph_gappy = nwk_parsed.graph;
    let mut branch_lengths_gappy = nwk_parsed.branch_lengths;
    let partitions_gappy = setup_sparse(&graph_gappy, &graph_gappy_names, &aln_gappy, &branch_lengths_gappy)?;
    {
      let total_length = total_sequence_length(&[], &partitions_gappy);
      let indel_counts = gather_edge_indel_counts(&graph_gappy, &[], &partitions_gappy);
      let sub_counts = gather_edge_sub_counts(&graph_gappy, &[], &partitions_gappy)?;
      let effective_lengths = gather_edge_effective_lengths(&graph_gappy, &[], &partitions_gappy)?;
      initial_guess_mixed(
        &graph_gappy,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        true,
        false,
        &mut branch_lengths_gappy,
      )?;
    }
    let bl_gappy = get_branch_lengths(&graph_gappy, &branch_lengths_gappy);

    // With 4 shared gap positions out of 16, the effective length is 12.
    // Substitutions at non-gap positions are the same, so the per-site rate
    // is higher by factor 16/12 = 4/3.
    let ratio = 16.0 / 12.0;
    for (clean, gappy) in bl_clean.iter().zip(bl_gappy.iter()) {
      if *clean > 0.0 {
        assert_abs_diff_eq!(gappy / clean, ratio, epsilon = 1e-10);
      } else {
        assert_abs_diff_eq!(*gappy, 0.0, epsilon = 1e-10);
      }
    }

    Ok(())
  }
}

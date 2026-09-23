#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

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
  use crate::optimize::params::ExistingBranchLengths;
  use crate::optimize::run_loop::{marginal_update_dense, marginal_update_sparse};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::optimize::contribution::OptimizationContribution;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use treetime_graph::graph::Graph;
  use treetime_primitives::AlignmentRecord;

  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  #[test]
  fn test_initial_guess_formula_sparse() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = divergent_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let partitions = setup_sparse(&graph, &names, &aln, &branch_lengths)?;

    {
      let total_length = total_sequence_length(&[], &partitions);
      let indel_counts = gather_edge_indel_counts(&graph, &[], &partitions);
      let sub_counts = gather_edge_sub_counts(&graph, &[], &partitions)?;
      let effective_lengths = gather_edge_effective_lengths(&graph, &[], &partitions)?;
      initial_guess_mixed(
        &graph,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        ExistingBranchLengths::Overwrite,
        false,
        &mut branch_lengths,
      )?;
    }

    let p = &partitions[0];
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let sub_count = p.edge_subs(edge_key)?.len();
      let effective_length = p.edge_effective_length(&graph, edge_key)?;
      let actual_bl = branch_lengths[&edge_key].unwrap_or(0.0);

      let expected_bl = if effective_length > 0 {
        sub_count as f64 / effective_length as f64
      } else {
        0.0
      };

      assert_abs_diff_eq!(expected_bl, actual_bl, epsilon = 1e-15);
    }
    Ok(())
  }

  #[test]
  fn test_initial_guess_formula_dense() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = divergent_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let partitions = setup_dense(&graph, &names, &aln, &branch_lengths)?;

    {
      let total_length = total_sequence_length(&partitions, &[]);
      let indel_counts = gather_edge_indel_counts(&graph, &partitions, &[]);
      let sub_counts = gather_edge_sub_counts(&graph, &partitions, &[])?;
      let effective_lengths = gather_edge_effective_lengths(&graph, &partitions, &[])?;
      initial_guess_mixed(
        &graph,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        ExistingBranchLengths::Overwrite,
        false,
        &mut branch_lengths,
      )?;
    }

    let p = &partitions[0];
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let sub_count = p.edge_subs(&graph, edge_key)?.len();
      let effective_length = p.edge_effective_length(&graph, edge_key)?;
      let actual_bl = branch_lengths[&edge_key].unwrap_or(0.0);

      let expected_bl = if effective_length > 0 {
        sub_count as f64 / effective_length as f64
      } else {
        0.0
      };

      assert_abs_diff_eq!(expected_bl, actual_bl, epsilon = 1e-15);
    }
    Ok(())
  }

  #[test]
  fn test_initial_guess_dense_sparse_ambiguous_r_reference_state_consistency() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = ambiguous_r_in_g_clade_alignment()?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_dense_names = nwk_parsed.names();
    let graph_dense = nwk_parsed.graph;
    let mut branch_lengths_dense = nwk_parsed.branch_lengths;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_sparse_names = nwk_parsed.names();
    let graph_sparse = nwk_parsed.graph;
    let mut branch_lengths_sparse = nwk_parsed.branch_lengths;

    let partitions_dense = setup_dense(&graph_dense, &graph_dense_names, &aln, &branch_lengths_dense)?;
    let partitions_sparse = setup_sparse(&graph_sparse, &graph_sparse_names, &aln, &branch_lengths_sparse)?;

    {
      let total_length = total_sequence_length(&partitions_dense, &[]);
      let indel_counts = gather_edge_indel_counts(&graph_dense, &partitions_dense, &[]);
      let sub_counts = gather_edge_sub_counts(&graph_dense, &partitions_dense, &[])?;
      let effective_lengths = gather_edge_effective_lengths(&graph_dense, &partitions_dense, &[])?;
      initial_guess_mixed(
        &graph_dense,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        ExistingBranchLengths::Overwrite,
        false,
        &mut branch_lengths_dense,
      )?;
    }
    {
      let total_length = total_sequence_length(&[], &partitions_sparse);
      let indel_counts = gather_edge_indel_counts(&graph_sparse, &[], &partitions_sparse);
      let sub_counts = gather_edge_sub_counts(&graph_sparse, &[], &partitions_sparse)?;
      let effective_lengths = gather_edge_effective_lengths(&graph_sparse, &[], &partitions_sparse)?;
      initial_guess_mixed(
        &graph_sparse,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        ExistingBranchLengths::Overwrite,
        false,
        &mut branch_lengths_sparse,
      )?;
    }

    let dense_branch_lengths = branch_lengths_by_child_name(&graph_dense, &graph_dense_names, &branch_lengths_dense)?;
    let sparse_branch_lengths =
      branch_lengths_by_child_name(&graph_sparse, &graph_sparse_names, &branch_lengths_sparse)?;

    assert_eq!(dense_branch_lengths, sparse_branch_lengths);
    Ok(())
  }

  #[test]
  fn test_optimize_contribution_dense_sparse_ambiguous_r_value_and_gradient_consistency() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = ambiguous_r_in_g_clade_alignment()?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_dense_names = nwk_parsed.names();
    let graph_dense = nwk_parsed.graph;
    let branch_lengths_dense = nwk_parsed.branch_lengths;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_sparse_names = nwk_parsed.names();
    let graph_sparse = nwk_parsed.graph;
    let branch_lengths_sparse = nwk_parsed.branch_lengths;

    let partitions_dense = setup_dense(&graph_dense, &graph_dense_names, &aln, &branch_lengths_dense)?;
    let partitions_sparse = setup_sparse(&graph_sparse, &graph_sparse_names, &aln, &branch_lengths_sparse)?;

    let dense_metrics = optimization_metrics_by_child_name(
      &graph_dense,
      &graph_dense_names,
      |key| Ok(partitions_dense[0].create_edge_contribution(key)),
      0.1,
    )?;
    let sparse_metrics = optimization_metrics_by_child_name(
      &graph_sparse,
      &graph_sparse_names,
      |key| partitions_sparse[0].create_edge_contribution(key),
      0.1,
    )?;

    assert_eq!(
      dense_metrics.keys().cloned().collect::<Vec<_>>(),
      sparse_metrics.keys().cloned().collect::<Vec<_>>()
    );
    for (edge_name, (dense_log_lh, dense_derivative, _dense_second_derivative)) in dense_metrics {
      let (sparse_log_lh, sparse_derivative, _sparse_second_derivative) = sparse_metrics[&edge_name];
      assert_abs_diff_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-12);
      assert_abs_diff_eq!(dense_derivative, sparse_derivative, epsilon = 1e-12);
    }
    Ok(())
  }

  fn divergent_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        AAAACCCCGGGGTTTT
        >B
        CCCCGGGGTTTTAAAA
        >C
        GGGGTTTTAAAACCCC
        >D
        TTTTAAAACCCCGGGG
      "#},
      &alphabet,
    )
  }

  fn ambiguous_r_in_g_clade_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        RCGTACGT
        >B
        GCGTACGT
        >C
        GCGTACGT
        >D
        GCGTACGT
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

  fn branch_lengths_by_child_name(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<BTreeMap<String, f64>, Report> {
    graph
      .get_edges()
      .map(|edge_ref| {
        let child_key = edge_ref.target();
        let child_name = names[&child_key].clone().unwrap();
        let branch_length = branch_lengths[&edge_ref.key()].unwrap_or(0.0);
        Ok((child_name, branch_length))
      })
      .collect()
  }

  fn optimization_metrics_by_child_name(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    contribution: impl Fn(GraphEdgeKey) -> Result<OptimizationContribution, Report>,
    branch_length: f64,
  ) -> Result<BTreeMap<String, (f64, f64, f64)>, Report> {
    graph
      .get_edges()
      .map(|edge_ref| {
        let child_key = edge_ref.target();
        let child_name = names[&child_key].clone().unwrap();
        let metrics = contribution(edge_ref.key())?
          .evaluate(branch_length)
          .expect("valid branch length");
        Ok((
          child_name,
          (metrics.log_lh.value(), metrics.derivative, metrics.second_derivative),
        ))
      })
      .collect()
  }
}

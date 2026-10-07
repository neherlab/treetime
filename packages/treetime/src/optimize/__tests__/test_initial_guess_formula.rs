#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::gather::{gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts};
  use crate::optimize::likelihood::evaluate_mixed;
  use crate::optimize::params::ExistingBranchLengths;
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction, SparseReconstruction};
  use crate::partition::optimize::contribution::OptimizationContribution;
  use crate::test_utils::leaf_seq_inputs;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::{FastaRecord, fasta_read};
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  #[test]
  fn test_initial_guess_formula_sparse() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = divergent_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_sparse(&graph, &names, &aln, &branch_lengths)?;

    {
      let total_length = reconstruction.sequence_length();
      let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
      let sub_counts = gather_edge_sub_counts(&graph, &reconstruction)?;
      let effective_lengths = gather_edge_effective_lengths(&graph, &reconstruction)?;
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

    let p = &reconstruction;
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
  fn test_initial_guess_formula_dense() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = divergent_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_dense(&graph, &names, &aln, &branch_lengths)?;

    {
      let total_length = reconstruction.sequence_length();
      let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
      let sub_counts = gather_edge_sub_counts(&graph, &reconstruction)?;
      let effective_lengths = gather_edge_effective_lengths(&graph, &reconstruction)?;
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

    let p = &reconstruction;
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

    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let graph_dense_names = nwk_parsed.names();
    let graph_dense = nwk_parsed.graph;
    let mut branch_lengths_dense = nwk_parsed.branch_lengths;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let graph_sparse_names = nwk_parsed.names();
    let graph_sparse = nwk_parsed.graph;
    let mut branch_lengths_sparse = nwk_parsed.branch_lengths;

    let reconstruction_dense = setup_dense(&graph_dense, &graph_dense_names, &aln, &branch_lengths_dense)?;
    let reconstruction_sparse = setup_sparse(&graph_sparse, &graph_sparse_names, &aln, &branch_lengths_sparse)?;

    {
      let total_length = reconstruction_dense.sequence_length();
      let indel_counts = gather_edge_indel_counts(&graph_dense, &reconstruction_dense);
      let sub_counts = gather_edge_sub_counts(&graph_dense, &reconstruction_dense)?;
      let effective_lengths = gather_edge_effective_lengths(&graph_dense, &reconstruction_dense)?;
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
      let total_length = reconstruction_sparse.sequence_length();
      let indel_counts = gather_edge_indel_counts(&graph_sparse, &reconstruction_sparse);
      let sub_counts = gather_edge_sub_counts(&graph_sparse, &reconstruction_sparse)?;
      let effective_lengths = gather_edge_effective_lengths(&graph_sparse, &reconstruction_sparse)?;
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

    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let graph_dense_names = nwk_parsed.names();
    let graph_dense = nwk_parsed.graph;
    let branch_lengths_dense = nwk_parsed.branch_lengths;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let graph_sparse_names = nwk_parsed.names();
    let graph_sparse = nwk_parsed.graph;
    let branch_lengths_sparse = nwk_parsed.branch_lengths;

    let reconstruction_dense = setup_dense(&graph_dense, &graph_dense_names, &aln, &branch_lengths_dense)?;
    let reconstruction_sparse = setup_sparse(&graph_sparse, &graph_sparse_names, &aln, &branch_lengths_sparse)?;

    let dense_metrics = optimization_metrics_by_child_name(
      &graph_dense,
      &graph_dense_names,
      |key| reconstruction_dense.create_edge_contribution(key),
      0.1,
    )?;
    let sparse_metrics = optimization_metrics_by_child_name(
      &graph_sparse,
      &graph_sparse_names,
      |key| reconstruction_sparse.create_edge_contribution(key),
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
    fasta_read(
      indoc! {r#"
        >A
        AAAACCCCGGGGTTTT
        >B
        CCCCGGGGTTTTAAAA
        >C
        GGGGTTTTAAAACCCC
        >D
        TTTTAAAACCCCGGGG
      "#}
      .as_bytes(),
      &alphabet,
    )
  }

  fn ambiguous_r_in_g_clade_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    fasta_read(
      indoc! {r#"
        >A
        RCGTACGT
        >B
        GCGTACGT
        >C
        GCGTACGT
        >D
        GCGTACGT
      "#}
      .as_bytes(),
      &alphabet,
    )
  }

  fn setup_sparse(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<MarginalReconstruction, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let fitch = create_fitch_partition(graph, alphabet, leaf_seq_inputs(graph, names, aln.to_vec()))?;
    let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
    let reconstruction = MarginalReconstruction::Sparse(SparseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
      node_states,
    ));
    let (reconstruction, _) = reconstruction.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;

    Ok(reconstruction)
  }

  fn setup_dense(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<MarginalReconstruction, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partition = PartitionMarginalDense::new(alphabet, graph, &leaf_seq_inputs(graph, names, aln.to_vec()))?;
    let reconstruction =
      MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, jc69(JC69Params::default())?));

    let (reconstruction, _) = reconstruction.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;

    Ok(reconstruction)
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
        let metrics = evaluate_mixed(&[contribution(edge_ref.key())?], branch_length).expect("valid branch length");
        Ok((
          child_name,
          (metrics.log_lh.value(), metrics.derivative, metrics.second_derivative),
        ))
      })
      .collect()
  }
}

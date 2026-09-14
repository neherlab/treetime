#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::gather::{
    gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts, total_sequence_length,
  };
  use crate::optimize::run_loop::{marginal_update_dense, marginal_update_sparse};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::traits::PartitionBranchOps;
  use crate::partition::traits::PartitionOptimizeOps;
  use crate::seq::alignment::get_common_length;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use treetime_graph::graph::Graph;

  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::{NwkParse, nwk_read_str};

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  /// Verify initial_guess_mixed sets branch_length = edge_subs().len() / edge_effective_length()
  /// for sparse partitions. This is the defining formula: the initial branch length
  /// is the fraction of non-gap positions with substitutions.
  #[test]
  fn test_initial_guess_formula_sparse() -> Result<(), Report> {
    let aln = divergent_alignment()?;
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
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
        true,
        false,
        &mut branch_lengths,
      )?;
    }

    let p = partitions[0].readout();
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
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

  /// Same formula verification for dense partitions. This test would have caught
  /// the soft Hamming override that was previously shadowing the sub count.
  #[test]
  fn test_initial_guess_formula_dense() -> Result<(), Report> {
    let aln = divergent_alignment()?;
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
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
        true,
        false,
        &mut branch_lengths,
      )?;
    }

    let p = partitions[0].readout();
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
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
    let aln = ambiguous_r_in_g_clade_alignment()?;

    let NwkParse {
      graph: graph_dense,
      names: graph_dense_names,
      branch_lengths: mut branch_lengths_dense,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let NwkParse {
      graph: graph_sparse,
      names: graph_sparse_names,
      branch_lengths: mut branch_lengths_sparse,
      ..
    } = nwk_read_str(TREE_NEWICK)?;

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
        true,
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
        true,
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
    let aln = ambiguous_r_in_g_clade_alignment()?;

    let NwkParse {
      graph: graph_dense,
      names: graph_dense_names,
      branch_lengths: branch_lengths_dense,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let NwkParse {
      graph: graph_sparse,
      names: graph_sparse_names,
      branch_lengths: branch_lengths_sparse,
      ..
    } = nwk_read_str(TREE_NEWICK)?;

    let partitions_dense = setup_dense(&graph_dense, &graph_dense_names, &aln, &branch_lengths_dense)?;
    let partitions_sparse = setup_sparse(&graph_sparse, &graph_sparse_names, &aln, &branch_lengths_sparse)?;

    let dense_metrics =
      optimization_metrics_by_child_name(&graph_dense, &graph_dense_names, &partitions_dense[0].readout(), 0.1)?;
    let sparse_metrics =
      optimization_metrics_by_child_name(&graph_sparse, &graph_sparse_names, &partitions_sparse[0].readout(), 0.1)?;

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
    aln: &[FastaRecord],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<Vec<SparseReconstruction>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let fitch = create_fitch_partition(graph, 0, alphabet, aln, names)?;
    let (partition, node_states) = fitch.into_marginal_sparse(jc69(JC69Params::default())?, graph)?;
    let partitions = vec![SparseReconstruction::seeded(partition, node_states)];
    let (partitions, _) = marginal_update_sparse(graph, &profile_branch_lengths(branch_lengths), partitions)?;

    Ok(partitions)
  }

  fn setup_dense(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[FastaRecord],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<Vec<DenseReconstruction>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partition = PartitionMarginalDense::new(0, jc69(JC69Params::default())?, alphabet, get_common_length(aln)?);
    let node_states = partition.attach_sequences(graph, aln, names)?;
    let partitions = vec![DenseReconstruction::seeded(partition, node_states)];

    let (partitions, _) = marginal_update_dense(graph, &profile_branch_lengths(branch_lengths), partitions)?;

    Ok(partitions)
  }

  fn branch_lengths_by_child_name(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<BTreeMap<String, f64>, Report> {
    graph
      .get_edges()
      .iter()
      .map(|edge_ref| {
        let edge_ref = edge_ref.read_arc();
        let child_key = edge_ref.target();
        let child_name = names[&child_key].clone().unwrap();
        let branch_length = branch_lengths[&edge_ref.key()].unwrap_or(0.0);
        Ok((child_name, branch_length))
      })
      .collect()
  }

  fn optimization_metrics_by_child_name<P: PartitionOptimizeOps>(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    partition: &P,
    branch_length: f64,
  ) -> Result<BTreeMap<String, (f64, f64, f64)>, Report> {
    graph
      .get_edges()
      .iter()
      .map(|edge_ref| {
        let edge_ref = edge_ref.read_arc();
        let child_key = edge_ref.target();
        let child_name = names[&child_key].clone().unwrap();
        let metrics = partition
          .create_edge_contribution(edge_ref.key())?
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

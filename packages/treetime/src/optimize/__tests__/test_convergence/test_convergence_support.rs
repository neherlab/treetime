#[cfg(test)]
pub mod tests {
  use treetime_io::nwk::nwk_fasta_node_inputs;
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
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use indoc::indoc;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;

  use std::sync::LazyLock;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};

  pub static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  // Small tree with 4 leaves
  pub const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  pub fn simple_alignment() -> Result<Vec<FastaRecord>, Report> {
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
      &*NUC_ALPHABET,
    )
  }

  pub fn setup_partitions(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[FastaRecord],
    branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<(Vec<DenseReconstruction>, Vec<SparseReconstruction>), Report> {
    let alphabet_dense = Alphabet::new(AlphabetName::Nuc)?;
    let alphabet_sparse = Alphabet::new(AlphabetName::Nuc)?;

    let dense_partition =
      PartitionMarginalDense::new(0, jc69(JC69Params::default())?, alphabet_dense, get_common_length(aln)?);
    let dense_node_states = dense_partition.attach_sequences(graph, &nwk_fasta_node_inputs(graph, names, aln.to_vec()))?;
    let dense_partitions = vec![DenseReconstruction::seeded(dense_partition, dense_node_states)];

    let fitch = create_fitch_partition(graph, 1, alphabet_sparse, &nwk_fasta_node_inputs(graph, names, aln.to_vec()))?;
    let (sparse_partition, sparse_node_states) = fitch.into_marginal_sparse(jc69(JC69Params::default())?, graph)?;
    let sparse_partitions = vec![SparseReconstruction::seeded(sparse_partition, sparse_node_states)];

    let (dense_partitions, _) =
      marginal_update_dense(graph, &profile_branch_lengths(branch_lengths), dense_partitions)?;
    let (sparse_partitions, _) =
      marginal_update_sparse(graph, &profile_branch_lengths(branch_lengths), sparse_partitions)?;

    {
      let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
      let indel_counts = gather_edge_indel_counts(graph, &dense_partitions, &sparse_partitions);
      let sub_counts = gather_edge_sub_counts(graph, &dense_partitions, &sparse_partitions)?;
      let effective_lengths = gather_edge_effective_lengths(graph, &dense_partitions, &sparse_partitions)?;
      initial_guess_mixed(
        graph,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        true,
        false,
        branch_lengths,
      )?;
    }

    Ok((dense_partitions, sparse_partitions))
  }

  pub fn compute_total_lh(
    graph: &Graph,
    dense_partitions: Vec<DenseReconstruction>,
    sparse_partitions: Vec<SparseReconstruction>,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<(Vec<DenseReconstruction>, Vec<SparseReconstruction>, f64), Report> {
    let (dense_partitions, dense_lh) =
      marginal_update_dense(graph, &profile_branch_lengths(branch_lengths), dense_partitions)?;
    let (sparse_partitions, sparse_lh) =
      marginal_update_sparse(graph, &profile_branch_lengths(branch_lengths), sparse_partitions)?;
    Ok((
      dense_partitions,
      sparse_partitions,
      dense_lh.value() + sparse_lh.value(),
    ))
  }
}

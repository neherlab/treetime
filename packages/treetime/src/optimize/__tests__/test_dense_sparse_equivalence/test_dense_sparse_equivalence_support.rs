#[cfg(test)]
pub mod tests {
  //! Tests for dense vs sparse optimization equivalence.
  //!
  //! Both dense and sparse partitions route through the unified optimizer
  //! (`run_optimize_mixed`), so the optimization loop (zero-branch shortcut,
  //! Newton iteration, grid search fallback) is identical. Differences in
  //! converged branch lengths arise from the coefficient representations:
  //! dense uses per-position probability vectors, sparse uses per-site
  //! multiplicity-weighted contributions. The tests verify:
  //! 1. Both produce finite, valid log-LH values
  //! 2. Both converge to stable values
  //! 3. Initial log-LH (before optimization) should be identical
  //! 4. Final log-LH difference should be bounded

  use std::collections::BTreeMap;

  use treetime_graph::node::GraphNodeKey;

  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::run_loop::{marginal_update_dense, marginal_update_sparse};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use indoc::indoc;
  use treetime_graph::graph::Graph;

  use std::sync::LazyLock;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};

  pub static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  pub const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  pub fn gap_free_alignment() -> Result<Vec<FastaRecord>, Report> {
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

  pub fn setup_dense_only(
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

  pub fn setup_sparse_only(
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

  pub fn get_branch_lengths(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Vec<f64> {
    graph
      .get_edges()
      .iter()
      .map(|edge| branch_lengths[&edge.read_arc().key()].unwrap_or(0.0))
      .collect()
  }
}

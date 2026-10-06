#[cfg(test)]
pub(crate) mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::gather::{gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts};
  use crate::optimize::params::ExistingBranchLengths;
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction, SparseReconstruction};
  use crate::test_utils::leaf_seq_inputs;
  use eyre::Report;
  use indoc::indoc;
  use std::collections::BTreeMap;
  use std::sync::LazyLock;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::fasta_read;
  use treetime_primitives::AlignmentRecord;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  pub(crate) const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  pub(crate) fn simple_alignment() -> Result<Vec<AlignmentRecord>, Report> {
    Ok(
      fasta_read(
        indoc! {r#"
      >A
      ACGTACGTACGTACGT
      >B
      ACGTACGTACGTACGA
      >C
      ACGTACGTACGTACGG
      >D
      ACGTACGTACGTACGC
    "#}
        .as_bytes(),
        &*NUC_ALPHABET,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect(),
    )
  }

  pub(crate) fn setup_reconstruction(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<MarginalReconstruction, Report> {
    let fitch = create_fitch_partition(
      graph,
      0,
      Alphabet::new(AlphabetName::Nuc)?,
      leaf_seq_inputs(graph, names, aln.to_vec()),
    )?;
    let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
    let reconstruction = MarginalReconstruction::Sparse(SparseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
      node_states,
    ));
    let (reconstruction, _) = reconstruction.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    apply_initial_guess(graph, &reconstruction, branch_lengths)?;
    Ok(reconstruction)
  }

  pub(crate) fn setup_dense_reconstruction(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<MarginalReconstruction, Report> {
    let partition = PartitionMarginalDense::new(
      0,
      Alphabet::new(AlphabetName::Nuc)?,
      graph,
      &leaf_seq_inputs(graph, names, aln.to_vec()),
    )?;
    let reconstruction =
      MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, jc69(JC69Params::default())?));
    let (reconstruction, _) = reconstruction.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    apply_initial_guess(graph, &reconstruction, branch_lengths)?;
    Ok(reconstruction)
  }

  pub(crate) fn compute_total_lh(
    graph: &Graph,
    reconstruction: MarginalReconstruction,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<(MarginalReconstruction, f64), Report> {
    let (reconstruction, log_lh) = reconstruction.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    Ok((reconstruction, log_lh.value()))
  }

  fn apply_initial_guess(
    graph: &Graph,
    reconstruction: &MarginalReconstruction,
    branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<(), Report> {
    initial_guess_mixed(
      graph,
      reconstruction.sequence_length(),
      &gather_edge_indel_counts(graph, reconstruction),
      &gather_edge_sub_counts(graph, reconstruction)?,
      &gather_edge_effective_lengths(graph, reconstruction)?,
      ExistingBranchLengths::Overwrite,
      false,
      branch_lengths,
    )
  }
}

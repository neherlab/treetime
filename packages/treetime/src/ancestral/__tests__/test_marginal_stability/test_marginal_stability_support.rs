#[cfg(test)]
pub(super) mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::gtr::gtr::GTR;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::storage::dense::DenseSeqDistribution;
  use crate::partition::storage::sparse::SparseSeqDistribution;
  use crate::pretty_assert_ulps_eq;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use eyre::Report;
  use treetime_graph::graph::Graph;
  use treetime_utils::{pretty_assert_array_finite, pretty_assert_array_nonneg};

  use std::sync::LazyLock;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  pub(crate) fn assert_dense_profile_stable(profile: &DenseSeqDistribution, max_ulps: u32) {
    assert!(
      profile.log_lh.value().is_finite(),
      "Profile log_lh is not finite: {}",
      profile.log_lh.value()
    );
    pretty_assert_array_finite!(profile.dis);
    pretty_assert_array_nonneg!(profile.dis, epsilon = 1e-15);
    for row in profile.dis.outer_iter() {
      pretty_assert_ulps_eq!(row.sum(), 1.0, max_ulps = max_ulps);
    }
  }

  pub(crate) fn assert_sparse_profile_stable(profile: &SparseSeqDistribution, max_ulps: u32) {
    assert!(
      profile.log_lh.value().is_finite(),
      "Profile log_lh is not finite: {}",
      profile.log_lh.value()
    );

    for var_pos in profile.variable.values() {
      pretty_assert_array_finite!(var_pos.dis);
      pretty_assert_array_nonneg!(var_pos.dis, epsilon = 1e-15);
      pretty_assert_ulps_eq!(var_pos.dis.sum(), 1.0, max_ulps = max_ulps);
    }

    for fixed_dis in profile.fixed.values() {
      pretty_assert_array_finite!(fixed_dis);
      pretty_assert_array_nonneg!(fixed_dis, epsilon = 1e-15);
      pretty_assert_ulps_eq!(fixed_dis.sum(), 1.0, max_ulps = max_ulps);
    }
  }

  pub(crate) fn run_dense_marginal_with_partitions(
    newick: &str,
    aln_str: &str,
    gtr: GTR,
  ) -> Result<(f64, DenseReconstruction), Report> {
    let nwk_parsed = nwk_read_str(newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let partition = PartitionMarginalDense::new(0, alphabet, get_common_length(&aln)?);
    let node_states = partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, aln))?;
    let recon = DenseReconstruction::seeded(partition, gtr, node_states);
    let (recon, log_lh) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh = log_lh.value();
    Ok((log_lh, recon))
  }

  pub(crate) fn run_sparse_marginal_with_partitions(
    newick: &str,
    aln_str: &str,
    gtr: GTR,
  ) -> Result<(f64, SparseReconstruction), Report> {
    let nwk_parsed = nwk_read_str(newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let fitch = create_fitch_partition(&graph, 0, alphabet, &node_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let recon = SparseReconstruction::seeded(partition, gtr, node_states);
    let (recon, log_lh) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh = log_lh.value();
    Ok((log_lh, recon))
  }
}

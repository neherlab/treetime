#[cfg(test)]
pub mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::gtr::GTR;
  use crate::partition::dense::DenseSeqDistribution;
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::partition::sparse::SparseSeqDistribution;
  use crate::payload::ancestral::GraphAncestral;
  use crate::pretty_assert_ulps_eq;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use treetime_utils::{pretty_assert_array_finite, pretty_assert_array_nonneg};

  use parking_lot::RwLock;
  use std::sync::{Arc, LazyLock};
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  /// Assert that a dense marginal profile is numerically stable.
  pub fn assert_dense_profile_stable(profile: &DenseSeqDistribution, max_ulps: u32) {
    assert!(
      profile.log_lh.is_finite(),
      "Profile log_lh is not finite: {}",
      profile.log_lh
    );
    pretty_assert_array_finite!(profile.dis);
    pretty_assert_array_nonneg!(profile.dis, epsilon = 1e-15);
    for row in profile.dis.outer_iter() {
      pretty_assert_ulps_eq!(row.sum(), 1.0, max_ulps = max_ulps);
    }
  }

  /// Assert that a sparse marginal profile is numerically stable.
  pub fn assert_sparse_profile_stable(profile: &SparseSeqDistribution, max_ulps: u32) {
    assert!(
      profile.log_lh.is_finite(),
      "Profile log_lh is not finite: {}",
      profile.log_lh
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

  /// Run dense marginal reconstruction with a custom GTR model and return the log-likelihood and partition array.
  pub fn run_dense_marginal_with_partitions(
    newick: &str,
    aln_str: &str,
    gtr: GTR,
  ) -> Result<(f64, [Arc<RwLock<PartitionMarginalDense>>; 1]), Report> {
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense::new(
      0,
      gtr,
      alphabet,
      get_common_length(&aln)?,
    )))];

    let log_lh = initialize_marginal(&graph, &partitions, &aln)?;
    Ok((log_lh, partitions))
  }

  /// Run sparse marginal reconstruction with a custom GTR model and return the log-likelihood and partition array.
  pub fn run_sparse_marginal_with_partitions(
    newick: &str,
    aln_str: &str,
    gtr: GTR,
  ) -> Result<(f64, [Arc<RwLock<PartitionMarginalSparse>>; 1]), Report> {
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;
    let partitions = [Arc::new(RwLock::new(fitch.into_marginal_sparse(gtr, &graph)?))];
    let log_lh = update_marginal(&graph, &partitions)?;
    Ok((log_lh, partitions))
  }
}

//! Tests for sparse GTR inference.

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{ancestral_reconstruction_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::infer_gtr::common::{InferGtrOptions, MutationCounts, infer_gtr_impl};
  use crate::gtr::infer_gtr::sparse::get_mutation_counts_sparse;
  use crate::pretty_assert_ulps_eq;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use ndarray::array;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
  }

  #[test]
  fn test_get_mutation_counts_sparse() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    compress_sequences(&graph, std::slice::from_ref(&partition), &aln)?;
    update_marginal(&graph, std::slice::from_ref(&partition))?;
    ancestral_reconstruction_marginal(&graph, true, std::slice::from_ref(&partition), |_, _| Ok(()))?;

    let counts_actual = get_mutation_counts_sparse(&graph, &partition)?;
    // Expected values reflect MAP-derived mutations from marginal posteriors.
    // Reconstruction updates seq.composition so Ti and root_state are consistent
    // with the post-marginal MAP states used by edge_subs_from_graph().
    let counts_expected = MutationCounts {
      nij: array![[0., 0., 1., 1.], [0., 0., 0., 1.], [1., 1., 0., 0.], [0., 0., 1., 0.]],
      Ti: array![1.47, 2.115, 3.495, 3.0],
      root_state: array![1.0, 3.0, 5.0, 5.0],
    };

    pretty_assert_ulps_eq!(counts_expected.nij, counts_actual.nij, epsilon = 1e-9);
    pretty_assert_ulps_eq!(counts_expected.Ti, counts_actual.Ti, epsilon = 1e-7);
    pretty_assert_ulps_eq!(counts_expected.root_state, counts_actual.root_state, epsilon = 1e-9);

    Ok(())
  }

  #[test]
  fn test_infer_gtr_sparse() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    compress_sequences(&graph, std::slice::from_ref(&partition), &aln)?;
    update_marginal(&graph, std::slice::from_ref(&partition))?;
    ancestral_reconstruction_marginal(&graph, true, std::slice::from_ref(&partition), |_, _| Ok(()))?;

    let counts = get_mutation_counts_sparse(&graph, &partition)?;
    let actual = infer_gtr_impl(
      &counts,
      &InferGtrOptions {
        pc: 0.1,
        ..InferGtrOptions::default()
      },
    )?;

    // Expected values reflect GTR inference from post-reconstruction MAP states
    // with consistent composition counts.
    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      array![
        [0.0, 0.42794649, 3.16989133, 1.92655196],
        [0.42794649, 0.0, 1.19560283, 1.34357392],
        [3.16989133, 1.19560283, 0.0, 0.85871142],
        [1.92655196, 1.34357392, 0.85871142, 0.0]
      ],
      &actual.W,
      epsilon = 1e-7
    );

    pretty_assert_ulps_eq!(
      array![0.12784332, 0.21149000, 0.34927102, 0.31139566],
      &actual.pi,
      epsilon = 1e-7
    );
    pretty_assert_ulps_eq!(0.5892584186085698, actual.mu, epsilon = 1e-7);

    Ok(())
  }
}

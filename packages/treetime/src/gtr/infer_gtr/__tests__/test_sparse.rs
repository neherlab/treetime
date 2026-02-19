//! Tests for sparse GTR inference.

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::update_marginal;
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

    let counts_actual = get_mutation_counts_sparse(&graph, &partition)?;
    let counts_expected = MutationCounts {
      nij: array![[0., 0., 0., 0.], [2., 0., 0., 1.], [3., 2., 0., 0.], [0., 1., 1., 0.]],
      Ti: array![1.98, 2.945, 2.515, 2.64],
      root_state: array![4.0, 3.0, 3.0, 4.0],
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

    let counts = get_mutation_counts_sparse(&graph, &partition)?;
    let actual = infer_gtr_impl(
      &counts,
      &InferGtrOptions {
        pc: 0.1,
        ..InferGtrOptions::default()
      },
    )?;

    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      array![
        [0.0, 2.1751124, 2.95601658, 0.18620301],
        [2.1751124, 0.0, 1.40528091, 1.41465696],
        [2.95601658, 1.40528091, 0.0, 0.74490315],
        [0.18620301, 1.41465696, 0.74490315, 0.0]
      ],
      &actual.W,
      epsilon = 1e-7
    );

    pretty_assert_ulps_eq!(
      array![0.14878846, 0.24051536, 0.31239203, 0.29830414],
      &actual.pi,
      epsilon = 1e-7
    );
    pretty_assert_ulps_eq!(0.9471364432348814, actual.mu, epsilon = 1e-7);

    Ok(())
  }
}

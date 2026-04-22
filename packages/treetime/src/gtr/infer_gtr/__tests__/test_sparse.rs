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
  use treetime_primitives::seq;

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
      root_sequence: seq![],
    }));

    compress_sequences(&graph, std::slice::from_ref(&partition), &aln)?;
    update_marginal(&graph, std::slice::from_ref(&partition))?;

    let counts_actual = get_mutation_counts_sparse(&graph, &partition)?;
    // Expected values reflect MAP-derived mutations from marginal posteriors
    // (not stale Fitch substitutions). After marginal inference, edge_subs()
    // compares MAP states between parent and child, which can differ from Fitch
    // assignments at ambiguous positions.
    let counts_expected = MutationCounts {
      nij: array![[0., 0., 1., 1.], [0., 0., 0., 1.], [1., 1., 0., 0.], [0., 0., 1., 0.]],
      Ti: array![1.72, 2.835, 2.775, 2.75],
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
      root_sequence: seq![],
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

    // Expected values reflect GTR inference from MAP-derived marginal mutations
    // (not stale Fitch substitutions).
    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      array![
        [0.0, 0.23187287, 2.51672889, 1.35482864],
        [0.23187287, 0.0, 1.29546764, 1.26410798],
        [2.51672889, 1.29546764, 0.0, 1.23030742],
        [1.35482864, 1.26410798, 1.23030742, 0.0]
      ],
      &actual.W,
      epsilon = 1e-7
    );

    pretty_assert_ulps_eq!(
      array![0.28554666, 0.21928785, 0.24010044, 0.25506505],
      &actual.pi,
      epsilon = 1e-7
    );
    pretty_assert_ulps_eq!(0.6041583893945609, actual.mu, epsilon = 1e-7);

    Ok(())
  }
}

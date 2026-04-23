//! Tests for sparse GTR inference.

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::update_marginal;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::infer_gtr::common::{InferGtrOptions, infer_gtr_impl};
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
  use std::slice::from_ref;
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

compress_sequences(&graph, from_ref(&partition), &aln)?;
    partition.write_arc().extract_root_sequence(&graph)?;
    update_marginal(&graph, from_ref(&partition))?;

    let counts_actual = get_mutation_counts_sparse(&graph, &partition)?;
    // Expected values reflect Fitch parsimony mutations. GTR inference reads
    // fitch_subs() directly because it runs before marginal inference.
    pretty_assert_ulps_eq!(
      counts_actual.nij,
      array![[0., 0., 0., 0.], [2., 0., 0., 1.], [3., 2., 0., 0.], [0., 1., 1., 0.]],
      epsilon = 1e-9
    );
    pretty_assert_ulps_eq!(counts_actual.root_state, array![4.0, 3.0, 3.0, 4.0], epsilon = 1e-9);
    pretty_assert_ulps_eq!(
      counts_actual.Ti,
      array![
        1.9800000227987766,
        2.9450000282377005,
        2.515000017359853,
        2.640000019222498
      ],
      epsilon = 1e-9
    );

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

compress_sequences(&graph, from_ref(&partition), &aln)?;
    partition.write_arc().extract_root_sequence(&graph)?;
    update_marginal(&graph, from_ref(&partition))?;

    let counts = get_mutation_counts_sparse(&graph, &partition)?;
    let actual = infer_gtr_impl(
      &counts,
      &InferGtrOptions {
        pc: 0.1,
        ..InferGtrOptions::default()
      },
    )?;

    // Expected values reflect GTR inference from Fitch parsimony mutations.
    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      array![
        [0.0,                  2.175112392280934,  2.956016575101041,  0.18620301184748264],
        [2.175112392280934,    0.0,                1.405280909422987,  1.4146569619049196 ],
        [2.956016575101041,    1.405280909422987,  0.0,                0.744903148823837  ],
        [0.18620301184748264,  1.4146569619049196, 0.744903148823837,  0.0                ]
      ],
      &actual.W,
      epsilon = 1e-7
    );

    pretty_assert_ulps_eq!(
      array![
        0.14878846342301447,
        0.2405153630584697,
        0.31239203299951585,
        0.29830414051899984
      ],
      &actual.pi,
      epsilon = 1e-7
    );
    pretty_assert_ulps_eq!(0.9471364348810554, actual.mu, epsilon = 1e-7);

    Ok(())
  }
}

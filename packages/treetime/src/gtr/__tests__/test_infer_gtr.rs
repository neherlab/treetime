#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::update_marginal;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::avg_transition;
  use crate::gtr::infer_gtr::{InferGtrOptions, MutationCounts, distance, get_mutation_counts, infer_gtr_impl};
  use crate::pretty_assert_ulps_eq;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
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
  fn test_infer_gtr_only() -> Result<(), Report> {
    let nij = array![
      [0.0, 1.0, 2.0, 1.0],
      [1.0, 0.0, 3.0, 2.0],
      [2.0, 3.0, 0.0, 1.0],
      [2.0, 3.0, 3.0, 0.0]
    ];
    let Ti = array![12.0, 20.0, 14.0, 12.4];
    let root_state = array![3.0, 2.0, 3.0, 4.0];

    let actual = infer_gtr_impl(
      &MutationCounts { nij, Ti, root_state },
      &InferGtrOptions {
        pc: 0.1,
        ..InferGtrOptions::default()
      },
    )?;

    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      array![
        [ 0.0000000000000000,  0.7358152606066288,  1.8250024541741396,  1.1480276310375928],
        [ 0.7358152606066288,  0.0000000000000000,  1.9426167902218199,  1.2779572118106166],
        [ 1.8250024541741396,  1.9426167902218199,  0.0000000000000000,  1.3712594037821175],
        [ 1.1480276310375928,  1.2779572118106166,  1.3712594037821175,  0.0000000000000000],
      ],
      &actual.W,
      epsilon = 1e-9
    );

    pretty_assert_ulps_eq!(
      array![
        0.209080146491563,
        0.245288114914305,
        0.209258588641852,
        0.336373149952278,
      ],
      &actual.pi,
      epsilon = 1e-9
    );

    pretty_assert_ulps_eq!(0.4004706866848004, actual.mu, epsilon = 1e-9);

    pretty_assert_ulps_eq!(1.0, avg_transition(&actual.W, &actual.pi)?, epsilon = 1e-5);

    Ok(())
  }

  #[test]
  fn test_infer_gtr_with_mutation_counts() -> Result<(), Report> {
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

    let counts_actual = get_mutation_counts(&graph, &partition)?;
    let counts_expected = MutationCounts {
      nij: array![[0., 0., 0., 0.], [2., 0., 0., 1.], [3., 2., 0., 0.], [0., 1., 1., 0.]],
      Ti: array![1.98, 2.945, 2.515, 2.64],
      root_state: array![4.0, 3.0, 3.0, 4.0],
    };

    pretty_assert_ulps_eq!(counts_expected.nij, counts_actual.nij, epsilon = 1e-9);
    pretty_assert_ulps_eq!(counts_expected.Ti, counts_actual.Ti, epsilon = 1e-7);
    pretty_assert_ulps_eq!(counts_expected.root_state, counts_actual.root_state, epsilon = 1e-9);

    let actual = infer_gtr_impl(
      &counts_actual,
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

  #[test]
  fn test_infer_gtr_distance_1() {
    let actual = distance(&array![0.0, 0.0, 0.0, 0.0], &array![0.25, 0.25, 0.25, 0.25]);
    let expected = 0.5;
    pretty_assert_ulps_eq!(actual, expected, epsilon = 1e-5);
  }

  #[test]
  fn test_infer_gtr_distance_2() {
    let actual = distance(
      &array![0.25, 0.25, 0.25, 0.25],
      &array![0.16031624, 0.24247873, 0.3087257, 0.28847933],
    );
    let expected = 0.11414514039292292;
    pretty_assert_ulps_eq!(actual, expected, epsilon = 1e-5);
  }

  #[test]
  fn test_infer_gtr_distance_3() {
    let actual = distance(
      &array![0.16031624, 0.24247873, 0.3087257, 0.28847933],
      &array![0.14968981, 0.24040983, 0.31181279, 0.29808757],
    );
    let expected = 0.014800329884495377;
    pretty_assert_ulps_eq!(actual, expected, epsilon = 1e-5);
  }

  #[test]
  fn test_infer_gtr_distance_4() {
    let actual = distance(
      &array![0.14878286, 0.24052002, 0.31238389, 0.29831323],
      &array![0.14878922, 0.24051699, 0.31239221, 0.29830159],
    );
    let expected = 1.59484629332816e-05;
    pretty_assert_ulps_eq!(actual, expected, epsilon = 1e-5);
  }

  #[test]
  fn test_infer_gtr_distance_5() {
    let same = &array![0.14878286, 0.24052002, 0.31238389, 0.29831323];
    let actual = distance(same, same);
    let expected = 0.0;
    pretty_assert_ulps_eq!(actual, expected, max_ulps = 3);
  }
}

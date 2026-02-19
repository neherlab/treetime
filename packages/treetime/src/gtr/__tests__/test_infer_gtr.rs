#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::avg_transition;
  use crate::gtr::infer_gtr::{
    InferGtrOptions, MutationCounts, distance, get_mutation_counts, get_mutation_counts_dense, infer_gtr_impl,
  };
  use crate::pretty_assert_ulps_eq;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use ndarray::array;
  use parking_lot::RwLock;
  use serde::Deserialize;
  use std::path::PathBuf;
  use std::sync::Arc;
  use treetime_io::fasta::{read_many_fasta, read_many_fasta_str};
  use treetime_io::nwk::{nwk_read_file, nwk_read_str};

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

  // ============================================================================
  // Golden master test for dense GTR inference
  // ============================================================================

  fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .map(PathBuf::from)
      .expect("project has workspace root")
  }

  #[derive(Debug, Deserialize)]
  struct GtrInferenceDenseFixture {
    mutation_counts: MutationCountsFixture,
    inferred_gtr: InferredGtrFixture,
  }

  #[derive(Debug, Deserialize)]
  struct MutationCountsFixture {
    nij: Vec<Vec<f64>>,
    #[serde(rename = "Ti")]
    ti: Vec<f64>,
    root_state: Vec<f64>,
  }

  #[derive(Debug, Deserialize)]
  struct InferredGtrFixture {
    #[serde(rename = "W")]
    w: Vec<Vec<f64>>,
    pi: Vec<f64>,
    mu: f64,
  }

  /// Golden master test: dense GTR inference matches Python v0 reference.
  ///
  /// Loads flu/h3n2/20 dataset, runs dense marginal reconstruction with JC69,
  /// extracts mutation counts, infers GTR, and compares against fixture from v0.
  #[test]
  fn test_infer_gtr_dense_golden_master() -> Result<(), Report> {
    let fixture_path =
      PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/gtr/__tests__/fixtures/gtr_inference_dense.json");
    let fixture_str = std::fs::read_to_string(&fixture_path)?;
    let fixture: GtrInferenceDenseFixture = serde_json::from_str(&fixture_str)?;

    let root = project_root();
    let tree_path = root.join("data/flu/h3n2/20/tree.nwk");
    let aln_path = root.join("data/flu/h3n2/20/aln.fasta.xz");

    let treat_gap_as_unknown = true;
    let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
    let aln = read_many_fasta(&[aln_path], &alphabet)?;

    let graph: GraphAncestral = nwk_read_file(&tree_path)?;

    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown,
      ..JC69Params::default()
    })?;

    let partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    initialize_marginal(&graph, std::slice::from_ref(&partition), &aln)?;

    let counts = get_mutation_counts_dense(&graph, &partition)?;

    let expected_nij = array![
      [
        fixture.mutation_counts.nij[0][0],
        fixture.mutation_counts.nij[0][1],
        fixture.mutation_counts.nij[0][2],
        fixture.mutation_counts.nij[0][3]
      ],
      [
        fixture.mutation_counts.nij[1][0],
        fixture.mutation_counts.nij[1][1],
        fixture.mutation_counts.nij[1][2],
        fixture.mutation_counts.nij[1][3]
      ],
      [
        fixture.mutation_counts.nij[2][0],
        fixture.mutation_counts.nij[2][1],
        fixture.mutation_counts.nij[2][2],
        fixture.mutation_counts.nij[2][3]
      ],
      [
        fixture.mutation_counts.nij[3][0],
        fixture.mutation_counts.nij[3][1],
        fixture.mutation_counts.nij[3][2],
        fixture.mutation_counts.nij[3][3]
      ],
    ];
    let expected_ti = array![
      fixture.mutation_counts.ti[0],
      fixture.mutation_counts.ti[1],
      fixture.mutation_counts.ti[2],
      fixture.mutation_counts.ti[3],
    ];
    let expected_root_state = array![
      fixture.mutation_counts.root_state[0],
      fixture.mutation_counts.root_state[1],
      fixture.mutation_counts.root_state[2],
      fixture.mutation_counts.root_state[3],
    ];

    pretty_assert_ulps_eq!(expected_nij, counts.nij, epsilon = 1e-7);
    pretty_assert_ulps_eq!(expected_ti, counts.Ti, epsilon = 1e-5);
    pretty_assert_ulps_eq!(expected_root_state, counts.root_state, epsilon = 1e-7);

    let result = infer_gtr_impl(&counts, &InferGtrOptions::default())?;

    let expected_w = array![
      [
        fixture.inferred_gtr.w[0][0],
        fixture.inferred_gtr.w[0][1],
        fixture.inferred_gtr.w[0][2],
        fixture.inferred_gtr.w[0][3]
      ],
      [
        fixture.inferred_gtr.w[1][0],
        fixture.inferred_gtr.w[1][1],
        fixture.inferred_gtr.w[1][2],
        fixture.inferred_gtr.w[1][3]
      ],
      [
        fixture.inferred_gtr.w[2][0],
        fixture.inferred_gtr.w[2][1],
        fixture.inferred_gtr.w[2][2],
        fixture.inferred_gtr.w[2][3]
      ],
      [
        fixture.inferred_gtr.w[3][0],
        fixture.inferred_gtr.w[3][1],
        fixture.inferred_gtr.w[3][2],
        fixture.inferred_gtr.w[3][3]
      ],
    ];
    let expected_pi = array![
      fixture.inferred_gtr.pi[0],
      fixture.inferred_gtr.pi[1],
      fixture.inferred_gtr.pi[2],
      fixture.inferred_gtr.pi[3],
    ];

    pretty_assert_ulps_eq!(expected_w, result.W, epsilon = 1e-7);
    pretty_assert_ulps_eq!(expected_pi, result.pi, epsilon = 1e-7);
    pretty_assert_ulps_eq!(fixture.inferred_gtr.mu, result.mu, epsilon = 1e-7);

    Ok(())
  }
}

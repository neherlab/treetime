#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::get_common_length;
  use crate::commands::ancestral::marginal::initialize_marginal;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::infer_gtr::common::{InferGtrOptions, infer_gtr_impl};
  use crate::gtr::infer_gtr::dense::get_mutation_counts_dense;
  use crate::pretty_assert_ulps_eq;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use ndarray::{Array1, Array2, array};
  use parking_lot::RwLock;
  use serde::Deserialize;
  use std::path::PathBuf;
  use std::sync::Arc;
  use treetime_io::fasta::{FastaRecord, read_many_fasta, read_many_fasta_str};
  use treetime_io::json::json_read_file;
  use treetime_io::nwk::{nwk_read_file, nwk_read_str};
  use treetime_utils::array::serde::{array1_from_vec, array2_from_vec};

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
  }

  fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .map(PathBuf::from)
      .expect("project has workspace root")
  }

  /// Helper to create a dense partition and run marginal reconstruction.
  fn setup_dense_partition(
    tree_nwk: &str,
    aln: &[FastaRecord],
    treat_gap_as_unknown: bool,
  ) -> Result<(GraphAncestral, Arc<RwLock<PartitionMarginalDense>>), Report> {
    let graph: GraphAncestral = nwk_read_str(tree_nwk)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown,
      ..JC69Params::default()
    })?;

    let partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    initialize_marginal(&graph, std::slice::from_ref(&partition), aln)?;
    Ok((graph, partition))
  }

  // ============================================================================
  // Golden master test
  // ============================================================================

  /// Golden master values for dense GTR inference test.
  ///
  /// Combines mutation counts (nij, Ti, root_state) and inferred GTR parameters (W, pi, mu).
  #[derive(Debug, Deserialize)]
  struct InferGtrDenseGoldenMaster {
    #[serde(deserialize_with = "array2_from_vec")]
    nij: Array2<f64>,
    #[serde(deserialize_with = "array1_from_vec")]
    Ti: Array1<f64>,
    #[serde(deserialize_with = "array1_from_vec")]
    root_state: Array1<f64>,
    #[serde(deserialize_with = "array2_from_vec")]
    W: Array2<f64>,
    #[serde(deserialize_with = "array1_from_vec")]
    pi: Array1<f64>,
    mu: f64,
  }

  /// Golden master test: dense GTR inference matches Python v0 reference.
  ///
  /// Dataset: flu/h3n2/20
  ///
  /// # Regenerate fixture
  ///
  /// Run from project root in container:
  /// ```bash
  /// ./dev/docker/python test_scripts/dump_gtr_inference_dense.py
  /// ```
  #[test]
  fn test_golden_master() -> Result<(), Report> {
    let fixture_path =
      PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/gtr/infer_gtr/__tests__/fixtures/gtr_inference_dense.json");
    let expected: InferGtrDenseGoldenMaster = json_read_file(&fixture_path)?;

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

    pretty_assert_ulps_eq!(expected.nij, counts.nij, epsilon = 1e-7);
    pretty_assert_ulps_eq!(expected.Ti, counts.Ti, epsilon = 1e-5);
    pretty_assert_ulps_eq!(expected.root_state, counts.root_state, epsilon = 1e-7);

    let result = infer_gtr_impl(&counts, &InferGtrOptions::default())?;

    pretty_assert_ulps_eq!(expected.W, result.W, epsilon = 1e-7);
    pretty_assert_ulps_eq!(expected.pi, result.pi, epsilon = 1e-7);
    pretty_assert_ulps_eq!(expected.mu, result.mu, epsilon = 1e-7);

    Ok(())
  }

  // ============================================================================
  // Unit tests
  // ============================================================================

  /// All sequences identical -> zero off-diagonal nij (no substitutions).
  #[test]
  fn test_uniform_sequences() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGT
      >B
      ACGT
      >C
      ACGT
      >D
      ACGT
      "#},
      &*NUC_ALPHABET,
    )?;

    let (graph, partition) =
      setup_dense_partition("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;", &aln, true)?;

    let counts = get_mutation_counts_dense(&graph, &partition)?;

    // No substitutions -> all off-diagonal elements are zero
    for i in 0..4 {
      for j in 0..4 {
        if i != j {
          pretty_assert_ulps_eq!(0.0, counts.nij[[i, j]], epsilon = 1e-9);
        }
      }
    }

    // Root composition: 1 of each nucleotide (ACGT)
    pretty_assert_ulps_eq!(array![1.0, 1.0, 1.0, 1.0], counts.root_state, epsilon = 1e-9);

    Ok(())
  }

  /// Single mutation from A to C at one position.
  #[test]
  fn test_single_mutation() -> Result<(), Report> {
    // Simple tree: root -> A, root -> B
    // A has "ACGT", B has "CCGT" (A->C at position 0)
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGT
      >B
      CCGT
      "#},
      &*NUC_ALPHABET,
    )?;

    let (graph, partition) = setup_dense_partition("(A:0.1,B:0.1)root:0.0;", &aln, true)?;

    let counts = get_mutation_counts_dense(&graph, &partition)?;

    // The root should be reconstructed with some state at position 0.
    // With marginal reconstruction on a symmetric tree, the root gets
    // probability split between A and C at position 0.
    // Due to argmax, it picks one (likely A or C based on alphabet order).

    // Check that nij has at least one non-zero off-diagonal element
    // (there's definitely a mutation from parent to child on one edge)
    let total_mutations: f64 = counts
      .nij
      .iter()
      .enumerate()
      .filter_map(|(idx, &val)| {
        let i = idx / 4;
        let j = idx % 4;
        (i != j).then_some(val)
      })
      .sum();

    // At least one mutation expected (A->C or C->A depending on root reconstruction)
    assert!(
      total_mutations >= 1.0,
      "Expected at least 1 mutation, got {total_mutations}"
    );

    Ok(())
  }

  /// Dense inference produces valid GTR model on realistic data.
  /// Verifies that the inferred model has valid properties (symmetric W, normalized pi, positive mu).
  #[test]
  fn test_produces_valid_model() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACATCGCCGTAGAC
      >B
      GCATCCCTGTAGGG
      >C
      CCGGCGATGTGTTG
      >D
      TCGGCCGTGTGTTG
      "#},
      &*NUC_ALPHABET,
    )?;

    let tree_nwk = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let (graph, partition) = setup_dense_partition(tree_nwk, &aln, true)?;

    let counts = get_mutation_counts_dense(&graph, &partition)?;
    let result = infer_gtr_impl(&counts, &InferGtrOptions::default())?;

    // W should be symmetric
    for i in 0..4 {
      for j in 0..4 {
        pretty_assert_ulps_eq!(result.W[[i, j]], result.W[[j, i]], epsilon = 1e-9);
      }
    }

    // pi should sum to 1
    let pi_sum: f64 = result.pi.iter().sum();
    pretty_assert_ulps_eq!(1.0, pi_sum, epsilon = 1e-9);

    // All pi values should be positive
    for &p in &result.pi {
      assert!(p > 0.0, "pi values should be positive, got {p}");
    }

    // mu should be positive
    assert!(result.mu > 0.0, "mu should be positive, got {}", result.mu);

    Ok(())
  }
}

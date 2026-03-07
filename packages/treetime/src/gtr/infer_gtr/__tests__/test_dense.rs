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
  use ndarray::array;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
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

  /// All sequences identical -> near-zero off-diagonal nij.
  ///
  /// With fractional expected counts, even identical sequences have small
  /// non-zero mutation probabilities due to the probabilistic nature of
  /// the joint distribution. These should be negligible compared to actual
  /// substitution counts.
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

    // Off-diagonal elements should be negligible (< 0.1 per site pair)
    // With fractional counts, there's small probability mass on off-diagonal
    // states even when sequences are identical
    for i in 0..4 {
      for j in 0..4 {
        if i != j {
          assert!(
            counts.nij[[i, j]] < 0.1,
            "Off-diagonal nij[{i},{j}] = {} should be negligible",
            counts.nij[[i, j]]
          );
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

  // TODO(dense-branch-clamp): this test originally verified clean mathematical properties of
  // zero branch lengths: `Ti == 0` (exact, epsilon 1e-15) and `off_diagonal_nij < 0.1`.
  // After adding `fix_branch_length` clamping to `get_mutation_counts_dense` (for consistency
  // with marginal passes and v0), zero-input BLs become `MIN_BRANCH_LENGTH_FRACTION / L`.
  //
  // Consequences:
  //  - Ti is small but non-zero (proportional to clamped BL ~2.5e-4 for L=4)
  //  - Off-diagonal nij ~2.0: the near-identity `expQt` has small off-diagonal entries, but
  //    profiles (computed with clamped BL during marginal reconstruction) carry actual mutation
  //    signal (A-C at position 0), which the joint distribution picks up
  //  - The `is_finite()` assertion is weaker than the original `< 0.1` bound
  //
  // If `fix_branch_length` is removed from `get_mutation_counts_dense`, restore original
  // assertions: `Ti == 0` (epsilon 1e-15) and `total_off_diagonal < 0.1`.
  /// Verifies that zero-input branch lengths (clamped internally) produce bounded results.
  #[test]
  fn test_zero_branch_lengths() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGT
      >B
      CCGT
      >C
      ACGT
      >D
      CCGT
      "#},
      &*NUC_ALPHABET,
    )?;

    let (graph, partition) = setup_dense_partition("((A:0.0,B:0.0)AB:0.0,(C:0.0,D:0.0)CD:0.0)root:0.0;", &aln, true)?;

    let counts = get_mutation_counts_dense(&graph, &partition)?;

    // Ti is small (proportional to clamped BL, not zero)
    for &ti in &counts.Ti {
      assert!(ti < 0.1, "Ti should be small with zero-input branch lengths, got {ti}");
    }

    // Off-diagonal nij is bounded (not zero due to clamping + profile mutation signal)
    let total_nij: f64 = counts.nij.sum();
    assert!(total_nij.is_finite(), "nij should be finite, got {total_nij}");

    // Root state should still be computed correctly
    assert!(counts.root_state.sum() > 0.0, "root_state should be populated");

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

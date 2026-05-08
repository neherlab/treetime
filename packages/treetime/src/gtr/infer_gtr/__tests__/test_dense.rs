#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::initialize_marginal;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::infer_gtr::common::{InferGtrOptions, infer_gtr_impl};
  use crate::gtr::infer_gtr::dense::{
    accumulate_mutation_counts, get_branch_mutation_matrix, get_mutation_counts_dense,
  };
  use crate::pretty_assert_abs_diff_eq;
  use crate::pretty_assert_ulps_eq;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use ndarray::{Array1, Array2, array};
  use parking_lot::RwLock;
  use std::slice::from_ref;
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
  ) -> Result<(GraphAncestral, Arc<RwLock<PartitionMarginalDense>>), Report> {
    let graph: GraphAncestral = nwk_read_str(tree_nwk)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
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

    initialize_marginal(&graph, from_ref(&partition), aln)?;
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

    let (graph, partition) = setup_dense_partition("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;", &aln)?;

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

    let (graph, partition) = setup_dense_partition("(A:0.1,B:0.1)root:0.0;", &aln)?;

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

  /// Verifies that zero-input branch lengths (clamped internally) produce bounded results.
  ///
  /// For a zero-length branch, the mathematically correct result is expQt(0) = I, Ti = 0,
  /// zero off-diagonal nij. But get_mutation_counts_dense applies fix_branch_length clamping
  /// for consistency with the marginal passes that computed the profiles it reads. Using
  /// expQt(0) when profiles were propagated through expQt(clamped_BL) would mix two different
  /// transition models on the same branch.
  ///
  /// Clamped BL = MIN_BRANCH_LENGTH_FRACTION / seq_length (here: 1e-3/4 = 2.5e-4).
  /// Consequences:
  ///  - Ti is small but non-zero (proportional to clamped BL)
  ///  - Off-diagonal nij ~2.0: the near-identity expQt has small off-diagonal entries, but
  ///    profiles (computed with clamped BL during marginal reconstruction) carry actual mutation
  ///    signal (A<->C at position 0), which the joint distribution picks up
  ///
  /// The mathematical invariant (Ti == 0, off-diagonal nij == 0) is tested separately in
  /// test_zero_branch_lengths_unclamped, which bypasses clamping.
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

    let (graph, partition) = setup_dense_partition("((A:0.0,B:0.0)AB:0.0,(C:0.0,D:0.0)CD:0.0)root:0.0;", &aln)?;

    let counts = get_mutation_counts_dense(&graph, &partition)?;

    // Ti proportional to clamped BL (~2.5e-4), bounded well below 1e-2
    let ti_max = counts.Ti.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    assert!(ti_max < 1e-2, "max(Ti) = {ti_max} should be < 1e-2 for clamped zero-BL");

    // Off-diagonal nij dominated by A<->C mutation signal from profiles (~2.0 total),
    // bounded well below 5.0
    let total_off_diagonal: f64 = counts
      .nij
      .indexed_iter()
      .filter(|&((i, j), _)| i != j)
      .map(|(_, &v)| v)
      .sum();
    assert!(
      total_off_diagonal < 5.0,
      "total off-diagonal nij = {total_off_diagonal} should be < 5.0 for clamped zero-BL"
    );

    // Root state should still be computed correctly
    assert!(counts.root_state.sum() > 0.0, "root_state should be populated");

    Ok(())
  }

  /// Mathematical invariant: with expQt = I and branch_length = 0, the accumulation
  /// functions produce exactly zero Ti and zero off-diagonal nij.
  ///
  /// Tests the accumulation logic directly, bypassing fix_branch_length clamping.
  /// When a branch has zero length, the transition matrix is the identity (no evolution),
  /// so child = parent with certainty. No mutations occur and no evolutionary time elapses.
  #[test]
  fn test_zero_branch_lengths_unclamped() -> Result<(), Report> {
    let n_states = 4;

    // Identity transition matrix: P(child=i | parent=j) = delta(i,j)
    let exp_qt = Array2::eye(n_states);

    // One-hot profiles: site 0 = A, site 1 = C, site 2 = G, site 3 = T
    let messages = Array2::eye(n_states);

    let mut_stack = get_branch_mutation_matrix(&messages, &messages, &exp_qt);

    let mut nij = Array2::zeros((n_states, n_states));
    let mut ti = Array1::zeros(n_states);

    accumulate_mutation_counts(&mut_stack, 0.0, &mut nij, &mut ti);

    // Ti == 0: no evolutionary time at zero branch length
    let expected_ti = Array1::zeros(n_states);
    pretty_assert_abs_diff_eq!(expected_ti, ti, epsilon = 1e-15);

    // nij == I: identity expQt means child = parent at every site, so each of the 4
    // sites contributes P(child=k, parent=k) = 1.0 on the diagonal, zero off-diagonal
    let expected_nij = Array2::eye(n_states);
    pretty_assert_abs_diff_eq!(expected_nij, nij, epsilon = 1e-15);

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
    let (graph, partition) = setup_dense_partition(tree_nwk, &aln)?;

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

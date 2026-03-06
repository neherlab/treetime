#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::get_common_length;
  use crate::commands::ancestral::marginal::initialize_marginal;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::pretty_assert_ulps_eq;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::prelude::*;
  use parking_lot::RwLock;
  use std::sync::{Arc, LazyLock};
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  /// Numerical tolerance for likelihood comparisons.
  /// The implementation performs multiple normalization steps that accumulate
  /// floating-point errors. A tolerance of 1e-7 allows for this while still
  /// catching any significant algorithmic errors.
  const LIKELIHOOD_EPSILON: f64 = 1e-7;

  /// Compute analytical likelihood for two-taxon tree: (A:t1,B:t2)root;
  ///
  /// Formula (Felsenstein pruning for two leaves):
  ///   L = sum_s pi[s] * P(obs_A|s,t1) * P(obs_B|s,t2)
  ///
  /// where P(obs|s,t) = expQt[obs,s] is the transition probability from state s
  /// to observed state obs over branch length t.
  fn analytical_two_taxon_likelihood(gtr: &GTR, obs_a: usize, obs_b: usize, t1: f64, t2: f64) -> f64 {
    let exp_qt1 = gtr.expQt(t1);
    let exp_qt2 = gtr.expQt(t2);

    let mut likelihood = 0.0;
    for s in 0..gtr.pi.len() {
      // P(obs_A|s,t1) * P(obs_B|s,t2) * pi[s]
      // expQt[i,j] = P(i|j,t) = probability of observing state i given ancestral state j
      likelihood += gtr.pi[s] * exp_qt1[[obs_a, s]] * exp_qt2[[obs_b, s]];
    }
    likelihood
  }

  /// Compute analytical likelihood for star tree: (A:t,B:t,C:t,D:t)root;
  ///
  /// Formula:
  ///   L = sum_s pi[s] * prod_i P(obs_i|s,t)
  fn analytical_star_tree_likelihood(gtr: &GTR, observations: &[usize], t: f64) -> f64 {
    let exp_qt = gtr.expQt(t);

    let mut likelihood = 0.0;
    for s in 0..gtr.pi.len() {
      let mut product = gtr.pi[s];
      for &obs in observations {
        product *= exp_qt[[obs, s]];
      }
      likelihood += product;
    }
    likelihood
  }

  /// Compute analytical likelihood for three-taxon tree: ((A:t_a,B:t_b)AB:t_ab,C:t_c)root;
  ///
  /// Formula (full Felsenstein pruning):
  ///   L = sum_{s_root} sum_{s_AB} pi[s_root]
  ///       * P(obs_C | s_root, t_C)
  ///       * P(s_AB | s_root, t_AB)
  ///       * P(obs_A | s_AB, t_A)
  ///       * P(obs_B | s_AB, t_B)
  fn analytical_three_taxon_likelihood(
    gtr: &GTR,
    obs_a: usize,
    obs_b: usize,
    obs_c: usize,
    t_a: f64,
    t_b: f64,
    t_ab: f64,
    t_c: f64,
  ) -> f64 {
    let exp_qt_a = gtr.expQt(t_a);
    let exp_qt_b = gtr.expQt(t_b);
    let exp_qt_ab = gtr.expQt(t_ab);
    let exp_qt_c = gtr.expQt(t_c);

    let mut likelihood = 0.0;
    for s_root in 0..4 {
      for s_ab in 0..4 {
        let msg_ab = exp_qt_a[[obs_a, s_ab]] * exp_qt_b[[obs_b, s_ab]];
        likelihood += gtr.pi[s_root] * exp_qt_c[[obs_c, s_root]] * exp_qt_ab[[s_ab, s_root]] * msg_ab;
      }
    }
    likelihood
  }

  /// Map a nucleotide character to its index in the state vector: A=0, C=1, G=2, T=3.
  fn state_index(c: char) -> usize {
    match c {
      'A' => 0,
      'C' => 1,
      'G' => 2,
      'T' => 3,
      _ => panic!("Invalid nucleotide: {c}"),
    }
  }

  /// Run dense marginal ancestral reconstruction and return only the total log-likelihood.
  ///
  /// Convenience wrapper: parses the Newick tree and FASTA alignment from strings,
  /// constructs a single dense partition with the given GTR model, and runs
  /// Felsenstein's pruning algorithm via `initialize_marginal`.
  fn run_dense_marginal_get_log_lh(newick: &str, aln_str: &str, gtr: GTR) -> Result<f64, Report> {
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    let log_lh = initialize_marginal(&graph, &partitions, &aln)?;
    Ok(log_lh)
  }

  // ============================================================================
  // T1: Two-taxon analytical likelihood tests
  // ============================================================================

  /// Two-taxon tree with identical leaf states under JC69 (Jukes-Cantor 1969).
  ///
  /// Tree: (A:0.1, B:0.2)root; both leaves observe state 'A'.
  ///
  /// Analytical: L = sum_s pi[s] * P(A|s,0.1) * P(A|s,0.2).
  /// Under JC69, pi = [0.25, 0.25, 0.25, 0.25] and P(i|j,t) = expQt[i,j].
  /// Matched states yield high likelihood because the diagonal of expQt dominates
  /// at short branch lengths.
  #[test]
  fn test_two_taxon_analytical_jc69_same_state() -> Result<(), Report> {
    // Tree: (A:0.1,B:0.2)root;
    // Both leaves have same nucleotide 'A' at single position
    // Analytical: L = sum_s pi[s] * expQt1[0,s] * expQt2[0,s]
    let gtr = jc69(JC69Params::default())?;
    let t1 = 0.1;
    let t2 = 0.2;

    let expected_lh = analytical_two_taxon_likelihood(&gtr, 0, 0, t1, t2);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.1,B:0.2)root;";
    let aln = ">A\nA\n>B\nA\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = LIKELIHOOD_EPSILON);
    Ok(())
  }

  /// Two-taxon tree with different leaf states under JC69.
  ///
  /// Tree: (A:0.1, B:0.2)root; leaf A observes 'A', leaf B observes 'T'.
  ///
  /// Analytical: L = sum_s pi[s] * P(A|s,0.1) * P(T|s,0.2).
  /// Mismatched states require at least one substitution event, so the likelihood is
  /// lower than for matched states at these branch lengths.
  #[test]
  fn test_two_taxon_analytical_jc69_different_states() -> Result<(), Report> {
    // Tree: (A:0.1,B:0.2)root;
    // A has 'A', B has 'T'
    let gtr = jc69(JC69Params::default())?;
    let t1 = 0.1;
    let t2 = 0.2;

    let expected_lh = analytical_two_taxon_likelihood(&gtr, state_index('A'), state_index('T'), t1, t2);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.1,B:0.2)root;";
    let aln = ">A\nA\n>B\nT\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = LIKELIHOOD_EPSILON);
    Ok(())
  }

  /// Two-taxon tree with non-uniform equilibrium frequencies (GTR model).
  ///
  /// Uses pi = [0.4, 0.1, 0.2, 0.3] instead of JC69's uniform [0.25, 0.25, 0.25, 0.25].
  /// Non-uniform pi tests that the implementation correctly weights ancestral states by
  /// their equilibrium probability, which is the root prior in Felsenstein's pruning.
  ///
  /// Leaf A observes 'A', leaf B observes 'C'. The A->C transition under non-uniform pi
  /// has different likelihood than under JC69 because pi[A]=0.4 makes 'A' a more probable
  /// root state.
  #[test]
  fn test_two_taxon_analytical_nonuniform_pi() -> Result<(), Report> {
    // Non-uniform equilibrium frequencies reveal bugs in pi weighting
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.4, 0.1, 0.2, 0.3],
    })?;

    let t1 = 0.15;
    let t2 = 0.25;

    // Test A->C transition
    let expected_lh = analytical_two_taxon_likelihood(&gtr, state_index('A'), state_index('C'), t1, t2);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.15,B:0.25)root;";
    let aln = ">A\nA\n>B\nC\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = LIKELIHOOD_EPSILON);
    Ok(())
  }

  /// Two-taxon tree with a 3-position alignment under JC69.
  ///
  /// Sequences: A="ACG", B="TCA". Each alignment position is independent under the
  /// Felsenstein model, so the total log-likelihood is the sum of per-position
  /// log-likelihoods:
  ///   ln(L_total) = sum_i ln(L_i)
  /// where L_i = sum_s pi[s] * P(obs_A_i|s,t1) * P(obs_B_i|s,t2).
  ///
  /// This verifies that multi-position likelihoods are correctly aggregated from
  /// independent single-position computations.
  #[test]
  fn test_two_taxon_analytical_multiple_positions() -> Result<(), Report> {
    // Multiple independent positions: total likelihood = product of per-position likelihoods
    let gtr = jc69(JC69Params::default())?;
    let t1 = 0.1;
    let t2 = 0.2;

    // Sequence: ACG for A, TCA for B
    let positions = [
      (state_index('A'), state_index('T')),
      (state_index('C'), state_index('C')),
      (state_index('G'), state_index('A')),
    ];

    let mut expected_log_lh = 0.0;
    for (obs_a, obs_b) in positions {
      let lh = analytical_two_taxon_likelihood(&gtr, obs_a, obs_b, t1, t2);
      expected_log_lh += lh.ln();
    }

    let newick = "(A:0.1,B:0.2)root;";
    let aln = ">A\nACG\n>B\nTCA\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = LIKELIHOOD_EPSILON);
    Ok(())
  }

  /// Two-taxon tree with highly asymmetric branch lengths under JC69.
  ///
  /// Tree: (A:0.01, B:1.0)root; both leaves observe 'G'.
  /// Branch A is very short (t=0.01, near-identity transition matrix) while branch B
  /// is long (t=1.0, significant probability of substitution).
  ///
  /// Tests numerical stability when expQt matrices have very different magnitudes.
  /// The short branch contributes near-1.0 diagonal elements while the long branch
  /// contributes values closer to equilibrium.
  #[test]
  fn test_two_taxon_analytical_asymmetric_branches() -> Result<(), Report> {
    // Highly asymmetric branch lengths
    let gtr = jc69(JC69Params::default())?;
    let t1 = 0.01; // Very short
    let t2 = 1.0; // Long

    let expected_lh = analytical_two_taxon_likelihood(&gtr, state_index('G'), state_index('G'), t1, t2);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.01,B:1.0)root;";
    let aln = ">A\nG\n>B\nG\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = LIKELIHOOD_EPSILON);
    Ok(())
  }

  // ============================================================================
  // T2: Star tree analytical likelihood tests
  // ============================================================================

  /// Star tree (4 leaves directly from root) with all identical states under JC69.
  ///
  /// Tree: (A:0.1, B:0.1, C:0.1, D:0.1)root; all leaves observe 'A'.
  ///
  /// Analytical: L = sum_s pi[s] * prod_i P(A|s,0.1).
  /// The star topology is a polytomy at the root with no internal nodes. This tests
  /// that Felsenstein's pruning correctly handles more than 2 children at the root.
  /// With all states identical, the diagonal of expQt dominates the product.
  #[test]
  fn test_star_tree_analytical_jc69_all_same() -> Result<(), Report> {
    // Star tree: (A:0.1,B:0.1,C:0.1,D:0.1)root;
    // All leaves have 'A'
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;
    let observations = [0, 0, 0, 0]; // All 'A'

    let expected_lh = analytical_star_tree_likelihood(&gtr, &observations, t);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.1,B:0.1,C:0.1,D:0.1)root;";
    let aln = ">A\nA\n>B\nA\n>C\nA\n>D\nA\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = LIKELIHOOD_EPSILON);
    Ok(())
  }

  /// Star tree with all four distinct nucleotide states under JC69.
  ///
  /// Tree: (A:0.2, B:0.2, C:0.2, D:0.2)root; leaves observe A, C, G, T respectively.
  ///
  /// Analytical: L = sum_s pi[s] * P(A|s,0.2) * P(C|s,0.2) * P(G|s,0.2) * P(T|s,0.2).
  /// With all four states present, no single ancestral state avoids at least 3 substitutions.
  /// Under JC69's symmetry, each root state contributes equally, yielding a low total
  /// likelihood that exercises the off-diagonal elements of expQt.
  #[test]
  fn test_star_tree_analytical_jc69_mixed_states() -> Result<(), Report> {
    // Star tree with different states at leaves
    let gtr = jc69(JC69Params::default())?;
    let t = 0.2;
    let observations = [state_index('A'), state_index('C'), state_index('G'), state_index('T')];

    let expected_lh = analytical_star_tree_likelihood(&gtr, &observations, t);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.2,B:0.2,C:0.2,D:0.2)root;";
    let aln = ">A\nA\n>B\nC\n>C\nG\n>D\nT\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = LIKELIHOOD_EPSILON);
    Ok(())
  }

  /// Star tree with non-uniform equilibrium frequencies (GTR model).
  ///
  /// Uses pi = [0.1, 0.2, 0.3, 0.4] with leaves observing T, T, G, C.
  ///
  /// Analytical: L = sum_s pi[s] * P(T|s,0.15) * P(T|s,0.15) * P(G|s,0.15) * P(C|s,0.15).
  /// Non-uniform pi breaks JC69's state symmetry: root state T (pi=0.4) contributes more
  /// than root state A (pi=0.1) to the sum, biasing the likelihood toward T-rich scenarios.
  #[test]
  fn test_star_tree_analytical_nonuniform_pi() -> Result<(), Report> {
    // Non-uniform equilibrium frequencies with star tree
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.1, 0.2, 0.3, 0.4],
    })?;

    let t = 0.15;
    let observations = [state_index('T'), state_index('T'), state_index('G'), state_index('C')];

    let expected_lh = analytical_star_tree_likelihood(&gtr, &observations, t);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.15,B:0.15,C:0.15,D:0.15)root;";
    let aln = ">A\nT\n>B\nT\n>C\nG\n>D\nC\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = LIKELIHOOD_EPSILON);
    Ok(())
  }

  // ============================================================================
  // T3: Single-position exhaustive verification
  // ============================================================================

  /// Three-taxon tree with a single position, verified against closed-form Felsenstein pruning.
  ///
  /// Tree: ((A:0.1, B:0.2)AB:0.15, C:0.25)root; leaves observe A, C, G.
  ///
  /// Analytical (full enumeration over internal node AB):
  ///   L = sum_{s_root} sum_{s_AB} pi[s_root]
  ///       * P(G|s_root, 0.25) * P(s_AB|s_root, 0.15)
  ///       * P(A|s_AB, 0.1) * P(C|s_AB, 0.2)
  ///
  /// This is the simplest tree topology that exercises the recursive message-passing
  /// structure of Felsenstein's algorithm: the AB subtree produces a partial likelihood
  /// vector that is combined with the C leaf message at the root.
  #[test]
  fn test_three_taxon_single_position_exhaustive() -> Result<(), Report> {
    // Tree: ((A:0.1,B:0.2)AB:0.15,C:0.25)root;
    let gtr = jc69(JC69Params::default())?;

    let t_a = 0.1;
    let t_b = 0.2;
    let t_ab = 0.15;
    let t_c = 0.25;

    let obs_a = state_index('A');
    let obs_b = state_index('C');
    let obs_c = state_index('G');

    let expected_lh = analytical_three_taxon_likelihood(&gtr, obs_a, obs_b, obs_c, t_a, t_b, t_ab, t_c);
    let expected_log_lh = expected_lh.ln();

    let newick = "((A:0.1,B:0.2)AB:0.15,C:0.25)root;";
    let aln = ">A\nA\n>B\nC\n>C\nG\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = LIKELIHOOD_EPSILON);
    Ok(())
  }

  /// Exhaustive verification of all 64 (4^3) single-position state combinations on a
  /// three-taxon tree under JC69.
  ///
  /// Tree: ((A:0.1, B:0.2)AB:0.15, C:0.25)root.
  ///
  /// For each triplet (state_A, state_B, state_C) in {A,C,G,T}^3, the implementation's
  /// log-likelihood is compared against the analytical formula that sums over all possible
  /// root and internal node states. This provides complete coverage of the state-space
  /// for this topology, catching any bugs in how specific state transitions are handled
  /// (e.g. transversions vs transitions under symmetric models).
  #[test]
  fn test_three_taxon_all_combinations() -> Result<(), Report> {
    // Verify likelihood computation for all 64 possible single-position triplets
    let gtr = jc69(JC69Params::default())?;

    let t_a = 0.1;
    let t_b = 0.2;
    let t_ab = 0.15;
    let t_c = 0.25;

    let states = ['A', 'C', 'G', 'T'];
    for &state_a in &states {
      for &state_b in &states {
        for &state_c in &states {
          let obs_a = state_index(state_a);
          let obs_b = state_index(state_b);
          let obs_c = state_index(state_c);

          let expected_lh = analytical_three_taxon_likelihood(&gtr, obs_a, obs_b, obs_c, t_a, t_b, t_ab, t_c);
          let expected_log_lh = expected_lh.ln();

          let newick = "((A:0.1,B:0.2)AB:0.15,C:0.25)root;";
          let aln = format!(">A\n{state_a}\n>B\n{state_b}\n>C\n{state_c}\n");
          let actual_log_lh = run_dense_marginal_get_log_lh(newick, &aln, gtr.clone())?;

          assert!(
            (expected_log_lh - actual_log_lh).abs() < LIKELIHOOD_EPSILON,
            "Mismatch for triplet ({state_a},{state_b},{state_c}): expected {expected_log_lh}, got {actual_log_lh}"
          );
        }
      }
    }
    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::alphabet::alphabet::AlphabetName;
  use crate::commands::ancestral::fitch::get_common_length;
  use crate::commands::ancestral::marginal::{ancestral_reconstruction_marginal, initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::pretty_assert_ulps_eq;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::payload::ancestral::GraphAncestral;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use ndarray::prelude::*;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::{Arc, LazyLock};
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::json::{JsonPretty, json_write_str};
  use treetime_io::nwk::nwk_read_str;

  /// Assert that every row of a dense probability matrix is a valid distribution.
  ///
  /// Checks three conditions per row:
  /// 1. Row sums to 1.0 within the given ULP tolerance (normalization).
  /// 2. All values are finite (no NaN or infinity from numerical overflow).
  /// 3. No significantly negative values (probabilities must be non-negative,
  ///    with a small tolerance of -1e-15 for floating-point rounding).
  fn assert_dense_rows_normalized(dis: &Array2<f64>, max_ulps: u32) {
    for (row_idx, row) in dis.rows().into_iter().enumerate() {
      let sum: f64 = row.sum();
      assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
      assert!(
        sum.is_finite(),
        "Row {row_idx} sum={sum} is not normalized to 1.0 within max_ulps={max_ulps}"
      );
      for (col_idx, &val) in row.iter().enumerate() {
        assert!(
          val.is_finite(),
          "Row {row_idx}, col {col_idx} has non-finite value: {val}"
        );
        assert!(val >= -1e-15, "Row {row_idx}, col {col_idx} has negative value: {val}");
      }
    }
  }

  /// Construct a GTR model with non-uniform equilibrium frequencies
  /// pi = [0.2, 0.3, 0.15, 0.35].
  ///
  /// The asymmetric stationary distribution (A=20%, C=30%, G=15%, T=35%) ensures
  /// the rate matrix Q = S * diag(pi) is not proportional to JC69, exercising
  /// the full GTR eigendecomposition. Uses default exchangeability (W = None),
  /// which produces equal off-diagonal rates before weighting by pi.
  fn make_nonuniform_gtr(treat_gap_as_unknown: bool) -> Result<GTR, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
    GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.2, 0.3, 0.15, 0.35],
    })
  }

  /// Run dense marginal ancestral reconstruction and return the log-likelihood
  /// and partition array.
  ///
  /// Initializes a `PartitionMarginalDense` with the given GTR model and
  /// alignment, then runs `initialize_marginal` which performs Felsenstein's
  /// pruning (upward pass) and the downward pass to compute marginal posterior
  /// distributions P(s|data) at every node and position.
  fn run_dense_marginal(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
    gtr: GTR,
    treat_gap_as_unknown: bool,
  ) -> Result<(f64, [Arc<RwLock<PartitionMarginalDense>>; 1]), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    let log_lh = initialize_marginal(graph, &partitions, aln)?;
    Ok((log_lh, partitions))
  }

  /// Parse a Newick tree string and compute the dense marginal log-likelihood.
  ///
  /// Convenience wrapper for root-invariance tests that need to evaluate the
  /// same alignment under different tree topologies (re-rootings). Returns only
  /// the scalar log-likelihood, discarding the partition data.
  fn run_dense_lh_for_newick(
    newick: &str,
    aln: &[FastaRecord],
    gtr: GTR,
    treat_gap_as_unknown: bool,
  ) -> Result<f64, Report> {
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let (log_lh, _) = run_dense_marginal(&graph, aln, gtr, treat_gap_as_unknown)?;
    Ok(log_lh)
  }

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  /// Verify that dense marginal ancestral reconstruction infers the expected
  /// ancestral sequences at internal nodes.
  ///
  /// Uses a 4-taxon tree with known leaf and internal-node sequences (provided
  /// as alignment input for leaves and as expected output for internal nodes).
  /// The alignment includes gaps, ambiguity codes (N, R), and all four
  /// nucleotides to exercise diverse code paths in the marginal reconstruction.
  ///
  /// The JC69 model with `treat_gap_as_unknown=true` treats gaps as fully
  /// ambiguous (equivalent to N). The test compares inferred ancestral sequences
  /// at root, AB, and CD against hand-verified expected sequences using
  /// `ancestral_reconstruction_marginal`, which picks the maximum a posteriori
  /// state argmax_s P(s|data) at each position.
  #[test]
  fn test_ancestral_reconstruction_marginal_dense() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let aln = read_many_fasta_str(
      indoc! {r#"
      >root
      TCAGCCATGTATTG--
      >AB
      ACATCCCTGTA-TG--
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >CD
      CCGGCCATGTATTG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
    "#},
      &*NUC_ALPHABET,
    )?;

    let expected = read_many_fasta_str(
      indoc! {r#"
      >root
      TCGGCGCTGTATTGAC
      >AB
      ACATCGCTGTA-TGAC
      >CD
      TCGGCGGTGTATTG--
    "#},
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let treat_gap_as_unknown = true;
    let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: true,
      ..JC69Params::default()
    })?;

    let partitions_marginal_dense = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    initialize_marginal(&graph, &partitions_marginal_dense, &aln)?;

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_marginal(&graph, false, &partitions_marginal_dense, |node, seq| {
      actual.insert(node.name.clone(), seq.to_string());
      Ok(())
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  /// Verify that all marginal posterior distributions are properly normalized
  /// after dense reconstruction.
  ///
  /// For every node and every edge message, checks that each row of the
  /// probability matrix sums to 1.0 (within 4 ULPs), contains only finite
  /// values, and has no significantly negative entries. This validates the
  /// probability axiom: sum_s P(s|data) = 1 for all positions.
  ///
  /// Also verifies the total log-likelihood against a known expected value
  /// for this specific tree/alignment/model combination.
  #[test]
  fn test_marginal_dense_probability_normalization() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >root
      TCAGCCATGTATTG--
      >AB
      ACATCCCTGTA-TG--
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >CD
      CCGGCCATGTATTG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
    "#},
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: true,
      ..JC69Params::default()
    })?;

    let (log_lh, partitions) = run_dense_marginal(&graph, &aln, gtr, true)?;

    // Verify log-likelihood is in expected range (tree with 7 nodes, 16 sites)
    pretty_assert_ulps_eq!(-57.712498930787206, log_lh, epsilon = 1e-6);

    let partition = partitions[0].read_arc();
    let max_ulps = 4;

    for node_data in partition.nodes.values() {
      if !node_data.profile.dis.is_empty() {
        assert_dense_rows_normalized(&node_data.profile.dis, max_ulps);
      }
    }

    for edge_data in partition.edges.values() {
      if !edge_data.msg_to_child.dis.is_empty() {
        assert_dense_rows_normalized(&edge_data.msg_to_child.dis, max_ulps);
      }
    }

    Ok(())
  }

  /// Verify that `update_marginal` is idempotent: running it twice produces
  /// the same log-likelihood.
  ///
  /// Marginal ancestral reconstruction via Felsenstein's algorithm converges in
  /// a single upward + downward pass on a tree (no iterative refinement needed).
  /// Calling `update_marginal` a second time must produce identical results
  /// because the message-passing equations have a unique fixed point on trees.
  ///
  /// This test first runs `initialize_marginal` (which includes one
  /// `update_marginal` call), then calls `update_marginal` twice more and
  /// verifies both calls return the same log-likelihood.
  #[test]
  fn test_marginal_dense_update_is_idempotent() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >root
      TCAGCCATGTATTG--
      >AB
      ACATCCCTGTA-TG--
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >CD
      CCGGCCATGTATTG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
    "#},
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: true,
      ..JC69Params::default()
    })?;

    let (_, partitions) = run_dense_marginal(&graph, &aln, gtr, true)?;

    let log_lh_first = update_marginal(&graph, &partitions)?;
    let log_lh_second = update_marginal(&graph, &partitions)?;

    // Verify log-likelihood value matches expected (same tree/alignment as normalization test)
    pretty_assert_ulps_eq!(-57.712498930787206, log_lh_first, epsilon = 1e-6);

    // Verify idempotency
    assert_ulps_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);

    Ok(())
  }

  /// Verify that the total log-likelihood is invariant under root placement
  /// for a time-reversible GTR model (Felsenstein's pulley principle).
  ///
  /// For time-reversible models where Q = S * diag(pi), detailed balance holds:
  /// pi[s] * P(s'|s,t) = pi[s'] * P(s|s',t). This implies the likelihood is
  /// independent of where the root is placed on the unrooted tree, because the
  /// root can be "pulled" along any edge without changing the product of
  /// transition probabilities (the "pulley principle").
  ///
  /// Tests three different rootings of the same 4-taxon unrooted tree topology
  /// {A,B,C,D} with a non-uniform GTR model (pi = [0.2, 0.3, 0.15, 0.35]):
  /// - tree1: root between (A,B) and (C,D) clades
  /// - tree2: root at the A-B ancestor (AB node becomes root)
  /// - tree3: root at the C-D ancestor (CD node becomes root)
  ///
  /// All three must produce the same log-likelihood.
  #[test]
  fn test_marginal_dense_log_lh_root_invariance_reversible_model() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGTACGTACGTACGT
      >B
      ACGTACGTACGTACGA
      >C
      ACGTACGTACGTACGG
      >D
      ACGTACGTACGTACGC
    "#},
      &*NUC_ALPHABET,
    )?;

    let treat_gap_as_unknown = true;
    let gtr1 = make_nonuniform_gtr(treat_gap_as_unknown)?;
    let gtr2 = make_nonuniform_gtr(treat_gap_as_unknown)?;
    let gtr3 = make_nonuniform_gtr(treat_gap_as_unknown)?;

    let tree1 = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let tree2 = "(A:0.1,B:0.2,(C:0.2,D:0.12)CD:0.15)AB:0.01;";
    let tree3 = "((A:0.1,B:0.2)AB:0.15,C:0.2,D:0.12)CD:0.01;";

    let log_lh1 = run_dense_lh_for_newick(tree1, &aln, gtr1, treat_gap_as_unknown)?;
    let log_lh2 = run_dense_lh_for_newick(tree2, &aln, gtr2, treat_gap_as_unknown)?;
    let log_lh3 = run_dense_lh_for_newick(tree3, &aln, gtr3, treat_gap_as_unknown)?;

    let epsilon = 1e-6;
    assert_ulps_eq!(log_lh1, log_lh2, epsilon = epsilon);
    assert_ulps_eq!(log_lh1, log_lh3, epsilon = epsilon);
    assert_ulps_eq!(log_lh2, log_lh3, epsilon = epsilon);

    Ok(())
  }

  /// Verify the probability axiom: the sum of likelihoods over all possible
  /// leaf state combinations equals 1.0.
  ///
  /// For a single alignment position on a 3-taxon tree, the total likelihood
  /// summed over all 4^3 = 64 possible nucleotide triplets (A, B, C) must equal
  /// 1.0 by the law of total probability:
  ///
  ///   sum_{s_A, s_B, s_C} L(s_A, s_B, s_C) = 1
  ///
  /// where L(s_A, s_B, s_C) = sum_{s_root} pi[s_root] * product over children c
  /// of (sum_{s_internal} P(s_internal|s_root, t_c) * P(s_leaf|s_internal, t_leaf)).
  ///
  /// Uses a skewed GTR model (pi = [0.9, 0.06, 0.02, 0.02]) with relatively
  /// long branches to stress-test numerical precision. Each of the 64 single-site
  /// alignments is evaluated independently, and the exponentiated log-likelihoods
  /// are summed.
  #[test]
  fn test_total_likelihood_marginal_dense_all_triplets() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let alphabet = Alphabet::default();
    let mut total_lh = 0.0;

    // Use asymmetric GTR model with non-uniform stationary distribution
    let mu = 1.0;
    let pi = array![0.9, 0.06, 0.02, 0.02];
    let gtr = GTR::new(GTRParams {
      alphabet: alphabet.clone(),
      W: None,
      pi,
      mu,
    })?;

    let graph: GraphAncestral = nwk_read_str("((A:0.6,B:0.3):0.1,C:0.2)root:0.001;")?;
    // Generate all possible triplets (4^3 = 64 combinations)
    let states = ['A', 'C', 'G', 'T'];
    for &state_a in &states {
      for &state_b in &states {
        for &state_c in &states {
          // Create alignment with single position containing this triplet
          let aln = read_many_fasta_str(format!(">A\n{state_a}\n>B\n{state_b}\n>C\n{state_c}\n"), &*NUC_ALPHABET)?;

          let partitions_marginal_dense = [Arc::new(RwLock::new(PartitionMarginalDense {
            index: 0,
            gtr: gtr.clone(),
            alphabet: alphabet.clone(),
            length: get_common_length(&aln)?,
            nodes: btreemap! {},
            edges: btreemap! {},
          }))];

          initialize_marginal(&graph, &partitions_marginal_dense, &aln)?;

          let log_lh = update_marginal(&graph, &partitions_marginal_dense)?;
          total_lh += log_lh.exp();
        }
      }
    }

    // since we test all possible triplets, the total likelihood should be 1
    pretty_assert_ulps_eq!(1.0, total_lh, epsilon = 1e-6);
    Ok(())
  }
}

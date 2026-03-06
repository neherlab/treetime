#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::test_utils::find_node_key_by_name;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use ndarray::array;
  use parking_lot::RwLock;
  use std::sync::{Arc, LazyLock};
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::make_report;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  /// 4-taxon balanced tree with moderate branch lengths (0.05-0.2 subs/site).
  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  /// 16bp gap-free nucleotide alignment for 4 taxa (A, B, C, D).
  ///
  /// All sequences share a 15bp prefix `ACGTACGTACGTACG` and differ only at
  /// position 16: A=T, B=A, C=G, D=C. This creates exactly one variable site,
  /// isolating the effect of the substitution model on ancestral state inference.
  fn gap_free_alignment() -> Result<Vec<FastaRecord>, Report> {
    read_many_fasta_str(
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
    )
  }

  /// Run dense marginal ancestral reconstruction on the given tree and alignment.
  ///
  /// Dense representation stores full N x K probability matrices where N is the
  /// alignment length and K is the alphabet size (4 for nucleotides). Computes
  /// Felsenstein's pruning algorithm (upward pass) followed by the downward pass
  /// to produce marginal posterior probabilities at every node.
  ///
  /// Returns the total log-likelihood and the populated partition.
  fn run_dense_marginal(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
    gtr: GTR,
  ) -> Result<(f64, Arc<RwLock<PartitionMarginalDense>>), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    let partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));
    let partitions = [Arc::clone(&partition)];

    let log_lh = initialize_marginal(graph, &partitions, aln)?;
    Ok((log_lh, partition))
  }

  /// Run sparse marginal ancestral reconstruction on the given tree and alignment.
  ///
  /// Sparse representation stores only "variable" positions where the marginal
  /// posterior differs from the fixed-character distribution. Positions where all
  /// descendants agree on a single state are compressed away, storing only the
  /// Fitch-compressed reference and per-character fixed distributions.
  ///
  /// The sparse path first compresses sequences via Fitch parsimony, then runs
  /// `update_marginal` which performs the upward and downward message-passing
  /// passes on the compressed representation.
  ///
  /// Returns the total log-likelihood and the populated partition.
  fn run_sparse_marginal(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
    gtr: GTR,
  ) -> Result<(f64, Arc<RwLock<PartitionMarginalSparse>>), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    let partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));
    let partitions = [Arc::clone(&partition)];

    compress_sequences(graph, &partitions, aln)?;
    let log_lh = update_marginal(graph, &partitions)?;
    Ok((log_lh, partition))
  }

  /// Verify that dense and sparse marginal implementations produce the same
  /// total log-likelihood on a gap-free alignment with JC69 model.
  ///
  /// Both representations compute the same quantity - the Felsenstein likelihood
  /// L = product over sites of sum over root states of pi[s] * L_root(s) - but via
  /// different code paths. Dense operates on full NxK matrices; sparse compresses
  /// invariant sites and operates only on variable positions. The log-likelihoods
  /// must agree to floating-point precision.
  #[test]
  fn test_marginal_dense_sparse_log_lh_consistency_gap_free() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let gtr_dense = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
      ..JC69Params::default()
    })?;
    let gtr_sparse = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
      ..JC69Params::default()
    })?;

    let (log_lh_dense, _) = run_dense_marginal(&graph, &aln, gtr_dense)?;
    let (log_lh_sparse, _) = run_sparse_marginal(&graph, &aln, gtr_sparse)?;

    assert_ulps_eq!(log_lh_dense, log_lh_sparse, epsilon = 1e-10);

    Ok(())
  }

  /// Verify that sparse variable-position distributions match the corresponding
  /// rows in the dense marginal probability matrix.
  ///
  /// At internal nodes (root, AB), for each position the sparse representation
  /// marks as "variable", the probability vector must equal the corresponding row
  /// of the dense NxK matrix. This checks profile-level consistency beyond the
  /// aggregate log-likelihood: the per-position marginal posteriors
  /// P(s|data) = pi[s] * L_node(s) / Z must agree between representations.
  #[test]
  fn test_marginal_sparse_varpos_matches_dense_profile_gap_free() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let gtr_dense = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
      ..JC69Params::default()
    })?;
    let gtr_sparse = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
      ..JC69Params::default()
    })?;

    let (_, dense_partition) = run_dense_marginal(&graph, &aln, gtr_dense)?;
    let (_, sparse_partition) = run_sparse_marginal(&graph, &aln, gtr_sparse)?;

    let root_key = find_node_key_by_name(&graph, "root").ok_or_else(|| make_report!("Root node not found"))?;
    let ab_key = find_node_key_by_name(&graph, "AB").ok_or_else(|| make_report!("AB node not found"))?;

    let dense = dense_partition.read_arc();
    let sparse = sparse_partition.read_arc();

    for node_key in [root_key, ab_key] {
      let dense_node = &dense.nodes[&node_key];
      let sparse_node = &sparse.nodes[&node_key];

      for (&pos, var_pos) in &sparse_node.profile.variable {
        let dense_row = dense_node.profile.dis.row(pos);
        let sparse_dis = &var_pos.dis;

        assert_eq!(
          dense_row.len(),
          sparse_dis.len(),
          "Profile dimensions mismatch at position {pos}"
        );

        for (idx, (&dense_val, &sparse_val)) in dense_row.iter().zip(sparse_dis.iter()).enumerate() {
          assert_ulps_eq!(dense_val, sparse_val, epsilon = 1e-6);
        }
      }
    }

    Ok(())
  }

  /// Verify dense/sparse consistency and probability normalization in the
  /// presence of IUPAC ambiguity codes (N = any nucleotide, R = A|G).
  ///
  /// Ambiguous characters assign equal probability to the states they represent
  /// (e.g., N -> uniform over {A,C,G,T}, R -> uniform over {A,G}). This test
  /// checks three properties:
  ///
  /// 1. Log-likelihood agreement: dense and sparse produce the same total
  ///    log-likelihood despite different handling of ambiguous tip states.
  /// 2. Dense normalization: every row of every node's marginal probability
  ///    matrix sums to 1.0 (valid probability distribution).
  /// 3. Sparse normalization: both variable-position distributions and
  ///    fixed-character distributions sum to 1.0, and all log-likelihoods
  ///    are finite.
  #[test]
  fn test_marginal_dense_sparse_ambiguous_character_expectations_documented() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGTACGTACGTACNT
      >B
      ACGTACGTACGTACRA
      >C
      ACGTACGTACGTACGG
      >D
      ACGTACGTACGTACGC
    "#},
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let gtr_dense = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
      ..JC69Params::default()
    })?;
    let gtr_sparse = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
      ..JC69Params::default()
    })?;

    let (log_lh_dense, dense_partition) = run_dense_marginal(&graph, &aln, gtr_dense)?;
    let (log_lh_sparse, sparse_partition) = run_sparse_marginal(&graph, &aln, gtr_sparse)?;

    assert_ulps_eq!(log_lh_dense, log_lh_sparse, epsilon = 1e-10);

    let dense = dense_partition.read_arc();
    let sparse = sparse_partition.read_arc();

    for node_data in dense.nodes.values() {
      if !node_data.profile.dis.is_empty() {
        for row in node_data.profile.dis.rows() {
          let sum: f64 = row.sum();
          assert!(sum.is_finite(), "Dense node profile row sum is not finite: {sum}");
          assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
          assert!(sum.is_finite(), "Dense node profile row not normalized: sum={sum}");
        }
      }
    }

    for node_data in sparse.nodes.values() {
      assert!(
        node_data.profile.log_lh.is_finite(),
        "Sparse node profile log_lh is not finite: {}",
        node_data.profile.log_lh
      );

      for (pos, var_pos) in &node_data.profile.variable {
        let sum: f64 = var_pos.dis.sum();
        assert!(
          sum.is_finite(),
          "Sparse variable position {pos} sum is not finite: {sum}"
        );
        assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
        assert!(
          sum.is_finite(),
          "Sparse variable position {pos} not normalized: sum={sum}"
        );
      }

      for (char_key, fixed_dis) in &node_data.profile.fixed {
        let sum: f64 = fixed_dis.sum();
        assert!(
          sum.is_finite(),
          "Sparse fixed distribution for char {char_key:?} sum is not finite: {sum}"
        );
        assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
        assert!(
          sum.is_finite(),
          "Sparse fixed distribution for char {char_key:?} not normalized: sum={sum}"
        );
      }
    }

    Ok(())
  }

  /// Verify marginal probability normalization under a highly skewed GTR model.
  ///
  /// Ported from Python v0 `test_treetime.py:137-155`. Uses equilibrium
  /// frequencies pi = [0.9, 0.06, 0.02, 0.02] with uniform exchangeability
  /// (W = ones with zero diagonal), creating extreme asymmetry that stresses
  /// numerical stability of the Felsenstein pruning algorithm.
  ///
  /// The 3-taxon tree with long branches (up to 0.601 subs/site) and the 64bp
  /// alignment covering all 4^3 = 64 state combinations ensure every code path
  /// through the transition probability matrix P(t) = exp(Q*t) is exercised.
  ///
  /// Asserts that every row of the marginal posterior matrix at every node sums
  /// to 1.0: sum_s P(s|data) = 1 for all positions and all nodes.
  #[test]
  fn test_marginal_likelihoods_sum_to_one_skewed_gtr() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;

    // Highly skewed equilibrium frequencies from Python test
    let pi = array![0.9, 0.06, 0.02, 0.02];

    // W = ones(4,4) with zero diagonal (symmetric rate matrix)
    let W = array![
      [0.0, 1.0, 1.0, 1.0],
      [1.0, 0.0, 1.0, 1.0],
      [1.0, 1.0, 0.0, 1.0],
      [1.0, 1.0, 1.0, 0.0],
    ];

    let gtr = GTR::new(GTRParams {
      alphabet: alphabet.clone(),
      mu: 1.0,
      W: Some(W),
      pi,
    })?;

    // Tree from Python: ((A:0.60100000009,B:0.3010000009):0.1,C:0.2):0.001
    let tree_newick = "((A:0.601,B:0.301):0.1,C:0.2):0.001;";
    let graph: GraphAncestral = nwk_read_str(tree_newick)?;

    // Alignment from Python test (64bp each)
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT
        >B
        AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
        >C
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
      "#},
      &alphabet,
    )?;

    let partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));
    let partitions = [Arc::clone(&partition)];

    initialize_marginal(&graph, &partitions, &aln)?;

    // Verify all marginal profile rows sum to 1.0
    let partition = partition.read_arc();
    for (node_key, node_data) in &partition.nodes {
      if node_data.profile.dis.is_empty() {
        continue;
      }
      for (pos, row) in node_data.profile.dis.rows().into_iter().enumerate() {
        let sum: f64 = row.sum();
        assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
        assert!(
          sum.is_finite(),
          "Node {node_key:?} position {pos}: marginal likelihood sum {sum} != 1.0"
        );
      }
    }

    Ok(())
  }
}

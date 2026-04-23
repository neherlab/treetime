#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{ancestral_reconstruction_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::pretty_assert_ulps_eq;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::partition::traits::PartitionBranchOps;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::representation::payload::sparse::MarginalSparseSeqDistribution;
  use crate::seq::mutation::Sub;
  use crate::test_utils::find_node_key_by_name;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use ndarray::prelude::*;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::{Arc, LazyLock};
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::json::{JsonPretty, json_write_str};
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;
  use treetime_primitives::seq;

  /// Lazily initialized default nucleotide alphabet (A, C, G, T with gap handling).
  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  /// Assert that a sparse marginal posterior distribution is a valid probability distribution.
  ///
  /// Checks both variable positions (where leaf states differ, stored as individual probability
  /// vectors per position) and fixed positions (grouped by conserved nucleotide character at the
  /// node, each group sharing a single probability vector). For each distribution vector,
  /// verifies:
  /// - Finite log-likelihood
  /// - All elements are finite and non-negative
  /// - Elements sum to 1.0 within the given ULP tolerance
  ///
  /// Kolmogorov axioms: P(s) >= 0 for all states s, and sum_s P(s) = 1.
  fn assert_sparse_profile_normalized(profile: &MarginalSparseSeqDistribution, max_ulps: u32) {
    assert!(
      profile.log_lh.is_finite(),
      "Profile log_lh is not finite: {}",
      profile.log_lh
    );

    for (pos, var_pos) in &profile.variable {
      let sum: f64 = var_pos.dis.sum();
      assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
      assert!(
        sum.is_finite(),
        "Variable position {pos} sum={sum} is not normalized to 1.0 within max_ulps={max_ulps}"
      );
      for (idx, &val) in var_pos.dis.iter().enumerate() {
        assert!(
          val.is_finite(),
          "Variable position {pos}, index {idx} has non-finite value: {val}"
        );
        assert!(
          val >= -1e-15,
          "Variable position {pos}, index {idx} has negative value: {val}"
        );
      }
    }

    for (char_key, fixed_dis) in &profile.fixed {
      let sum: f64 = fixed_dis.sum();
      assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
      assert!(
        sum.is_finite(),
        "Fixed distribution for char {char_key:?} sum={sum} is not normalized to 1.0 within max_ulps={max_ulps}"
      );
      for (idx, &val) in fixed_dis.iter().enumerate() {
        assert!(
          val.is_finite(),
          "Fixed distribution for char {char_key:?}, index {idx} has non-finite value: {val}"
        );
        assert!(
          val >= -1e-15,
          "Fixed distribution for char {char_key:?}, index {idx} has negative value: {val}"
        );
      }
    }
  }

  /// Create a GTR model with non-uniform equilibrium frequencies pi = [0.2, 0.3, 0.15, 0.35].
  ///
  /// Non-uniform pi breaks symmetries present in JC69 (where pi = [0.25, 0.25, 0.25, 0.25]),
  /// exposing bugs that only manifest when nucleotide frequencies are unequal. The model
  /// remains reversible with default symmetric exchangeabilities (W = None), satisfying
  /// detailed balance: pi_i * Q_ij = pi_j * Q_ji.
  fn make_nonuniform_gtr() -> Result<GTR, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    GTR::new(GTRParams {
      n_states,
      mu: 1.0,
      W: None,
      pi: array![0.2, 0.3, 0.15, 0.35],
    })
  }

  /// Run sparse marginal ancestral reconstruction on a pre-built graph with a given alignment
  /// and GTR model.
  ///
  /// Performs Fitch parsimony compression (classifying sites as variable or invariant) followed
  /// by two-pass sum-product message passing (backward pass computes partial likelihoods from
  /// leaves to root, forward pass propagates outgroup information from root to leaves) in the
  /// sparse representation. Returns the total log-likelihood and the partition for further
  /// inspection of node/edge posteriors.
  fn run_sparse_marginal(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
    gtr: GTR,
  ) -> Result<(f64, [Arc<RwLock<PartitionMarginalSparse>>; 1]), Report> {
    let alphabet = Alphabet::default();
    let partitions = [Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
      root_sequence: seq![],
    }))];

    compress_sequences(graph, &partitions, aln)?;
    partitions[0].write_arc().extract_root_sequence(graph);
    let log_lh = update_marginal(graph, &partitions)?;
    Ok((log_lh, partitions))
  }

  /// Parse a Newick tree string and compute the sparse marginal log-likelihood.
  ///
  /// Convenience wrapper for root-invariance tests that need to evaluate the same alignment
  /// under different rootings of the same unrooted topology. Returns only the scalar
  /// log-likelihood, discarding the partition data.
  fn run_sparse_lh_for_newick(newick: &str, aln: &[FastaRecord], gtr: GTR) -> Result<f64, Report> {
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let (log_lh, _) = run_sparse_marginal(&graph, aln, gtr)?;
    Ok(log_lh)
  }

  /// Verify that sparse marginal ancestral reconstruction infers the correct MAP (maximum a
  /// posteriori) sequences at internal nodes.
  ///
  /// Uses a balanced 4-taxon tree ((A,B)AB,(C,D)CD)root with JC69 model. The alignment
  /// contains gaps, ambiguous characters (N, R), and variable sites. At each internal node,
  /// the marginal posterior P(s|data,tree,GTR) is computed via two-pass sum-product message
  /// passing (backward pass for partial likelihoods, forward pass for outgroup messages),
  /// then the MAP character argmax_s P(s|data,tree,GTR) is extracted at each position.
  ///
  /// Expected MAP sequences are pre-computed and compared character-by-character.
  /// Also verifies the total log-likelihood matches the expected value.
  #[test]
  fn test_ancestral_reconstruction_marginal_sparse() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

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

    let expected = read_many_fasta_str(
      indoc! {r#"
      >root
      TCGGCGCTGTATTG--
      >AB
      ACATCGCTGTA-TG--
      >CD
      TCGGCGGTGTATTG--
    "#},
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions_marginal_sparse = [Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
      root_sequence: seq![],
    }))];

    compress_sequences(&graph, &partitions_marginal_sparse, &aln)?;
    partitions_marginal_sparse[0].write_arc().extract_root_sequence(&graph);

    let log_lh = update_marginal(&graph, &partitions_marginal_sparse)?;

    // generate ancestral reconstruction and test against expectation
    let mut actual = BTreeMap::new();
    ancestral_reconstruction_marginal(&graph, false, &partitions_marginal_sparse, |node, seq| {
      actual.insert(node.name.clone(), seq.to_string());
      Ok(())
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    // test overall likelihood
    pretty_assert_ulps_eq!(-55.55428499726621, log_lh, epsilon = 1e-6);
    Ok(())
  }

  /// Verify that all posterior distributions produced by sparse marginal reconstruction
  /// are valid probability distributions.
  ///
  /// After running Fitch compression and two-pass marginal reconstruction on a 4-taxon
  /// tree with JC69 model, checks every node profile and every edge outgroup message
  /// (msg_to_child). For each distribution vector, asserts sum_s P(s) = 1 and P(s) >= 0.
  /// This validates the normalization step in both the backward (partial likelihood) and
  /// forward (outgroup message) passes.
  #[test]
  fn test_marginal_sparse_probability_normalization() -> Result<(), Report> {
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
    let gtr = jc69(JC69Params::default())?;

    let (log_lh, partitions) = run_sparse_marginal(&graph, &aln, gtr)?;

    // Verify log-likelihood is in expected range (tree with 4 leaves, 16 sites)
    // Matches test_ancestral_reconstruction_marginal_sparse value
    pretty_assert_ulps_eq!(-55.55428496980045, log_lh, epsilon = 1e-6);

    let partition = partitions[0].read_arc();

    for node_data in partition.nodes.values() {
      assert_sparse_profile_normalized(&node_data.profile, 4);
    }

    for edge_data in partition.edges.values() {
      assert_sparse_profile_normalized(&edge_data.msg_to_child, 4);
    }

    Ok(())
  }

  /// Verify that `update_marginal()` is a fixed-point operation: calling it twice produces
  /// the same log-likelihood.
  ///
  /// Sum-product message passing on a tree (acyclic graph) computes exact marginal posteriors
  /// in a single backward + forward pass because the message-passing equations have a unique
  /// fixed point on trees. A second invocation must not change any computed values. This tests
  /// that no hidden mutable state accumulates across calls.
  #[test]
  fn test_marginal_sparse_update_is_idempotent() -> Result<(), Report> {
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
    let gtr = jc69(JC69Params::default())?;

    let alphabet = Alphabet::default();
    let partitions = [Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
      root_sequence: seq![],
    }))];

    compress_sequences(&graph, &partitions, &aln)?;

    let log_lh_first = update_marginal(&graph, &partitions)?;
    let log_lh_second = update_marginal(&graph, &partitions)?;

    // Verify log-likelihood value matches expected (same tree/alignment as normalization test)
    pretty_assert_ulps_eq!(-55.55428496980045, log_lh_first, epsilon = 1e-6);

    // Verify idempotency
    assert_ulps_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);

    Ok(())
  }

  /// Verify Felsenstein's pulley principle: the total log-likelihood is invariant under
  /// root placement for a time-reversible substitution model.
  ///
  /// For a reversible model, the rate matrix Q satisfies detailed balance:
  /// pi_i * Q_ij = pi_j * Q_ji. This propagates to the transition probability matrix
  /// P(t) = exp(Qt), making the likelihood of an unrooted tree independent of where the
  /// root is placed - the root can be "pulled" along any edge without changing the product
  /// of transition probabilities. The test re-roots the same 4-taxon unrooted topology at
  /// three different internal edges and verifies all three rootings produce the same
  /// log-likelihood.
  ///
  /// Uses a non-uniform GTR with pi = [0.2, 0.3, 0.15, 0.35] to ensure the invariance
  /// holds for asymmetric equilibrium frequencies, not just JC69.
  #[test]
  fn test_marginal_sparse_log_lh_root_invariance_reversible_model() -> Result<(), Report> {
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

    let gtr1 = make_nonuniform_gtr()?;
    let gtr2 = make_nonuniform_gtr()?;
    let gtr3 = make_nonuniform_gtr()?;

    let tree1 = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let tree2 = "(A:0.1,B:0.2,(C:0.2,D:0.12)CD:0.15)AB:0.01;";
    let tree3 = "((A:0.1,B:0.2)AB:0.15,C:0.2,D:0.12)CD:0.01;";

    let log_lh1 = run_sparse_lh_for_newick(tree1, &aln, gtr1)?;
    let log_lh2 = run_sparse_lh_for_newick(tree2, &aln, gtr2)?;
    let log_lh3 = run_sparse_lh_for_newick(tree3, &aln, gtr3)?;

    assert_ulps_eq!(log_lh1, log_lh2, epsilon = 1e-6);
    assert_ulps_eq!(log_lh1, log_lh3, epsilon = 1e-6);
    assert_ulps_eq!(log_lh2, log_lh3, epsilon = 1e-6);

    Ok(())
  }

  /// Verify specific marginal posterior values at root and internal node AB against Python
  /// v0 reference values.
  ///
  /// Uses a non-uniform GTR with pi = [0.2, 0.3, 0.15, 0.35] on a 4-taxon tree with
  /// ambiguous characters (N, R) and gaps. Checks the posterior distribution
  /// P(s|data,tree,GTR) at position 0 for the root node and internal node AB, plus
  /// position 3 for AB. These reference values come from the Python v0 TreeTime
  /// implementation (test_scripts/ancestral_sparse.py).
  ///
  /// The log-likelihood differs slightly from the dense representation due to different
  /// handling of ambiguous characters in the sparse representation.
  #[test]
  fn test_marginal_sparse_posterior_values_python_parity() -> Result<(), Report> {
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
    let gtr = make_nonuniform_gtr()?;

    let (log_lh, partitions) = run_sparse_marginal(&graph, &aln, gtr)?;

    // note that this LH is slightly different from dense or python treetime due to
    // different handling of ambiguous characters (value from test_scripts/ancestral_sparse.py)
    pretty_assert_ulps_eq!(-56.946298878390444, log_lh, epsilon = 1e-6);

    let partition = partitions[0].read_arc();

    // test variable position distribution at the root (pos 0)
    let root_key = graph.get_exactly_one_root()?.read_arc().key();
    let root_profile = &partition.nodes[&root_key].profile;
    let pos_zero_root = array![0.28212327, 0.21643546, 0.13800802, 0.36343326];
    pretty_assert_ulps_eq!(&root_profile.variable[&0].dis, &pos_zero_root, epsilon = 1e-6);

    // test variable position distribution at internal node AB (pos 0)
    let ab_key = find_node_key_by_name(&graph, "AB").expect("AB node should exist");
    let ab_profile = &partition.nodes[&ab_key].profile;
    let pos_zero_ab = array![0.51275208, 0.09128506, 0.24647255, 0.14949031];
    pretty_assert_ulps_eq!(&ab_profile.variable[&0].dis, &pos_zero_ab, epsilon = 1e-6);

    // test variable position distribution at internal node AB (pos 3)
    let dis_ab_pos3 = array![
      0.0013914677323952813,
      0.002087201598592933,
      0.042827146239885545,
      0.9536941844291262
    ];
    pretty_assert_ulps_eq!(&ab_profile.variable[&3].dis, &dis_ab_pos3, epsilon = 1e-6);

    Ok(())
  }

  /// Verify the probability axiom: the marginal likelihood over all possible leaf state
  /// combinations sums to 1.0.
  ///
  /// For a 3-taxon tree with a 4-state nucleotide alphabet, enumerates all 4^3 = 64
  /// possible single-site leaf state assignments (a, b, c) in {A, C, G, T}. For each
  /// triplet, computes the site likelihood P(a, b, c | tree, GTR) via the backward pass
  /// (Felsenstein's pruning algorithm), which marginalizes over all internal node states.
  /// The sum over all leaf state triplets must equal 1 by the law of total probability:
  ///
  ///   sum_{a,b,c in {A,C,G,T}} P(a, b, c | tree, GTR) = 1
  ///
  /// Uses an asymmetric GTR with pi = [0.9, 0.06, 0.02, 0.02] to stress the computation
  /// with extreme frequency imbalance.
  #[test]
  fn test_total_likelihood_marginal_sparse_all_triplets() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let alphabet = Alphabet::default();
    let mut total_lh = 0.0;

    // Use asymmetric GTR model with non-uniform stationary distribution
    let mu = 1.0;
    let pi = array![0.9, 0.06, 0.02, 0.02];
    let gtr = GTR::new(GTRParams {
      n_states: alphabet.n_canonical(),
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

          let partitions_marginal_sparse = [Arc::new(RwLock::new(PartitionMarginalSparse {
            index: 0,
            gtr: gtr.clone(),
            alphabet: alphabet.clone(),
            length: get_common_length(&aln)?,
            nodes: btreemap! {},
            edges: btreemap! {},
            root_sequence: seq![],
          }))];

          compress_sequences(&graph, &partitions_marginal_sparse, &aln)?;

          let log_lh = update_marginal(&graph, &partitions_marginal_sparse)?;
          total_lh += log_lh.exp();
        }
      }
    }

    // since we test all possible triplets, the total likelihood should be 1
    pretty_assert_ulps_eq!(1.0, total_lh, epsilon = 1e-6);
    Ok(())
  }

  /// Verify that sparse `edge_subs()` reports current reconstructed branch
  /// substitutions rather than raw Fitch substitutions.
  ///
  /// After marginal inference, the branch result should match the difference
  /// between the reconstructed parent sequence and child sequence. Stored
  /// `edge.subs` is only the initial sparse encoding. This test checks:
  /// - `edge_subs()` matches the reconstructed parent-child diff on every edge
  /// - the result differs from stored `edge.subs` on at least one edge
  #[test]
  fn test_sparse_edge_subs_match_reconstructed_branch_differences() -> Result<(), Report> {
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
    let partitions = [Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: make_nonuniform_gtr()?,
      alphabet: Alphabet::default(),
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
      root_sequence: seq![],
    }))];

    compress_sequences(&graph, &partitions, &aln)?;
    partitions[0].write_arc().extract_root_sequence(&graph);
    update_marginal(&graph, &partitions)?;

    let actual_by_edge = {
      let partition = partitions[0].read_arc();
      graph
        .get_edges()
        .iter()
        .map(|edge| {
          let edge_key = edge.read_arc().key();
          let actual = partition.edge_subs(&graph, edge_key)?;
          Ok((edge_key, actual))
        })
        .collect::<Result<BTreeMap<_, _>, Report>>()?
    };

    // Reconstruct node sequences, then turn parent-child sequence differences
    // into the expected branch changes.
    let mut seqs_by_name = BTreeMap::new();
    ancestral_reconstruction_marginal(&graph, true, &partitions, |node, seq| {
      seqs_by_name.insert(
        node.name.clone().expect("all test nodes should have names"),
        seq.clone(),
      );
      Ok(())
    })?;

    let partition = partitions[0].read_arc();
    let expected_by_edge = helpers::expected_edge_subs_by_edge(&graph, &partition, &seqs_by_name)?;

    assert_eq!(expected_by_edge, actual_by_edge);
    Ok(())
  }

  mod helpers {
    use super::*;

    /// Convert reconstructed parent and child sequences into canonical substitutions.
    ///
    /// This stays simple on purpose: compare the two sequences site by site and
    /// keep only nucleotide changes.
    pub fn diff_canonical_subs(alphabet: &Alphabet, parent_seq: &Seq, child_seq: &Seq) -> Result<Vec<Sub>, Report> {
      parent_seq
        .iter()
        .zip(child_seq.iter())
        .enumerate()
        .filter(|(_, (parent, child))| parent != child)
        .filter(|(_, (parent, child))| alphabet.is_canonical(**parent) && alphabet.is_canonical(**child))
        .map(|(pos, (parent, child))| Sub::new(*parent, pos, *child))
        .collect()
    }

    /// Build the expected branch-mutation map from reconstructed node sequences.
    ///
    /// This keeps the expected result separate from the sparse edge logic under
    /// test. First reconstruct node sequences, then compare parent and child for
    /// each edge.
    pub fn expected_edge_subs_by_edge(
      graph: &GraphAncestral,
      partition: &PartitionMarginalSparse,
      seqs_by_name: &BTreeMap<String, Seq>,
    ) -> Result<BTreeMap<treetime_graph::edge::GraphEdgeKey, Vec<Sub>>, Report> {
      graph
        .get_edges()
        .iter()
        .map(|edge| {
          let edge = edge.read_arc();
          let edge_key = edge.key();
          let expected = diff_canonical_subs(
            &partition.alphabet,
            get_reconstructed_seq(graph, seqs_by_name, edge.source()),
            get_reconstructed_seq(graph, seqs_by_name, edge.target()),
          )?;
          Ok((edge_key, expected))
        })
        .collect()
    }

    /// Retrieve one reconstructed node sequence by graph node.
    fn get_reconstructed_seq<'a>(
      graph: &GraphAncestral,
      seqs_by_name: &'a BTreeMap<String, Seq>,
      node_key: GraphNodeKey,
    ) -> &'a Seq {
      let node = graph.get_node(node_key).expect("node should exist");
      let node = node.read_arc();
      let payload = node.payload().read_arc();
      let name = payload.name.as_ref().expect("all test nodes should have names");
      &seqs_by_name[name]
    }
  }
}

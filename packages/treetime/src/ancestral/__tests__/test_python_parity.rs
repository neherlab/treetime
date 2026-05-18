#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::{ancestral_reconstruction_marginal, initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::pretty_assert_ulps_eq;
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::seq::alignment::get_common_length;

  use crate::partition::payload::ancestral::GraphAncestral;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  
  use ndarray::array;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::path::PathBuf;
  use std::slice::from_ref;
  use std::sync::{Arc, LazyLock};
  use treetime_io::fasta::{read_many_fasta, read_many_fasta_str};
  use treetime_io::nwk::{nwk_read_file, nwk_read_str};

  use treetime_utils::make_report;

  /// Resolve the workspace root directory by walking up from the crate's manifest directory.
  ///
  /// Used to locate test data files (Newick trees, FASTA alignments) under `data/`.
  fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .map(PathBuf::from)
      .expect("project has workspace root")
  }

  /// Parity test: root sequence from two-pass marginal ancestral reconstruction under
  /// JC69 (Jukes-Cantor 1969) matches Python v0.
  ///
  /// The reconstruction uses Felsenstein's peeling algorithm (backward pass, leaves to
  /// root) followed by a forward pass (root to leaves) to compute marginal posterior
  /// distributions P(state | data) at each internal node. The root sequence is the
  /// argmax of the marginal posterior at each alignment position.
  ///
  /// Dataset: H3N2 influenza neuraminidase, 20 sequences.
  /// Model: JC69 with equal equilibrium frequencies pi_i = 0.25 and equal
  /// substitution rates.
  ///
  /// Python test reference: packages/legacy/treetime/test/test_treetime.py:116-135
  ///
  /// Python command:
  /// ```
  /// t = TreeAnc(gtr='Jukes-Cantor', tree=nwk, aln=fasta, rng_seed=1234)
  /// t.reconstruct_anc(method='ml', marginal=True)
  /// ```
  #[test]
  fn test_root_sequence_matches_python_h3n2_na_20() -> Result<(), Report> {
    let root = project_root();
    let tree_path = root.join("data/flu/h3n2/20/tree.nwk");
    let aln_path = root.join("data/flu/h3n2/20/aln.fasta.xz");

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta(&[aln_path], &alphabet)?;

    let graph: GraphAncestral = nwk_read_file(&tree_path)?;

    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense::new(0, gtr, alphabet, get_common_length(&aln)?)))];

    initialize_marginal(&graph, &partitions, &aln)?;

    let mut root_seq = String::new();
    ancestral_reconstruction_marginal(&graph, false, &partitions, |node, seq| {
      if node.name.as_deref() == Some("NODE_0000000") {
        root_seq = seq.to_string();
      }
      Ok(())
    })?;

    // Expected root sequence from Python v0 (packages/legacy/treetime/test/test_treetime.py:134)
    let expected = "ATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTTCCACAATATGCTTCTTCATGCAAATTGCCATCTTGATAACTACTGTAACATTGCATTTCAAGCAATATGAATTCAACTCCCCCCCAAACAACCAAGTGATGCTGTGTGAACCAACAATAATAGAAAGAAACATAACAGAGATAGTGTATCTGACCAACACCACCATAGAGAAGGAAATATGCCCCAAACCAGCAGAATACAGAAATTGGTCAAAACCGCAATGTGGCATTACAGGATTTGCACCTTTCTCTAAGGACAATTCGATTAGGCTTTCCGCTGGTGGGGACATCTGGGTGACAAGAGAACCTTATGTGTCATGCGATCCTGACAAGTGTTATCAATTTGCCCTTGGACAGGGAACAACACTAAACAACGTGCATTCAAATAACACAGTACGTGATAGGACCCCTTATCGGACTCTATTGATGAATGAGTTGGGTGTTCCTTTTCATCTGGGGACCAAGCAAGTGTGCATAGCATGGTCCAGCTCAAGTTGTCACGATGGAAAAGCATGGCTGCATGTTTGTATAACGGGGGATGATAAAAATGCAACTGCTAGCTTCATTTACAATGGGAGGCTTGTAGATAGTGTTGTTTCATGGTCCAAAGAAATTCTCAGGACCCAGGAGTCAGAATGCGTTTGTATCAATGGAACTTGTACAGTAGTAATGACTGATGGAAGTGCTTCAGGAAAAGCTGATACTAAAATACTATTCATTGAGGAGGGGAAAATCGTTCATACTAGCACATTGTCAGGAAGTGCTCAGCATGTCGAAGAGTGCTCTTGCTATCCTCGATATCCTGGTGTCAGATGTGTCTGCAGAGACAACTGGAAAGGCTCCAATCGGCCCATCGTAGATATAAACATAAAGGATCATAGCATTGTTTCCAGTTATGTGTGTTCAGGACTTGTTGGAGACACACCCAGAAAAAACGACAGCTCCAGCAGTAGCCATTGTTTGGATCCTAACAATGAAGAAGGTGGTCATGGAGTGAAAGGCTGGGCCTTTGATGATGGAAATGACGTGTGGATGGGAAGAACAATCAACGAGACGTCACGCTTAGGGTATGAAACCTTCAAAGTCATTGAAGGCTGGTCCAACCCTAAGTCCAAATTGCAGATAAATAGGCAAGTCATAGTTGACAGAGGTGATAGGTCCGGTTATTCTGGTATTTTCTCTGTTGAAGGCAAAAGCTGCATCAATCGGTGCTTTTATGTGGAGTTGATTAGGGGAAGAAAAGAGGAAACTGAAGTCTTGTGGACCTCAAACAGTATTGTTGTGTTTTGTGGCACCTCAGGTACATATGGAACAGGCTCATGGCCTGATGGGGCGGACCTCAATCTCATGCCTATA";

    assert_eq!(expected, root_seq, "Root sequence mismatch with Python v0");

    Ok(())
  }

  /// Lazily initialized default nucleotide alphabet (A, C, G, T with gaps treated as
  /// missing data, i.e. uniform likelihood across all states).
  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  /// Create the primary GTR model used in Python v0 parity tests, with non-uniform
  /// equilibrium frequencies pi = [0.2, 0.3, 0.15, 0.35].
  ///
  /// All exchangeability parameters are equal (W = None uses default uniform
  /// exchangeabilities). Combined with non-uniform pi, this is equivalent to the F81
  /// model (Felsenstein 1981): equal rates of exchange, unequal equilibrium frequencies.
  /// The non-uniform pi breaks the symmetry present in JC69, producing distinct posterior
  /// distributions for each state at internal nodes.
  ///
  /// Python reference: test_scripts/ancestral_sparse.py:235
  fn make_python_reference_gtr() -> Result<GTR, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    GTR::new(GTRParams {
      n_states,
      mu: 1.0,
      W: None,
      pi: array![0.2, 0.3, 0.15, 0.35],
    })
  }

  /// Create the secondary F81 model for multi-partition parity tests, with equilibrium
  /// frequencies pi = [0.4, 0.15, 0.12, 0.33].
  ///
  /// Used alongside `make_python_reference_gtr()` to verify that two partitions with
  /// different substitution models on the same tree compute independently. The different
  /// pi values produce distinct posterior distributions at each node, confirming that
  /// partition-level model parameters do not leak across partitions.
  ///
  /// Python reference: test_scripts/ancestral_sparse.py:255
  fn make_python_reference_gtr2() -> Result<GTR, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    GTR::new(GTRParams {
      n_states,
      mu: 1.0,
      W: None,
      pi: array![0.4, 0.15, 0.12, 0.33],
    })
  }

  /// Newick tree shared across Python parity tests: topologically balanced (symmetric
  /// bifurcating) 4-taxon tree with asymmetric branch lengths. The branch length
  /// asymmetry (A:0.1 vs B:0.2, C:0.2 vs D:0.12, AB:0.1 vs CD:0.05) produces unequal
  /// transition probability matrices on sister branches, testing that the algorithm
  /// correctly weights contributions from closer vs. more distant leaves.
  const PYTHON_TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  /// FASTA alignment shared across Python parity tests. Contains IUPAC ambiguity codes
  /// N (any base) and R (purine: A or G), plus gap characters (-). Ambiguous positions
  /// produce non-delta leaf likelihood vectors (uniform for N, partial for R), and gaps
  /// treated as missing data produce uniform vectors. These exercise the general case of
  /// Felsenstein's peeling where leaf contributions are not simple indicator vectors.
  const PYTHON_ALN: &str = ">A\nACATCGCCNNA--GAC\n>B\nGCATCCCTGTA-NG--\n>C\nCCGGCGATGTRTTG--\n>D\nTCGGCCGTGTRTTG--\n";

  /// Verify marginal posterior distribution at internal node AB (position 0) matches
  /// Python v0 reference values.
  ///
  /// The marginal posterior P(s|D) at an internal node is the product of the upward
  /// partial likelihood L_down(s) (Felsenstein peeling from leaves below) and the
  /// downward partial likelihood L_up(s) (from the rest of the tree above), normalized
  /// to sum to 1. At node AB position 0, leaf A has character A (branch length 0.1) and
  /// leaf B has character G (branch length 0.2). The posterior is biased toward A
  /// (P(A)=0.51) primarily because A's shorter branch produces a sharper transition
  /// probability vector (less diffusion), and secondarily because the downward message
  /// incorporates the equilibrium prior pi and the CD subtree's contribution.
  ///
  /// Expected: [0.51275208, 0.09128506, 0.24647255, 0.14949031] (A, C, G, T order).
  /// Python reference: test_scripts/ancestral_sparse.py:249-250
  #[test]
  fn test_internal_node_ab_profile_matches_python() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?;
    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense::new(0, gtr, alphabet, get_common_length(&aln)?)))];

    initialize_marginal(&graph, &partitions, &aln)?;

    // Find node AB and check profile at position 0
    let ab_key = find_node_key_by_name(&graph, "AB").ok_or_else(|| make_report!("Node AB not found"))?;
    let partition = partitions[0].read_arc();
    let ab_profile = &partition.data.nodes[&ab_key].profile.dis;

    // Position 0 is variable in Python (first variable position)
    // Python: [0.51275208, 0.09128506, 0.24647255, 0.14949031]
    // Order: A, C, G, T
    let pos0_profile = ab_profile.row(0);

    pretty_assert_ulps_eq!(pos0_profile[0], 0.51275208, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_profile[1], 0.09128506, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_profile[2], 0.24647255, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_profile[3], 0.14949031, epsilon = 1e-6);

    Ok(())
  }

  /// Verify marginal posterior distribution at the root node (position 0) matches Python
  /// v0 reference values.
  ///
  /// The root posterior is P(s|D) proportional to L_root(s) * pi(s), where L_root(s) is
  /// the product of upward partial likelihoods from the AB and CD subtrees, and pi is the
  /// equilibrium frequency vector. At position 0 the four leaves contribute characters
  /// A, G, C, T respectively - no single state dominates the likelihood. The prior
  /// pi = [0.2, 0.3, 0.15, 0.35] gives T the highest equilibrium weight, pulling
  /// P(T|D) = 0.36 above what the likelihood alone would produce.
  ///
  /// Expected: [0.28212327, 0.21643546, 0.13800802, 0.36343326] (A, C, G, T order).
  /// Python reference: test_scripts/ancestral_sparse.py:245
  #[test]
  fn test_root_profile_matches_python() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?;
    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense::new(0, gtr, alphabet, get_common_length(&aln)?)))];

    initialize_marginal(&graph, &partitions, &aln)?;

    // Find root node
    let root_key = find_node_key_by_name(&graph, "root").ok_or_else(|| make_report!("Node root not found"))?;
    let partition = partitions[0].read_arc();
    let root_profile = &partition.data.nodes[&root_key].profile.dis;

    // Position 0 profile
    // Python: [0.28212327, 0.21643546, 0.13800802, 0.36343326]
    let pos0_profile = root_profile.row(0);

    pretty_assert_ulps_eq!(pos0_profile[0], 0.28212327, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_profile[1], 0.21643546, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_profile[2], 0.13800802, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_profile[3], 0.36343326, epsilon = 1e-6);

    Ok(())
  }

  /// Verify that internal node CD has valid normalized posterior distributions at every
  /// position.
  ///
  /// Complements `test_internal_node_ab_profile_matches_python()` by checking the other
  /// major subtree (C,D)CD. While AB-side posteriors are checked against exact Python
  /// reference values, the CD-side is validated structurally: each row of the posterior
  /// matrix must be a valid probability distribution (non-negative, finite, sums to 1).
  /// This catches asymmetric bugs where one subtree computes correctly but the other
  /// does not.
  #[test]
  fn test_internal_node_cd_profile_valid() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?;
    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense::new(0, gtr, alphabet, get_common_length(&aln)?)))];

    initialize_marginal(&graph, &partitions, &aln)?;

    let cd_key = find_node_key_by_name(&graph, "CD").ok_or_else(|| make_report!("Node CD not found"))?;
    let partition = partitions[0].read_arc();
    let cd_profile = &partition.data.nodes[&cd_key].profile.dis;

    // Verify all positions are normalized and valid
    for (pos, row) in cd_profile.rows().into_iter().enumerate() {
      let sum: f64 = row.sum();
      assert!(
        (sum - 1.0).abs() < 1e-10,
        "CD profile row {pos} not normalized: sum={sum}"
      );
      for (state, &val) in row.iter().enumerate() {
        assert!(val >= 0.0, "CD profile row {pos} state {state} negative: {val}");
        assert!(val.is_finite(), "CD profile row {pos} state {state} not finite: {val}");
      }
    }

    Ok(())
  }

  /// Verify that every internal node in the tree has normalized posterior distributions
  /// at all alignment positions.
  ///
  /// Iterates over all node profiles in the partition and checks sum_s P(s|D) = 1 at
  /// every position. This is the simplex constraint: marginal posteriors from the
  /// two-pass peeling algorithm (backward + forward) must be valid probability
  /// distributions everywhere in the tree, not just at spot-checked nodes.
  #[test]
  fn test_all_internal_nodes_normalized() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?;
    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense::new(0, gtr, alphabet, get_common_length(&aln)?)))];

    initialize_marginal(&graph, &partitions, &aln)?;

    let partition = partitions[0].read_arc();

    for (key, node_data) in &partition.data.nodes {
      let profile = &node_data.profile.dis;
      for (pos, row) in profile.rows().into_iter().enumerate() {
        let sum: f64 = row.sum();
        assert!(
          (sum - 1.0).abs() < 1e-10,
          "Node {key:?} row {pos} not normalized: sum={sum}"
        );
      }
    }

    Ok(())
  }

  /// Verify that two partitions with different F81 models on the same tree compute
  /// independently and match their respective Python v0 reference values.
  ///
  /// The partition system allows multiple data types (or the same data with different
  /// models) to share a single tree topology. Each partition maintains its own node/edge
  /// data and runs the two-pass marginal peeling algorithm independently. This test uses
  /// two F81 models with pi1 = [0.2, 0.3, 0.15, 0.35] and pi2 = [0.4, 0.15, 0.12, 0.33]
  /// on the same alignment and tree, verifying that root posteriors differ according to
  /// each model's equilibrium frequencies and match their respective Python v0 values.
  ///
  /// Python reference: test_scripts/ancestral_sparse.py:252-274
  #[test]
  fn test_multi_partition_independent_computation() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?;

    let gtr1 = make_python_reference_gtr()?;
    let gtr2 = make_python_reference_gtr2()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let length = get_common_length(&aln)?;

    // Create two partitions with different F81 models
    let partition1 = Arc::new(RwLock::new(PartitionMarginalDense::new(0, gtr1, alphabet.clone(), length)));

    let partition2 = Arc::new(RwLock::new(PartitionMarginalDense::new(1, gtr2, alphabet, length)));

    let partitions = [Arc::clone(&partition1), Arc::clone(&partition2)];

    initialize_marginal(&graph, &partitions, &aln)?;
    let root_key = find_node_key_by_name(&graph, "root").ok_or_else(|| make_report!("Node root not found"))?;

    let p1 = partition1.read_arc();
    let p2 = partition2.read_arc();

    let root1 = &p1.data.nodes[&root_key].profile.dis;
    let root2 = &p2.data.nodes[&root_key].profile.dis;

    // Position 0 profiles should differ due to different pi values
    // Python F81 (pi1): [0.28212327, 0.21643546, 0.13800802, 0.36343326]
    // Python F81 (pi2): [0.29664652, 0.20554302, 0.13582953, 0.36198093]
    let pos0_p1 = root1.row(0);
    let pos0_p2 = root2.row(0);

    // Verify partition 1 matches Python F81 (pi1) reference
    pretty_assert_ulps_eq!(pos0_p1[0], 0.28212327, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p1[1], 0.21643546, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p1[2], 0.13800802, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p1[3], 0.36343326, epsilon = 1e-6);

    // Verify partition 2 matches Python F81 (pi2) reference
    pretty_assert_ulps_eq!(pos0_p2[0], 0.29664652, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p2[1], 0.20554302, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p2[2], 0.13582953, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p2[3], 0.36198093, epsilon = 1e-6);

    Ok(())
  }

  /// Verify internal node AB posteriors in both partitions match their respective Python
  /// v0 reference values.
  ///
  /// Complements `test_multi_partition_independent_computation()` by checking an internal
  /// node rather than the root. At internal nodes, the downward message L_up(s) carries
  /// the equilibrium prior through the root, so different pi values between partitions
  /// alter the AB posterior even though both partitions share the same tree and alignment.
  /// Partition 1 (pi1) yields [0.51275208, 0.09128506, 0.24647255, 0.14949031] and
  /// partition 2 (pi2) yields [0.52331521, 0.08336271, 0.24488808, 0.148434].
  ///
  /// Python reference: test_scripts/ancestral_sparse.py:273-274
  #[test]
  fn test_multi_partition_internal_node_ab() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?;

    let gtr1 = make_python_reference_gtr()?;
    let gtr2 = make_python_reference_gtr2()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let length = get_common_length(&aln)?;

    let partition1 = Arc::new(RwLock::new(PartitionMarginalDense::new(0, gtr1, alphabet.clone(), length)));

    let partition2 = Arc::new(RwLock::new(PartitionMarginalDense::new(1, gtr2, alphabet, length)));

    let partitions = [Arc::clone(&partition1), Arc::clone(&partition2)];

    initialize_marginal(&graph, &partitions, &aln)?;

    let ab_key = find_node_key_by_name(&graph, "AB").ok_or_else(|| make_report!("Node AB not found"))?;

    let p1 = partition1.read_arc();
    let ab1 = &p1.data.nodes[&ab_key].profile.dis;
    let pos0_ab1 = ab1.row(0);

    pretty_assert_ulps_eq!(pos0_ab1[0], 0.51275208, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab1[1], 0.09128506, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab1[2], 0.24647255, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab1[3], 0.14949031, epsilon = 1e-6);

    // Partition 2: [0.52331521, 0.08336271, 0.24488808, 0.148434]
    let p2 = partition2.read_arc();
    let ab2 = &p2.data.nodes[&ab_key].profile.dis;
    let pos0_ab2 = ab2.row(0);

    pretty_assert_ulps_eq!(pos0_ab2[0], 0.52331521, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab2[1], 0.08336271, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab2[2], 0.24488808, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab2[3], 0.148434, epsilon = 1e-6);

    Ok(())
  }

  /// Verify that dense and sparse representations produce identical log-likelihoods on
  /// clean sequences (no gaps, no ambiguous characters).
  ///
  /// The sparse representation groups invariant positions by observed character and stores
  /// only variable positions as full probability vectors. The dense representation stores
  /// full vectors at every position. For unambiguous ACGT-only sequences, leaf
  /// observations are Kronecker delta indicator vectors in both representations, so the
  /// Felsenstein site likelihoods must be identical. Any discrepancy would indicate a bug
  /// in the sparse compression or aggregation of invariant site contributions. This
  /// cross-validates the sparse optimization against the dense implementation.
  #[test]
  fn test_multi_partition_sparse_dense_consistency() -> Result<(), Report> {
    // Use simple alignment without gaps or ambiguous characters
    let simple_aln = ">A\nACATCGCCTTACGGAC\n>B\nGCATCCCTGTACTGAC\n>C\nCCGGCGATGTATTGAC\n>D\nTCGGCCGTGTATTGAC\n";

    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(simple_aln, &*NUC_ALPHABET)?;

    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let length = get_common_length(&aln)?;

    // Dense partition
    let dense_partition = Arc::new(RwLock::new(PartitionMarginalDense::new(0, gtr.clone(), alphabet.clone(), length)));

    let dense_log_lh = initialize_marginal(&graph, from_ref(&dense_partition), &aln)?;

    // Sparse partition
    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;
    let sparse_partition = Arc::new(RwLock::new(fitch.into_marginal_sparse(gtr, &graph)?));
    let sparse_log_lh = update_marginal(&graph, from_ref(&sparse_partition))?;

    // Log-likelihoods should match for clean sequences
    pretty_assert_ulps_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-10);

    Ok(())
  }
}

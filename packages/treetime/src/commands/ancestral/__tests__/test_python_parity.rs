#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{ancestral_reconstruction_marginal, initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::pretty_assert_ulps_eq;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::array;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::path::PathBuf;
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

  /// Parity test: root sequence from marginal reconstruction matches Python v0.
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

    let treat_gap_as_unknown = true;
    let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
    let aln = read_many_fasta(&[aln_path], &alphabet)?;

    let graph: GraphAncestral = nwk_read_file(&tree_path)?;

    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown,
      ..JC69Params::default()
    })?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

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

  // ============================================================================
  // T12: Comprehensive internal node validation
  // ============================================================================

  /// Lazily initialized default nucleotide alphabet (A, C, G, T with gap handling).
  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  /// Create the primary GTR model used in Python v0 parity tests, with non-uniform
  /// equilibrium frequencies pi = [0.2, 0.3, 0.15, 0.35].
  ///
  /// Symmetric rate matrix W (default) with these frequencies produces a reversible model
  /// satisfying detailed balance. The non-uniform pi ensures posteriors are not trivially
  /// symmetric across states.
  ///
  /// Python reference: test_scripts/ancestral_sparse.py:241-250
  fn make_python_reference_gtr() -> Result<GTR, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.2, 0.3, 0.15, 0.35],
    })
  }

  /// Create the secondary GTR model for multi-partition parity tests, with equilibrium
  /// frequencies pi = [0.4, 0.15, 0.12, 0.33].
  ///
  /// Used alongside `make_python_reference_gtr()` to verify that two partitions with
  /// different GTR parameters on the same tree compute independently. The different pi
  /// values produce distinct posterior distributions at each node.
  ///
  /// Python reference: test_scripts/ancestral_sparse.py:255
  fn make_python_reference_gtr2() -> Result<GTR, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.4, 0.15, 0.12, 0.33],
    })
  }

  /// Newick tree shared across Python parity tests: balanced 4-taxon tree
  /// ((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01.
  const PYTHON_TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  /// FASTA alignment shared across Python parity tests. Contains ambiguous characters
  /// (N, R) and gaps (-) to exercise non-trivial likelihood handling at those positions.
  const PYTHON_ALN: &str = ">A\nACATCGCCNNA--GAC\n>B\nGCATCCCTGTA-NG--\n>C\nCCGGCGATGTRTTG--\n>D\nTCGGCCGTGTRTTG--\n";

  /// Verify marginal posterior distribution at internal node AB (position 0) matches
  /// Python v0 reference values.
  ///
  /// The marginal posterior P(s|data) at an internal node is computed by combining the
  /// upward partial likelihood (from leaf data below) with the downward message (from the
  /// rest of the tree above). At node AB position 0, leaves A and B have characters A and
  /// G respectively, producing a posterior biased toward A (0.51) due to A's shorter
  /// branch length (0.1 vs 0.2).
  ///
  /// Expected: [0.51275208, 0.09128506, 0.24647255, 0.14949031] (A, C, G, T order).
  /// Python reference: test_scripts/ancestral_sparse.py:249-250
  #[test]
  fn test_internal_node_ab_profile_matches_python() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?;
    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    initialize_marginal(&graph, &partitions, &aln)?;

    // Find node AB and check profile at position 0
    let ab_key = find_node_key_by_name(&graph, "AB").ok_or_else(|| make_report!("Node AB not found"))?;
    let partition = partitions[0].read();
    let ab_profile = &partition.nodes[&ab_key].profile.dis;

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
  /// The root posterior combines upward partial likelihoods from both subtrees (AB and CD)
  /// with the equilibrium distribution pi as prior. With pi = [0.2, 0.3, 0.15, 0.35],
  /// T has the highest prior weight, which shifts the posterior toward T (0.36) despite
  /// the leaf data favoring other states.
  ///
  /// Expected: [0.28212327, 0.21643546, 0.13800802, 0.36343326] (A, C, G, T order).
  /// Python reference: test_scripts/ancestral_sparse.py:245
  #[test]
  fn test_root_profile_matches_python() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?;
    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    initialize_marginal(&graph, &partitions, &aln)?;

    // Find root node
    let root_key = find_node_key_by_name(&graph, "root").ok_or_else(|| make_report!("Node root not found"))?;
    let partition = partitions[0].read();
    let root_profile = &partition.nodes[&root_key].profile.dis;

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
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    initialize_marginal(&graph, &partitions, &aln)?;

    let cd_key = find_node_key_by_name(&graph, "CD").ok_or_else(|| make_report!("Node CD not found"))?;
    let partition = partitions[0].read();
    let cd_profile = &partition.nodes[&cd_key].profile.dis;

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
  /// Iterates over all node profiles in the partition and checks sum_s P(s|data) = 1 at
  /// every position. This is a global invariant of Felsenstein's marginal reconstruction:
  /// the upward-downward message passing with proper normalization must produce valid
  /// probability distributions everywhere in the tree, not just at spot-checked nodes.
  #[test]
  fn test_all_internal_nodes_normalized() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?;
    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    initialize_marginal(&graph, &partitions, &aln)?;

    let partition = partitions[0].read();

    for (key, node_data) in &partition.nodes {
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

  // ============================================================================
  // T13: Multi-partition test
  // ============================================================================

  /// Verify that two partitions with different GTR models on the same tree compute
  /// independently and match their respective Python v0 reference values.
  ///
  /// The partition system allows multiple data types (or the same data with different
  /// models) to share a single tree topology. Each partition maintains its own node/edge
  /// data and runs Felsenstein's pruning independently. This test uses
  /// pi1 = [0.2, 0.3, 0.15, 0.35] and pi2 = [0.4, 0.15, 0.12, 0.33] on the same
  /// alignment and tree, verifying that root posteriors differ according to each model's
  /// equilibrium frequencies and match their respective Python v0 values.
  ///
  /// Python reference: test_scripts/ancestral_sparse.py:252-274
  #[test]
  fn test_multi_partition_independent_computation() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?;

    let gtr1 = make_python_reference_gtr()?;
    let gtr2 = make_python_reference_gtr2()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let length = get_common_length(&aln)?;

    // Create two partitions with different GTRs
    let partition1 = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: gtr1,
      alphabet: alphabet.clone(),
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    let partition2 = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 1,
      gtr: gtr2,
      alphabet,
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    let partitions = [Arc::clone(&partition1), Arc::clone(&partition2)];

    initialize_marginal(&graph, &partitions, &aln)?;

    // Check root profiles differ between partitions
    let root_key = find_node_key_by_name(&graph, "root").ok_or_else(|| make_report!("Node root not found"))?;

    let p1 = partition1.read();
    let p2 = partition2.read();

    let root1 = &p1.nodes[&root_key].profile.dis;
    let root2 = &p2.nodes[&root_key].profile.dis;

    // Position 0 profiles should differ due to different pi values
    // Python GTR1: [0.28212327, 0.21643546, 0.13800802, 0.36343326]
    // Python GTR2: [0.29664652, 0.20554302, 0.13582953, 0.36198093]
    let pos0_p1 = root1.row(0);
    let pos0_p2 = root2.row(0);

    // Verify partition 1 matches Python GTR1 reference
    pretty_assert_ulps_eq!(pos0_p1[0], 0.28212327, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p1[1], 0.21643546, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p1[2], 0.13800802, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p1[3], 0.36343326, epsilon = 1e-6);

    // Verify partition 2 matches Python GTR2 reference
    pretty_assert_ulps_eq!(pos0_p2[0], 0.29664652, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p2[1], 0.20554302, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p2[2], 0.13582953, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p2[3], 0.36198093, epsilon = 1e-6);

    Ok(())
  }

  /// Verify internal node AB posteriors in both partitions match their respective Python
  /// v0 reference values.
  ///
  /// Extends `test_multi_partition_independent_computation()` by checking an internal node
  /// rather than the root. The two partitions use different equilibrium frequencies, so
  /// the AB posterior at position 0 differs between them. Partition 1 (pi1) yields
  /// [0.51275208, 0.09128506, 0.24647255, 0.14949031] and partition 2 (pi2) yields
  /// [0.52331521, 0.08336271, 0.24488808, 0.148434].
  ///
  /// Python reference: test_scripts/ancestral_sparse.py:273-274
  #[test]
  fn test_multi_partition_internal_node_ab() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?;

    let gtr1 = make_python_reference_gtr()?;
    let gtr2 = make_python_reference_gtr2()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let length = get_common_length(&aln)?;

    let partition1 = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: gtr1,
      alphabet: alphabet.clone(),
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    let partition2 = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 1,
      gtr: gtr2,
      alphabet,
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    let partitions = [Arc::clone(&partition1), Arc::clone(&partition2)];

    initialize_marginal(&graph, &partitions, &aln)?;

    let ab_key = find_node_key_by_name(&graph, "AB").ok_or_else(|| make_report!("Node AB not found"))?;

    // Partition 1: [0.51275208, 0.09128506, 0.24647255, 0.14949031]
    let p1 = partition1.read();
    let ab1 = &p1.nodes[&ab_key].profile.dis;
    let pos0_ab1 = ab1.row(0);

    pretty_assert_ulps_eq!(pos0_ab1[0], 0.51275208, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab1[1], 0.09128506, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab1[2], 0.24647255, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab1[3], 0.14949031, epsilon = 1e-6);

    // Partition 2: [0.52331521, 0.08336271, 0.24488808, 0.148434]
    let p2 = partition2.read();
    let ab2 = &p2.nodes[&ab_key].profile.dis;
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
  /// full vectors at every position. For unambiguous ACGT-only sequences, both approaches
  /// must compute the same Felsenstein likelihood because the leaf observations are
  /// identical delta distributions in both cases. This cross-validates the sparse
  /// optimization against the simpler dense implementation.
  #[test]
  fn test_multi_partition_sparse_dense_consistency() -> Result<(), Report> {
    // Use simple alignment without gaps or ambiguous characters
    let simple_aln = ">A\nACATCGCCTTACGGAC\n>B\nGCATCCCTGTACTGAC\n>C\nCCGGCGATGTATTGAC\n>D\nTCGGCCGTGTATTGAC\n";

    let graph: GraphAncestral = nwk_read_str(PYTHON_TREE)?;
    let aln = read_many_fasta_str(simple_aln, &*NUC_ALPHABET)?;

    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let length = get_common_length(&aln)?;

    // Dense partition
    let dense_partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: gtr.clone(),
      alphabet: alphabet.clone(),
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    let dense_log_lh = initialize_marginal(&graph, std::slice::from_ref(&dense_partition), &aln)?;

    // Sparse partition
    let sparse_partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr,
      alphabet,
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    compress_sequences(&graph, std::slice::from_ref(&sparse_partition), &aln)?;
    let sparse_log_lh = update_marginal(&graph, std::slice::from_ref(&sparse_partition))?;

    // Log-likelihoods should match for clean sequences
    pretty_assert_ulps_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-10);

    Ok(())
  }
}

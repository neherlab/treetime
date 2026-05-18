#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::{ancestral_reconstruction_marginal, initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::pretty_assert_ulps_eq;
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::partition::traits::PartitionBranchOps;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use crate::seq::mutation::Sub;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use indoc::indoc;
  
  use ndarray::{Array1, array};
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::slice::from_ref;
  use std::sync::{Arc, LazyLock};
  use treetime_graph::node::Named;
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

  /// 4-taxon alignment shaped like the sparse ambiguity bug.
  ///
  /// At position 1, leaf A is `R = {A,G}` while all other leaves are canonical `G`.
  /// Fitch forward can resolve the A leaf canonically to `G` from its parent context,
  /// but sparse marginal leaf encoding must not rebuild a false canonical `A` from
  /// the ambiguous state set.
  fn ambiguous_r_in_g_clade_alignment() -> Result<Vec<FastaRecord>, Report> {
    read_many_fasta_str(
      indoc! {r#"
      >A
      RCGTACGT
      >B
      GCGTACGT
      >C
      GCGTACGT
      >D
      GCGTACGT
    "#},
      &*NUC_ALPHABET,
    )
  }

  /// Run dense marginal ancestral reconstruction on the given tree and alignment.
  ///
  /// Dense representation stores full N x K probability matrices where N is the
  /// alignment length and K is the alphabet size (4 for nucleotides). Runs
  /// Felsenstein's pruning algorithm (backward pass, leaves to root) to compute
  /// partial likelihoods and the total log-likelihood, then the forward pass
  /// (root to leaves) to produce marginal posterior probabilities at every node.
  ///
  /// Returns the total log-likelihood and the populated partition.
  fn run_dense_marginal(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
    gtr: GTR,
  ) -> Result<(f64, Arc<RwLock<PartitionMarginalDense>>), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partition = Arc::new(RwLock::new(PartitionMarginalDense::new(0, gtr, alphabet, get_common_length(aln)?)));
    let partitions = [Arc::clone(&partition)];

    let log_lh = initialize_marginal(graph, &partitions, aln)?;
    Ok((log_lh, partition))
  }

  /// Run sparse marginal ancestral reconstruction on the given tree and alignment.
  ///
  /// Sparse representation stores only "variable" positions where the marginal
  /// posterior differs from the fixed-character distribution. Positions where all
  /// descendants agree on a single state are compressed away via Fitch parsimony,
  /// storing only the Fitch-compressed reference and per-character fixed
  /// distributions.
  ///
  /// The sparse path first runs `compress_sequences` (Fitch parsimony), then
  /// `update_marginal` which performs the backward pass (leaves to root) and
  /// forward pass (root to leaves) on the compressed representation.
  ///
  /// Returns the total log-likelihood and the populated partition.
  fn run_sparse_marginal(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
    gtr: GTR,
  ) -> Result<(f64, Arc<RwLock<PartitionMarginalSparse>>), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let fitch = create_fitch_partition(graph, 0, alphabet, aln)?;
    let partition = Arc::new(RwLock::new(fitch.into_marginal_sparse(gtr, graph)?));
    let partitions = [Arc::clone(&partition)];

    let log_lh = update_marginal(graph, &partitions)?;
    Ok((log_lh, partition))
  }

  /// Verify that dense and sparse marginal implementations produce the same
  /// total log-likelihood on a gap-free alignment with JC69 model.
  ///
  /// Both representations compute the same quantity - the Felsenstein tree
  /// likelihood evaluated at the root after the backward pass:
  /// log L = sum over sites of log(sum over states s of pi[s] * L_root(s)),
  /// where L_root(s) is the partial likelihood at the root for state s. Dense
  /// operates on full NxK matrices; sparse compresses invariant sites and
  /// operates only on variable positions. The log-likelihoods must agree to
  /// floating-point precision.
  #[test]
  fn test_marginal_dense_sparse_log_lh_consistency_gap_free() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let gtr_dense = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let gtr_sparse = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let (log_lh_dense, _) = run_dense_marginal(&graph, &aln, gtr_dense)?;
    let (log_lh_sparse, _) = run_sparse_marginal(&graph, &aln, gtr_sparse)?;

    pretty_assert_ulps_eq!(log_lh_dense, log_lh_sparse, epsilon = 1e-10);

    Ok(())
  }

  /// Verify that sparse variable-position distributions match the corresponding
  /// rows in the dense marginal posterior matrix.
  ///
  /// At internal nodes (root, AB), for each position the sparse representation
  /// marks as "variable", the probability vector must equal the corresponding row
  /// of the dense NxK matrix. This checks profile-level consistency beyond the
  /// aggregate log-likelihood: the per-position marginal posteriors must agree
  /// between representations. At each node, the posterior combines the subtree
  /// partial likelihood (from the backward pass) with the outgroup message (from
  /// the forward pass), normalized to sum to 1.
  #[test]
  fn test_marginal_sparse_varpos_matches_dense_profile_gap_free() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let gtr_dense = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let gtr_sparse = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let (_, dense_partition) = run_dense_marginal(&graph, &aln, gtr_dense)?;
    let (_, sparse_partition) = run_sparse_marginal(&graph, &aln, gtr_sparse)?;

    let root_key = find_node_key_by_name(&graph, "root").ok_or_else(|| make_report!("Root node not found"))?;
    let ab_key = find_node_key_by_name(&graph, "AB").ok_or_else(|| make_report!("AB node not found"))?;

    let dense = dense_partition.read_arc();
    let sparse = sparse_partition.read_arc();

    for node_key in [root_key, ab_key] {
      let dense_node = &dense.data.nodes[&node_key];
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
          pretty_assert_ulps_eq!(dense_val, sparse_val, epsilon = 1e-6);
        }
      }
    }

    Ok(())
  }

  /// Verify dense/sparse consistency and probability normalization in the
  /// presence of IUPAC ambiguity codes (N = any nucleotide, R = A|G).
  ///
  /// Ambiguous characters assign equal probability to the states they represent
  /// (e.g., N -> uniform over {A,C,G,T}, R -> uniform over {A,G}). Both dense
  /// and sparse use the same tip-state encoding, but differ in how invariant
  /// vs variable positions are stored and processed. This test checks three
  /// properties:
  ///
  /// 1. Log-likelihood agreement: dense and sparse produce the same total
  ///    log-likelihood despite different internal representations.
  /// 2. Dense normalization: every row of every node's marginal posterior
  ///    matrix sums to 1.0 (valid probability distribution).
  /// 3. Sparse normalization: variable-position distributions, fixed-character
  ///    distributions, and per-node log-likelihoods are all finite and
  ///    distributions sum to 1.0.
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
      ..JC69Params::default()
    })?;
    let gtr_sparse = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let (log_lh_dense, dense_partition) = run_dense_marginal(&graph, &aln, gtr_dense)?;
    let (log_lh_sparse, sparse_partition) = run_sparse_marginal(&graph, &aln, gtr_sparse)?;

    pretty_assert_ulps_eq!(log_lh_dense, log_lh_sparse, epsilon = 1e-10);

    let dense = dense_partition.read_arc();
    let sparse = sparse_partition.read_arc();

    for node_data in dense.data.nodes.values() {
      if !node_data.profile.dis.is_empty() {
        for row in node_data.profile.dis.rows() {
          let sum: f64 = row.sum();
          assert!(sum.is_finite(), "Dense node profile row sum is not finite: {sum}");
          pretty_assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
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
        pretty_assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
      }

      for (char_key, fixed_dis) in &node_data.profile.fixed {
        let sum: f64 = fixed_dis.sum();
        assert!(
          sum.is_finite(),
          "Sparse fixed distribution for char {char_key:?} sum is not finite: {sum}"
        );
        pretty_assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
      }
    }

    Ok(())
  }

  #[test]
  fn test_marginal_dense_sparse_ambiguous_r_reference_state_consistency() -> Result<(), Report> {
    let aln = ambiguous_r_in_g_clade_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let gtr_dense = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let gtr_sparse = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let (log_lh_dense, dense_partition) = run_dense_marginal(&graph, &aln, gtr_dense)?;
    let (log_lh_sparse, sparse_partition) = run_sparse_marginal(&graph, &aln, gtr_sparse)?;

    pretty_assert_ulps_eq!(log_lh_dense, log_lh_sparse, epsilon = 1e-10);

    let dense_sequences = reconstruct_named_sequences(&graph, from_ref(&dense_partition))?;
    let sparse_sequences = reconstruct_named_sequences(&graph, from_ref(&sparse_partition))?;
    assert_eq!(dense_sequences, sparse_sequences);

    let dense_branch_subs = edge_subs_by_edge_name(&graph, &*dense_partition.read_arc())?;
    let sparse_branch_subs = edge_subs_by_edge_name(&graph, &*sparse_partition.read_arc())?;
    assert_eq!(dense_branch_subs, sparse_branch_subs);

    Ok(())
  }

  fn reconstruct_named_sequences<P>(
    graph: &GraphAncestral,
    partitions: &[Arc<RwLock<P>>],
  ) -> Result<BTreeMap<String, String>, Report>
  where
    P: crate::partition::traits::PartitionMarginalOps<
        crate::payload::ancestral::NodeAncestral,
        crate::payload::ancestral::EdgeAncestral,
      > + crate::partition::traits::HasLogLh,
  {
    let mut actual = BTreeMap::new();
    ancestral_reconstruction_marginal(graph, false, partitions, |node, seq| {
      actual.insert(node.name.clone().expect("all test nodes are named"), seq.to_string());
      Ok(())
    })?;
    Ok(actual)
  }

  fn edge_subs_by_edge_name<P>(graph: &GraphAncestral, partition: &P) -> Result<BTreeMap<String, Vec<Sub>>, Report>
  where
    P: PartitionBranchOps,
  {
    graph
      .get_edges()
      .iter()
      .map(|edge_ref| {
        let edge = edge_ref.read_arc();
        let parent_name = graph
          .get_node(edge.source())
          .expect("parent node exists")
          .read_arc()
          .payload()
          .read_arc()
          .name()
          .map(|name| name.as_ref().to_owned())
          .expect("named parent");
        let child_name = graph
          .get_node(edge.target())
          .expect("child node exists")
          .read_arc()
          .payload()
          .read_arc()
          .name()
          .map(|name| name.as_ref().to_owned())
          .expect("named child");
        let edge_name = format!("{parent_name}->{child_name}");
        let subs = partition.edge_subs(graph, edge.key())?;
        Ok((edge_name, subs))
      })
      .collect()
  }

  /// Verify marginal posterior normalization under a highly skewed GTR model.
  ///
  /// Uses the tree, alignment, and GTR parameters from Python v0
  /// `test_treetime.py:137-155`, but tests a different invariant. The Python
  /// test checks that per-site likelihoods summed over all 4^3 = 64 data
  /// patterns equal 1.0 (likelihood normalization over observable data). This
  /// test checks that the marginal posterior at every node is a valid
  /// probability distribution: sum_s P(s|data) = 1 for each position and node.
  ///
  /// The skewed equilibrium frequencies pi = [0.9, 0.06, 0.02, 0.02] with
  /// uniform exchangeability (W_ij = 1 for i != j) create extreme asymmetry in
  /// P(t) = exp(Q*t). Combined with long branches (A: 0.601, B: 0.301
  /// subs/site) and a 64bp alignment covering all 4^3 three-taxon state
  /// combinations, this stresses numerical stability of the backward and
  /// forward passes.
  #[test]
  fn test_marginal_posteriors_sum_to_one_skewed_gtr() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    // Highly skewed equilibrium frequencies: A dominates at 90%
    let pi = array![0.9, 0.06, 0.02, 0.02];

    // Uniform exchangeability (W_ij = 1 for i != j, zero diagonal)
    let W = array![
      [0.0, 1.0, 1.0, 1.0],
      [1.0, 0.0, 1.0, 1.0],
      [1.0, 1.0, 0.0, 1.0],
      [1.0, 1.0, 1.0, 0.0],
    ];

    let gtr = GTR::new(GTRParams {
      n_states: alphabet.n_canonical(),
      mu: 1.0,
      W: Some(W),
      pi,
    })?;

    // 3-taxon tree with long branches (from v0 test, branch lengths rounded)
    let tree_newick = "((A:0.601,B:0.301):0.1,C:0.2):0.001;";
    let graph: GraphAncestral = nwk_read_str(tree_newick)?;

    // 64bp alignment: all 4^3 = 64 three-taxon state combinations (A x B x C)
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

    let partition = Arc::new(RwLock::new(PartitionMarginalDense::new(0, gtr, alphabet, get_common_length(&aln)?)));
    let partitions = [Arc::clone(&partition)];

    initialize_marginal(&graph, &partitions, &aln)?;

    // Verify all marginal posterior rows sum to 1.0
    let partition = partition.read_arc();
    for (node_key, node_data) in &partition.data.nodes {
      if node_data.profile.dis.is_empty() {
        continue;
      }
      for (pos, row) in node_data.profile.dis.rows().into_iter().enumerate() {
        let sum: f64 = row.sum();
        assert!(
          sum.is_finite(),
          "Node {node_key:?} position {pos}: posterior sum is not finite: {sum}"
        );
        pretty_assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
      }
    }

    Ok(())
  }

  /// G3: Uniform per-site rates produce the same sparse log-likelihood as
  /// scalar mu. Verifies the per-site rate code path does not introduce new
  /// divergence beyond the existing dense-sparse gap.
  #[test]
  fn test_marginal_sparse_uniform_site_rates_matches_scalar() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let seq_len = get_common_length(&aln)?;

    // Run without site_rates (scalar mu)
    let gtr_scalar = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let (log_lh_scalar, _) = run_sparse_marginal(&graph, &aln, gtr_scalar)?;

    // Run with uniform site_rates (all 1.0) - must produce identical results
    let mut gtr_uniform = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    gtr_uniform.set_site_rates(Array1::ones(seq_len));
    let (log_lh_uniform, _) = run_sparse_marginal(&graph, &aln, gtr_uniform)?;

    pretty_assert_ulps_eq!(log_lh_scalar, log_lh_uniform, epsilon = 1e-10);

    Ok(())
  }

  /// G3b: Dense with uniform per-site rates matches dense with scalar mu.
  #[test]
  fn test_marginal_dense_uniform_site_rates_matches_scalar() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let seq_len = get_common_length(&aln)?;

    let gtr_scalar = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let (log_lh_scalar, _) = run_dense_marginal(&graph, &aln, gtr_scalar)?;

    let mut gtr_uniform = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    gtr_uniform.set_site_rates(Array1::ones(seq_len));
    let (log_lh_uniform, _) = run_dense_marginal(&graph, &aln, gtr_uniform)?;

    pretty_assert_ulps_eq!(log_lh_scalar, log_lh_uniform, epsilon = 1e-10);

    Ok(())
  }
}

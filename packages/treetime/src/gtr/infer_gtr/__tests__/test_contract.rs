//! Contract tests for `MutationCounts`: orientation, scaling, consistency, root state.
//!
//! These tests verify invariants that golden master tests do not isolate:
//! - nij orientation (child-row, parent-column)
//! - Ti linear scaling with branch length
//! - Dense-sparse agreement on well-resolved inputs
//! - Root state correctness for known ancestral sequences

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::infer_gtr::dense::get_mutation_counts_dense;
  use crate::gtr::infer_gtr::sparse::get_mutation_counts_sparse;
  use crate::pretty_assert_ulps_eq;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
  }

  // State index constants: A=0, C=1, G=2, T=3 (from canonical order ['A','C','G','T'])
  const IDX_A: usize = 0;
  const IDX_C: usize = 1;
  const IDX_G: usize = 2;
  const IDX_T: usize = 3;

  fn setup_dense(
    tree_nwk: &str,
    aln: &[FastaRecord],
  ) -> Result<(GraphAncestral, Arc<RwLock<PartitionMarginalDense>>), Report> {
    let graph: GraphAncestral = nwk_read_str(tree_nwk)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: true,
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

  fn setup_sparse(
    tree_nwk: &str,
    aln: &[FastaRecord],
  ) -> Result<(GraphAncestral, Arc<RwLock<PartitionMarginalSparse>>), Report> {
    let graph: GraphAncestral = nwk_read_str(tree_nwk)?;
    let alphabet = Alphabet::default();
    let partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));
    compress_sequences(&graph, std::slice::from_ref(&partition), aln)?;
    update_marginal(&graph, std::slice::from_ref(&partition))?;
    Ok((graph, partition))
  }

  /// Dense nij orientation: a single A->C mutation produces nij[C,A] >> nij[A,C].
  ///
  /// 4-leaf tree with 3 leaves having A at position 0, 1 leaf having C. The 3:1
  /// majority strongly anchors the root as A at position 0, making the A->C mutation
  /// unambiguously located on the branch leading to the mutant leaf.
  #[test]
  fn test_nij_orientation_dense() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >ref1
      AAAAAAAA
      >ref2
      AAAAAAAA
      >ref3
      AAAAAAAA
      >mut1
      CAAAAAAA
      "#},
      &*NUC_ALPHABET,
    )?;

    let (graph, partition) = setup_dense(
      "((ref1:0.1,ref2:0.1)R12:0.05,(ref3:0.1,mut1:0.1)R3M:0.05)root:0.0;",
      &aln,
    )?;
    let counts = get_mutation_counts_dense(&graph, &partition)?;

    // nij[C, A] should dominate: mutation from parent=A to child=C
    assert!(
      counts.nij[[IDX_C, IDX_A]] > 0.1,
      "nij[C,A] should carry A->C mutation signal, got {}",
      counts.nij[[IDX_C, IDX_A]]
    );
    // nij[A, C] should be near-zero: no mutation from parent=C to child=A
    assert!(
      counts.nij[[IDX_C, IDX_A]] > 5.0 * counts.nij[[IDX_A, IDX_C]],
      "nij[C,A] ({}) should be >> nij[A,C] ({}): child-row, parent-column convention",
      counts.nij[[IDX_C, IDX_A]],
      counts.nij[[IDX_A, IDX_C]]
    );

    Ok(())
  }

  /// Sparse nij orientation: a single A->C mutation produces nij[C,A] == 1 and nij[A,C] == 0.
  ///
  /// Same tree as dense test. Fitch parsimony reconstructs root as A at position 0
  /// (majority rule among leaves). The substitution list has exactly one A->C entry
  /// on the root->leaf_c edge.
  #[test]
  fn test_nij_orientation_sparse() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >leaf_a
      AAAAAAAA
      >leaf_c
      CAAAAAAA
      "#},
      &*NUC_ALPHABET,
    )?;

    let (graph, partition) = setup_sparse("(leaf_a:0.1,leaf_c:0.1)root:0.0;", &aln)?;
    let counts = get_mutation_counts_sparse(&graph, &partition)?;

    // Sparse gives exact integer counts: one A->C substitution
    pretty_assert_ulps_eq!(1.0, counts.nij[[IDX_C, IDX_A]], epsilon = 1e-9);
    pretty_assert_ulps_eq!(0.0, counts.nij[[IDX_A, IDX_C]], epsilon = 1e-9);

    Ok(())
  }

  /// Ti scales linearly with branch length: doubling BL doubles Ti.
  ///
  /// Star topology with 4 identical sequences (all ACGTACGT). No mutations exist,
  /// so Ti is purely branch_length * composition_count per state. Doubling BL
  /// doubles every Ti entry.
  #[rustfmt::skip]
  #[rstest]
  #[case::baseline_vs_double((0.1, 0.2))]
  #[case::small_vs_large(    (0.05, 0.5))]
  #[trace]
  fn test_ti_scaling_sparse(#[case] (bl1, bl2): (f64, f64)) -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGTACGT
      >B
      ACGTACGT
      >C
      ACGTACGT
      >D
      ACGTACGT
      "#},
      &*NUC_ALPHABET,
    )?;

    let tree1 = format!("((A:{bl1},B:{bl1})AB:{bl1},(C:{bl1},D:{bl1})CD:{bl1})root:0.0;");
    let tree2 = format!("((A:{bl2},B:{bl2})AB:{bl2},(C:{bl2},D:{bl2})CD:{bl2})root:0.0;");

    let (graph1, partition1) = setup_sparse(&tree1, &aln)?;
    let (graph2, partition2) = setup_sparse(&tree2, &aln)?;

    let counts1 = get_mutation_counts_sparse(&graph1, &partition1)?;
    let counts2 = get_mutation_counts_sparse(&graph2, &partition2)?;

    let ratio = bl2 / bl1;
    for k in 0..4 {
      pretty_assert_ulps_eq!(
        counts2.Ti[k],
        counts1.Ti[k] * ratio,
        epsilon = 1e-7
      );
    }

    Ok(())
  }

  /// Ti scales linearly with branch length for dense path.
  ///
  /// Same setup as sparse: star topology, identical sequences, no mutations.
  /// The midpoint approximation Ti += 0.5 * BL * (P_parent + P_child) also scales
  /// linearly because the posterior probabilities are BL-independent for identical
  /// sequences (profile is always [1,0,0,0] at A-positions etc).
  #[rustfmt::skip]
  #[rstest]
  #[case::baseline_vs_double((0.1, 0.2))]
  #[case::small_vs_large(    (0.05, 0.5))]
  #[trace]
  fn test_ti_scaling_dense(#[case] (bl1, bl2): (f64, f64)) -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGTACGT
      >B
      ACGTACGT
      >C
      ACGTACGT
      >D
      ACGTACGT
      "#},
      &*NUC_ALPHABET,
    )?;

    let tree1 = format!("((A:{bl1},B:{bl1})AB:{bl1},(C:{bl1},D:{bl1})CD:{bl1})root:0.0;");
    let tree2 = format!("((A:{bl2},B:{bl2})AB:{bl2},(C:{bl2},D:{bl2})CD:{bl2})root:0.0;");

    let (graph1, partition1) = setup_dense(&tree1, &aln)?;
    let (graph2, partition2) = setup_dense(&tree2, &aln)?;

    let counts1 = get_mutation_counts_dense(&graph1, &partition1)?;
    let counts2 = get_mutation_counts_dense(&graph2, &partition2)?;

    let ratio = bl2 / bl1;
    for k in 0..4 {
      // Measured max diff: 8.94e-8 (case small_vs_large, Ti[0])
      pretty_assert_ulps_eq!(
        counts2.Ti[k],
        counts1.Ti[k] * ratio,
        epsilon = 1e-7
      );
    }

    Ok(())
  }

  /// Dense and sparse mutation counts agree for unambiguous sequences.
  ///
  /// 4-taxon tree with a well-anchored root and few leaf-level mutations.
  /// Dense fractional counts should approximate sparse integer counts because
  /// the posterior concentrates on a single (child, parent) state pair per site.
  #[test]
  fn test_dense_sparse_consistency() -> Result<(), Report> {
    // Mostly-identical sequences with 1-2 mutations on leaf branches only.
    // The 3:1 majority at each position anchors the root state, avoiding
    // the ambiguity that causes dense-sparse divergence.
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGTACGT
      >B
      ACGTACGT
      >C
      ACGTACGT
      >D
      GCGTACGT
      "#},
      &*NUC_ALPHABET,
    )?;

    let tree_nwk = "((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;";

    let (graph_d, partition_d) = setup_dense(tree_nwk, &aln)?;
    let (graph_s, partition_s) = setup_sparse(tree_nwk, &aln)?;

    let dense = get_mutation_counts_dense(&graph_d, &partition_d)?;
    let sparse = get_mutation_counts_sparse(&graph_s, &partition_s)?;

    // nij: dense fractional counts should approximate sparse integer counts.
    // Measured max nij_diff: 6.68e-2
    let nij_diff = (&dense.nij - &sparse.nij).mapv(f64::abs).sum();
    assert!(
      nij_diff < 1e-1,
      "Dense-sparse nij total absolute difference should be small, got {nij_diff}"
    );

    // Ti: same direction, proportional magnitudes
    // Measured max rel_diff: 4.36e-3 (Ti[0], state A absorbs the single mutation)
    for k in 0..4 {
      let d = dense.Ti[k];
      let s = sparse.Ti[k];
      assert!(d > 0.0, "Dense Ti[{k}] should be positive, got {d}");
      assert!(s > 0.0, "Sparse Ti[{k}] should be positive, got {s}");
      let rel_diff = (d - s).abs() / s.max(1e-10);
      assert!(
        rel_diff < 1e-2,
        "Dense-sparse Ti[{k}] relative difference too large: dense={d}, sparse={s}, rel_diff={rel_diff}"
      );
    }

    // root_state: both should agree on total and dominant state
    let dense_total = dense.root_state.sum();
    let sparse_total = sparse.root_state.sum();
    pretty_assert_ulps_eq!(dense_total, sparse_total, epsilon = 1e-9);

    let dense_argmax = dense
      .root_state
      .iter()
      .copied()
      .enumerate()
      .max_by(|a, b| a.1.total_cmp(&b.1))
      .map(|(i, _)| i);
    let sparse_argmax = sparse
      .root_state
      .iter()
      .copied()
      .enumerate()
      .max_by(|a, b| a.1.total_cmp(&b.1))
      .map(|(i, _)| i);
    assert_eq!(
      dense_argmax, sparse_argmax,
      "Dense and sparse should agree on dominant root state"
    );

    Ok(())
  }

  /// Dense root_state reflects the ancestral composition.
  ///
  /// All leaves are mostly A (7/8 positions). The root reconstructs as predominantly A.
  /// root_state[A] should be the largest entry.
  #[test]
  fn test_root_state_dense() -> Result<(), Report> {
    // Leaves differ only at position 0: one has C, rest have A.
    // Root should reconstruct as all-A (position 0 resolves to A by majority).
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      AAAAAAAA
      >B
      AAAAAAAA
      >C
      AAAAAAAA
      >D
      CAAAAAAA
      "#},
      &*NUC_ALPHABET,
    )?;

    let (graph, partition) = setup_dense("((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;", &aln)?;
    let counts = get_mutation_counts_dense(&graph, &partition)?;

    // root_state[A] should dominate: 8 positions, root reconstructed as mostly A
    assert!(
      counts.root_state[IDX_A] >= 7.0,
      "root_state[A] should be >= 7 (out of 8 positions), got {}",
      counts.root_state[IDX_A]
    );
    // Total should equal alignment length
    pretty_assert_ulps_eq!(8.0, counts.root_state.sum(), epsilon = 1e-9);

    Ok(())
  }

  /// Sparse root_state reflects the Fitch consensus composition.
  ///
  /// Same alignment as dense test. Fitch consensus at the root should be all-A
  /// (3/4 leaves have A at position 0, majority wins).
  #[test]
  fn test_root_state_sparse() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      AAAAAAAA
      >B
      AAAAAAAA
      >C
      AAAAAAAA
      >D
      CAAAAAAA
      "#},
      &*NUC_ALPHABET,
    )?;

    let (graph, partition) = setup_sparse("((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;", &aln)?;
    let counts = get_mutation_counts_sparse(&graph, &partition)?;

    // Sparse root_state comes from Fitch composition counts.
    // All 8 positions should be A at root (3/4 leaves have A at pos 0).
    pretty_assert_ulps_eq!(8.0, counts.root_state[IDX_A], epsilon = 1e-9);
    pretty_assert_ulps_eq!(0.0, counts.root_state[IDX_C], epsilon = 1e-9);
    pretty_assert_ulps_eq!(0.0, counts.root_state[IDX_G], epsilon = 1e-9);
    pretty_assert_ulps_eq!(0.0, counts.root_state[IDX_T], epsilon = 1e-9);

    Ok(())
  }

  /// nij orientation verified across multiple mutation types.
  ///
  /// Alignment with three distinct mutations at separate positions:
  /// A->G (position 0), C->T (position 1), G->A (position 2).
  /// Each mutation should appear in the correct nij cell.
  #[test]
  fn test_nij_orientation_multiple_mutations_sparse() -> Result<(), Report> {
    // Two leaves: leaf_ref has the ancestral states, leaf_mut has mutations.
    // With equal branch lengths, Fitch places the root at leaf_ref states
    // (2 leaves, ties broken by parsimony with outgroup = majority).
    // Using 3 extra leaves to anchor the root at the "ref" state.
    let aln = read_many_fasta_str(
      indoc! {r#"
      >ref1
      ACGAAAAA
      >ref2
      ACGAAAAA
      >ref3
      ACGAAAAA
      >mut1
      GTAAAAAA
      "#},
      &*NUC_ALPHABET,
    )?;

    let tree_nwk = "((ref1:0.1,ref2:0.1)R12:0.05,(ref3:0.1,mut1:0.1)R3M:0.05)root:0.0;";
    let (graph, partition) = setup_sparse(tree_nwk, &aln)?;
    let counts = get_mutation_counts_sparse(&graph, &partition)?;

    // Position 0: A->G on the mut1 branch. nij[G, A] += 1
    assert!(
      counts.nij[[IDX_G, IDX_A]] >= 1.0,
      "Expected A->G mutation in nij[G,A], got {}",
      counts.nij[[IDX_G, IDX_A]]
    );

    // Position 1: C->T on the mut1 branch. nij[T, C] += 1
    assert!(
      counts.nij[[IDX_T, IDX_C]] >= 1.0,
      "Expected C->T mutation in nij[T,C], got {}",
      counts.nij[[IDX_T, IDX_C]]
    );

    // Position 2: G->A on the mut1 branch. nij[A, G] += 1
    assert!(
      counts.nij[[IDX_A, IDX_G]] >= 1.0,
      "Expected G->A mutation in nij[A,G], got {}",
      counts.nij[[IDX_A, IDX_G]]
    );

    // C->A and T->C reverses are zero (no such mutations in input).
    // A->G and G->A are mutual reverses: both nij[G,A] and nij[A,G] are non-zero.
    pretty_assert_ulps_eq!(0.0, counts.nij[[IDX_A, IDX_C]], epsilon = 1e-9); // nij[A, C] = 0
    pretty_assert_ulps_eq!(0.0, counts.nij[[IDX_C, IDX_T]], epsilon = 1e-9); // nij[C, T] = 0

    Ok(())
  }

  /// Dense nij accumulates across edges correctly.
  ///
  /// 6-leaf tree with 4 ancestral (A) and 2 mutant (T) leaves, one mutant per
  /// subtree. Root is unambiguously A at position 0 (4:2 majority). Two
  /// independent A->T mutations contribute to nij[T, A] from separate edges.
  #[test]
  fn test_nij_accumulation_dense() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >a1
      AAAAAAAA
      >a2
      AAAAAAAA
      >t1
      TAAAAAAA
      >a3
      AAAAAAAA
      >a4
      AAAAAAAA
      >t2
      TAAAAAAA
      "#},
      &*NUC_ALPHABET,
    )?;

    let tree_nwk = "((a1:0.1,a2:0.1,t1:0.1)left:0.05,(a3:0.1,a4:0.1,t2:0.1)right:0.05)root:0.0;";
    let (graph, partition) = setup_dense(tree_nwk, &aln)?;
    let counts = get_mutation_counts_dense(&graph, &partition)?;

    // nij[T, A] should reflect two A->T mutations (one per subtree)
    assert!(
      counts.nij[[IDX_T, IDX_A]] > 1.0,
      "nij[T,A] should reflect two A->T mutations, got {}",
      counts.nij[[IDX_T, IDX_A]]
    );

    // nij[A, T] should be much smaller (no T->A mutations)
    assert!(
      counts.nij[[IDX_T, IDX_A]] > 3.0 * counts.nij[[IDX_A, IDX_T]],
      "nij[T,A] ({}) should be >> nij[A,T] ({})",
      counts.nij[[IDX_T, IDX_A]],
      counts.nij[[IDX_A, IDX_T]]
    );

    Ok(())
  }

  /// Ti is proportional across states for uniform composition.
  ///
  /// Alignment with equal counts of each nucleotide. With identical sequences
  /// (no mutations), Ti should be approximately equal for all states since
  /// each state occupies the same fraction of positions.
  #[test]
  fn test_ti_proportional_to_composition_sparse() -> Result<(), Report> {
    // Each state appears exactly twice across 8 positions
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      AACCGGTT
      >B
      AACCGGTT
      >C
      AACCGGTT
      >D
      AACCGGTT
      "#},
      &*NUC_ALPHABET,
    )?;

    let tree_nwk = "((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;";
    let (graph, partition) = setup_sparse(tree_nwk, &aln)?;
    let counts = get_mutation_counts_sparse(&graph, &partition)?;

    // All Ti values should be equal (uniform composition, no mutations)
    pretty_assert_ulps_eq!(counts.Ti[IDX_A], counts.Ti[IDX_C], epsilon = 1e-9);
    pretty_assert_ulps_eq!(counts.Ti[IDX_C], counts.Ti[IDX_G], epsilon = 1e-9);
    pretty_assert_ulps_eq!(counts.Ti[IDX_G], counts.Ti[IDX_T], epsilon = 1e-9);

    Ok(())
  }

  /// Dense-sparse nij agreement on asymmetric mutations.
  ///
  /// Tree with a single clear mutation direction (all leaves on one side
  /// have C at position 0, all on the other have A). Both paths should
  /// agree on which nij cell receives the mutation count.
  #[test]
  fn test_dense_sparse_nij_direction_agreement() -> Result<(), Report> {
    // Clear asymmetry: clade AB has A at pos 0, clade CD has C at pos 0.
    // Root is ambiguous but parsimony and marginal should agree on mutation direction.
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGTACGT
      >B
      ACGTACGT
      >C
      CCGTACGT
      >D
      CCGTACGT
      "#},
      &*NUC_ALPHABET,
    )?;

    let tree_nwk = "((A:0.05,B:0.05)AB:0.1,(C:0.05,D:0.05)CD:0.1)root:0.0;";

    let (graph_d, partition_d) = setup_dense(tree_nwk, &aln)?;
    let (graph_s, partition_s) = setup_sparse(tree_nwk, &aln)?;

    let dense = get_mutation_counts_dense(&graph_d, &partition_d)?;
    let sparse = get_mutation_counts_sparse(&graph_s, &partition_s)?;

    // Both should have their largest off-diagonal nij entry in the same cell
    let dense_max_cell = dense
      .nij
      .indexed_iter()
      .filter(|((i, j), _)| i != j)
      .max_by(|a, b| a.1.total_cmp(b.1))
      .map(|((i, j), _)| (i, j));

    let sparse_max_cell = sparse
      .nij
      .indexed_iter()
      .filter(|((i, j), _)| i != j)
      .max_by(|a, b| a.1.total_cmp(b.1))
      .map(|((i, j), _)| (i, j));

    assert_eq!(
      dense_max_cell, sparse_max_cell,
      "Dense and sparse should agree on dominant mutation direction"
    );

    Ok(())
  }

  /// Root state total equals alignment length for dense path.
  #[test]
  fn test_root_state_total_equals_alignment_length_dense() -> Result<(), Report> {
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

    let (graph, partition) = setup_dense("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;", &aln)?;
    let counts = get_mutation_counts_dense(&graph, &partition)?;

    // root_state sums to alignment length (one count per position)
    pretty_assert_ulps_eq!(14.0, counts.root_state.sum(), epsilon = 1e-9);

    Ok(())
  }

  /// nij diagonal is zero after accumulation.
  ///
  /// No-change events (parent == child) are excluded from mutation counts.
  /// Both dense and sparse should zero the diagonal.
  #[test]
  fn test_nij_diagonal_zero() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACATCGCC
      >B
      GCATCCCT
      >C
      CCGGCGAT
      >D
      TCGGCCGT
      "#},
      &*NUC_ALPHABET,
    )?;

    let tree_nwk = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    let (graph_d, partition_d) = setup_dense(tree_nwk, &aln)?;
    let (graph_s, partition_s) = setup_sparse(tree_nwk, &aln)?;

    let dense = get_mutation_counts_dense(&graph_d, &partition_d)?;
    let sparse = get_mutation_counts_sparse(&graph_s, &partition_s)?;

    for k in 0..4 {
      pretty_assert_ulps_eq!(0.0, dense.nij[[k, k]], epsilon = 1e-15);
      pretty_assert_ulps_eq!(0.0, sparse.nij[[k, k]], epsilon = 1e-15);
    }

    Ok(())
  }

  /// nij is non-negative everywhere.
  #[test]
  fn test_nij_non_negative() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACATCGCC
      >B
      GCATCCCT
      >C
      CCGGCGAT
      >D
      TCGGCCGT
      "#},
      &*NUC_ALPHABET,
    )?;

    let tree_nwk = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    let (graph_d, partition_d) = setup_dense(tree_nwk, &aln)?;
    let (graph_s, partition_s) = setup_sparse(tree_nwk, &aln)?;

    let dense = get_mutation_counts_dense(&graph_d, &partition_d)?;
    let sparse = get_mutation_counts_sparse(&graph_s, &partition_s)?;

    for ((i, j), &v) in dense.nij.indexed_iter() {
      assert!(v >= 0.0, "Dense nij[{i},{j}] should be non-negative, got {v}");
    }
    for ((i, j), &v) in sparse.nij.indexed_iter() {
      assert!(v >= 0.0, "Sparse nij[{i},{j}] should be non-negative, got {v}");
    }

    Ok(())
  }

  /// Ti is non-negative everywhere.
  #[test]
  fn test_ti_non_negative() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACATCGCC
      >B
      GCATCCCT
      >C
      CCGGCGAT
      >D
      TCGGCCGT
      "#},
      &*NUC_ALPHABET,
    )?;

    let tree_nwk = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    let (graph_d, partition_d) = setup_dense(tree_nwk, &aln)?;
    let (graph_s, partition_s) = setup_sparse(tree_nwk, &aln)?;

    let dense = get_mutation_counts_dense(&graph_d, &partition_d)?;
    let sparse = get_mutation_counts_sparse(&graph_s, &partition_s)?;

    for (k, &v) in dense.Ti.iter().enumerate() {
      assert!(v >= 0.0, "Dense Ti[{k}] should be non-negative, got {v}");
    }
    for (k, &v) in sparse.Ti.iter().enumerate() {
      assert!(v >= 0.0, "Sparse Ti[{k}] should be non-negative, got {v}");
    }

    Ok(())
  }
}

#[cfg(test)]
#[allow(
  clippy::float_cmp,
  reason = "integer-valued f64s from += 1.0 accumulation and explicit = 0.0 assignment"
)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::gtr_inference::get_mutation_counts_fitch;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::fitch::partition::PartitionFitch;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::shared::update::MarginalPasses;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use pretty_assertions::assert_eq;

  use eyre::Report;
  use indoc::indoc;
  use std::sync::LazyLock;
  use treetime_graph::graph::Graph;
  use treetime_utils::{
    pretty_assert_array_diag_abs, pretty_assert_array_nonneg, pretty_assert_array_positive, pretty_assert_ulps_eq,
  };

  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  const IDX_A: usize = 0;
  const IDX_C: usize = 1;
  const IDX_G: usize = 2;
  const IDX_T: usize = 3;

  fn setup_dense(
    tree_nwk: &str,
    aln: &[AlignmentRecord],
  ) -> Result<(Graph, DenseReconstruction, BTreeMap<GraphEdgeKey, Option<f64>>), Report> {
    let nwk_parsed = nwk_read_str(tree_nwk)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let partition = PartitionMarginalDense::new(0, alphabet, get_common_length(aln)?);
    let node_states = partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, aln.to_vec()))?;
    let recon = DenseReconstruction::seeded(partition, gtr, node_states);
    let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    Ok((graph, recon, branch_lengths))
  }

  fn setup_sparse(
    tree_nwk: &str,
    aln: &[AlignmentRecord],
  ) -> Result<(Graph, PartitionFitch, BTreeMap<GraphEdgeKey, Option<f64>>), Report> {
    let nwk_parsed = nwk_read_str(tree_nwk)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let alphabet = Alphabet::default();
    let fitch = create_fitch_partition(&graph, 0, alphabet, &node_seq_inputs(&graph, &names, aln.to_vec()))?;
    Ok((graph, fitch, branch_lengths))
  }

  #[test]
  fn test_nij_orientation_dense() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let (graph, partition, branch_lengths) = setup_dense(
      "((ref1:0.1,ref2:0.1)R12:0.05,(ref3:0.1,mut1:0.1)R3M:0.05)root:0.0;",
      &aln,
    )?;
    let counts = partition.partition.count_transitions(
      &partition.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &partition.node_states,
      &partition.edges.backward,
      &partition.edges.forward,
    )?;

    assert!(
      counts.nij[[IDX_C, IDX_A]] > 0.1,
      "nij[C,A] should carry A->C mutation signal, got {}",
      counts.nij[[IDX_C, IDX_A]]
    );
    assert!(
      counts.nij[[IDX_C, IDX_A]] > 5.0 * counts.nij[[IDX_A, IDX_C]],
      "nij[C,A] ({}) should be >> nij[A,C] ({}): child-row, parent-column convention",
      counts.nij[[IDX_C, IDX_A]],
      counts.nij[[IDX_A, IDX_C]]
    );

    Ok(())
  }

  #[test]
  fn test_nij_orientation_sparse() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
      indoc! {r#"
      >leaf_a
      AAAAAAAA
      >leaf_c
      CAAAAAAA
      "#},
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let (graph, fitch, branch_lengths) = setup_sparse("(leaf_a:0.1,leaf_c:0.1)root:0.0;", &aln)?;
    let counts = get_mutation_counts_fitch(&graph, &fitch, &branch_lengths_or_zero(&branch_lengths))?;

    assert_eq!(1.0, counts.nij[[IDX_C, IDX_A]]);
    assert_eq!(0.0, counts.nij[[IDX_A, IDX_C]]);

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::baseline_vs_double((0.1, 0.2))]
  #[case::small_vs_large(    (0.05, 0.5))]
  #[trace]
  fn test_ti_scaling_sparse(#[case] (bl1, bl2): (f64, f64)) -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let tree1 = format!("((A:{bl1},B:{bl1})AB:{bl1},(C:{bl1},D:{bl1})CD:{bl1})root:0.0;");
    let tree2 = format!("((A:{bl2},B:{bl2})AB:{bl2},(C:{bl2},D:{bl2})CD:{bl2})root:0.0;");

    let (graph1, fitch1, branch_lengths1) = setup_sparse(&tree1, &aln)?;
    let (graph2, fitch2, branch_lengths2) = setup_sparse(&tree2, &aln)?;

    let counts1 = get_mutation_counts_fitch(&graph1, &fitch1, &branch_lengths_or_zero(&branch_lengths1))?;
    let counts2 = get_mutation_counts_fitch(&graph2, &fitch2, &branch_lengths_or_zero(&branch_lengths2))?;

    let ratio = bl2 / bl1;
    pretty_assert_ulps_eq!(counts2.Ti, &counts1.Ti * ratio, epsilon = 1e-7);

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::baseline_vs_double((0.1, 0.2))]
  #[case::small_vs_large(    (0.05, 0.5))]
  #[trace]
  fn test_ti_scaling_dense(#[case] (bl1, bl2): (f64, f64)) -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let tree1 = format!("((A:{bl1},B:{bl1})AB:{bl1},(C:{bl1},D:{bl1})CD:{bl1})root:0.0;");
    let tree2 = format!("((A:{bl2},B:{bl2})AB:{bl2},(C:{bl2},D:{bl2})CD:{bl2})root:0.0;");

    let (graph1, partition1, branch_lengths1) = setup_dense(&tree1, &aln)?;
    let (graph2, partition2, branch_lengths2) = setup_dense(&tree2, &aln)?;

    let counts1 = partition1.partition.count_transitions(&partition1.gtr, &graph1, &branch_lengths_or_zero(&branch_lengths1), &partition1.node_states, &partition1.edges.backward, &partition1.edges.forward)?;
    let counts2 = partition2.partition.count_transitions(&partition2.gtr, &graph2, &branch_lengths_or_zero(&branch_lengths2), &partition2.node_states, &partition2.edges.backward, &partition2.edges.forward)?;

    let ratio = bl2 / bl1;
    pretty_assert_ulps_eq!(counts2.Ti, &counts1.Ti * ratio, epsilon = 1e-7);

    Ok(())
  }

  #[test]
  fn test_dense_sparse_consistency() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let tree_nwk = "((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;";

    let (graph_d, partition_d, branch_lengths_d) = setup_dense(tree_nwk, &aln)?;
    let (graph_s, fitch_s, branch_lengths_s) = setup_sparse(tree_nwk, &aln)?;

    let dense = partition_d.partition.count_transitions(
      &partition_d.gtr,
      &graph_d,
      &branch_lengths_or_zero(&branch_lengths_d),
      &partition_d.node_states,
      &partition_d.edges.backward,
      &partition_d.edges.forward,
    )?;
    let sparse = get_mutation_counts_fitch(&graph_s, &fitch_s, &branch_lengths_or_zero(&branch_lengths_s))?;

    let nij_diff = (&dense.nij - &sparse.nij).mapv(f64::abs).sum();
    assert!(
      nij_diff < 1e-1,
      "Dense-sparse nij total absolute difference should be small, got {nij_diff}"
    );

    pretty_assert_array_positive!(dense.Ti);
    pretty_assert_array_positive!(sparse.Ti);
    approx::assert_relative_eq!(dense.Ti, sparse.Ti, max_relative = 1e-2);

    let dense_total = dense.root_state.sum();
    let sparse_total = sparse.root_state.sum();
    assert_eq!(dense_total, sparse_total);

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

  #[test]
  fn test_root_state_dense() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let (graph, partition, branch_lengths) = setup_dense("((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;", &aln)?;
    let counts = partition.partition.count_transitions(
      &partition.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &partition.node_states,
      &partition.edges.backward,
      &partition.edges.forward,
    )?;

    assert!(
      counts.root_state[IDX_A] >= 7.0,
      "root_state[A] should be >= 7 (out of 8 positions), got {}",
      counts.root_state[IDX_A]
    );
    assert_eq!(8.0, counts.root_state.sum());

    Ok(())
  }

  #[test]
  fn test_root_state_sparse() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let (graph, fitch, branch_lengths) = setup_sparse("((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;", &aln)?;
    let counts = get_mutation_counts_fitch(&graph, &fitch, &branch_lengths_or_zero(&branch_lengths))?;

    assert_eq!(8.0, counts.root_state[IDX_A]);
    assert_eq!(0.0, counts.root_state[IDX_C]);
    assert_eq!(0.0, counts.root_state[IDX_G]);
    assert_eq!(0.0, counts.root_state[IDX_T]);

    Ok(())
  }

  #[test]
  fn test_nij_orientation_multiple_mutations_sparse() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let tree_nwk = "((ref1:0.1,ref2:0.1)R12:0.05,(ref3:0.1,mut1:0.1)R3M:0.05)root:0.0;";
    let (graph, fitch, branch_lengths) = setup_sparse(tree_nwk, &aln)?;
    let counts = get_mutation_counts_fitch(&graph, &fitch, &branch_lengths_or_zero(&branch_lengths))?;

    assert!(
      counts.nij[[IDX_G, IDX_A]] >= 1.0,
      "Expected A->G mutation in nij[G,A], got {}",
      counts.nij[[IDX_G, IDX_A]]
    );

    assert!(
      counts.nij[[IDX_T, IDX_C]] >= 1.0,
      "Expected C->T mutation in nij[T,C], got {}",
      counts.nij[[IDX_T, IDX_C]]
    );

    assert!(
      counts.nij[[IDX_A, IDX_G]] >= 1.0,
      "Expected G->A mutation in nij[A,G], got {}",
      counts.nij[[IDX_A, IDX_G]]
    );

    assert_eq!(0.0, counts.nij[[IDX_A, IDX_C]]);
    assert_eq!(0.0, counts.nij[[IDX_C, IDX_T]]);

    Ok(())
  }

  #[test]
  fn test_nij_accumulation_dense() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let tree_nwk = "((a1:0.1,a2:0.1,t1:0.1)left:0.05,(a3:0.1,a4:0.1,t2:0.1)right:0.05)root:0.0;";
    let (graph, partition, branch_lengths) = setup_dense(tree_nwk, &aln)?;
    let counts = partition.partition.count_transitions(
      &partition.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &partition.node_states,
      &partition.edges.backward,
      &partition.edges.forward,
    )?;

    assert!(
      counts.nij[[IDX_T, IDX_A]] > 1.0,
      "nij[T,A] should reflect two A->T mutations, got {}",
      counts.nij[[IDX_T, IDX_A]]
    );

    assert!(
      counts.nij[[IDX_T, IDX_A]] > 3.0 * counts.nij[[IDX_A, IDX_T]],
      "nij[T,A] ({}) should be >> nij[A,T] ({})",
      counts.nij[[IDX_T, IDX_A]],
      counts.nij[[IDX_A, IDX_T]]
    );

    Ok(())
  }

  #[test]
  fn test_ti_proportional_to_composition_sparse() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let tree_nwk = "((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;";
    let (graph, fitch, branch_lengths) = setup_sparse(tree_nwk, &aln)?;
    let counts = get_mutation_counts_fitch(&graph, &fitch, &branch_lengths_or_zero(&branch_lengths))?;

    assert_eq!(counts.Ti[IDX_A], counts.Ti[IDX_C]);
    assert_eq!(counts.Ti[IDX_C], counts.Ti[IDX_G]);
    assert_eq!(counts.Ti[IDX_G], counts.Ti[IDX_T]);

    Ok(())
  }

  #[test]
  fn test_dense_sparse_nij_direction_agreement() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let tree_nwk = "((A:0.05,B:0.05)AB:0.1,(C:0.05,D:0.05)CD:0.1)root:0.0;";

    let (graph_d, partition_d, branch_lengths_d) = setup_dense(tree_nwk, &aln)?;
    let (graph_s, fitch_s, branch_lengths_s) = setup_sparse(tree_nwk, &aln)?;

    let dense = partition_d.partition.count_transitions(
      &partition_d.gtr,
      &graph_d,
      &branch_lengths_or_zero(&branch_lengths_d),
      &partition_d.node_states,
      &partition_d.edges.backward,
      &partition_d.edges.forward,
    )?;
    let sparse = get_mutation_counts_fitch(&graph_s, &fitch_s, &branch_lengths_or_zero(&branch_lengths_s))?;

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

  #[test]
  fn test_root_state_total_equals_alignment_length_dense() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let (graph, partition, branch_lengths) =
      setup_dense("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;", &aln)?;
    let counts = partition.partition.count_transitions(
      &partition.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &partition.node_states,
      &partition.edges.backward,
      &partition.edges.forward,
    )?;

    assert_eq!(14.0, counts.root_state.sum());

    Ok(())
  }

  #[test]
  fn test_nij_diagonal_zero() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let tree_nwk = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    let (graph_d, partition_d, branch_lengths_d) = setup_dense(tree_nwk, &aln)?;
    let (graph_s, fitch_s, branch_lengths_s) = setup_sparse(tree_nwk, &aln)?;

    let dense = partition_d.partition.count_transitions(
      &partition_d.gtr,
      &graph_d,
      &branch_lengths_or_zero(&branch_lengths_d),
      &partition_d.node_states,
      &partition_d.edges.backward,
      &partition_d.edges.forward,
    )?;
    let sparse = get_mutation_counts_fitch(&graph_s, &fitch_s, &branch_lengths_or_zero(&branch_lengths_s))?;

    pretty_assert_array_diag_abs!(dense.nij, epsilon = 1e-15);
    pretty_assert_array_diag_abs!(sparse.nij, epsilon = 1e-15);

    Ok(())
  }

  #[test]
  fn test_nij_non_negative() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let tree_nwk = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    let (graph_d, partition_d, branch_lengths_d) = setup_dense(tree_nwk, &aln)?;
    let (graph_s, fitch_s, branch_lengths_s) = setup_sparse(tree_nwk, &aln)?;

    let dense = partition_d.partition.count_transitions(
      &partition_d.gtr,
      &graph_d,
      &branch_lengths_or_zero(&branch_lengths_d),
      &partition_d.node_states,
      &partition_d.edges.backward,
      &partition_d.edges.forward,
    )?;
    let sparse = get_mutation_counts_fitch(&graph_s, &fitch_s, &branch_lengths_or_zero(&branch_lengths_s))?;

    pretty_assert_array_nonneg!(dense.nij);
    pretty_assert_array_nonneg!(sparse.nij);

    Ok(())
  }

  #[test]
  fn test_ti_non_negative() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let tree_nwk = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    let (graph_d, partition_d, branch_lengths_d) = setup_dense(tree_nwk, &aln)?;
    let (graph_s, fitch_s, branch_lengths_s) = setup_sparse(tree_nwk, &aln)?;

    let dense = partition_d.partition.count_transitions(
      &partition_d.gtr,
      &graph_d,
      &branch_lengths_or_zero(&branch_lengths_d),
      &partition_d.node_states,
      &partition_d.edges.backward,
      &partition_d.edges.forward,
    )?;
    let sparse = get_mutation_counts_fitch(&graph_s, &fitch_s, &branch_lengths_or_zero(&branch_lengths_s))?;

    pretty_assert_array_nonneg!(dense.Ti);
    pretty_assert_array_nonneg!(sparse.Ti);

    Ok(())
  }
}
